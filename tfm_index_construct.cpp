#include <algorithm>
#include <cstdint>
#include <ctime>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>
#include <sstream>
#include <stdexcept>
#include <stdlib.h>
#include <string>
#include <sys/types.h>
#include <utility>
#include <vector>

#include <stdint.h>
#include <stdio.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <semaphore.h>
#include <errno.h>
#include <assert.h>

#include <sdsl/util.hpp>
#include "tfm_index.hpp"
#include "dbg_algorithms.hpp"

extern "C" {
#include "gsacak.c"
#include "utils.c"
}

using namespace std;
using namespace sdsl;
using namespace __gnu_cxx;

// =============== algorithm limits ===================
// maximum number of distinct words
#define MAX_DISTINCT_WORDS (INT32_MAX - 1)
// typedef uint32_t word_int_t;
// maximum number of occurrences of a single word
#define MAX_WORD_OCC (UINT32_MAX)
// typedef uint32_t occ_int_t;

// clang-format off
uint8_t asc2dnacat[] = {
    /*   0 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  16 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  32 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0,
           /*                                        -     */
    /*  48 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  64 */ 0, 1, 2, 1, 2, 0, 0, 1, 2, 0, 0, 2, 0, 2, 2, 0,
           /*    A  B  C  D        G  H        K     M  N  */
    /*  80 */ 0, 0, 2, 2, 1, 0, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0,
           /*       R  S  T     V  W  X  Y */
    /*  96 */ 0, 1, 2, 1, 2, 0, 0, 1, 2, 0, 0, 2, 0, 2, 2, 0,
           /*    a  b  c  d        g  h        k     m  n  */
    /* 112 */ 0, 0, 2, 2, 1, 0, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0,
           /*       r  s  t     v  w  x  y                 */
    /* 128 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 144 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 160 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 176 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 192 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 208 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 224 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 240 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};
// clang-format on

struct word_stats {
    string str;             // word
    uint32_t occ;           // its number of occurences
    uint32_t rank = 0;      // its rank
};

struct KR_window {
    int wsize;
    int *window;
    int asize;
    const uint64_t prime = 1999999973;
    uint64_t hash;
    uint64_t tot_char;
    uint64_t asize_pot; // asize^(wsize-1) mod prime

    KR_window(int w) : wsize(w) {
        asize = 256;
        asize_pot = 1;
        for (int i = 1; i < wsize; i++)
            asize_pot =
                (asize_pot * asize) % prime; // ugly linear-time power algorithm
        // alloc and clear window
        window = new int[wsize];
        reset();
    }

    void reset() {
        for (int i = 0; i < wsize; i++)
            window[i] = 0;
        hash = tot_char = 0;
    }

    uint64_t addchar(int c) {
        int k = tot_char++ % wsize;
        // complex expression to avoid negative numbers
        hash += (prime - (window[k] * asize_pot) % prime);  // remove window[k] contribution
        hash = (asize * hash + c) % prime;                  //  add char i
        window[k] = c;
        return hash;
    }

    ~KR_window() { delete[] window; }
};
// -----------------------------------------------------------

uint64_t kr_hash(string s) {
    uint64_t hash = 0;
    const uint64_t prime = 27162335252586509; // next prime (2**54 + 2**53 + 2**47 + 2**13)
    for (size_t k = 0; k < s.size(); k++) {
        int c = (unsigned char)s[k];
        assert(c >= 0 && c < 256);
        hash = (256 * hash + c) % prime; //  add char k
    }
    return hash;
}

static void save_update_word(string &w, unsigned int minsize, map<uint64_t, word_stats> &freq, vector<uint64_t> &parse, vector<char> &last, uint64_t &pos) {
    assert(pos == 0 || w.size() > minsize);
    if (w.size() <= minsize)
        return;
    // get the hash value and write it to the temporary parse file
    uint64_t hash = kr_hash(w);
    parse.push_back(hash);

    // update frequency table for current hash
    if (freq.find(hash) == freq.end()) {
        freq[hash].occ = 1; // new hash
        freq[hash].str = w;
    } else {
        freq[hash].occ += 1; // known hash
        if (freq[hash].occ <= 0) {
            cerr << "Emergency exit! Maximum # of occurence of dictionary word "
                    "(";
            cerr << MAX_WORD_OCC << ") exceeded\n";
            exit(1);
        }
        if (freq[hash].str != w) {
            cerr << "Emergency exit! Hash collision for strings:\n";
            cerr << freq[hash].str << "\n  vs\n" << w << endl;
            exit(1);
        }
    }
    // output char w+1 from the end
    last.push_back(w[w.size() - minsize - 1]);
    // compute ending position +1 of current word and write it to sa file
    // pos is the ending position+1 of the previous word and is updated here
    if (pos == 0)
        pos = w.size() - 1; // -1 is for the initial $ of the first word
    else
        pos += w.size() - minsize;
    // keep only the overlapping part of the window
    w.erase(0, w.size() - minsize);
}

uint64_t process_file(string &filename, size_t w, size_t p, map<uint64_t, word_stats> &wordFreq, vector<uint64_t> &g_vec, vector<char> &last_vec) {
    uint64_t pos = 0;
    assert( IBYTES <= sizeof(pos) );
    string word("");
    word.append(1, Dollar);
    KR_window krw(w);
    std::string line;
    ifstream f(filename);
    if (!f.rdbuf()->is_open()) { // is_open does not work on igzstreams
        perror(__func__);
        throw new std::runtime_error("Cannot open input file " + filename);
    }
    int c;
    while ((c = f.get()) != EOF) {
        if (c <= Dollar) {
            cerr << "Invalid char found in input file: no additional chars "
                    "will be read\n";
            break;
        }
        word.append(1, c);
        uint64_t hash = krw.addchar(c);
        if (hash % p == 0) {
            save_update_word(word, w, wordFreq, g_vec, last_vec, pos);
        }
    }
    f.close();
    word.append(w, Dollar);
    save_update_word(word, w, wordFreq, g_vec, last_vec, pos);

    return krw.tot_char;
}

bool pstringCompare(const string *a, const string *b) { return *a <= *b; }

void writeDictOcc(map<uint64_t, word_stats> &wfreq, vector<const string *> &sortedDict, vector<char> &dict) {
    assert(sortedDict.size() == wfreq.size());
    vector<uint32_t> vocc{};

    uint32_t wrank = 1; // current word rank (1 based)
    for (auto x : sortedDict) {
        const char *word = (*x).data(); // current dictionary word
        size_t len = (*x).size(); // offset and length of word
        // assert(len > (size_t)arg.w);
        for (size_t i = 0; i < len; i++) {
            dict.push_back(word[i]);
        }
        dict.push_back(EndOfWord);

        uint64_t hash = kr_hash(*x);
        struct word_stats &wf = wfreq.at(hash);
        assert(wf.occ > 0);
        vocc.push_back(wf.occ);

        assert(wf.rank == 0);
        wf.rank = wrank++;
    }
    dict.push_back(EndOfDict);
}

void calculate_word_frequencies(string &filename, size_t w, size_t p, map<uint64_t, word_stats> &wordFreq, vector<uint64_t> &parse, vector<char> &last) {
    try {
        process_file(filename, w, p, wordFreq, parse, last);
    } catch (const std::bad_alloc &) {
        cout << "Out of memory (parsing phase)... emergency exit\n";
        die("bad alloc exception");
    }


    if (wordFreq.size() > MAX_DISTINCT_WORDS) {
        cerr << "Emergency exit! The number of distinc words (" << wordFreq.size()
             << ")\n";
        cerr << "is larger than the current limit (" << MAX_DISTINCT_WORDS
             << ")\n";
        exit(1);
    }
}

struct Dict {
    uint8_t *d; // pointer to the dictionary
    uint64_t *end; // end[i] is the index of the ending symbol of the i-th phrase
    uint64_t dsize;  // dicionary size in symbols
    uint64_t dwords; // the number of phrases of the dicionary
};

// binary search for x in an array a[0..n-1] that doesn't contain x
// return the lowest position that is larger than x
static long binsearch(uint_t x, uint_t a[], long n) {
    long lo = 0;
    long hi = n - 1;
    while (hi > lo) {
        assert(((lo == 0) || x > a[lo - 1]) && x < a[hi]);
        int mid = (lo + hi) / 2;
        assert(x != a[mid]); // x is not in a[]
        if (x < a[mid])
            hi = mid;
        else
            lo = mid + 1;
    }
    assert(((hi == 0) || x > a[hi - 1]) && x < a[hi]);
    return hi;
}

// return the length of the suffix starting in position p.
// also write to seqid the id of the sequence containing that suffix
// n is the # of distinct words in the dictionary, hence the length of eos[]
static int_t getlen(uint_t p, uint_t eos[], long n, uint32_t *seqid) {
    assert(p < eos[n - 1]);
    *seqid = binsearch(p, eos, n);
    assert(eos[*seqid] > p); // distance between position p and the next $
    return eos[*seqid] - p;
}

static void compute_dict_bwt_lcp(uint8_t *d, long dsize, long dwords, int w, uint_t **sap, int_t **lcpp) {
    uint_t *sa = new uint_t[dsize];
    int_t *lcp = new int_t[dsize];
    (void)dwords;
    (void)w;

    gsacak(d, sa, lcp, NULL, dsize);
    assert(d[dsize - 1] == EndOfDict);
    assert(sa[0] == (unsigned long)dsize - 1); // sa[0] is the EndOfDict symbol
    for (long i = 0; i < dwords; i++)
        assert(d[sa[i + 1]] == EndOfWord); // there are dwords EndOfWord symbols
    // EndOfWord symbols are in position order, so the last is d[dsize-2]
    assert(sa[dwords] == (unsigned long)dsize - 2);
    // there are wsize+1 $ symbols:
    // one at the beginning of the first word, wsize at the end of the last word
    for (long i = 0; i <= w; i++)
        assert(d[sa[i + dwords + 1]] == Dollar);
    // in sa[dwords+w+1] we have the first word in the parsing since that $ is
    // the lex. larger
    assert(d[0] == Dollar);
    assert(sa[dwords + w + 1] == 0);
    assert(
        d[sa[dwords + w + 2]] > Dollar
    ); // end of Dollar chars in the first column
    assert(lcp[dwords + w + 2] == 0);
    // copy sa and lcp address
    *sap = sa;
    *lcpp = lcp;
}

struct SeqId {
    uint32_t id;   // lex. id of the dictionary word to which the suffix belongs
    int remaining; // remaining copies of the suffix to be considered
    uint32_t *bwtpos;   // list of bwt positions of this dictionary word
    uint8_t char2write; // char to be written (is the one preceeding the suffix)

    // constructor
    SeqId(uint32_t i, int r, uint32_t *b, int8_t c)
        : id(i), remaining(r), bwtpos(b) {
        char2write = c;
    }

    // advance to the next bwt position, return false if there are no more
    // positions
    bool next() {
        remaining--;
        bwtpos += 1;
        return remaining > 0;
    }
    bool operator<(const SeqId &a);
};

bool SeqId::operator<(const SeqId &a) { return *bwtpos > *(a.bwtpos); }

inline uint8_t get_prev(int w, uint8_t *d, uint64_t *end, uint32_t seqid) {
    return d[end[seqid] - w - 1];
}

void write_bitvector(FILE *f, bool bit, uint8_t &cnt, uint8_t &buffer, bool hard_write = false) {
    buffer |= (bit << (7 - cnt++));
    if (hard_write || (cnt == 8)) {
        if (fputc(buffer, f) == EOF)
            die("Din/Dout write error 0");
        cnt = 0;
        buffer = 0;
    }
}

vector<char> store_bwt(size_t w, uint8_t *d, long dsize, uint64_t *end_to_phrase, uint32_t *ilist, tfm_index &tfmp, long dwords, uint_t *sa, int_t *lcp) {
    // starting point in ilist for each word and # words
    // set d[0]==0 as this is the EOF char in the final BWT
    assert(d[0] == Dollar);
    d[0] = 0;

    // derive eos from sa. for i=0...dwords-1, eos[i] is the eos position of
    // string i in d
    uint_t *eos = sa + 1;
    for (int i = 0; i < dwords - 1; i++)
        assert(eos[i] < eos[i + 1]);

    vector<char> out{};

    // main loop: consider each entry in the SA of dict
    long full_words = 0, easy_bwts = 0, hard_bwts = 0, next;
    uint32_t seqid;
    for (long i = dwords + w + 1; i < dsize; i = next) {
        // we are considering d[sa[i]....]
        next = i + 1; // prepare for next iteration
        // compute length of this suffix and sequence it belongs
        int_t suffixLen = getlen(sa[i], eos, dwords, &seqid);
        //  ignore suffixes of lenght <= w
        if (suffixLen <= (int_t)w)
            continue;
        // ----- simple case: the suffix is a full word
        if (sa[i] == 0 || d[sa[i] - 1] == EndOfWord) {
            full_words++;
            uint32_t start = tfmp.C[seqid + 1], end = tfmp.C[seqid + 2];
            assert(tfmp.din[start] == 1);
            for (uint32_t j = start; j < end; j++) {
                if (tfmp.din[j] == 1) {
                    uint32_t pos = tfmp.dout_select(tfmp.din_rank(j + 1));
                    while (1) {
                        if (tfmp.L[pos] == 0)
                            pos = 0;
                        uint32_t act_phrase = tfmp.L[pos] - 1;
                        uint8_t char_to_write = get_prev(w, d, end_to_phrase, act_phrase);
                        easy_bwts++;
                        out.push_back(char_to_write);
                        // if (fputc(char_to_write, fbwt) == EOF)
                            // die("L write error 0");
                        if (tfmp.dout[++pos] == 1)
                            break;
                    }
                }
            }
            continue; // proceed with next i
        } else {
            // ----- hard case: there can be a group of equal suffixes starting
            // at i save seqid and the corresponding char
            vector<uint32_t> id2merge(1, seqid);
            vector<uint8_t> char2write(1, d[sa[i] - 1]);
            while (next < dsize && lcp[next] >= suffixLen) {
                assert(
                    lcp[next] == suffixLen
                ); // the lcp cannot be greater than suffixLen
                assert(
                    sa[next] > 0 && d[sa[next] - 1] != EndOfWord
                ); // sa[next] cannot be a full word
                int_t nextsuffixLen = getlen(sa[next], eos, dwords, &seqid);
                assert(nextsuffixLen >= suffixLen);
                if (nextsuffixLen == suffixLen) {
                    id2merge.push_back(seqid); // sequence to consider
                    char2write.push_back(d[sa[next] - 1]); // corresponding char
                    next++;
                } else
                    break;
            }
            // output to fbwt the bwt chars corresponding to the current
            // dictionary suffix, and, if requested, some SA values
            // fwrite_chars_same_suffix(id2merge, char2write, tfmp, ilist, fbwt, easy_bwts, hard_bwts);
            // static void fwrite_chars_same_suffix(vector<uint32_t> &id2merge, vector<uint8_t> &char2write, tfm_index &tfmp, uint32_t *ilist, FILE *fbwt, long &easy_bwts, long &hard_bwts) {
            size_t numwords = id2merge.size(); // numwords dictionary words contain the same suffix
            bool samechar = true;
            for (size_t i = 1; (i < numwords) && samechar; i++) {
                samechar = (char2write[i - 1] == char2write[i]);
            }

            if (samechar) {
                for (size_t i = 0; i < numwords; i++) {
                    uint32_t s = id2merge[i] + 1;
                    for (uint64_t j = tfmp.C[s]; j < tfmp.C[s + 1]; j++) {
                        out.push_back(char2write[0]);
                    }
                    easy_bwts += tfmp.C[s + 1] - tfmp.C[s];
                }
            } else {
                // many words, many chars...
                vector<SeqId> heap; // create heap
                for (size_t i = 0; i < numwords; i++) {
                    uint32_t s = id2merge[i] + 1;
                    heap.push_back(SeqId(
                        s, tfmp.C[s + 1] - tfmp.C[s], ilist + (tfmp.C[s] - 1),
                        char2write[i]
                    ));
                }
                std::make_heap(heap.begin(), heap.end());
                while (heap.size() > 0) {
                    // output char for the top of the heap
                    SeqId s = heap.front();
                    out.push_back(s.char2write);
                    hard_bwts += 1;
                    // remove top
                    pop_heap(heap.begin(), heap.end());
                    heap.pop_back();
                    // if remaining positions, reinsert to heap
                    if (s.next()) {
                        heap.push_back(s);
                        push_heap(heap.begin(), heap.end());
                    }
                }
            }
        }
    }
    assert(full_words == dwords);

    return out;
}

vector<bool> store_din(size_t w, uint8_t *d, long dsize, tfm_index &tfmp, long dwords, uint_t *sa, int_t *lcp) {
    // starting point in ilist for each word and # words

    // derive eos from sa. for i=0...dwords-1, eos[i] is the eos position of
    // string i in d
    uint_t *eos = sa + 1;
    for (int i = 0; i < dwords - 1; i++)
        assert(eos[i] < eos[i + 1]);

    vector<bool> din{};

    // main loop: consider each entry in the SA of dict
    long next;
    uint32_t seqid;
    for (long i = dwords + w + 1; i < dsize; i = next) {
        // we are considering d[sa[i]....]
        next = i + 1; // prepare for next iteration
        // compute length of this suffix and sequence it belongs
        int_t suffixLen = getlen(sa[i], eos, dwords, &seqid);
        // cout << suffixLen << " " << seqid << endl;
        //  ignore suffixes of lenght <= w
        if (suffixLen <= (int_t)w)
            continue;
        // ----- simple case: the suffix is a full word
        if (sa[i] == 0 || d[sa[i] - 1] == EndOfWord) {
            uint32_t start = tfmp.C[seqid + 1], end = tfmp.C[seqid + 2];
            assert(tfmp.din[start] == 1);
            for (uint32_t j = start; j < end; j++) {
                din.push_back(tfmp.din[j]);
                // write_bitvector(fdin, tfmp.din[j], cnt, buffer);
                // din.resize(din.capacity()*2);
            }
            continue; // proceed with next i
        } else {
            // ----- hard case: there can be a group of equal suffixes starting
            // at i save seqid and the corresponding char
            int bits_to_write = tfmp.C[seqid + 2] - tfmp.C[seqid + 1];
            while (next < dsize && lcp[next] >= suffixLen) {
                assert(
                    lcp[next] == suffixLen
                ); // the lcp cannot be greater than suffixLen
                assert(
                    sa[next] > 0 && d[sa[next] - 1] != EndOfWord
                ); // sa[next] cannot be a full word
                int_t nextsuffixLen = getlen(sa[next], eos, dwords, &seqid);
                assert(nextsuffixLen >= suffixLen);
                if (nextsuffixLen == suffixLen) {
                    bits_to_write += tfmp.C[seqid + 2] - tfmp.C[seqid + 1];
                    next++;
                } else
                    break;
            }
            for (int k = 0; k < bits_to_write; k++)
                din.push_back(1);
        }
    }
    din.push_back(1);
    return din;
}

vector<bool> store_dout(size_t w, uint8_t *d, long dsize, tfm_index &tfmp, long dwords, uint_t *sa, int_t *lcp) {
    // derive eos from sa. for i=0...dwords-1, eos[i] is the eos position of
    // string i in d
    uint_t *eos = sa + 1;
    for (int i = 0; i < dwords - 1; i++)
        assert(eos[i] < eos[i + 1]);

    vector<bool> dout{};

    // main loop: consider each entry in the SA of dict
    long next;
    uint32_t seqid;
    for (long i = dwords + w + 1; i < dsize; i = next) {
        // we are considering d[sa[i]....]
        next = i + 1; // prepare for next iteration
        // compute length of this suffix and sequence it belongs
        int_t suffixLen = getlen(sa[i], eos, dwords, &seqid);
        // cout << suffixLen << " " << seqid << endl;
        //  ignore suffixes of lenght <= w
        if (suffixLen <= (int_t)w)
            continue;
        // ----- simple case: the suffix is a full word
        if (sa[i] == 0 || d[sa[i] - 1] == EndOfWord) {
            uint32_t start = tfmp.C[seqid + 1], end = tfmp.C[seqid + 2];
            assert(tfmp.din[start] == 1);
            for (uint32_t j = start; j < end; j++) {
                if (tfmp.din[j] == 1) {
                    uint32_t pos = tfmp.dout_select(tfmp.din_rank(j + 1));
                    if (tfmp.L[pos] == 0)
                        pos = 0;
                    while (1) {
                        dout.push_back(tfmp.dout[pos]);
                        if (tfmp.dout[++pos] == 1)
                            break;
                    }
                }
            }
            continue; // proceed with next i
        } else {
            // ----- hard case: there can be a group of equal suffixes starting
            // at i save seqid and the corresponding char
            int bits_to_write = tfmp.C[seqid + 2] - tfmp.C[seqid + 1];
            while (next < dsize && lcp[next] >= suffixLen) {
                assert(
                    lcp[next] == suffixLen
                ); // the lcp cannot be greater than suffixLen
                assert(
                    sa[next] > 0 && d[sa[next] - 1] != EndOfWord
                ); // sa[next] cannot be a full word
                int_t nextsuffixLen = getlen(sa[next], eos, dwords, &seqid);
                assert(nextsuffixLen >= suffixLen);
                if (nextsuffixLen == suffixLen) {
                    bits_to_write += tfmp.C[seqid + 2] - tfmp.C[seqid + 1];
                    next++;
                } else
                    break;
            }
            for (int k = 0; k < bits_to_write; k++) {
                dout.push_back(1);
            }
        }
    }
    dout.push_back(1);
    return dout;
}

Dict read_dictionary(vector<char> &dict) {
    size_t dsize = dict.size();
    uint8_t *d = new uint8_t[dsize];
    for (size_t i = 0; i < dsize; i++) {
        d[i] = dict[i];
    }

    uint64_t dwords = 0;
    for (size_t i = 0; i < dsize; i++) {
        if (d[i] == EndOfWord)
            dwords++;
    }

    uint64_t *end = new uint64_t[dwords];
    int cnt = 0;
    for (size_t i = 0; i < dsize; i++) {
        if (d[i] == EndOfWord)
            end[cnt++] = i;
    }

    Dict res = {d, end, (uint64_t)dsize, dwords};
    return res;
}

void generate_ilist(uint32_t *ilist, tfm_index &tfmp, uint64_t dwords) {
    vector<vector<uint32_t>> phrase_sources(dwords);
    for (uint64_t i = 0; i < tfmp.L.size(); i++) {
        uint32_t act_char = tfmp.L[i];
        if (act_char == 0)
            continue;
        phrase_sources[act_char - 1].push_back(i);
    }
    uint64_t cnt = 0;
    for (uint64_t i = 0; i < phrase_sources.size(); i++) {
        for (int j = 0; j < (int)phrase_sources[i].size(); j++)
            ilist[cnt++] = phrase_sources[i][j];
    }
}

tfm_index create_tfm(size_t size, sdsl::int_vector_buffer<> &L_buf, sdsl::bit_vector &din, sdsl::bit_vector &dout) {
    tfm_index tfm;
    tfm.text_len = size;
    tfm.m_L = tfm_index::wt_type(L_buf, L_buf.size());
    tfm.m_C = vector<uint64_t>(tfm.m_L.sigma + 1, 0);
    for (uint64_t i = 0; i < L_buf.size(); i++) {
        // cout << L_buf[i] + 1 << endl;
        tfm.m_C[L_buf[i] + 1] += 1;
    }
    for (uint64_t i = 0; i < tfm.m_L.sigma; i++) tfm.m_C[i + 1] += tfm.m_C[i];
    tfm.m_dout = tfm_index::bit_vector_type(std::move(dout));
    sdsl::util::init_support(tfm.m_dout_rank, &tfm.m_dout);
    sdsl::util::init_support(tfm.m_dout_select, &tfm.m_dout);
    tfm.m_din = tfm_index::bit_vector_type(std::move(din));
    sdsl::util::init_support(tfm.m_din_rank, &tfm.m_din);
    sdsl::util::init_support(tfm.m_din_select, &tfm.m_din);

    return tfm;
}

void load_bitvector(sdsl::int_vector<1> &B, const std::string filename, const uint64_t n) {
    FILE *fin = fopen(filename.c_str(), "rb");
    uint64_t cnt = 0;
    uint8_t buffer = 0;
    for (uint64_t i = 0; i < (n + 7) / 8; i++) {
        int e = fread(&buffer, sizeof(uint8_t), 1, fin);
        if (e < 0)
            std::cout << "ERROR during bitvector loading!" << std::endl;
        // std::cout << (int) buffer << std::endl;
        for (int j = 0; j < 8; j++) {
            bool bit = 1 & (buffer >> (7 - j));
            B[cnt++] = bit;
            if (cnt == n) {
                fclose(fin);
                return;
            }
        }
    }
}

bit_vector create_from_boolvec(vector<bool> &v) {
    bit_vector b(v.size(), 0);
    for (size_t i=0; i < v.size(); i++) {
        b[i] = v[i];
    }
    return b;
}

tfm_index unparse(tfm_index &wg_parse, Dict &dict, size_t w, size_t size) {
    uint32_t *ilist = new uint32_t[wg_parse.L.size() - 1];
    generate_ilist(ilist, wg_parse, dict.dwords);

    uint_t *sa;
    int_t *lcp;
    compute_dict_bwt_lcp(dict.d, dict.dsize, dict.dwords, w, &sa, &lcp);

    vector<char> l_array = store_bwt(w, dict.d, dict.dsize, dict.end, ilist, wg_parse, dict.dwords, sa, lcp);
    vector<bool> bdin = store_din(w, dict.d, dict.dsize, wg_parse, dict.dwords, sa, lcp);
    bit_vector din = create_from_boolvec(bdin);
    vector<bool> bdout = store_dout(w, dict.d, dict.dsize, wg_parse, dict.dwords, sa, lcp);
    bit_vector dout = create_from_boolvec(bdout);

    string filename = "tmp.L";
    FILE *fbwt = fopen(filename.c_str(), "wb");
    for (char c: l_array) { if (fputc(c, fbwt) == EOF) die("L write error 0"); }
    fclose(fbwt);
    int_vector_buffer<> L_buf(filename, std::ios::in, 1024*1024, 8, true);

    tfm_index tfm = create_tfm(size, L_buf, din, dout);
    remove(filename);
    return tfm;
}

void write_tfm(tfm_index &unparsed, string &output) {
    string filename = output;
    filename += ".L";
    FILE *fbwt = fopen(filename.c_str(), "wb");
    for (char c: unparsed.L) { if (fputc(c, fbwt) == EOF) die("L write error 0"); }
    fclose(fbwt);

    string din_file = output + ".din";
    FILE *fdin = fopen(din_file.c_str(), "wb");
    uint8_t cnt = 0, buffer = 0;
    for (size_t i = 0; i < unparsed.din.size() - 1; i++) {
        write_bitvector(fdin, unparsed.din[i], cnt, buffer);
    }
    write_bitvector(fdin, 1, cnt, buffer, true);
    fclose(fdin);

    string dout_file = output + ".dout";
    FILE *fdout = fopen(dout_file.c_str(), "wb");
    cnt = 0, buffer = 0;
    for (size_t i = 0; i < unparsed.dout.size() - 1; i++) {
        write_bitvector(fdout, unparsed.dout[i], cnt, buffer);
    }
    write_bitvector(fdout, 1, cnt, buffer, true);
    fclose(fdout);
}

//------------------------------------------------------------------------------

struct Args {
    string input;
    string output;
    size_t w;       // sliding window size and its default
    size_t p;       // modulus for establishing stopping w-tuples
};

void print_help(char **argv) {
    cout << "Usage: " << argv[0] << " -w 10 -p 100 -i input.txt -o output.wg" << endl
         << "\tOptions: " << endl
         << "\t-w W\tsliding window size" << endl
         << "\t-p M\tmodulo for defining phrases" << endl
         << "\t-i I\tinput file (text)" << endl
         << "\t-o O\toutput file (binary representation of WG)" << endl
         << "\t-h  \tshow help and exit" << endl;
}

Args parse_args(int argc, char **argv) {
    extern char *optarg;

    Args arg;
    int c;
    string sarg;

    while ((c = getopt(argc, argv, "p:w:i:o:h")) != -1) {
        switch (c) {
            case 'i':
                arg.input.assign(optarg);
                break;
            case 'o':
                arg.output.assign(optarg);
                break;
            case 'w':
                sarg.assign(optarg);
                arg.w = stoi(sarg);
                break;
            case 'p':
                sarg.assign(optarg);
                arg.p = stoi(sarg);
                break;
            case 'h':
                print_help(argv);
                exit(1);
            case '?':
                cout << "Unknown option. Use -h for help." << endl;
                exit(1);
        }
    }
    return arg;
}

vector<uint64_t> remapParse(map<uint64_t, word_stats> &wfreq, vector<uint64_t> &parse) {
    vector<uint64_t> new_parse{};

    vector<uint32_t> occ(wfreq.size() + 1, 0); // ranks are zero based
    for (uint64_t hash : parse) {
        uint32_t rank = wfreq.at(hash).rank;
        occ[rank]++;
        new_parse.push_back(rank);
    }
    new_parse.push_back(0);
    return new_parse;
}

void pf_parse(string &input, size_t w, size_t p, vector<uint64_t> &parse, Dict &dict) {
    map<uint64_t, word_stats> wordFreq;
    vector<char> last{}; // this is maybe not needed
    calculate_word_frequencies(input, w, p, wordFreq, parse, last);

    // create array of dictionary words
    vector<const string *> dictArray;
    uint64_t totDWord = wordFreq.size();
    dictArray.reserve(totDWord);
    for (auto &x : wordFreq) { dictArray.push_back(&x.second.str); }
    assert(dictArray.size() == totDWord);
    sort(dictArray.begin(), dictArray.end(), pstringCompare);
    // write plain dictionary, also compute rank for each hash
    vector<char> dictionary{};
    writeDictOcc(wordFreq, dictArray, dictionary);
    dictArray.clear(); // reclaim memory

    dict = read_dictionary(dictionary);
    parse = remapParse(wordFreq, parse);
}

vector<uint64_t> compute_bwt(vector<uint64_t> &text) {
    uint64_t sigma = 0; // = 183416 + 1 + 2;
    for (size_t i = 0; i < text.size(); i++) {
        if (sigma < text[i])
            sigma = text[i];
    }
    sigma += 1 + 2;
    // +1 because {0,1,2} => 3, +2 because gsacak reserves 0 and 1
    // so all number needs to be shifted

    size_t n = text.size() + 2;
    // appending 1 ends "words", appending 0 ends input to gsacak

    uint32_t *t = (uint32_t *)malloc(n * sizeof(*t));
    for (size_t i = 0; i < text.size(); i++) t[i] = text[i]+2;
    t[n-2] = 1; t[n-1] = 0;

    uint32_t *sa = (uint32_t *)malloc(n * sizeof(*sa));
    gsacak_int(t, sa, NULL, NULL, n, sigma);

    vector<uint64_t> bwt{};
    for (size_t i = 0; i < n; i++) {
        if (sa[i] == 0) { bwt.push_back(0); continue; }
        if (t[sa[i]-1] == 1) continue;
        if (t[sa[i]-1] == 2) continue;
        bwt.push_back(t[sa[i]-1]-2);
    }
    free(sa);
    free(t);
    return bwt;
}

tfm_index construct_tfm_index(vector<uint64_t> &bwt) {
    string bwt_filename = "bwt.tmp";
    FILE *fout = fopen(bwt_filename.c_str(), "wb");
    fwrite(bwt.data(), sizeof(bwt[0]), bwt.size(), fout);
    fclose(fout);

    int_vector_buffer<> L(bwt_filename, std::ios::in, bwt.size(), 64, true);
    wt_blcd_int<> wt_L = wt_blcd_int<>(L, bwt.size());
    vector<uint64_t> C = vector<uint64_t>(wt_L.sigma + 1, 0);
    for (uint64_t i = 0; i < bwt.size(); i++) C[L[i] + 1] += 1;
    for (uint64_t i = 0; i < wt_L.sigma; i++) C[i + 1] += C[i];

    bit_vector B;
    pair<tfm_index::size_type, tfm_index::size_type> dbg_res;
    dbg_res = dbg_algorithms::find_min_dbg(wt_L, C, B);

    sdsl::bit_vector dout = B;
    sdsl::bit_vector din;
    std::swap(din, B);
    dbg_algorithms::mark_prefix_intervals(wt_L, C, dout, din);

    sdsl::remove(bwt_filename);
    string tmp_file_name = "construct_tfm_index.tmp";

    sdsl::int_vector_buffer<> L_buf(tmp_file_name, std::ios::out);

    tfm_index::size_type p = 0;
    tfm_index::size_type q = 0;
    for (tfm_index::size_type i = 0; i < wt_L.size(); i++) {
        if (din[i] == 1) {
            L_buf.push_back(wt_L[i]);
            dout[p++] = dout[i];
        }
        if (dout[i] == 1) {
            din[q++] = din[i];
        }
    }
    dout[p++] = 1;
    din[q++] = 1;
    dout.resize(p);
    din.resize(q);

    tfm_index tfm_index = create_tfm(bwt.size(), L_buf, din, dout);

    sdsl::remove(tmp_file_name);
    return tfm_index;
}

void untunnel(tfm_index &tfm, string &filename) {
    char *original = new char[tfm.size()];
    auto p = tfm.end();
    for (sdsl::int_vector<>::size_type i = 0; i < tfm.size(); i++) {
        char c = (char)tfm.backwardstep(p);
        original[tfm.size() - i - 1] = c;
    }

    FILE *fout = fopen(filename.c_str(), "w");
    fwrite(original, sizeof(char), tfm.size(), fout);
    fclose(fout);
}

int main(int argc, char **argv) {
    Args arg = parse_args(argc, argv);

    vector<uint64_t> parse{};
    Dict dict;
    size_t size = 12156306; // yeast.raw
    pf_parse(arg.input, arg.w, arg.p, parse, dict);
    vector<uint64_t> bwt = compute_bwt(parse);
    tfm_index tfm = construct_tfm_index(bwt);
    tfm_index unparsed = unparse(tfm, dict, arg.w, size);

    write_tfm(unparsed, arg.input);
    store_to_file(unparsed, arg.output);
    tfm_index loaded;
    load_from_file(loaded, arg.output);
    string filename = arg.input + ".untunneled";
    write_tfm(loaded, filename);
    // untunnel(loaded, filename);

    return 0;
}
