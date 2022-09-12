#include <algorithm>
#include <cstdint>
#include <ctime>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
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
// #include <zlib.h>

#include "tfm_index.hpp"

extern "C" {
#include "gsacak.h"
#include "utils.h"
}

using namespace std;
using namespace sdsl;
using namespace __gnu_cxx;

// =============== algorithm limits ===================
// maximum number of distinct words
#define MAX_DISTINCT_WORDS (INT32_MAX - 1)
typedef uint32_t word_int_t;
// maximum number of occurrences of a single word
#define MAX_WORD_OCC (UINT32_MAX)
typedef uint32_t occ_int_t;

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

// values of the wordFreq map: word, its number of occurrences, and its rank
struct word_stats {
    string str;
    occ_int_t occ;
    word_int_t rank = 0;
};

// -------------------------------------------------------------
// struct containing command line parameters and other globals
struct Args {
    string inputFileName = "";
    size_t w = 10;         // sliding window size and its default
    size_t p = 100;        // modulus for establishing stopping w-tuples
};

// -----------------------------------------------------------------
// class to maintain a window in a string and its KR fingerprint
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

    // init window, hash, and tot_char
    void reset() {
        for (int i = 0; i < wsize; i++)
            window[i] = 0;
        // init hash value and related values
        hash = tot_char = 0;
    }

    uint64_t addchar(int c) {
        int k = tot_char++ % wsize;
        // complex expression to avoid negative numbers
        hash +=
            (prime - (window[k] * asize_pot) % prime
            );                             // remove window[k] contribution
        hash = (asize * hash + c) % prime; //  add char i
        window[k] = c;
        // cerr << get_window() << " ~~ " << window << " --> " << hash << endl;
        return hash;
    }
    // debug only
    string get_window() {
        string w = "";
        int k = (tot_char - 1) % wsize;
        for (int i = k + 1; i < k + 1 + wsize; i++)
            w.append(1, window[i % wsize]);
        return w;
    }

    ~KR_window() { delete[] window; }
};
// -----------------------------------------------------------

// compute 64-bit KR hash of a string
// to avoid overflows in 64 bit aritmethic the prime is taken < 2**55
uint64_t kr_hash(string s) {
    uint64_t hash = 0;
    // const uint64_t prime = 3355443229;     // next prime(2**31+2**30+2**27)
    const uint64_t prime =
        27162335252586509; // next prime (2**54 + 2**53 + 2**47 + 2**13)
    for (size_t k = 0; k < s.size(); k++) {
        int c = (unsigned char)s[k];
        assert(c >= 0 && c < 256);
        hash = (256 * hash + c) % prime; //  add char k
    }
    return hash;
}

// save current word in the freq map and update it leaving only the
// last minsize chars which is the overlap with next word
static void save_update_word(
    string &w, unsigned int minsize, map<uint64_t, word_stats> &freq,
    FILE *tmp_parse_file, FILE *last, FILE *sa, uint64_t &pos
) {
    assert(pos == 0 || w.size() > minsize);
    if (w.size() <= minsize)
        return;
    // get the hash value and write it to the temporary parse file
    uint64_t hash = kr_hash(w);
    if (fwrite(&hash, sizeof(hash), 1, tmp_parse_file) != 1)
        die("parse write error");

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
    if (fputc(w[w.size() - minsize - 1], last) == EOF)
        die("Error writing to .last file");
    // compute ending position +1 of current word and write it to sa file
    // pos is the ending position+1 of the previous word and is updated here
    if (pos == 0)
        pos = w.size() - 1; // -1 is for the initial $ of the first word
    else
        pos += w.size() - minsize;
    if (sa)
        if (fwrite(&pos, IBYTES, 1, sa) != 1)
            die("Error writing to sa info file");
    // keep only the overlapping part of the window
    w.erase(0, w.size() - minsize);
}

// prefix free parse of file fnam. w is the window size, p is the modulus
// use a KR-hash as the word ID that is immediately written to the parse file
uint64_t process_file(Args &arg, map<uint64_t, word_stats> &wordFreq) {
    // open a, possibly compressed, input file
    string fnam = arg.inputFileName;
    // open the 1st pass parsing file
    FILE *g = open_aux_file(arg.inputFileName.c_str(), EXTPARS0, "wb");
    // open output file containing the char at position -(w+1) of each word
    FILE *last_file = open_aux_file(arg.inputFileName.c_str(), EXTLST, "wb");

    // main loop on the chars of the input file
    uint64_t pos = 0; // ending position +1 of previous word in the original
                      // text, used for computing sa_info
    assert( IBYTES <= sizeof(pos) );
    // init first word in the parsing with a Dollar char
    string word("");
    word.append(1, Dollar);
    KR_window krw(arg.w);
    std::string line;
    ifstream f(fnam);
    if (!f.rdbuf()->is_open()) { // is_open does not work on igzstreams
        perror(__func__);
        throw new std::runtime_error("Cannot open input file " + fnam);
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
        if (hash % arg.p == 0) {
            // end of word, save it and write its full hash to the output
            // file cerr << "~"<< c << "~ " << hash << " ~~ <" << word << ">
            // ~~ <" << krw.get_window() << ">" <<  endl;
            save_update_word(
                word, arg.w, wordFreq, g, last_file, NULL, pos
            );
        }
    }
    f.close();
    // virtually add w null chars at the end of the file and add the last word
    // in the dict
    word.append(arg.w, Dollar);
    save_update_word(word, arg.w, wordFreq, g, last_file, NULL, pos);
    // close input and output files
    if (fclose(last_file) != 0) die("Error closing last file");
    if (fclose(g) != 0) die("Error closing parse file");
    if (pos != krw.tot_char + arg.w)
        cerr << "Pos: " << pos << " tot " << krw.tot_char << endl;
    return krw.tot_char;
}

// function used to compare two string pointers
bool pstringCompare(const string *a, const string *b) { return *a <= *b; }

// given the sorted dictionary and the frequency map write the dictionary and
// occ files also compute the 1-based rank for each hash
void writeDictOcc(
    Args &arg, map<uint64_t, word_stats> &wfreq,
    vector<const string *> &sortedDict
) {
    assert(sortedDict.size() == wfreq.size());
    // open dictionary and occ files
    FILE *fdict = open_aux_file(arg.inputFileName.c_str(), EXTDICT, "wb");
    FILE *focc = open_aux_file(arg.inputFileName.c_str(), EXTOCC, "wb");

    word_int_t wrank = 1; // current word rank (1 based)
    for (auto x : sortedDict) {
        const char *word = (*x).data(); // current dictionary word
        int offset = 0;
        size_t len = (*x).size(); // offset and length of word
        assert(len > (size_t)arg.w);
        size_t s = fwrite(word + offset, 1, len, fdict);
        if (s != len)
            die("Error writing to DICT file");
        if (fputc(EndOfWord, fdict) == EOF)
            die("Error writing EndOfWord to DICT file");
        uint64_t hash = kr_hash(*x);
        auto &wf = wfreq.at(hash);
        assert(wf.occ > 0);
        s = fwrite(&wf.occ, sizeof(wf.occ), 1, focc);
        if (s != 1)
            die("Error writing to OCC file");
        assert(wf.rank == 0);
        wf.rank = wrank++;
    }
    if (fputc(EndOfDict, fdict) == EOF)
        die("Error writing EndOfDict to DICT file");
    if (fclose(focc) != 0)
        die("Error closing OCC file");
    if (fclose(fdict) != 0)
        die("Error closing DICT file");
}

void remapParse(Args &arg, map<uint64_t, word_stats> &wfreq) {
    // open parse files. the old parse can be stored in a single file or in
    // multiple files
    mFile *moldp = mopen_aux_file(arg.inputFileName.c_str(), EXTPARS0, 0);
    FILE *newp = open_aux_file(arg.inputFileName.c_str(), EXTPARSE, "wb");

    // recompute occ as an extra check
    vector<occ_int_t> occ(wfreq.size() + 1, 0); // ranks are zero based
    uint64_t hash;
    while (true) {
        size_t s = mfread(&hash, sizeof(hash), 1, moldp);
        if (s == 0) break;
        if (s != 1) die("Unexpected parse EOF");
        word_int_t rank = wfreq.at(hash).rank;
        occ[rank]++;
        s = fwrite(&rank, sizeof(rank), 1, newp);
        if (s != 1) die("Error writing to new parse file");
    }
    if (fclose(newp) != 0)
        die("Error closing new parse file");
    if (mfclose(moldp) != 0)
        die("Error closing old parse segment");
}

void print_help(char **argv, Args &args) {
    cout << "Usage: " << argv[0] << " <input filename> [options]" << endl;
    cout << "  Options: " << endl
         << "\t-w W\tsliding window size, def. " << args.w << endl
         << "\t-p M\tmodulo for defining phrases, def. " << args.p << endl
         << "\t-h  \tshow help and exit" << endl;
    exit(1);
}

void parseArgs(int argc, char **argv, Args &arg) {
    int c;
    extern char *optarg;
    extern int optind;

    string sarg;
    while ((c = getopt(argc, argv, "p:w:h")) != -1) {
        switch (c) {
        case 'w':
            sarg.assign(optarg);
            arg.w = stoi(sarg);
            break;
        case 'p':
            sarg.assign(optarg);
            arg.p = stoi(sarg);
            break;
        case 'h':
            print_help(argv, arg);
            exit(1);
        case '?':
            cout << "Unknown option. Use -h for help." << endl;
            exit(1);
        }
    }
    // the only input parameter is the file name
    if (argc == optind + 1) {
        arg.inputFileName.assign(argv[optind]);
    } else {
        cout << "Invalid number of arguments" << endl;
        print_help(argv, arg);
    }
    // check algorithm parameters
    if (arg.w < 4) {
        cout << "Windows size must be at least 4\n";
        exit(1);
    }
    if (arg.p < 10) {
        cout << "Modulus must be at leas 10\n";
        exit(1);
    }
}

// -------------------------------------------------------
// type used to represent an entry in the SA
// this is currently 32 bit for gsacak and 64 bit for gsacak-64
// note that here we use sacak (SA computation for a single string of 32 bit
// symbols)
typedef uint_t sa_index_t;

void printUsage(char **argv) {
    cerr << "USAGE: " << argv[0] << " INFILE TFMOUTFILE" << endl;
    cerr << "INFILE:" << endl;
    cerr << "  Parse file to construct tunneled FM index from, nullbytes are "
            "permitted"
         << endl;
    cerr << "TFMOUTFILE:" << endl;
    cerr
        << "  File where to store the serialized tunneled FM index of the parse"
        << endl;
};

void compute_BWT(uint32_t *Text, long n, long k, string filename) {
    sa_index_t *SA = (sa_index_t *)malloc(n * sizeof(*SA));
    sacak_int(Text, SA, n, k);

    sa_index_t *BWTsa = SA; // BWT overlapping SA
    assert(n > 1);
    // first BWT symbol
    assert(SA[0] == n);
    // 2nd, 3rd etc BWT symbols
    for (long i = 0; i < n; i++) {
        if (SA[i] == 0) { BWTsa[i] = 0; }
        else { BWTsa[i] = Text[SA[i] - 1]; }
    }

    FILE *fout = fopen(filename.c_str(), "wb");
    fwrite(BWTsa, sizeof(BWTsa[0]), n, fout);
    fclose(fout);
}

uint32_t *load_parse(const string &infile, size_t &psize) {
    FILE *fp = fopen(infile.c_str(), "r");

    fseek(fp, 0, SEEK_END);
    size_t n;
    n = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    n = n / 4;
    psize = n;

    uint32_t *parse = new uint32_t[n + 1];
    size_t s = fread(parse, sizeof(*parse), n, fp);
    if (s != n) {
        printf("Parse loading error!");
        exit(1);
    }
    parse[n] = 0; // sacak needs this
    return parse;
}

size_t compute_sigma(const uint32_t *parse, const size_t psize) {
    uint32_t max = 0;
    for (size_t i = 0; i < psize; i++) {
        if (max < parse[i])
            max = parse[i];
    }
    return (max + 1);
}

void calculate_word_frequencies(Args &arg, map<uint64_t, word_stats> &wordFreq) {
    try {
        process_file(arg, wordFreq);
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
    uint64_t
        *end; // end[i] is the index of the ending symbol of the i-th phrase
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

// compute the SA and LCP array for the set of (unique) dictionary words
// using gSACA-K. Also do some checking based on the number and order of the
// special symbols d[0..dsize-1] is the dictionary consisting of the
// concatenation of dictionary words in lex order with EndOfWord (0x1) at the
// end of each word and d[size-1] = EndOfDict (0x0) at the very end. It is
// d[0]=Dollar (0x2) since the first words starts with $. There is another word
// somewhere ending with Dollar^wEndOfWord (it is the last word in the parsing,
// but its lex rank is unknown).
static void compute_dict_bwt_lcp(
    uint8_t *d, long dsize, long dwords, int w, uint_t **sap, int_t **lcpp
) // output parameters
{
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

// class representing the suffix of a dictionary word
// instances of this class are stored to a heap to handle the hard bwts
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

// write to the bwt all the characters preceding a given suffix
// doing a merge operation if necessary
static void fwrite_chars_same_suffix(
    vector<uint32_t> &id2merge, vector<uint8_t> &char2write, tfm_index<> &tfmp,
    uint32_t *ilist, FILE *fbwt, long &easy_bwts, long &hard_bwts
) {
    size_t numwords =
        id2merge.size(); // numwords dictionary words contain the same suffix
    bool samechar = true;
    for (size_t i = 1; (i < numwords) && samechar; i++) {
        samechar = (char2write[i - 1] == char2write[i]);
    }

    if (samechar) {
        for (size_t i = 0; i < numwords; i++) {
            uint32_t s = id2merge[i] + 1;
            for (uint64_t j = tfmp.C[s]; j < tfmp.C[s + 1]; j++) {
                if (fputc(char2write[0], fbwt) == EOF)
                    die("L write error 1");
            }
            easy_bwts += tfmp.C[s + 1] - tfmp.C[s];
        }
    } else {
        // many words, many chars...
        vector<SeqId> heap; // create heap
        for (size_t i = 0; i < numwords; i++) {
            uint32_t s = id2merge[i] + 1;
            // cout << "phrase: " << s << " pos: " << ilist[tfmp.C[s]-1] <<
            // endl;
            heap.push_back(SeqId(
                s, tfmp.C[s + 1] - tfmp.C[s], ilist + (tfmp.C[s] - 1),
                char2write[i]
            ));
        }
        std::make_heap(heap.begin(), heap.end());
        while (heap.size() > 0) {
            // output char for the top of the heap
            SeqId s = heap.front();
            if (fputc(s.char2write, fbwt) == EOF)
                die("L write error 2");
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

inline uint8_t get_prev(int w, uint8_t *d, uint64_t *end, uint32_t seqid) {
    return d[end[seqid] - w - 1];
}

void write_bitvector(
    FILE *f, bool bit, uint8_t &cnt, uint8_t &buffer, bool hard_write = false
) {
    buffer |= (bit << (7 - cnt++));
    if (hard_write || (cnt == 8)) {
        if (fputc(buffer, f) == EOF)
            die("Din/Dout write error 0");
        cnt = 0;
        buffer = 0;
    }
}

/* *******************************************************************
 * Computation of the final BWT
 *
 * istart[] and islist[] are used together. For each dictionary word i
 * (considered in lexicographic order) for k=istart[i]...istart[i+1]-1
 * ilist[k] contains the ordered positions in BWT(P) containing word i
 * ******************************************************************* */
void bwt(
    Args &arg, uint8_t *d, long dsize,
    uint64_t *end_to_phrase,            // dictionary and its size
    uint32_t *ilist, tfm_index<> &tfmp, // ilist, last and their size
    long dwords, uint_t *sa, int_t *lcp
) { // starting point in ilist for each word and # words

    // set d[0]==0 as this is the EOF char in the final BWT
    assert(d[0] == Dollar);
    d[0] = 0;

    // derive eos from sa. for i=0...dwords-1, eos[i] is the eos position of
    // string i in d
    uint_t *eos = sa + 1;
    for (int i = 0; i < dwords - 1; i++)
        assert(eos[i] < eos[i + 1]);

    // open output file
    FILE *fbwt = open_aux_file(arg.inputFileName.c_str(), "L", "wb");

    // main loop: consider each entry in the SA of dict
    long full_words = 0, easy_bwts = 0, hard_bwts = 0, next;
    uint32_t seqid;
    for (long i = dwords + arg.w + 1; i < dsize; i = next) {
        // we are considering d[sa[i]....]
        next = i + 1; // prepare for next iteration
        // compute length of this suffix and sequence it belongs
        int_t suffixLen = getlen(sa[i], eos, dwords, &seqid);
        // cout << suffixLen << " " << seqid << endl;
        //  ignore suffixes of lenght <= w
        if (suffixLen <= (int_t)arg.w)
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
                        uint8_t char_to_write =
                            get_prev(arg.w, d, end_to_phrase, act_phrase);
                        easy_bwts++;
                        // cout << easy_bwts + hard_bwts << " " << seqid << " "
                        // << char_to_write << endl;
                        if (fputc(char_to_write, fbwt) == EOF)
                            die("L write error 0");
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
            fwrite_chars_same_suffix(
                id2merge, char2write, tfmp, ilist, fbwt, easy_bwts, hard_bwts
            );
        }
    }
    assert(full_words == dwords);
    fclose(fbwt);
}

void din(
    Args &arg, uint8_t *d, long dsize, // dictionary and its size
    tfm_index<> &tfmp, long dwords, uint_t *sa, int_t *lcp
) { // starting point in ilist for each word and # words

    // derive eos from sa. for i=0...dwords-1, eos[i] is the eos position of
    // string i in d
    uint_t *eos = sa + 1;
    for (int i = 0; i < dwords - 1; i++)
        assert(eos[i] < eos[i + 1]);

    // open output file
    FILE *fdin = open_aux_file(arg.inputFileName.c_str(), "din", "wb");

    // main loop: consider each entry in the SA of dict
    long next;
    uint32_t seqid;
    uint8_t cnt = 0, buffer = 0;
    for (long i = dwords + arg.w + 1; i < dsize; i = next) {
        // we are considering d[sa[i]....]
        next = i + 1; // prepare for next iteration
        // compute length of this suffix and sequence it belongs
        int_t suffixLen = getlen(sa[i], eos, dwords, &seqid);
        // cout << suffixLen << " " << seqid << endl;
        //  ignore suffixes of lenght <= w
        if (suffixLen <= (int_t)arg.w)
            continue;
        // ----- simple case: the suffix is a full word
        if (sa[i] == 0 || d[sa[i] - 1] == EndOfWord) {
            uint32_t start = tfmp.C[seqid + 1], end = tfmp.C[seqid + 2];
            assert(tfmp.din[start] == 1);
            for (uint32_t j = start; j < end; j++) {
                write_bitvector(fdin, tfmp.din[j], cnt, buffer);
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
                write_bitvector(fdin, 1, cnt, buffer);
        }
    }
    write_bitvector(fdin, 1, cnt, buffer, true);
    fclose(fdin);
}

void dout(
    Args &arg, uint8_t *d, long dsize, // dictionary and its size
    tfm_index<> &tfmp, long dwords, uint_t *sa, int_t *lcp
) { // starting point in ilist for each word and # words

    // derive eos from sa. for i=0...dwords-1, eos[i] is the eos position of
    // string i in d
    uint_t *eos = sa + 1;
    for (int i = 0; i < dwords - 1; i++)
        assert(eos[i] < eos[i + 1]);

    // open output file
    FILE *fdout = open_aux_file(arg.inputFileName.c_str(), "dout", "wb");

    // main loop: consider each entry in the SA of dict
    long next;
    uint32_t seqid;
    uint8_t cnt = 0, buffer = 0;
    for (long i = dwords + arg.w + 1; i < dsize; i = next) {
        // we are considering d[sa[i]....]
        next = i + 1; // prepare for next iteration
        // compute length of this suffix and sequence it belongs
        int_t suffixLen = getlen(sa[i], eos, dwords, &seqid);
        // cout << suffixLen << " " << seqid << endl;
        //  ignore suffixes of lenght <= w
        if (suffixLen <= (int_t)arg.w)
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
                        write_bitvector(fdout, tfmp.dout[pos], cnt, buffer);
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
            for (int k = 0; k < bits_to_write; k++)
                write_bitvector(fdout, 1, cnt, buffer);
        }
    }
    write_bitvector(fdout, 1, cnt, buffer, true);
    fclose(fdout);
}

Dict read_dictionary(const char *filename) {
    FILE *g = open_aux_file(filename, EXTDICT, "rb");
    fseek(g, 0, SEEK_END);
    long dsize = ftell(g);
    if (dsize < 0)
        die("ftell dictionary");
    if (dsize <= 1 + 4)
        die("invalid dictionary file");
    if (dsize > 0x7FFFFFFE) {
        printf("Dictionary size greater than  2^31-2!\n");
        printf("Please use 64 bit version\n");
        exit(1);
    }

    uint8_t *d = new uint8_t[dsize];
    rewind(g);
    uint64_t e = fread(d, 1, dsize, g);
    if (e != (uint64_t)dsize)
        die("Dictionary fread errror!");
    fclose(g);

    uint64_t dwords = 0;
    for (int i = 0; i < dsize; i++) {
        if (d[i] == EndOfWord)
            dwords++;
    }

    uint64_t *end = new uint64_t[dwords];
    int cnt = 0;
    for (int i = 0; i < dsize; i++) {
        if (d[i] == EndOfWord)
            end[cnt++] = i;
    }

    Dict res = {d, end, e, dwords};
    return res;
}

void generate_ilist(uint32_t *ilist, tfm_index<> &tfmp, uint64_t dwords) {
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


int main(int argc, char **argv) {
    Args arg;
    parseArgs(argc, argv, arg);

    map<uint64_t, word_stats> wordFreq;
    calculate_word_frequencies(arg, wordFreq); // + <fn>.last <fn>.parse_old

    // create array of dictionary words
    vector<const string *> dictArray;
    uint64_t totDWord = wordFreq.size();
    dictArray.reserve(totDWord);
    for (auto &x : wordFreq) {
        dictArray.push_back(&x.second.str);
    }
    assert(dictArray.size() == totDWord);
    sort(dictArray.begin(), dictArray.end(), pstringCompare);
    // write plain dictionary and occ file, also compute rank for each hash
    writeDictOcc(arg, wordFreq, dictArray); // + <fn>.dict <fn>.occ
    dictArray.clear(); // reclaim memory

    remapParse(arg, wordFreq); // + <fn>.parse

    // construct tunneled fm index
    tfm_index<> tfm;

    cache_config config(true, "./", util::basename(arg.inputFileName));
    size_t psize;
    uint32_t *parse = load_parse(arg.inputFileName + ".parse", psize);
    size_t sigma = compute_sigma(parse, psize);
    compute_BWT(parse, psize + 1, sigma, arg.inputFileName + ".bwt"); // <fn>.bwt
    delete parse;

    construct_tfm_index(tfm, arg.inputFileName + ".bwt", psize + 1, config);

//-------------------------------------------------------------------------------

    struct Dict dict = read_dictionary(arg.inputFileName.c_str());

    uint32_t *ilist = new uint32_t[tfm.L.size() - 1];
    generate_ilist(ilist, tfm, dict.dwords);

    // compute SA and BWT of D and do some checking on them
    uint_t *sa;
    int_t *lcp;
    compute_dict_bwt_lcp(dict.d, dict.dsize, dict.dwords, arg.w, &sa, &lcp);

    bwt(arg, dict.d, dict.dsize, dict.end, ilist, tfm, dict.dwords, sa, lcp);
    din(arg, dict.d, dict.dsize, tfm, dict.dwords, sa, lcp);
    dout(arg, dict.d, dict.dsize, tfm, dict.dwords, sa, lcp);
    delete[] ilist;
    delete[] dict.d;
    delete[] dict.end;
    delete[] lcp;
    delete[] sa;

    return 0;
}
