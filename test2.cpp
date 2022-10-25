#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sdsl/int_vector.hpp>
#include <vector>

#include <sdsl/io.hpp>
#include <sdsl/util.hpp>
#include "tfm_index.hpp"

extern "C" {
#include "gsacak.c"
#include "utils.c"
}

using namespace std;
using namespace sdsl;

struct Dict {
    uint8_t *d; // pointer to the dictionary
    uint64_t *end; // end[i] is the index of the ending symbol of the i-th phrase
    uint64_t dsize;  // dicionary size in symbols
    uint64_t dwords; // the number of phrases of the dicionary
};

string untunnel(tfm_index &tfm) {
    string original(tfm.size(), ' ');

    auto p = tfm.end();
    for (tfm_index::size_type i = 0; i < tfm.size(); i++) {
        char c = (char)tfm.backwardstep(p);
        original[tfm.size() - i - 1] = c;
    }

    return original;
}

//------------------------------------------------------------------------------
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

uint8_t get_prev(int w, uint8_t *d, uint64_t *end, uint32_t seqid) {
    return d[end[seqid] - w - 1];
}

long binsearch(uint_t x, uint_t a[], long n) {
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

int_t getlen(uint_t p, uint_t eos[], long n, uint32_t *seqid) {
    assert(p < eos[n - 1]);
    *seqid = binsearch(p, eos, n);
    assert(eos[*seqid] > p); // distance between position p and the next $
    return eos[*seqid] - p;
}

void tmp(tfm_index &wg, Dict &dict, size_t w, uint32_t *sa, int32_t *lcp, uint32_t *ilist) {
    for (size_t i=0; i<dict.dsize; i++) {
        if (dict.d[i] == '\001') dict.d[i] = '$';
        if (dict.d[i] == '\002') dict.d[i] = '#';
    }

    for (size_t i=0; i<dict.dsize; i++) {
        // uint32_t seqid = binsearch(sa_d[i], sa_d + 1, dict.dwords);
        uint32_t seqid = 0;
        // int32_t suffixLen = getlen(sa_d[i], sa_d + 1, dict.dwords, &seqid);
        cout << i << "\t"
             << sa[i] << "\t"
             << seqid << "\t"
             << dict.d + sa[i] << endl;
    }

    for (size_t i=0; i<dict.dsize; i++) {
        if (dict.d[i] == '$') dict.d[i] = '\001';
        if (dict.d[i] == '#') dict.d[i] = '\002';
    }

    // seqid = binsearch(sa[i], eos, dwords);
    //     int_t suffixLen = eos[seqid] - sa[i];

    // uint32_t seqid;
    // cout << dict.dwords + w + 1;
    // for (uint64_t i = 0; i < dict.dsize; i++) {
        // cout << dict.d + sa[i] << endl;
        // int32_t suffixLen = getlen(sa[i], sa + 1, dict.dwords, &seqid);
        // if (suffixLen <= (int32_t)w) continue;
        // char prev = dict.d[sa[i] - 1];
        // uint64_t parse_occ = wg.C[seqid + 1] - wg.C[seqid];
        // uint32_t *rank = ilist + (wg.C[seqid] - 1);

        // printf("%.*s\n", suffixLen, dict.d + sa[i]);
        // cout << prev << "\t" << nextsuffixLen << "\n";
    // }
}

size_t get_untunneled_size(tfm_index &tfmp, Dict &dict, size_t w, uint32_t *sa, int32_t *lcp, uint32_t *ilist) {
    size_t size = 0;

    uint64_t next;
    uint32_t seqid;
    for (uint64_t i = dict.dwords + w + 1; i < dict.dsize; i = next) {
        next = i + 1;
        int32_t suffixLen = getlen(sa[i], sa + 1, dict.dwords, &seqid);
        if (suffixLen <= (int32_t)w) continue;

        if (sa[i] == 0 || dict.d[sa[i] - 1] == EndOfWord) {
            // ----- simple case: the suffix is a full word
            uint32_t start = tfmp.C[seqid + 1];
            uint32_t end = tfmp.C[seqid + 2];
            for (uint32_t j = start; j < end; j++) {
                if (tfmp.din[j] == 1) {
                    uint32_t pos = tfmp.dout_select(tfmp.din_rank(j + 1));
                    size++; pos++;
                    while (tfmp.dout[pos] != 1) { size++; pos++; }
                }
            }
        } else {
            // ----- hard case: there can be a group of equal suffixes starting
            // at i save seqid and the corresponding char
            vector<uint32_t> id2merge(1, seqid);
            vector<uint8_t> char2write(1, dict.d[sa[i] - 1]);
            while (next < dict.dsize && lcp[next] >= suffixLen) {
                int32_t nextsuffixLen = getlen(sa[next], sa + 1, dict.dwords, &seqid);
                if (nextsuffixLen != suffixLen) break;
                id2merge.push_back(seqid); // sequence to consider
                char2write.push_back(dict.d[sa[next] - 1]); // corresponding char
                next++;
            }

            size_t numwords = id2merge.size(); // numwords dictionary words contain the same suffix
            bool samechar = true;
            for (size_t i = 1; (i < numwords) && samechar; i++) {
                samechar = (char2write[i - 1] == char2write[i]);
            }

            if (samechar) {
                for (size_t i = 0; i < numwords; i++) {
                    uint32_t s = id2merge[i] + 1;
                    size += tfmp.C[s + 1] - tfmp.C[s];
                }
            } else {
                // many words, many chars...
                vector<SeqId> heap; // create heap
                for (size_t i = 0; i < numwords; i++) {
                    uint32_t s = id2merge[i] + 1;
                    heap.push_back(SeqId(
                        s,                          // letter from parse
                        tfmp.C[s + 1] - tfmp.C[s],  // num of occ in parse
                        ilist + (tfmp.C[s] - 1),    // ???
                        char2write[i]               // previous chars
                    ));
                }
                std::make_heap(heap.begin(), heap.end());
                while (heap.size() > 0) {
                    // output char for the top of the heap
                    SeqId s = heap.front();
                    pop_heap(heap.begin(), heap.end());
                    heap.pop_back();

                    size++;
                    // if remaining positions, reinsert to heap
                    if (s.next()) {
                        heap.push_back(s);
                        push_heap(heap.begin(), heap.end());
                    }
                }
            }
        }
    }

    return size;
}

int_vector<> compute_L(size_t w, uint8_t *d, long dsize, uint64_t *end_to_phrase, uint32_t *ilist, tfm_index &tfmp, long dwords, uint_t *sa, int_t *lcp) {
    uint_t *eos = sa + 1;
    vector<char> out{};

    long next;
    uint32_t seqid;
    for (long i = dwords + w + 1; i < dsize; i = next) {
        next = i + 1;
        int_t suffixLen = getlen(sa[i], eos, dwords, &seqid);
        if (suffixLen <= (int_t)w) continue;

        if (sa[i] == 0 || d[sa[i] - 1] == EndOfWord) {
            // ----- simple case: the suffix is a full word
            uint32_t start = tfmp.C[seqid + 1];
            uint32_t end = tfmp.C[seqid + 2];
            for (uint32_t j = start; j < end; j++) {
                if (tfmp.din[j] == 1) {
                    uint32_t pos = tfmp.dout_select(tfmp.din_rank(j + 1));
                    do {
                        if (tfmp.L[pos] == 0) pos = 0;
                        uint32_t act_phrase = tfmp.L[pos] - 1;
                        uint8_t char_to_write = get_prev(w, d, end_to_phrase, act_phrase);
                        out.push_back(char_to_write);
                    } while (tfmp.dout[++pos] != 1);
                }
            }
        } else {
            // ----- hard case: there can be a group of equal suffixes starting
            // at i save seqid and the corresponding char
            vector<uint32_t> id2merge(1, seqid);
            vector<uint8_t> char2write(1, d[sa[i] - 1]);
            while (next < dsize && lcp[next] >= suffixLen) {
                int_t nextsuffixLen = getlen(sa[next], eos, dwords, &seqid);
                if (nextsuffixLen != suffixLen) break;
                id2merge.push_back(seqid); // sequence to consider
                char2write.push_back(d[sa[next] - 1]); // corresponding char
                next++;
            }

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

    int_vector<> L(out.size(), 0);
    for (size_t i=0; i<L.size(); i++) { L[i] = out[i]; }

    return L;
}

void compute_degrees(
    tfm_index &tfmp, Dict &dict, size_t w, uint32_t *sa, int32_t *lcp,
    bit_vector &din, bit_vector &dout
) {
    uint8_t *d = dict.d;
    long dsize = dict.dsize;
    long dwords = dict.dwords;

    uint32_t *eos = sa + 1;
    size_t p = 0;
    size_t q = 0;

    long next;
    uint32_t seqid;
    for (long i = dwords + w + 1; i < dsize; i = next) {
        next = i + 1;
        int32_t suffixLen = getlen(sa[i], eos, dwords, &seqid);
        if (suffixLen <= (int32_t)w) continue;

        if (sa[i] == 0 || d[sa[i] - 1] == EndOfWord) {
            // ----- simple case: the suffix is a full word
            uint32_t start = tfmp.C[seqid + 1];
            uint32_t end = tfmp.C[seqid + 2];
            for (uint32_t j = start; j < end; j++) {
                din[p++] = tfmp.din[j];
                if (tfmp.din[j] == 1) {
                    uint32_t pos = tfmp.dout_select(tfmp.din_rank(j + 1));
                    if (tfmp.L[pos] == 0) pos = 0;
                    do { dout[q++] = tfmp.dout[pos]; }
                    while (tfmp.dout[++pos] != 1);
                }
            }
        } else {
            // ----- hard case: there can be a group of equal suffixes starting
            // at i save seqid and the corresponding char
            int bits_to_write = tfmp.C[seqid + 2] - tfmp.C[seqid + 1];
            while (next < dsize && lcp[next] >= suffixLen) {
                int32_t nextsuffixLen = getlen(sa[next], eos, dwords, &seqid);
                if (nextsuffixLen != suffixLen) break;
                bits_to_write += tfmp.C[seqid + 2] - tfmp.C[seqid + 1];
                next++;
            }
            for (int k = 0; k < bits_to_write; k++) {
                din[p++] = 1;
                dout[q++] = 1;
            }
        }
    }
    din[p++] = 1;
    dout[q++] = 1;
    din.resize(p);
    dout.resize(q);
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

tfm_index unparse(tfm_index &wg_parse, Dict &dict, size_t w, size_t size) {
    uint32_t *inverted_list = new uint32_t[wg_parse.L.size() - 1];
    generate_ilist(inverted_list, wg_parse, dict.dwords);

    uint32_t *sa_d = new uint32_t[dict.dsize];
    int32_t *lcp_d = new int32_t[dict.dsize];
    // separators s[i]=1 and with s[n-1]=0
    // cout << dict.d << "\n" << dict.dsize << endl;;
    gsacak(dict.d, sa_d, lcp_d, NULL, dict.dsize);

    tmp(wg_parse, dict, w, sa_d, lcp_d, inverted_list);

    size_t s = get_untunneled_size(wg_parse, dict, w, sa_d, lcp_d, inverted_list);
    int_vector<> L = compute_L(w, dict.d, dict.dsize, dict.end, inverted_list, wg_parse, dict.dwords, sa_d, lcp_d);
    cout << s << " " << L.size() << endl;
    bit_vector din(L.size() + 1, 1);
    bit_vector dout(L.size() + 1, 1);

    compute_degrees(wg_parse, dict, w, sa_d, lcp_d, din, dout);

    tfm_index tfm(size, L, din, dout);
    return tfm;
}

int main() {
    string input = "GTAGGTGGGTTGGTAGGTGGGTTGGTTT";
    size_t orig_size = input.size();
    size_t w = 2;
    // triggers = [GT, TC]

    Dict dict;
    dict.dsize = 33;
    dict.dwords = 5;
    dict.d = (uint8_t *)malloc(dict.dsize);                         // {2GT, GTAGGT, GTGGGT, GTTGGT, GTTT22}
    dict.end = (uint64_t *)malloc(dict.dwords * sizeof(uint64_t));  // {3, 10, 17, 24, 31}

    uint8_t d[] = "\002GT\001GTAGGT\001GTGGGT\001GTTGGT\001GTTT\002\002\001\000";
    for (size_t i=0; i<dict.dsize; i++) { dict.d[i] = d[i]; }
    uint64_t end[] = {3, 10, 17, 24, 31};
    for (size_t i=0; i<dict.dwords; i++) { dict.end[i] = end[i]; }

    // parse = [1, 2, 3, 4, 2, 3, 4, 5, 0]
    // bwt   = [5, 0, 1, 4, 2, 2, 3, 3, 4]
    // L     = [5, 0, 1, 4, 2,    3,    4]

    size_t size = 9;
    int_vector<> l = {5, 0, 1, 4, 2, 3,    4};
    bit_vector    out{1, 1, 1, 0, 1, 1,    1, 1};
    bit_vector     in{1, 1, 1,    1, 1, 0, 1, 1};

    tfm_index tfm(size, l, in, out);

    tfm_index unparsed = unparse(tfm, dict, w, orig_size);

    cout << input << endl;
    cout << untunnel(unparsed) << endl;
}
