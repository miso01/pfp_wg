#include <cstddef>
#include <cstdint>
#include <cstdio>
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
    uint8_t  *d;        // pointer to the dictionary
    uint64_t *end;      // end[i] is the index of the ending symbol of the i-th phrase
    uint64_t dsize;     // dicionary size in symbols
    uint64_t dwords;    // the number of phrases of the dicionary
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
    SeqId(uint32_t i, int r, uint32_t *b, int8_t c) {
        id = i;
        remaining = r;
        bwtpos = b;
        char2write = c;
    }

    // advance to the next bwt position, return false if there are no more
    // positions
    bool next() {
        remaining--;
        bwtpos += 1;
        return remaining > 0;
    };

    bool operator<(const SeqId &a) { return *bwtpos > *(a.bwtpos); };
};

long binsearch(uint_t x, uint_t a[], long n) {
    long lo = 0;
    long hi = n - 1;
    while (hi > lo) {
        int mid = (lo + hi) / 2;
        if (x < a[mid]) hi = mid;
        else lo = mid + 1;
    }
    return hi;
}

struct row_t {
    bool is_phrase;
    uint32_t id;
    uint32_t len;
    uint32_t id_in;
    uint32_t id_out;
    uint32_t n_occs;
    uint32_t *ranks;
    uint32_t n_prevs;
    uint32_t *prevs;
};

row_t get_row(uint32_t i, tfm_index &wg, Dict &dict, size_t w, uint32_t *sa, int32_t *lcp, uint32_t *ilist) {
    row_t row;
    row.id = binsearch(sa[i], sa + 1, dict.dwords) + 1; // +1 because seqid 0 is parse end symbol
    row.len = sa[row.id] - sa[i];

    uint32_t start = wg.C[row.id];

    row.id_in = wg.din_select(row.id + 2) - wg.din_select(row.id + 1);
    row.id_out = wg.dout_select(row.id + 2) - wg.dout_select(row.id + 1);
    row.is_phrase = false;

    string prev = "";
    vector<int32_t> ranks;

    if (sa[i] == 0 || dict.d[sa[i] - 1] == '$') {
        row.is_phrase = true;

        for (uint32_t j = start; j < start + row.n_occs; j++) {
            if (wg.din[j] == 1) {
                uint32_t pos = wg.dout_select(wg.din_rank(j + 1));
                do {
                    if (wg.L[pos] == 0) pos = 0;
                    uint32_t act_phrase = wg.L[pos];    // F -> L

                    prev += dict.d[sa[act_phrase] - w - 1];
                    pos++;
                } while (wg.dout[pos] == 0);
            }
        }
    } else {
        prev += dict.d[sa[i] - 1];

        ranks.clear();
        uint32_t start = wg.C[row.id];
        for (uint32_t j = start; j < start + row.n_occs; j++) {
            ranks.push_back(ilist[j - 1]);
        }
    }
    return row;
}

// void untunnel_2(tfm_index &wg, Dict &dict, size_t w, uint32_t *sa, int32_t *lcp, uint32_t *ilist) {
//     uint32_t n_special = dict.dwords + w + 1;
//     for (size_t i = n_special; i < dict.dsize; i++) {
//         row row_i = get_row(i, wg, ilist, dict, sa, lcp);
//         if (row_i.len <= w) continue;
//         if (row_i.is_phrase) {
//             for (size_t j = 0; j < row_i.n_prevs; j++) { L += row_i.prevs[j]; }
//             dout += 1;
//             for (size_t j = 1; j < row_i.n_prevs; j++) { dout += 0; }
//             din += 1;
//             for (size_t j = 1; j < row_i.p_occ; j++) { din += 0; }
//         } else {
//             if (lcp[i+1] >= row_i.len) {
//                 for (size_t j = 0; j < row_i.p_occ; j++) {
//                     ...
//                 }
//             }
//         }
//     }
// }

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

void print_table(tfm_index &wg, Dict &dict, size_t w, uint32_t *sa, int32_t *lcp) {
    uint32_t *ilist = new uint32_t[wg.L.size() - 1];
    generate_ilist(ilist, wg, dict.dwords);

    // make special symbols readable
    for (size_t i=0; i<dict.dsize; i++) {
        if (dict.d[i] == '\001') dict.d[i] = '$';
        if (dict.d[i] == '\002') dict.d[i] = '#';
    }

    cout << "i\tsa[i]\tlcp[i]\tlen\tphrase\tid\tid_in\tid_out\tranks\tprevs\tsuffix" << endl;
    uint32_t n_special = dict.dwords + w + 1; // dwords*$, w*# and 1*\000
    for (size_t i = n_special; i < dict.dsize; i++) {
        row_t row = get_row(i, wg, dict, w, sa, lcp, ilist);
        if (row.len <= (uint32_t)w) continue;

        string prev = "";
        vector<int32_t> ranks;

        if (row.is_phrase) {
            uint32_t start = wg.C[row.id];
            for (uint32_t j = start; j < start + row.id_in; j++) {
                if (wg.din[j] == 1) {
                    uint32_t pos = wg.dout_select(wg.din_rank(j + 1));
                    do {
                        if (wg.L[pos] == 0) pos = 0;
                        uint32_t act_phrase = wg.L[pos];    // F -> L

                        prev += dict.d[sa[act_phrase] - w - 1];
                        pos++;
                    } while (wg.dout[pos] == 0);
                }
            }
        } else {
            prev += dict.d[sa[i] - 1];

            ranks.clear();
            uint32_t start = wg.C[row.id];
            for (uint32_t j = start; j < start + row.id_in; j++) {
                ranks.push_back(ilist[j - 1]);
            }
        }

        // cout << "i\tsa[i]\tlcp[i]\tlen\tis_phrase\tid\tid_in\tid_out\tranks\tprevs\tsuffix" << endl;
        cout << i << "\t"
             << sa[i] << "\t"
             << lcp[i] << "\t"
             << row.len << "\t"
             << row.is_phrase << "\t"
             << row.id << "\t"
             << row.id_in << "\t"
             << row.id_out << "\t"
             << "[";

        for (int32_t x: ranks) cout << x;

        cout << "]\t"
             << prev << "\t"
             << dict.d + sa[i] << endl;
    }

    // make special symbols special
    for (size_t i=0; i<dict.dsize; i++) {
        if (dict.d[i] == '$') dict.d[i] = '\001';
        if (dict.d[i] == '#') dict.d[i] = '\002';
    }
}

size_t get_untunneled_size(tfm_index &wg, Dict &dict, size_t w, uint32_t *sa) {
    size_t size = 0;

    int32_t seqid = -1;
    uint32_t len = 0;
    uint64_t parse_occ = 0;

    for (uint64_t i = dict.dwords + w + 1; i < dict.dsize; i++) {
        // sa + 1 should be dict.end
        seqid = binsearch(sa[i], sa + 1, dict.dwords) + 1;

        // sa[seqid] is the length of the phrase seqid, check table
        len = sa[seqid] - sa[i];
        if (len <= (uint32_t)w) continue;

        parse_occ = wg.C[seqid + 1] - wg.C[seqid];
        size += parse_occ;
    }

    return size;
}

// int_vector<> compute_L(size_t w, uint8_t *d, long dsize, uint64_t *end_to_phrase, uint32_t *ilist, tfm_index &tfmp, long dwords, uint_t *sa, int_t *lcp) {
void compute_L(tfm_index &tfm, Dict dict, size_t w, uint32_t *sa, int32_t *lcp, int_vector<> &L) {
    uint32_t *ilist = new uint32_t[tfm.L.size() - 1];
    generate_ilist(ilist, tfm, dict.dwords);

    uint_t *eos = sa + 1;
    size_t p = 0;

    long next;
    uint32_t seqid;
    for (long i = dict.dwords + w + 1; i < dict.dsize; i = next) {
        next = i + 1;
        seqid = binsearch(sa[i], eos, dict.dwords);
        int_t suffixLen = eos[seqid] - sa[i];
        if (suffixLen <= (int_t)w) continue;

        if (sa[i] == 0 || dict.d[sa[i] - 1] == EndOfWord) {
            // ----- simple case: the suffix is a full word
            uint32_t start = tfm.C[seqid + 1];
            uint32_t end = tfm.C[seqid + 2];
            for (uint32_t j = start; j < end; j++) {
                if (tfm.din[j] == 1) {
                    uint32_t pos = tfm.dout_select(tfm.din_rank(j + 1));
                    do {
                        if (tfm.L[pos] == 0) pos = 0;
                        uint32_t act_phrase = tfm.L[pos] - 1;
                        uint8_t char_to_write = dict.d[dict.end[act_phrase] - w - 1];
                        // out.push_back(char_to_write);
                        L[p] = char_to_write;
                        p++;
                    } while (tfm.dout[++pos] != 1);
                }
            }
        } else {
            // ----- hard case: there can be a group of equal suffixes starting
            // at i save seqid and the corresponding char
            vector<uint32_t> id2merge(1, seqid);
            vector<uint8_t> char2write(1, dict.d[sa[i] - 1]);
            while (next < dict.dsize && lcp[next] >= suffixLen) {
                seqid = binsearch(sa[next], eos, dict.dwords);
                int_t nextsuffixLen = eos[seqid] - sa[next];
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
                    for (uint64_t j = tfm.C[s]; j < tfm.C[s + 1]; j++) {
                        // out.push_back(char2write[0]);
                        L[p] = char2write[0];
                        p++;
                    }
                }
            } else {
                // many words, many chars...
                vector<SeqId> heap; // create heap
                for (size_t i = 0; i < numwords; i++) {
                    uint32_t s = id2merge[i] + 1;
                    heap.push_back(SeqId(
                        s, tfm.C[s + 1] - tfm.C[s], ilist + (tfm.C[s] - 1),
                        char2write[i]
                    ));
                }
                std::make_heap(heap.begin(), heap.end());
                while (heap.size() > 0) {
                    // output char for the top of the heap
                    SeqId s = heap.front();
                    // out.push_back(s.char2write);
                    L[p] = s.char2write;
                    p++;
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
    return;
}

void compute_degrees(tfm_index &tfmp, Dict &dict, size_t w, uint32_t *sa, int32_t *lcp, bit_vector &din, bit_vector &dout) {
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
        seqid = binsearch(sa[i], eos, dwords);
        int32_t suffixLen = eos[seqid] - sa[i];
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
                seqid = binsearch(sa[next], eos, dwords);
                int32_t nextsuffixLen = eos[seqid] - sa[next];
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

tfm_index unparse(tfm_index &wg_parse, Dict &dict, size_t w, size_t size) {

    uint32_t *sa_d = new uint32_t[dict.dsize];
    int32_t *lcp_d = new int32_t[dict.dsize];
    // separators s[i]=1 and with s[n-1]=0
    // cout << dict.d << "\n" << dict.dsize << endl;;
    gsacak(dict.d, sa_d, lcp_d, NULL, dict.dsize);

    size_t s = get_untunneled_size(wg_parse, dict, w, sa_d);
    int_vector<> L(s);
    bit_vector din(s + 1, 1);
    bit_vector dout(s + 1, 1);

    print_table(    wg_parse, dict, w, sa_d, lcp_d);
    compute_L(      wg_parse, dict, w, sa_d, lcp_d, L);
    compute_degrees(wg_parse, dict, w, sa_d, lcp_d, din, dout);

    tfm_index tfm(size, L, din, dout);
    return tfm;
}

void print_wg(tfm_index &wg) {
    for (uint i=0; i < wg.L.size(); i++)
        cout << (char)wg.L[i];
    cout << "\n";

    cout << wg.dout << endl;
    cout << wg.din << endl << endl;
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

    print_wg(tfm);
    print_wg(unparsed);

    cout << input << endl;
    cout << untunnel(unparsed) << endl;
}
