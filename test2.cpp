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
    for (tfm_index::size_type i = 1; i < tfm.size(); i++) {
        char c = (char)tfm.backwardstep(p);
        original[tfm.size() - i - 1] = c;
    }
    original[tfm.size() - 1] = (char)tfm.preceding_char(p);

    return original;
}

void tmp(tfm_index &wg, Dict &dict, size_t w, uint32_t *sa, int32_t *lcp, uint32_t *ilist) {
    uint32_t seqid;
    cout << dict.dwords + w + 1;
    for (uint64_t i = 0; i < dict.dsize; i++) {
        // cout << dict.d + sa[i] << endl;
        // int32_t suffixLen = getlen(sa[i], sa + 1, dict.dwords, &seqid);
        // if (suffixLen <= (int32_t)w) continue;
        // char prev = dict.d[sa[i] - 1];
        // uint64_t parse_occ = wg.C[seqid + 1] - wg.C[seqid];
        // uint32_t *rank = ilist + (wg.C[seqid] - 1);

        // printf("%.*s\n", suffixLen, dict.d + sa[i]);
        // cout << prev << "\t" << nextsuffixLen << "\n";
    }
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

tfm_index unparse(tfm_index &wg_parse, Dict &dict, size_t w, size_t size) {
    uint32_t *inverted_list = new uint32_t[wg_parse.L.size() - 1];
    generate_ilist(inverted_list, wg_parse, dict.dwords);

    uint32_t *sa_d = new uint32_t[dict.dsize];
    int32_t *lcp_d = new int32_t[dict.dsize];
    // separators s[i]=1 and with s[n-1]=0
    // cout << dict.d << "\n" << dict.dsize << endl;;
    gsacak(dict.d, sa_d, lcp_d, NULL, dict.dsize);
    dict.d[0] = 0;

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
    dict.d = (uint8_t *)malloc(dict.dsize);     // {2GT, GTAGGT, GTGGGT, GTTGGT, GTTT22}
    dict.end = (uint64_t *)malloc(dict.dwords); // {3, 10, 17, 24, 31}

    uint8_t d[] = "\002GT\001GTAGGT\001GTGGGT\001GTTGGT\001GTTT\002\002\001\000";
    uint64_t end[] = {3, 10, 17, 24, 31};

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