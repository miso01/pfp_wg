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

struct Dict {
    uint8_t *d; // pointer to the dictionary
    uint64_t *end; // end[i] is the index of the ending symbol of the i-th phrase
    uint64_t dsize;  // dicionary size in symbols
    uint64_t dwords; // the number of phrases of the dicionary
};

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

int main() {
    size_t size = 9;
    int_vector<> L = {5, 0, 1, 4, 2, 3, 4};
    bit_vector in{1, 1, 1, 1, 1, 0, 1, 1};
    bit_vector out{1, 1, 1, 0, 1, 1, 1, 1};

    // cout << "L size: " << L.size() << "\tL width: " << (size_t)L.width() << endl;

    tfm_index tfm(size, L, in, out);
    // cout << untunnel(tfm) << endl;

    string T = "#ABDACDAEDABDACDAEDAZ$";
    // {#A, ABDA, ACDA, AEDA, AZ$}
    // parse = [1, 2, 3, 4, 2, 3, 4, 5]
    int w = 1;

    Dict dict;
    dict.dsize = 22;

    uint8_t *tmp_var = (uint8_t *)malloc(10);
    uint8_t tmp2[] = "#A\001ABDA\0ACDA\0AEDA\0AZ$\0";
    strcpy((char *)tmp_var, (char *)tmp2);
    dict.d = tmp_var;

    dict.dwords = 5;

    dict.end = (uint64_t *)malloc(10);
    dict.end[0] = 2;
    dict.end[1] = 7;

    uint32_t *sa;
    int32_t *lcp;
    uint32_t *ilist;

    tmp(tfm, dict, w, sa, lcp, ilist);
}