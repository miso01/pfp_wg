#include <cstdio>
#include <deque>
#include <iostream>
#include <map>
#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>
#include <utility>
#include <vector>

#include "tfm_index.hpp"

using namespace std;
using namespace sdsl;

typedef typename sdsl::int_vector<>::size_type size_type;

void printUsage(char **argv) {
    cerr << "USAGE: " << argv[0] << " TFMFILE" << endl;
    cerr << "TFMFILE:" << endl;
    cerr << "  File where to store the serialized trie" << endl;
};

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

tfm_index create_tfm(size_t size, int_vector<8> &L, bit_vector &din, bit_vector &dout) {
    tfm_index tfm;
    tfm.text_len = size;

    string tmp = "tmp2.L";
    FILE *fbwt = fopen(tmp.c_str(), "wb");
    for (char c: L) { fputc(c, fbwt); }
    fclose(fbwt);
    int_vector_buffer<> buf(tmp, std::ios::in, L.size(), 8, true);
    tfm.m_L = tfm_index::wt_type(buf, L.size());
    remove(tmp);

    tfm.m_C = vector<uint64_t>(255, 0);
    for (uint64_t i = 0; i < L.size(); i++) tfm.m_C[L[i] + 1] += 1;
    for (uint64_t i = 0; i < tfm.m_C.size() - 1; i++) tfm.m_C[i + 1] += tfm.m_C[i];

    tfm.m_dout = tfm_index::bit_vector_type(std::move(dout));
    sdsl::util::init_support(tfm.m_dout_rank, &tfm.m_dout);
    sdsl::util::init_support(tfm.m_dout_select, &tfm.m_dout);

    tfm.m_din = tfm_index::bit_vector_type(std::move(din));
    sdsl::util::init_support(tfm.m_din_rank, &tfm.m_din);
    sdsl::util::init_support(tfm.m_din_select, &tfm.m_din);

    return tfm;
};

tfm_index construct_from_pfwg(const string basename) {

    size_t orig_size = 12156306;

    int_vector<8> L;
    load_vector_from_file(L, basename + ".L", 1);
    // uint64_t size = ;
    // sdsl::int_vector<1> din, dout;

    bit_vector din;
    din.resize(L.size() + 1);
    load_bitvector(din, basename + ".din", L.size() + 1);

    bit_vector dout;
    dout.resize(L.size() + 1);
    load_bitvector(dout, basename + ".dout", L.size() + 1);

    tfm_index tfm = create_tfm(orig_size, L, din, dout);
    return tfm;
}

void untunnel(tfm_index &tfm, string &filename) {
    char *original = new char[tfm.size()];
    auto p = tfm.end();
    for (size_type i = 0; i < tfm.size(); i++) {
        char c = (char)tfm.backwardstep(p);
        original[tfm.size() - i - 1] = c;
    }

    FILE *fout = fopen(filename.c_str(), "w");
    fwrite(original, sizeof(char), tfm.size(), fout);
    fclose(fout);
}

int main(int argc, char **argv) {
    // check parameters
    if (argc < 2) {
        printUsage(argv);
        cerr << "At least 1 parameter expected" << endl;
        return 1;
    }

    // find original string size
    sdsl::int_vector<8> original;
    load_vector_from_file(original, argv[1], 1);
    cout << original.size() << endl;

    tfm_index tfm = construct_from_pfwg(argv[1]);

    string filename = argv[1];
    filename += ".untunneled";
    untunnel(tfm, filename);
}
