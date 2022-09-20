#include <cstdio>
#include <deque>
#include <iostream>
#include <map>
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

void symbol_frequencies(std::vector<uint64_t> &C, sdsl::int_vector<8> &L) {
    C = std::vector<uint64_t>(255, 0); // lol I hope it's enough :D
    for (uint64_t i = 0; i < L.size(); i++)
        C[L[i] + 1] += 1;
    for (uint64_t i = 0; i < C.size() - 1; i++) {
        C[i + 1] += C[i];
    }
}

void construct_from_pfwg(tfm_index &tfm_index, const std::string basename) {
    // find original string size
    sdsl::int_vector<8> original;
    load_vector_from_file(original, basename, 1);

    sdsl::int_vector<8> L;
    load_vector_from_file(L, basename + ".L", 1);
    uint64_t size = L.size();
    sdsl::int_vector<1> din, dout;
    din.resize(size + 1);
    dout.resize(size + 1);
    // one additional bit at the end
    load_bitvector(din, basename + ".din", size + 1);
    load_bitvector(dout, basename + ".dout", size + 1);

    typedef ::tfm_index::wt_type wt_type;
    typedef ::tfm_index::bit_vector_type bv_type;

    tfm_index.text_len = original.size();
    sdsl::int_vector_buffer<> buf(basename + ".L", std::ios::in, size, 8, true);
    tfm_index.m_L = wt_type(buf, size);

    symbol_frequencies(tfm_index.m_C, L);
    tfm_index.m_dout = bv_type(std::move(dout));
    sdsl::util::init_support(tfm_index.m_dout_rank, &tfm_index.m_dout);
    sdsl::util::init_support(tfm_index.m_dout_select, &tfm_index.m_dout);

    tfm_index.m_din = bv_type(std::move(din));
    sdsl::util::init_support(tfm_index.m_din_rank, &tfm_index.m_din);
    sdsl::util::init_support(tfm_index.m_din_select, &tfm_index.m_din);
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

    // load tunneled fm index
    tfm_index tfm;
    construct_from_pfwg(tfm, argv[1]);

    string filename = argv[1];
    filename += ".untunneled";
    untunnel(tfm, filename);
}
