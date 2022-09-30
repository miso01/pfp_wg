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
    if (argc < 3) {
        printUsage(argv);
        cerr << "At least 2 parameter expected" << endl;
        return 1;
    }

    tfm_index loaded;
    load_from_file(loaded, argv[1]);
    string filename = argv[2];
    untunnel(loaded, filename);
}
