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

    char *S = new char[tfm.size()];
    auto p = tfm.end();
    for (size_type i = 0; i < tfm.size(); i++) {
        char c = (char)tfm.backwardstep(p);
        S[tfm.size() - i - 1] = c;
        // cout << "\t pos:" << p.first << " c: " << c << endl;
    }

    char *name;
    int e = asprintf(&name, "%s.%s", argv[1], "untunneled");
    if (e == -1) {
        cout << "ERROR while creating the output file name!" << endl;
        return 1;
    }
    FILE *fout = fopen(name, "w");
    fwrite(S, sizeof(char), tfm.size(), fout);
    fclose(fout);
}
