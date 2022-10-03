#include <cstddef>
#include <cstdint>
#include <cstdlib>
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
    for (tfm_index::size_type i = 0; i < tfm.size(); i++) {
        char c = (char)tfm.backwardstep(p);
        original[tfm.size() - i - 1] = c;
    }

    return original;
}

int main() {
    size_t size = 12;
    int_vector<> L = {'i', 'p', 's', 's', 'm', '$', 'p', 'i', 's', 's', 'i', 'i'};
    bit_vector in{1,1,1,1,1,1,1,1,1,1,1,1,1};
    bit_vector out{1,1,1,1,1,1,1,1,1,1,1,1,1};

    cout << "L size: " << L.size() << "\tL width: " << (size_t)L.width() << endl;

    tfm_index tfm(size, L, in, out);
    cout << untunnel(tfm) << endl;
}