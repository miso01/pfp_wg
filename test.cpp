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

bit_vector create_from_boolvec(vector<bool> &v) {
    bit_vector b(v.size(), 0);
    for (size_t i=0; i < v.size(); i++) {
        b[i] = v[i];
    }
    return b;
}

int main() {
    bit_vector b(10, 0);
    b[1] = 1;

    store_to_file(b, "tmp");

    bit_vector b2;
    load_from_file(b2, "tmp");

    vector<bool> b3{false, true, false, false, false, false, false, false, false, false};
    bit_vector b4 = create_from_boolvec(b3);
    store_to_file(b, "tmp2");


}