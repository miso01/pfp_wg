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

vector<uint64_t> get_C(int_vector<> &L, size_t sigma) {
    // 255 handles char alphabet, sigma + 1 handles int alphabets
    vector<uint64_t> v(max((size_t)255, sigma + 1), 0);
    for (uint64_t i = 0; i < L.size(); i++) v[L[i] + 1] += 1;
    for (uint64_t i = 0; i < v.size() - 1; i++) v[i + 1] += v[i];
    return v;
}

tfm_index construct(size_t size, int_vector<> &L, bit_vector &din, bit_vector &dout) {
    tfm_index tfm;

    tfm.text_len = size;
    construct_im(tfm.m_L, L);
    tfm.m_C = get_C(L, tfm.m_L.sigma);
    tfm.m_dout = dout;
    tfm.m_din = din;

    sdsl::util::init_support(tfm.m_dout_rank,   &tfm.m_dout);
    sdsl::util::init_support(tfm.m_dout_select, &tfm.m_dout);
    sdsl::util::init_support(tfm.m_din_rank,    &tfm.m_din);
    sdsl::util::init_support(tfm.m_din_select,  &tfm.m_din);

    return tfm;
}

int main() {
    size_t size = 12;
    int_vector<> L = {'i', 'p', 's', 's', 'm', '$', 'p', 'i', 's', 's', 'i', 'i'};
    bit_vector in{1,1,1,1,1,1,1,1,1,1,1,1,1};
    bit_vector out{1,1,1,1,1,1,1,1,1,1,1,1,1};

    cout << "L size: " << L.size() << "\tL width: " << (size_t)L.width() << endl;

    tfm_index tfm = construct(size, L, in, out);
    cout << untunnel(tfm) << endl;
}