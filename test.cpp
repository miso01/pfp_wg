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

void visualize_wt_rec(const wt_blcd<>& wt, typename wt_blcd<>::node_type v, size_t level, vector<string>& out)
{
    if (!wt.is_leaf(v)) {
        if (out.size() < level+4) {
            out.push_back("");
            out.push_back("");
        }
        while (out[level+2].size() < out[level].size()) {
            out[level+2] += " ";
            out[level+3] += " ";
        }

        auto vs = wt.expand(v);
        if (!wt.empty(vs[0])) {
            visualize_wt_rec(wt, vs[0], level+2, out);
        }
        if (!wt.empty(vs[0]) and !wt.empty(vs[1])) {
            out[level+2] += " ";
            out[level+3] += " ";
        }
        if (!wt.empty(vs[1])) {
            visualize_wt_rec(wt, vs[1], level+2, out);
        }

        size_t begin = out[level].size();
        size_t end   = out[level+2].size();
        size_t size  = wt.size(v);
        size_t delta = (end-begin)-size;

        for (size_t i=0; i < delta/2; ++i) {
            out[level] += " ";
            out[level+1] += " ";
        }
        auto seq_vec = wt.seq(v);
        auto bit_it = wt.bit_vec(v).begin();
        for (auto it = seq_vec.begin(); it!=seq_vec.end(); ++it, ++bit_it) {
            out[level]   += *it;
            out[level+1] += *bit_it ? "1" : "0";
        }

        for (size_t i=0; i < (delta+1)/2; ++i) {
            out[level] += " ";
            out[level+1] += " ";
        }
    } else {
        auto seq = wt.seq(v);
        for (auto it = seq.begin(); it!=seq.end(); ++it) {
            out[level] += *it;
            out[level+1] += " ";
        }
    }
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

    // int_vector<> L{1, 4, 5};
    vector<bool> bdin{false, true, false, false, false, false, false, false, false, false};
    bit_vector din = create_from_boolvec(bdin);
    vector<bool> bdout{false, true, false, false, false, false, false, false, false, false};
    bit_vector ddout = create_from_boolvec(bdout);

    int_vector<8> L{1, 2, 3, 5, 9, 4};

    cout << L.size() << "\n" << (uint64_t)L.width() << endl;

    wt_blcd<> wt;
    construct_im(wt, L);

    vector<string> vs(2,"");
    visualize_wt_rec(wt, wt.root(), 0, vs);
    for (size_t i=0; i<vs.size(); ++i)
        cout<<vs[i]<<endl;

    tfm_index tfm;
    // tfm.L = wt;
    // tfm.C =
    // tfm.din =
    // tfm.dout =

}