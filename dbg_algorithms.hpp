#ifndef DEBRUIJNGRAPH_ALGORITHMS_HPP
#define DEBRUIJNGRAPH_ALGORITHMS_HPP

#include <sdsl/csa_wt.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/wavelet_trees.hpp>

#include <deque>
#include <limits>
#include <utility>
#include <vector>

class dbg_algorithms {
  public:
    typedef sdsl::int_vector<>::size_type size_type;

  private:
    typedef typename std::deque<std::pair<size_type, size_type>> kmer_queue;

    // edge minimization algorithm
    //  stop_on_min: wether algorithm should stop if minimum is certain or
    //  proceed to kmax csa: fm-index max_k: maximal order k to which
    //  minimization should be performed B: a bitvector of size csa.size()
    //  filled with ones nb_buffer: a buffer used to write the node bound
    //  indices occuring during the execution of the algorithm
    //  order_edgecnt_output: function that is called for every order k with the
    //  current edgecount and the order
    template <bool stop_on_min = true, class f_order_edgecnt_output>
    static std::pair<size_type, size_type> minimize_dbg_edges(
        const sdsl::wt_blcd_int<> &csa, const std::vector<uint64_t> &C,
        size_type kmax, sdsl::bit_vector &B,
        sdsl::int_vector_buffer<> &nb_buffer,
        f_order_edgecnt_output order_edgecnt_output
    );

  public:
    // function behaves as calling dbg_edgecount for each values from 1 to
    // (including) k,
    //  but in a more efficient manner.
    //  csa: a compressed suffix array of the string using a wavelet tree and
    //  the BWT k: maximal k from which the number of edges of the DBG shall be
    //  determined function returns a vector of size k containing the edge count
    //  of each order k, where number of edges of an order k edge-reduced graph
    //  can be queried at index (k-1).
    static std::vector<size_type> dbg_edgespectrum(
        const sdsl::wt_blcd_int<> &csa, const std::vector<uint64_t> &C,
        size_type k
    ) {
        assert(k > 0 && k <= csa.size());

        // initialize B
        sdsl::bit_vector B(csa.size() + 1, 0);

        // edge spectrum vector
        std::vector<size_type> ES(k);

        // send node bounds to nothing
        sdsl::int_vector_buffer<> nb_buffer("/dev/null", std::ios::out);

        // compute spectrum
        minimize_dbg_edges<false>(
            csa, C, k, B, nb_buffer,
            [&ES](size_type k, size_type m) {
                if (k > 0u)
                    ES[k - 1] = m;
            }
        );

        return ES;
    };

    // function computes the order k such that the number of edges of
    // the edge-reduced DeBruijn graph of order k are minimal under all possible
    // values of k csa: a compressed suffix array of the string using a wavelet
    // tree and the BWT B: a bitvector where the upper bounds of the k-mer
    // intervals of the smallest DBG will be marked.
    //    After function call, B has size csa.size() + 1, and B[csa.size()] will
    //    always be 1
    // config: a config indicating where temporary files should be saved
    // function returns a pair of two integers, where the first integer
    // corresponds to the order k whilst the second integer indicates the
    // minimal amount of edges.
    static std::pair<size_type, size_type> find_min_dbg(
        const sdsl::wt_blcd_int<> &csa, const std::vector<uint64_t> &C,
        sdsl::bit_vector &B
    ) {
        // prepare bitvector B
        // initialize B
        B.resize(csa.size() + 1);
        sdsl::util::set_to_value(B, 0);

        // prepare buffer with node bounds
        std::string tmp_file_name = "find_min_dbg.tmp";
        sdsl::int_vector_buffer<> nb_buffer(tmp_file_name, std::ios::out);

        auto result = minimize_dbg_edges<true>(
            csa, C, csa.size(), B, nb_buffer, [](size_type, size_type) {}
        );

        // retain node bounds of optimal solution
        for (auto nb : nb_buffer)
            B[nb] = 0;

        // delete node bound buffer
        nb_buffer.close();
        sdsl::remove(tmp_file_name);

        return result;
    };

    // function computes a prefix interval marking.
    // starts of a possible prefix interval are indicated by a 10...0 sequence
    // in din, ends of a possible prefix interval are indicated by a 10...0
    // sequence in dout. in case that dout and din are copies of bitvector B
    // from find_min_dbg function, this computes a marking of k-mer prefix
    // intervals.
    static void mark_prefix_intervals(
        const sdsl::wt_blcd_int<> &csa, const std::vector<uint64_t> &C,
        sdsl::bit_vector &dout, sdsl::bit_vector &din
    ) {

        // variables needed for interval_symbols function
        typename sdsl::wt_blcd_int<>::size_type iv_chars;
        std::vector<long long unsigned int> cs(csa.sigma);
        std::vector<size_type> rank_c_l(csa.sigma);
        std::vector<size_type> rank_c_r(csa.sigma);

        size_type i = 0;
        for (size_type j = 1; j < din.size(); j++) {
            if (din[j] == 1) { //[i,j) is a start of a possible prefix interval

                sdsl::interval_symbols(
                    csa, i, j, iv_chars, cs, rank_c_l, rank_c_r
                );

                // check if start has multiple predecessors (iv_chars > 1)
                // or [LF[i],LF[j-1]] is no possible end of a tunnel
                size_type lb = C[cs[0]] + rank_c_l[0];
                size_type rb = C[cs[0]] + rank_c_r[0];

                bool pi = (iv_chars == 1) && (dout[lb] == 1) &&
                          (dout[rb] == 1); // pi: possible prefix interval
                while (pi && (++lb < rb)) {
                    pi = pi && (dout[lb] == 0);
                }

                // remove markings in case that no prefix interval is possible
                if (!pi) {
                    // clean zeros, no column of a prefix interval
                    for (size_type c = 0; c < iv_chars; c++) {
                        auto Cc = C[cs[c]];
                        do {
                            din[i++] = 1;
                            dout[Cc + rank_c_l[c]++] = 1;
                        } while (rank_c_l[c] < rank_c_r[c]);
                    }
                }
                i = j; // process next interval
            }
        }
    };
};

///////////////////////////////////////////////////////////////////////////////

//// EDGE MINIMIZATION ALGORITHM ////
template <bool stop_on_min, class f_order_edgecnt_output>
std::pair<dbg_algorithms::size_type, dbg_algorithms::size_type>
dbg_algorithms::minimize_dbg_edges(
    const sdsl::wt_blcd_int<> &csa, const std::vector<uint64_t> &C,
    size_type kmax, sdsl::bit_vector &B, sdsl::int_vector_buffer<> &nb_buffer,
    f_order_edgecnt_output order_edgecnt_output
) {
    assert(csa.size() < std::numeric_limits<size_type>::max());
    // ensure some input is given
    if (csa.size() == 0) {
        return std::make_pair((size_type)0, (size_type)0);
    }
    size_type n = csa.size(); // initialize
    sdsl::bit_vector F(n);
    std::vector<std::deque<std::pair<size_type, size_type>>> Q(csa.sigma);
    std::vector<size_type> qsize(csa.sigma);

    B[0] = B[n] = 1;          // boundaries of the root node
    size_type nG = 1, m = n;  // node counter and edge counter
    size_type ks = 1, ms = n; // order with minimum number of edges ms
    size_type fusible = 0;    // inherited fusions

    // variables needed for interval_symbols function
    size_type iv_chars;
    // std::vector<unsigned char> cs ( csa.sigma );
    std::vector<long long unsigned int> cs(csa.sigma);
    std::vector<size_type> rank_c_i(csa.sigma);
    std::vector<size_type> rank_c_j(csa.sigma);

    // queue initialization
    Q.front().emplace_back(0, csa.size() - 1);

    for (size_type k = 0; k <= kmax; k++) {
        fusible = 0;
        for (size_type c = 0; c < csa.sigma; c++) {
            qsize[c] = Q[c].size();
        }
        for (size_type c = 0; c < csa.sigma; c++) { // in alphabetical order
            for (size_type q = 0; q < qsize[c]; q++) {
                auto iv = Q[c].front();
                Q[c].pop_front();
                size_type lb = iv.first, rb = iv.second;

                sdsl::interval_symbols(
                    csa, lb, rb + 1, iv_chars, cs, rank_c_i, rank_c_j
                );
                if (iv_chars == 1) { // reclassify edges
                    auto c = cs[0];
                    auto i = C[c] + rank_c_i[0];
                    auto j = C[c] + rank_c_j[0] - 1;
                    if (B[i] == 1 && B[j + 1] == 1) { // node has no siblings
                        m -= (rb - lb);
                    } else {
                        fusible += (rb - lb); // possible fusion in next order
                    }
                    F[rb] = 1;
                }
                for (size_type c_i = 0; c_i < iv_chars; c_i++) {
                    auto c = cs[c_i];
                    auto i = C[c] + rank_c_i[c_i];
                    auto j = C[c] + rank_c_j[c_i] - 1;
                    if (B[i] == 0 || B[j + 1] == 0) {
                        if (B[j + 1] == 0) {
                            nG++;
                            if (stop_on_min && nG >= ms) {
                                return std::make_pair(ks, ms);
                            }
                        }
                        Q[c].emplace_back(i, j);
                    }
                }
            }
        }
        order_edgecnt_output(k, m);
        if (m < ms) {
            ks = k;
            ms = m;
            nb_buffer.reset(); // EXTERNAL CLEAR
        }
        m -= fusible; // establish new inherited fusions
        size_type last = std::numeric_limits<size_type>::max();
        for (size_type c = 0; c < csa.sigma; c++) {
            for (auto iv : Q[c]) {
                size_type lb = iv.first, rb = iv.second;
                if (B[rb + 1] == 0) {
                    B[rb + 1] = 1;
                    nb_buffer.push_back(rb + 1); // EXTERNAL WRITE
                    if (last == std::numeric_limits<size_type>::max()) {
                        last = lb;
                    }
                } else {
                    if (F[rb] == 1) {
                        m += (rb - last);
                        F[rb] = 0;
                    }
                    last = std::numeric_limits<size_type>::max();
                }
            }
        }
    }
    return std::make_pair(ks, ms);
}

#endif
