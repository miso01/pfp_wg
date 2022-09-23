#ifndef TFM_INDEX_HPP
#define TFM_INDEX_HPP

#include <assert.h>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/construct_lcp_helper.hpp>
#include <sdsl/csa_wt.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/int_vector_buffer.hpp>
#include <sdsl/io.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/sdsl_concepts.hpp>
#include <sdsl/select_support.hpp>
#include <sdsl/util.hpp>
#include <sdsl/wavelet_trees.hpp>

#include <algorithm>
#include <iostream>
#include <limits>
#include <string>
#include <utility>


//! a class representing a tunneled fm-index
class tfm_index {
  public:
    typedef sdsl::int_vector<> text_type;
    typedef sdsl::int_vector<>::size_type size_type;

    typedef sdsl::wt_blcd_int<> wt_type;
    typedef sdsl::wt_blcd_int<>::value_type value_type;
    typedef sdsl::wt_blcd_int<>::bit_vector_type bit_vector_type;
    typedef sdsl::wt_blcd_int<>::bit_vector_type::rank_1_type rank_type;
    typedef sdsl::wt_blcd_int<>::bit_vector_type::select_1_type select_type;

    // first index is next outgoing edge, second index is tunnel entry offset
    typedef std::pair<size_type, size_type> nav_type;

  private:
    friend tfm_index create_tfm(size_t size, sdsl::int_vector_buffer<> &L_buf, sdsl::bit_vector &din, sdsl::bit_vector &dout);
    // friend void construct_from_pfwg(tfm_index &tfm_index, const std::string filename);
    // friend tfm_index construct_from_pfwg(const std::string filename);
    friend tfm_index create_tfm(size_t size, sdsl::int_vector<8> L, sdsl::bit_vector din, sdsl::bit_vector dout);

    size_type text_len; // original textlen
    sdsl::wt_blcd_int<> m_L;
    std::vector<size_type> m_C;
    sdsl::wt_blcd_int<>::bit_vector_type m_dout;
    rank_type m_dout_rank;
    select_type m_dout_select;
    sdsl::wt_blcd_int<>::bit_vector_type m_din;
    rank_type m_din_rank;
    select_type m_din_select;

  public:
    const wt_type &L = m_L;
    const std::vector<size_type> &C = m_C;
    const bit_vector_type &dout = m_dout;
    const rank_type &dout_rank = m_dout_rank;
    const select_type &dout_select = m_dout_select;
    const bit_vector_type &din = m_din;
    const rank_type &din_rank = m_din_rank;
    const select_type &din_select = m_din_select;

    //! returns the size of the original string
    size_type size() const { return text_len; };

    //! returns the end, i.e. the position in L where the string ends
    nav_type end() const { return std::make_pair((size_type)0, (size_type)0); }

    //! returns the character preceding the current position
    value_type preceding_char(const nav_type &pos) const {
        return L[pos.first];
    }

    //! Operation performs an backward step from current position.
    //! function sets posm to the new value and returns the result
    //! of preceding_char( pos ) before the backward step was performed
    value_type backwardstep(nav_type &pos) const {
        size_type &i = pos.first; // create references into position pair
        size_type &o = pos.second;

        // navigate to next entry
        auto is = L.inverse_select(i);
        auto c = is.second;
        i = C[c] + is.first;

        // check for the start of a tunnel
        auto din_rank_ip1 = din_rank(i + 1);
        if (din[i] == 0) {
            o = i -
                din_select(din_rank_ip1); // save offset to uppermost entry edge
        }
        // navigate to outedges of current node
        i = dout_select(din_rank_ip1);

        // check for end of a tunnel
        if (dout[i + 1] == 0) {
            i += o; // jump back offset
            o = 0;
        }
        return c;
    };

    //! serializes opbject
    size_type serialize(
        std::ostream &out, sdsl::structure_tree_node *v, std::string name
    ) const {

        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(
            v, name, sdsl::util::class_name(*this)
        );
        size_type written_bytes = 0;
        written_bytes += sdsl::write_member(text_len, out, child, "text_len");
        written_bytes += m_L.serialize(out, child, "L");
        written_bytes += sdsl::serialize(m_C, out, child, "C");
        written_bytes += m_dout.serialize(out, child, "dout");
        written_bytes += m_dout_rank.serialize(out, child, "dout_rank");
        written_bytes += m_dout_select.serialize(out, child, "dout_select");
        written_bytes += m_din.serialize(out, child, "din");
        written_bytes += m_din_rank.serialize(out, child, "din_rank");
        written_bytes += m_din_select.serialize(out, child, "din_select");
        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    };

    //! loads a serialized object
    void load(std::istream &in) {
        sdsl::read_member(text_len, in);
        m_L.load(in);
        sdsl::load(m_C, in);
        m_dout.load(in);
        m_dout_rank.load(in, &m_dout);
        m_dout_select.load(in, &m_dout);
        m_din.load(in);
        m_din_rank.load(in, &m_din);
        m_din_select.load(in, &m_din);
    };
};

#endif
