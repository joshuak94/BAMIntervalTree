#pragma once
#include <seqan3/io/sam_file/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "Record.hpp"

class IntervalNode
{
using my_fields = seqan3::fields<seqan3::field::id,
                                 seqan3::field::ref_id,
                                 seqan3::field::ref_offset,
                                 seqan3::field::cigar>;
using sam_file_input_t = seqan3::sam_file_input<seqan3::sam_file_input_default_traits<>,
                                                my_fields,
                                                seqan3::type_list<seqan3::format_sam, seqan3::format_bam>>;
private:
    bool isLeaf{};
    std::vector<Record> records{};
    IntervalNode * lNode{};
    IntervalNode * rNode{};
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr IntervalNode()                        = default; //!< Defaulted.
    IntervalNode(IntervalNode const &)              = default; //!< Defaulted.
    IntervalNode(IntervalNode &&)                   = default; //!< Defaulted.
    IntervalNode & operator=(IntervalNode const &)  = default; //!< Defaulted.
    IntervalNode & operator=(IntervalNode &&)       = default; //!< Defaulted.
    ~IntervalNode()                                 = default; //!< Defaulted.

    IntervalNode(IntervalNode * lNode_i, IntervalNode * rNode_i) :
                 lNode{std::move(lNode_i)}, rNode{std::move(rNode_i)} {}
     //!\}

     void set_left_node(IntervalNode * new_node);

     void set_right_node(IntervalNode * new_node);

     IntervalNode * get_left_node();

     IntervalNode * get_right_node();

     std::vector<Record> & get_records();

     IntervalNode * construct_tree(std::vector<Record> & records_i);

     void print(int32_t level);

     int32_t calculate_median(std::vector<Record> & records_i);
};
