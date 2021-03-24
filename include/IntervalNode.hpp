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
    bool isLeaf{false};
    int32_t median{};
    std::vector<Record> records{};
    std::unique_ptr<IntervalNode> lNode{nullptr};
    std::unique_ptr<IntervalNode> rNode{nullptr};
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
     //!\}

     std::unique_ptr<IntervalNode> & get_left_node();

     std::unique_ptr<IntervalNode> & get_right_node();

     int32_t const & get_median();

     std::vector<Record> & get_records();

     void set_is_leaf();

     void set_median(int32_t const & m);

     void print(int32_t level);

};
int32_t calculate_median(std::vector<Record> const & records_i);

void construct_tree(std::unique_ptr<IntervalNode> & node, std::vector<Record> const & records_i);

std::vector<Record> overlap(std::unique_ptr<IntervalNode> & root, int32_t start, int32_t end);
