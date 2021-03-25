#include <gtest/gtest.h>
#include <math.h>

#include "IntervalNode.hpp"

// Recursive function to go through the tree and check the medians.
void check_medians(std::unique_ptr<IntervalNode> & root, int level, int pos)
{
    std::vector<int32_t> expected_medians{944, 467, 1397, 257, 676, 1132, 1664, 140, 363, 571, 802, 1032, 1257, 1531,
                                          1822, 77, 197, 0, 414, 519, 624, 730, 886, 0, 0, 1194, 1323, 1464, 1593, 1742,
                                          1928};
    if (expected_medians[pow(2, level) - 1 + pos] > 0)
    {
        EXPECT_EQ(root->get_median(), expected_medians[pow(2, level) - 1 + pos]);
        if(root->get_left_node())
        {
            check_medians(root->get_left_node(), level + 1, (pos * 2));
        }
        if(root->get_right_node())
        {
            check_medians(root->get_right_node(), level + 1, (pos * 2) + 1);
        }
    }

}

TEST(tree_construct, simulated_chr1_small_golden)
{
    // Construct Tree
    using my_fields = seqan3::fields<seqan3::field::id,
                                     seqan3::field::ref_id,
                                     seqan3::field::ref_offset,
                                     seqan3::field::cigar,
                                     seqan3::field::seq>;

    seqan3::sam_file_input input{DATADIR"simulated_chr1_small_golden.bam", my_fields{}};

    std::vector<Record> records{};
    for (auto & r : input)
    {
        int32_t start = std::get<2>(r).value();
        int32_t end = std::get<2>(r).value() + std::get<4>(r).size();
        Record rec{start, end};
        records.push_back(rec);
    }

    std::unique_ptr<IntervalNode> root(nullptr);
    construct_tree(root, records);

    // Test Tree
    check_medians(root, 0, 0);
}
