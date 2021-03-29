#include <gtest/gtest.h>
#include <math.h>

#include "IntervalNode.hpp"
#include "Record.hpp"

// Recursive function to go through the tree and check the medians.
void check_medians(std::unique_ptr<bamit::IntervalNode> const & root, int level, int pos,
                   std::vector<int32_t> const & expected_medians)
{
    if (expected_medians[pow(2, level) - 1 + pos] > 0)
    {
        EXPECT_EQ(root->get_median(), expected_medians[pow(2, level) - 1 + pos]);
        if(root->get_left_node())
        {
            check_medians(root->get_left_node(), level + 1, (pos * 2), expected_medians);
        }
        if(root->get_right_node())
        {
            check_medians(root->get_right_node(), level + 1, (pos * 2) + 1, expected_medians);
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

    std::vector<bamit::Record> records{};
    std::vector<uint32_t> cumulative_length{};
    std::filesystem::path input{DATADIR"simulated_chr1_small_golden.bam"};
    bamit::parse_file(input, cumulative_length, records);

    std::unique_ptr<bamit::IntervalNode> root(nullptr);
    bamit::construct_tree(root, records);

    // Test Tree
    std::vector<int32_t> expected_medians{944, 467, 1397, 257, 676, 1132, 1664, 140, 363, 571, 802, 1032, 1257,
                                          1531, 1822, 77, 197, 0, 414, 519, 624, 730, 886, 0, 0, 1194, 1323,
                                          1464, 1593, 1742, 1928};
    check_medians(root, 0, 0, expected_medians);
}

TEST(tree_construct, simulated_mult_chr_small_golden)
{
    // Construct Tree
    using my_fields = seqan3::fields<seqan3::field::id,
                                     seqan3::field::ref_id,
                                     seqan3::field::ref_offset,
                                     seqan3::field::cigar,
                                     seqan3::field::seq>;

    std::vector<bamit::Record> records{};
    std::vector<uint32_t> cumulative_length{};
    std::filesystem::path input{DATADIR"simulated_mult_chr_small_golden.bam"};
    bamit::parse_file(input, cumulative_length, records);

    std::unique_ptr<bamit::IntervalNode> root(nullptr);
    bamit::construct_tree(root, records);


    std::vector<int32_t> expected_medians{991, 502, 1611, 283, 824, 1276, 1821, 130, 388, 623, 909, 1142, 1505,
                                          1714, 1938, 71, 196, 336, 448, 563, 759, 0, 0, 1064, 1209, 1450, 1558,
                                          0, 0, 1885, 2003};
    // Test Tree
    check_medians(root, 0, 0, expected_medians);
}

TEST(tree_construct, unsorted)
{
    std::filesystem::path const tmp_dir = std::filesystem::temp_directory_path();     // get the temp directory
    // Create a blank SAM file without a sorting indicator.
    std::filesystem::path const unsorted_sam_path{tmp_dir/"unsorted.sam"};
    std::ofstream unsorted_sam{unsorted_sam_path.c_str()};
    unsorted_sam << "@HD\tVN:1.6\n" <<
                    "@SQ\tSN:testchr\tLN:1000\n" <<
                    "test1\t16\ttestchr\t1\t60\t10M\t=\t1\t0\tGCGCGCGCGC\tFFFFFFFFFF\n";
    unsorted_sam.close();

    std::vector<bamit::Record> records{};
    std::vector<uint32_t> cumulative_length{};
    EXPECT_THROW(bamit::parse_file(unsorted_sam_path, cumulative_length, records), seqan3::format_error);

    std::filesystem::remove(unsorted_sam_path);
}
