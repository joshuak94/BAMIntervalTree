#include <gtest/gtest.h>
#include <math.h>

#include <bamit/all.hpp>

// Recursive function to go through the tree and check the medians.
void check_medians(std::unique_ptr<bamit::IntervalNode> const & root, int level, int pos,
                   std::vector<bamit::Position> const & expected_medians)
{
    if (std::get<0>(expected_medians[pow(2, level) - 1 + pos]) != -1)
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
    std::filesystem::path input{DATADIR"simulated_chr1_small_golden.bam"};
    bamit::parse_file(input, records);

    std::unique_ptr<bamit::IntervalNode> root(nullptr);
    bamit::construct_tree(root, records);

    // Test Tree
    std::vector<bamit::Position> expected_medians{std::make_pair(0, 944), std::make_pair(0, 467), std::make_pair(0, 1397),
                                                 std::make_pair(0, 257), std::make_pair(0, 676), std::make_pair(0, 1132),
                                                 std::make_pair(0, 1664), std::make_pair(0, 140), std::make_pair(0, 363),
                                                 std::make_pair(0, 571), std::make_pair(0, 802), std::make_pair(0, 1032),
                                                 std::make_pair(0, 1257), std::make_pair(0, 1531), std::make_pair(0, 1822),
                                                 std::make_pair(0, 77), std::make_pair(0, 197), std::make_pair(-1, 0),
                                                 std::make_pair(0, 414), std::make_pair(0, 519), std::make_pair(0, 624),
                                                 std::make_pair(0, 730), std::make_pair(0, 886), std::make_pair(-1, 0),
                                                 std::make_pair(-1, 0), std::make_pair(0, 1194), std::make_pair(0, 1323),
                                                 std::make_pair(0, 1464), std::make_pair(0, 1593), std::make_pair(0, 1742),
                                                 std::make_pair(0, 1928)};
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
    std::filesystem::path input{DATADIR"simulated_mult_chr_small_golden.bam"};
    bamit::parse_file(input, records);

    std::unique_ptr<bamit::IntervalNode> root(nullptr);
    bamit::construct_tree(root, records);


    std::vector<bamit::Position> expected_medians{std::make_pair(1, 291), std::make_pair(0, 502), std::make_pair(2, 211),
                                                  std::make_pair(0, 283), std::make_pair(1, 124), std::make_pair(1, 576),
                                                  std::make_pair(2, 421), std::make_pair(0, 130), std::make_pair(0, 388),
                                                  std::make_pair(0, 623), std::make_pair(1, 209), std::make_pair(1, 442),
                                                  std::make_pair(2, 105), std::make_pair(2, 314), std::make_pair(2, 538),
                                                  std::make_pair(0, 71), std::make_pair(0, 196), std::make_pair(0, 336),
                                                  std::make_pair(0, 448), std::make_pair(0, 563), std::make_pair(1, 59),
                                                  std::make_pair(-1, 0), std::make_pair(-1, 0), std::make_pair(1, 364),
                                                  std::make_pair(1, 509), std::make_pair(2, 50), std::make_pair(2, 158),
                                                  std::make_pair(-1, 0), std::make_pair(-1, 0), std::make_pair(2, 485),
                                                  std::make_pair(2, 603)};
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
    EXPECT_THROW(bamit::parse_file(unsorted_sam_path, records), seqan3::format_error);

    std::filesystem::remove(unsorted_sam_path);
}
