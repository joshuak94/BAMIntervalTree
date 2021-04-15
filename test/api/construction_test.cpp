#include <gtest/gtest.h>
#include <math.h>

#include <bamit/all.hpp>

// Recursive function to go through the tree and check the medians.
void check_tree(std::unique_ptr<bamit::IntervalNode> const & root, int level, int pos,
                   std::vector<std::tuple<uint32_t, uint32_t>> const & expected_values)
{
    EXPECT_EQ(std::make_tuple(root->get_start(), root->get_end()), expected_values[pow(2, level) - 1 + pos]);
    if(root->get_left_node())
    {
        check_tree(root->get_left_node(), level + 1, (pos * 2), expected_values);
    }
    if(root->get_right_node())
    {
        check_tree(root->get_right_node(), level + 1, (pos * 2) + 1, expected_values);
    }

}

TEST(tree_construct, simulated_chr1_small_golden)
{
    std::vector<std::vector<bamit::Record>> records{};
    std::filesystem::path input{DATADIR"simulated_chr1_small_golden.bam"};
    bamit::sam_file_input_type input_file{input};

    std::vector<std::unique_ptr<bamit::IntervalNode>> node_list = bamit::construct_tree(input_file);

    // Test Tree
    // Compare chr, start, and end for each node.
    std::vector<std::tuple<uint32_t, uint32_t>> expected_values{
        std::make_tuple(849, 1016), std::make_tuple(368, 561), std::make_tuple(1310, 1495),
        std::make_tuple(157, 357), std::make_tuple(575, 775), std::make_tuple(1031, 1232),
        std::make_tuple(1603, 1760), std::make_tuple(43, 231), std::make_tuple(262, 453),
        std::make_tuple(472, 665), std::make_tuple(702, 902), std::make_tuple(946, 1119),
        std::make_tuple(1156, 1354), std::make_tuple(1443, 1617), std::make_tuple(1725, 1920),
        std::make_tuple(18, 137), std::make_tuple(146, 248), std::make_tuple(0, 0),
        std::make_tuple(364, 465), std::make_tuple(469, 570), std::make_tuple(574, 675),
        std::make_tuple(677, 783), std::make_tuple(830, 942), std::make_tuple(0, 0),
        std::make_tuple(0, 0), std::make_tuple(1135, 1254), std::make_tuple(1258, 1388),
        std::make_tuple(1400, 1529), std::make_tuple(1534, 1653), std::make_tuple(1669, 1816),
        std::make_tuple(1842, 2014)};
    check_tree(node_list[0], 0, 0, expected_values);

    std::filesystem::remove(input.replace_extension("bam.bit"));
}

TEST(tree_construct, simulated_mult_chr_small_golden)
{
    std::vector<std::vector<bamit::Record>> records{};
    std::filesystem::path input{DATADIR"simulated_mult_chr_small_golden.bam"};
    bamit::sam_file_input_type input_file{input};

    std::vector<std::unique_ptr<bamit::IntervalNode>> node_list = bamit::construct_tree(input_file);


    std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> expected_values{
        std::make_tuple(1, 193, 381), std::make_tuple(0, 404, 601), std::make_tuple(2, 113, 288),
        std::make_tuple(0, 182, 384), std::make_tuple(1, 38, 220), std::make_tuple(1, 476, 677),
        std::make_tuple(2, 321, 522), std::make_tuple(0, 31, 205), std::make_tuple(0, 287, 489),
        std::make_tuple(0, 505, 655), std::make_tuple(1, 130, 289), std::make_tuple(1, 385, 539),
        std::make_tuple(2, 24, 203), std::make_tuple(2, 220, 408), std::make_tuple(2, 437, 614),
        std::make_tuple(0, 14, 128), std::make_tuple(0, 134, 259), std::make_tuple(0, 286, 387),
        std::make_tuple(0, 398, 499), std::make_tuple(0, 0, 0), std::make_tuple(1, 6, 113),
        std::make_tuple(0, 0, 0), std::make_tuple(0, 0, 0), std::make_tuple(1, 295, 434),
        std::make_tuple(1, 446, 572), std::make_tuple(2, 0, 101), std::make_tuple(2, 108, 209),
        std::make_tuple(0, 0, 0), std::make_tuple(0, 0, 0), std::make_tuple(2, 434, 537),
        std::make_tuple(2, 539, 667)};
    // Test Tree
    // check_tree(root, 0, 0, expected_values);
    for (auto const & node : node_list)
    {
        node->print(0);
    }

    std::filesystem::remove(input.replace_extension("bam.bit"));
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

    bamit::sam_file_input_type input_file{unsorted_sam_path};
    EXPECT_THROW(bamit::construct_tree(input_file), seqan3::format_error);

    std::filesystem::remove(unsorted_sam_path);
}
