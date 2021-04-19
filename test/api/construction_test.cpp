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

    std::vector<std::unique_ptr<bamit::IntervalNode>> node_list = bamit::index(input_file);

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

    std::vector<std::unique_ptr<bamit::IntervalNode>> node_list = bamit::index(input_file);


    std::vector<std::vector<std::tuple<uint32_t, uint32_t>>> expected_values{
        {std::make_tuple(286, 469), std::make_tuple(134,322), std::make_tuple(440,625), std::make_tuple(14, 205),
         std::make_tuple(228, 384), std::make_tuple(388, 510), std::make_tuple(529, 655)},
        {std::make_tuple(174, 369), std::make_tuple(62, 251), std::make_tuple(408, 595), std::make_tuple(6, 146),
         std::make_tuple(152, 266), std::make_tuple(295, 434), std::make_tuple(497, 677), std::make_tuple(0, 0),
         std::make_tuple(0, 0), std::make_tuple(0, 0), std::make_tuple(0, 0), std::make_tuple(280, 381),
         std::make_tuple(385, 486), std::make_tuple(0, 0), std::make_tuple(0, 0)},
        {std::make_tuple(231, 433), std::make_tuple(54, 248), std::make_tuple(402, 603), std::make_tuple(0, 153),
         std::make_tuple(156, 328), std::make_tuple(334, 478), std::make_tuple(511, 667)}};
    // Test Tree
    for (size_t i = 0; i < node_list.size(); ++i)
    {
        check_tree(node_list[i], 0, 0, expected_values[i]);
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
    EXPECT_THROW(bamit::index(input_file), seqan3::format_error);

    std::filesystem::remove(unsorted_sam_path);
}
