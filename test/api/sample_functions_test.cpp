#include <gtest/gtest.h>

#include <bamit/sample_functions.hpp>

TEST(sample_functions_test, sample_read_depth_test)
{
    std::filesystem::path input{DATADIR"simulated_mult_chr_small_golden.bam"};
    seqan3::sam_file_input input_file{input};

    std::vector<std::unique_ptr<bamit::IntervalNode>> node_list = bamit::index(input_file);
    bamit::EstimationResult result{};
    EXPECT_NO_THROW(result = bamit::sample_read_depth(input_file, node_list, 10));
    EXPECT_NO_THROW(result = bamit::sample_read_depth(input_file, node_list, 9));
    EXPECT_THROW(result = bamit::sample_read_depth(input_file, node_list, 1), std::invalid_argument);
}
