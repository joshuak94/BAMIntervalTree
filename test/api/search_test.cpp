#include <gtest/gtest.h>
#include <math.h>

#include <seqan3/test/expect_range_eq.hpp>

#include <bamit/all.hpp>

TEST(get_overlap_records, simulated_chr1_small_golden)
{
    // Construct Tree
    std::filesystem::path input{DATADIR"simulated_mult_chr_small_golden.bam"};
    bamit::sam_file_input_type input_file{input};
    std::vector<std::unique_ptr<bamit::IntervalNode>> node_list = bamit::index(input_file);

    // Search Tree
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path();     // get the temp directory
    std::filesystem::path result_sam_path{tmp_dir/"result.sam"};
    bamit::Position start{1, 100};
    bamit::Position end{1, 110};
    // Potentially undefined behavior in our seqan branch: input_file is parsed in the tree construction,
    // so using the same input_file variable here leads to incorrect results (probably because the iterator no longer
    // points to the beginning of the file).
    // Temporary fix: use new object.
    // Permanent fix: seqan3 is changing how you get/use the file position.
    bamit::sam_file_input_type input_file_2{input};
    get_overlap_records(input_file_2, node_list, start, end, result_sam_path);
    seqan3::sam_file_input expected{DATADIR"samtools_result.sam",  bamit::sam_file_output_fields{}};
    seqan3::sam_file_input result{result_sam_path,  bamit::sam_file_output_fields{}};

    EXPECT_RANGE_EQ(expected, result);

     std::filesystem::remove(result_sam_path);
     std::filesystem::remove(input.replace_extension("bam.bit"));
}
