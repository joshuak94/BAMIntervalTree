#include <gtest/gtest.h>
#include <math.h>

#include <seqan3/test/expect_range_eq.hpp>

#include <bamit/all.hpp>

TEST(get_records, simulated_chr1_small_golden)
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

    // Search Tree


    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path();     // get the temp directory
    std::filesystem::path result_sam_path{tmp_dir/"result.sam"};
    bamit::sam_file_input_type input_bam{input};
    bamit::Position start = std::make_pair(1, 100);
    bamit::Position end{1, 110};
    {get_records(input_bam, root, start, end, result_sam_path);}
    seqan3::sam_file_input expected{DATADIR"samtools_result.sam",  bamit::sam_file_output_fields{}};
    seqan3::sam_file_input result{result_sam_path,  bamit::sam_file_output_fields{}};

    EXPECT_RANGE_EQ(expected, result);


    std::filesystem::remove(result_sam_path);
}