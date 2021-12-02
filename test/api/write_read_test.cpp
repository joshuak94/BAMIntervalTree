#include <gtest/gtest.h>
#include <math.h>

#include <cereal/archives/binary.hpp>

#include <bamit/all.hpp>

TEST(write_read_test, write_read_test)
{
    std::filesystem::path tmp = std::filesystem::temp_directory_path()/"intervaltree";
    bamit::Position start{1, 100};
    bamit::Position end{1, 150};
    std::filesystem::path input{DATADIR"simulated_mult_chr_small_golden.bam"};

    // Initialize variables for writing.
    seqan3::sam_file_input sam_in{input};
    std::ofstream out{tmp, std::ios_base::binary | std::ios_base::out};
    std::vector<std::unique_ptr<bamit::IntervalNode>> node_list{};
    std::vector<std::vector<bamit::Record>> records{};
    cereal::BinaryOutputArchive ar_out(out);

    node_list = bamit::index(sam_in);
    auto result = bamit::get_overlap_records(sam_in, node_list, start, end);

    // Write tree to output.
    bamit::write(node_list, ar_out);
    out.close();

    // Initialize variables for reading.
    std::ifstream in{tmp, std::ios_base::binary | std::ios_base::in};
    cereal::BinaryInputArchive ar_in(in);

    // Read tree from input.
    node_list.clear();
    bamit::read(node_list, ar_in);
    auto result_after_reading = bamit::get_overlap_records(sam_in, node_list, start, end);
    in.close();

    // Compare search reults from tree constructed directly and tree loaded from file.
    EXPECT_EQ(result.size(), result_after_reading.size());

    for (size_t i = 0; i < result.size(); ++i)
    {
        EXPECT_EQ(result[i].reference_id().value(), result_after_reading[i].reference_id().value());
        EXPECT_EQ(result[i].reference_position().value(), result_after_reading[i].reference_position().value());
        EXPECT_EQ(result[i].id(), result_after_reading[i].id());
    }
}
