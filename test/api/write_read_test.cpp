#include <gtest/gtest.h>
#include <math.h>

#include <cereal/archives/binary.hpp>

#include <bamit/all.hpp>

TEST(write_read_test, write_read_test)
{
    std::filesystem::path tmp = std::filesystem::temp_directory_path()/"intervaltree";
    std::streamoff result{-1};
    std::streamoff result_after_reading{-1};
    bamit::Position start{1, 100};
    bamit::Position end{1, 150};
    std::filesystem::path input{DATADIR"simulated_mult_chr_small_golden.bam"};
    bamit::sam_file_input_type sam_in{input};
    {
        // Initialize variables for writing.
        std::ofstream out{tmp, std::ios_base::binary | std::ios_base::out};
        std::vector<std::unique_ptr<bamit::IntervalNode>> node_list{};
        std::vector<std::vector<bamit::Record>> records{};
        cereal::BinaryOutputArchive ar(out);

        node_list = bamit::index(sam_in);
        result = bamit::get_overlap_records(sam_in, node_list, start, end);

        // Write tree to output.
        bamit::write(node_list, ar);
        out.close();
    }
    {
        // Initialize variables for reading.
        std::ifstream in{tmp, std::ios_base::binary | std::ios_base::in};
        std::vector<std::unique_ptr<bamit::IntervalNode>> node_list{};
        cereal::BinaryInputArchive ar(in);

        // Read tree from input.
        bamit::read(node_list, ar);
        result_after_reading = bamit::get_overlap_records(sam_in, node_list, start, end);
        in.close();
    }

    // Compare search reults from tree constructed directly and tree loaded from file.
    EXPECT_EQ(result, result_after_reading);
}
