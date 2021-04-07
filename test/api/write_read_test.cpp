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
        std::unique_ptr<bamit::IntervalNode> root(nullptr);
        std::vector<bamit::Record> records{};
        cereal::BinaryOutputArchive ar(out);

        bamit::parse_file(input, records);
        bamit::construct_tree(root, records);
        bamit::overlap(sam_in, root, start, end, result);

        // Write tree to output.
        bamit::write(root, ar);
        out.close();
    }
    {
        // Initialize variables for reading.
        std::ifstream in{tmp, std::ios_base::binary | std::ios_base::in};
        std::unique_ptr<bamit::IntervalNode> in_node = std::make_unique<bamit::IntervalNode>();
        cereal::BinaryInputArchive ar(in);

        // Read tree from input.
        bamit::read(in_node, ar);
        bamit::overlap(sam_in, in_node, start, end, result_after_reading);
        in.close();
    }

    // Compare search reults from tree constructed directly and tree loaded from file.
    EXPECT_EQ(result, result_after_reading);
}
