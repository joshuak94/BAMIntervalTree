#include <gtest/gtest.h>
#include <math.h>
#include <filesystem>
#include <ios>

#include <cereal/archives/binary.hpp>

#include <bamit/all.hpp>

TEST(write_read_test, write_read_test)
{
    std::filesystem::path tmp = std::filesystem::temp_directory_path()/"intervaltree";
    std::vector<bamit::Record> results{};
    std::vector<bamit::Record> results_after_reading{};
    bamit::Position start{1, 100};
    bamit::Position end{1, 150};
    {
        // Initialize variables for writing.
        std::filesystem::path input{DATADIR"simulated_mult_chr_small_golden.bam"};
        std::ofstream out{tmp, std::ios_base::binary | std::ios_base::out};
        std::unique_ptr<bamit::IntervalNode> root(nullptr);
        std::vector<bamit::Record> records{};
        cereal::BinaryOutputArchive ar(out);

        // Construct Tree
        using my_fields = seqan3::fields<seqan3::field::id,
                                         seqan3::field::ref_id,
                                         seqan3::field::ref_offset,
                                         seqan3::field::cigar,
                                         seqan3::field::seq>;

        bamit::parse_file(input, records);
        bamit::construct_tree(root, records);
        bamit::overlap(root, start, end, results);

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
        bamit::overlap(in_node, start, end, results_after_reading);
        in.close();
    }

    // Compare search reults from tree constructed directly and tree loaded from file.
    EXPECT_EQ(results, results_after_reading);
}
