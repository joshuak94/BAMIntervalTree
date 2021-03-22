#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "IntervalNode.hpp"

int main()
{
    using my_fields = seqan3::fields<seqan3::field::id,
                                     seqan3::field::ref_id,
                                     seqan3::field::ref_offset,
                                     seqan3::field::cigar>;

    seqan3::debug_stream << "Reading SAM file.\n";
    std::filesystem::path input_path{"/home/kim_j/testFiles/simulated_reads/simulated_chr1_small_golden.bam"};
    seqan3::sam_file_input input{input_path, my_fields{}};

    seqan3::debug_stream << "Putting reads into vector.\n";
    std::vector<std::string> records{};
    for (auto & r : input)
    {
        std::string id = std::get<0>(r);
        records.push_back(id);
    }

    seqan3::debug_stream << "Creating Node.\n";
    IntervalNode root{NULL, NULL, records};

    std::vector<std::string> result = root.get_records();

    seqan3::debug_stream << result << std::endl;

    root.construct_tree(input);
    return 0;
}
