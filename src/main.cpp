#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "IntervalNode.hpp"

int main()
{
    using my_fields = seqan3::fields<seqan3::field::id,
                                     seqan3::field::ref_id,
                                     seqan3::field::ref_offset,
                                     seqan3::field::cigar,
                                     seqan3::field::seq>;

    seqan3::debug_stream << "Reading SAM file.\n";
    std::filesystem::path input_path{"/home/kim_j/development/BAMIntervalTree/test/data/simulated_chr1_small_golden.bam"};
    seqan3::sam_file_input input{input_path, my_fields{}};

    seqan3::debug_stream << "Extracting info from reads\n";
    // Calculate total length of genome.
    // size_t length{0};
    // for (auto c : input.header().ref_id_info)
    // {
    //     length += std::get<0>(c);
    // }
    std::vector<Record> records{};
    for (auto & r : input)
    {
        int32_t start = std::get<2>(r).value();
        int32_t end = std::get<2>(r).value() + std::get<4>(r).size();
        Record rec{start, end};
        records.push_back(rec);
    }

    seqan3::debug_stream << "Creating Node.\n";
    IntervalNode root{NULL, NULL};


    root.construct_tree(records);

    // root.print(0);
    // std::vector<Record> result = root.get_records();
    //
    // for (auto r : result)
    // {
    //     seqan3::debug_stream << "[" << r.start << ", " << r.end << "]" << std::endl;
    // }
    return 0;
}
