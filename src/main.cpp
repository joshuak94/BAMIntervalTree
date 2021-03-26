#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "IntervalNode.hpp"

struct CmdOptions
{
    std::filesystem::path input_path{};
    int32_t start{};
    int32_t end{};
};

void initialize_argument_parser(seqan3::argument_parser & parser, CmdOptions & options)
{
    parser.info.author = "Joshua Kim, Mitra Darvish";
    parser.info.app_name = "BAMIntervalTree";
    parser.info.man_page_title = "An Interval Tree indexer for BAM/SAM files.";
    parser.info.short_description = "Create an Interval Tree over an aligment file for quick range queries.";
    parser.info.version = "0.0.1";
    parser.info.date = "24-03-2021";    // last update
    parser.info.email = "kim_j@molgen.mpg.de";
    parser.info.long_copyright = "long_copyright";
    parser.info.short_copyright = "short_copyright";
    parser.info.url = "https://github.com/joshuak94/BAMIntervalTree/";

    parser.add_option(options.input_path, 'i', "input_bam",
                      "Input a sorted BAM/SAM file.", seqan3::option_spec::standard,
                      seqan3::input_file_validator{{"sam", "bam"}});
    parser.add_option(options.start, 's', "start", "Start of queried interval");
    parser.add_option(options.end, 'e', "end", "End of queried interval");
}

int main(int argc, char ** argv)
{
    seqan3::argument_parser my_parser{"BAMIntervalTree", argc, argv};
    CmdOptions options{};

    initialize_argument_parser(my_parser, options);

    // Parse the given arguments and catch possible errors.
    try
    {
        my_parser.parse();                                               // trigger command line parsing
    }
    catch (seqan3::argument_parser_error const & ext)                   // catch user errors
    {
        seqan3::debug_stream << "[Error] " << ext.what() << '\n';       // customise your error message
        return -1;
    }

    using my_fields = seqan3::fields<seqan3::field::id,
                                     seqan3::field::ref_id,
                                     seqan3::field::ref_offset,
                                     seqan3::field::cigar,
                                     seqan3::field::seq>;

    seqan3::debug_stream << "Reading SAM file.\n";
    seqan3::sam_file_input input{options.input_path, my_fields{}};

    seqan3::debug_stream << "Extracting info from reads\n";

    std::vector<bamit::Record> records{};
    std::vector<uint32_t> cumulative_length{};
    bamit::parse_file(options.input_path, cumulative_length, records);

    seqan3::debug_stream << "Creating Node.\n";
    std::unique_ptr<bamit::IntervalNode> root(nullptr);

    bamit::construct_tree(root, records);

    root->print(0);

    if (options.start && options.end)
    {
        seqan3::debug_stream << "Search: " << options.start << " " << options.end << "\n";
        std::vector<bamit::Record> results{};
        bamit::overlap(root, options.start, options.end, results);
        for (auto & r: results)
            seqan3::debug_stream << r.start << " " << r.end << "\n";
    }
    // std::vector<Record> result = root.get_records();
    //
    // for (auto r : result)
    // {
    //     seqan3::debug_stream << "[" << r.start << ", " << r.end << "]" << std::endl;
    // }
    return 0;
}
