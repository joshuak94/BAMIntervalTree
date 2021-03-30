#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/all.hpp>

#include <bamit/all.hpp>

struct CmdOptions
{
    std::filesystem::path input_path{};
    std::string start{};
    std::string end{};
};

/*!
   \brief Parse the start and end strings from the given options and put them into the given bamit::Position containers.
   \param start Where the resulting start bamit::Position will go.
   \param end Where the resulting end bamit::Position will go.
   \param options The CmdOptions object storing the start and end strings from the user.
   \param ref_ids The reference chromosome names, stored in a deque by seqan3.
   \return 0 if the parsing was successful, -1 otherwise.
*/
int32_t const parse_overlap_query(bamit::Position & start,
                                  bamit::Position & end,
                                  CmdOptions const & options,
                                  std::deque<std::string> const & ref_ids)
{
    // Get the location of the colon in the start and end.
    size_t const start_split = options.start.find(',');
    size_t const end_split = options.end.find(',');
    if (start_split == options.start.size() || end_split == options.start.size())
    {
        seqan3::debug_stream << "[ERROR] Start and end positions must be in the format chr_name:position!\n";
        return -1;
    }

    // Get the iterator position of the chromosomes given.
    auto start_ref_id_it = std::find(ref_ids.begin(), ref_ids.end(), options.start.substr(0, start_split));
    auto end_ref_id_it = std::find(ref_ids.begin(), ref_ids.end(), options.end.substr(0, end_split));
    if (start_ref_id_it == ref_ids.end() || end_ref_id_it == ref_ids.end())
    {
        seqan3::debug_stream << "[ERROR] One or both of the chromosome names supplied could not be found.\n";
        return -1;
    }

    // Convert the start and end strings to integers.
    int32_t start_position{};
    int32_t end_position{};
    try
    {
        start_position = std::stoi(options.start.substr(start_split + 1));
        end_position = std::stoi(options.end.substr(end_split + 1));
    }
    catch (...)
    {
        seqan3::debug_stream << "[ERROR] There was a formatting error with the overlap query!\n";
        return -1;
    }

    // Convert query string into two Positions.
    int32_t start_ref_id = start_ref_id_it - ref_ids.begin();
    int32_t end_ref_id = end_ref_id_it - ref_ids.begin();
    start = std::make_pair(start_ref_id, start_position);
    end = std::make_pair(end_ref_id, end_position);
    return 0;
}

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

    // Ref ID regex taken from the SAM format manual.
    std::string ref_id_match{"[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*"};
    seqan3::regex_validator query_validator{ref_id_match + ",{1}[0-9]+"};

    parser.add_option(options.input_path, 'i', "input_bam",
                      "Input a sorted BAM/SAM file.", seqan3::option_spec::standard,
                      seqan3::input_file_validator{{"sam", "bam"}});
    parser.add_option(options.start, 's', "start",
                      "The start of the interval to query, in the format chrA,posA."
                      " Note that when start and end are the same, this queries for reads overlapping a point.",
                      seqan3::option_spec::standard,
                      query_validator);
    parser.add_option(options.end, 'e', "end",
                      "The end of the interval to query, in the format chrB,posB."
                      " Note that when start and end are the same, this queries for reads overlapping a point.",
                      seqan3::option_spec::standard,
                      query_validator);
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

    seqan3::debug_stream << "Reading SAM file.\n";

    seqan3::debug_stream << "Extracting info from reads\n";

    std::vector<bamit::Record> records{};
    bamit::parse_file(options.input_path, records);

    seqan3::debug_stream << "Creating Interval Tree.\n";
    std::unique_ptr<bamit::IntervalNode> root(nullptr);
    bamit::construct_tree(root, records);

    root->print(0);

    if (!options.start.empty() && !options.end.empty())
    {
        using my_fields = seqan3::fields<seqan3::field::id,
                                         seqan3::field::ref_id,
                                         seqan3::field::ref_offset,
                                         seqan3::field::cigar,
                                         seqan3::field::seq>;
        seqan3::sam_file_input input{options.input_path, my_fields{}};
        bamit::Position start, end;
        if (parse_overlap_query(start, end, options, input.header().ref_ids()) == -1)
        {
            return -1;
        }
        seqan3::debug_stream << "Search: " << start << " " << end << "\n";
        std::vector<bamit::Record> results{};
        bamit::overlap(root, start, end, results);
        for (auto & r: results)
            seqan3::debug_stream << r.start << " " << r.end << "\n";
    }
    return 0;
}
