#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/all.hpp>

#include <cereal/archives/binary.hpp>

#include <bamit/all.hpp>

struct IndexOptions
{
    std::filesystem::path input_path{};
};

struct OverlapOptions
{
    std::filesystem::path input_path{};
    std::filesystem::path out_file{};
    std::string start{};
    std::string end{};
};

void initialize_top_parser(seqan3::argument_parser & parser)
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
}

void initialize_index_parser(seqan3::argument_parser & parser, IndexOptions & options)
{
    parser.add_option(options.input_path, 'i', "input_bam",
                      "Input a sorted BAM/SAM file to construct an index over.", seqan3::option_spec::required,
                      seqan3::input_file_validator{{"sam", "bam"}});
}

void initialize_overlap_parser(seqan3::argument_parser & parser, OverlapOptions & options)
{
    // Ref ID regex taken from the SAM format manual.
    std::string ref_id_match{"[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*"};
    seqan3::regex_validator query_validator{ref_id_match + ",{1}[0-9]+"};

    parser.add_option(options.input_path, 'i', "input_bam",
                      "The name of the SAM/BAM file to query.", seqan3::option_spec::required,
                      seqan3::input_file_validator{{"sam", "bam"}});
    parser.add_option(options.out_file, 'o', "output_sam",
                      "The SAM file, where the results should be stored.", seqan3::option_spec::standard,
                      seqan3::input_file_validator{{"sam"}});
    parser.add_option(options.start, 's', "start",
                      "The start of the interval to query, in the format chrA,posA."
                      " Note that when start and end are the same, this queries for reads overlapping a point.",
                      seqan3::option_spec::required,
                      query_validator);
    parser.add_option(options.end, 'e', "end",
                      "The end of the interval to query, in the format chrB,posB."
                      " Note that when start and end are the same, this queries for reads overlapping a point.",
                      seqan3::option_spec::required,
                      query_validator);
}

/*!
   \brief Parse the start and end strings from the given options and put them into the given bamit::Position containers.
   \param start Where the resulting start bamit::Position will go.
   \param end Where the resulting end bamit::Position will go.
   \param options The CmdOptions object storing the start and end strings from the user.
   \param ref_ids The reference chromosome names, stored in a deque by seqan3.
   \return 0 if the parsing was successful, -1 otherwise.
*/
int const parse_overlap_query(bamit::Position & start,
                              bamit::Position & end,
                              OverlapOptions const & options,
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

int run_index(std::unique_ptr<bamit::IntervalNode> & root, std::filesystem::path const & input)
{
    std::vector<std::vector<bamit::Record>> records{};
    bamit::parse_file(input, records);

    seqan3::debug_stream << "Creating Interval Tree.\n";
    bamit::construct_tree(root, records);
    seqan3::debug_stream << "Writing to file.\n";
    {
        std::filesystem::path index_path{input};
        index_path.replace_extension("bit");
        std::ofstream out_file(index_path,
                               std::ios_base::binary | std::ios_base::out);
        cereal::BinaryOutputArchive archive(out_file);
        bamit::write(root, archive);
        out_file.close();
    }
    return 0;
}

int parse_index(seqan3::argument_parser & parser)
{
    IndexOptions options{};

    initialize_index_parser(parser, options);

    // Parse the given arguments and catch possible errors.
    try
    {
        parser.parse();                                               // trigger command line parsing
    }
    catch (seqan3::argument_parser_error const & ext)                   // catch user errors
    {
        seqan3::debug_stream << "[Error] " << ext.what() << '\n';       // customise your error message
        return -1;
    }
    std::unique_ptr<bamit::IntervalNode> root{nullptr};
    run_index(root, options.input_path);

    return 0;
}

int parse_overlap(seqan3::argument_parser & parser)
{
    OverlapOptions options{};

    initialize_overlap_parser(parser, options);

    // Parse the given arguments and catch possible errors.
    try
    {
      parser.parse();                                               // trigger command line parsing
    }
    catch (seqan3::argument_parser_error const & ext)                   // catch user errors
    {
      seqan3::debug_stream << "[Error] " << ext.what() << '\n';       // customise your error message
      return -1;
    }

    bamit::sam_file_input_type input{options.input_path};
    std::filesystem::path index_path{options.input_path};
    std::unique_ptr<bamit::IntervalNode> root(nullptr);
    index_path.replace_extension("bit");
    if (!std::filesystem::exists(index_path))
    {
        run_index(root, options.input_path);
    }
    else
    {
        seqan3::debug_stream << "Reading index file...\n";
        root = std::make_unique<bamit::IntervalNode>();
        {
            std::ifstream in_file{index_path, std::ios_base::binary | std::ios_base::in};
            cereal::BinaryInputArchive iarchive(in_file);
            bamit::read(root, iarchive);
            in_file.close();
        }
    }
    seqan3::debug_stream << "Searching...\n";
    bamit::Position start, end;
    if (parse_overlap_query(start, end, options, input.header().ref_ids()) == -1)
    {
        return -1;
    }
    seqan3::debug_stream << "Search: " << start << " " << end << "\n";
    bamit::get_overlap_records(input, root, start, end, options.out_file);

    return 0;
}

int main(int argc, char ** argv)
{
    seqan3::argument_parser top_level_parser{"BAMIntervalTree", argc, argv,
                                             seqan3::update_notifications::on,
                                             {"index", "overlap"}};

    initialize_top_parser(top_level_parser);

    // Parse the given arguments and catch possible errors.
    try
    {
        top_level_parser.parse();                                               // trigger command line parsing
    }
    catch (seqan3::argument_parser_error const & ext)                   // catch user errors
    {
        seqan3::debug_stream << "[Error] " << ext.what() << '\n';       // customise your error message
        return -1;
    }

    seqan3::argument_parser & sub_parser = top_level_parser.get_sub_parser();
    if (sub_parser.info.app_name == std::string_view{"BAMIntervalTree-index"})
        return parse_index(sub_parser);
    else if (sub_parser.info.app_name == std::string_view{"BAMIntervalTree-overlap"})
        return parse_overlap(sub_parser);
    else
        seqan3::debug_stream << "Unhandled subparser named " << sub_parser.info.app_name << '\n';

    return 0;
}
