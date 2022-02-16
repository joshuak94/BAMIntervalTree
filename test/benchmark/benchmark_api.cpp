#include <gtest/gtest.h>
#include <bamit/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <random>
#include <chrono>
#include <ctime>
#include <cstring>

#include "test_htslib.hpp"

#define RUN(x, y) {startTimeMessage(y);x;endTimeMessage(y);}

time_t _t1, _t2;
std::chrono::time_point<std::chrono::system_clock> _m1, _m2;

void printTimeMessage(std::string msg)
{
    time_t t = time(0);
    struct tm* now = localtime(&t);
    seqan3::debug_stream << "[ " << std::put_time(now, "%d-%m-%Y %H:%M:%S") << "] " << msg << '\n';
}

void startTimeMessage(std::string msg)
{
    time(&_t1);
    _m1 = std::chrono::time_point_cast<std::chrono::microseconds>(std::chrono::system_clock::now());
    printTimeMessage("[START] " + msg);
}

void endTimeMessage(std::string msg)
{
    using std::chrono::operator""ms;
    using std::chrono::operator""us;
    time(&_t2);
    _m2 = std::chrono::time_point_cast<std::chrono::microseconds>(std::chrono::system_clock::now());
    auto difference = std::chrono::duration_cast<std::chrono::microseconds>(_m2 - _m1);
    if (difference < 1000us)
    {
        printTimeMessage("[END] " + msg + " (" + std::to_string(difference.count()) + " microseconds.)");
    }
    else if (difference >= 1000us && difference < 1000ms)
    {
        printTimeMessage("[END] " + msg + " (" +
                         std::to_string((std::chrono::duration_cast<std::chrono::milliseconds>(difference)).count()) +
                         " milliseconds.)");
    }
    else
    {
        printTimeMessage("[END] " + msg + " (" + std::to_string((int)difftime(_t2,_t1)) + " seconds.)");
    }
}

void get_random_position(bamit::Position & start, bamit::Position & end,
                         seqan3::sam_file_header<std::deque<std::string>> & header)
{
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::uniform_int_distribution<> distr_start(0, header.ref_ids().size() - 1); // define the range
    uint32_t rand_chr_start = distr_start(gen); // Get start random chromosome.
    std::uniform_int_distribution<> distr_end(rand_chr_start, std::min((rand_chr_start+1), static_cast<uint32_t>(header.ref_ids().size() - 1))); // define the range
    uint32_t rand_chr_end = distr_end(gen); // Get end random chromosome, guaranteed to be >= start.

    std::uniform_int_distribution<> distr_pos_start(0, std::get<0>(header.ref_id_info[rand_chr_start]) - 1);
    uint32_t rand_pos_start = distr_pos_start(gen);

    // Get end position which is greater than or equal to start position.
    uint32_t end_start = (rand_chr_start == rand_chr_end) ? rand_pos_start : 0;
    std::uniform_int_distribution<> distr_pos_end(end_start, std::get<0>(header.ref_id_info[rand_chr_end]) - 1);
    uint32_t rand_pos_end = distr_pos_end(gen);

    start = std::make_tuple(rand_chr_start, rand_pos_start);
    end = std::make_tuple(rand_chr_end, rand_pos_end);
}

TEST(benchmark, construct_and_search)
{
    using std::chrono::operator""us;
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path();     // get the temp directory
    std::filesystem::path large_file{DATADIR"large_file.bam"};
    if (!std::filesystem::exists(large_file))
    {
        seqan3::debug_stream << "large_file.bam does not exist in the data directory.";
        return;
    }
    seqan3::contrib::bgzf_thread_count = 2;
    seqan3::sam_file_input<seqan3::sam_file_input_default_traits<>,
                           seqan3::fields<seqan3::field::ref_id,
                                          seqan3::field::ref_offset,
                                          seqan3::field::cigar,
                                          seqan3::field::flag>,
                           seqan3::type_list<seqan3::format_bam,
                                             seqan3::format_sam>> input_bam{large_file};

    // Construct tree.
    std::vector<std::unique_ptr<bamit::IntervalNode>> node_list{};
    RUN((node_list = bamit::index(input_bam)), "Construction");

    // Index via HTSlib
    RUN(htslib_index(large_file.c_str()), "HTSlib indexing");
    hts_idx_t * index;
    htsFile * in = hts_open(large_file.c_str(), "r");
    index = sam_index_load(in, large_file.c_str());

    // Load header via HTSlib
    sam_hdr_t * header = sam_hdr_read(in);

    // Collect average times.
    auto avg_bit_offset{0us}, avg_bit_overlap{0us}, avg_hts_offset{0us}, avg_hts_overlap{0us};
    // Generate 100 overlaps.
    bamit::Position start, end;
    std::streamoff result{-1};
    for (int i = 0; i < 100; i++)
    {
        seqan3::sam_file_input <seqan3::sam_file_input_default_traits<>,
                               seqan3::fields<seqan3::field::ref_id,
                                              seqan3::field::ref_offset,
                                              seqan3::field::cigar,
                                              seqan3::field::flag>,
                               seqan3::type_list<seqan3::format_bam,
                                                 seqan3::format_sam>> input_bam_write{large_file};
        seqan3::sam_file_input <seqan3::sam_file_input_default_traits<>,
                               seqan3::fields<seqan3::field::ref_id,
                                              seqan3::field::ref_offset,
                                              seqan3::field::cigar,
                                              seqan3::field::flag>,
                               seqan3::type_list<seqan3::format_bam,
                                                 seqan3::format_sam>> input_bam_offset{large_file};
        get_random_position(start, end, input_bam.header());
        std::string query{"[" + input_bam.header().ref_ids()[std::get<0>(start)] + ", " +
                                std::to_string(std::get<1>(start)) + "] - [" +
                                input_bam.header().ref_ids()[std::get<0>(end)] + ", " +
                                std::to_string(std::get<1>(end)) + "]\n"};
        std::filesystem::path result_sam_path = tmp_dir/(std::to_string(i) + "_bit.bam");
        std::filesystem::path result_htslib_path = tmp_dir/(std::to_string(i) + "_hts.bam");
        htsFile * out = hts_open(result_htslib_path.c_str(), "w");

        // Convert start and end to a string of regions for htslib.
        std::vector<std::string> region_list{};
        // char ** region = (char**) malloc( (1 + std::get<0>(end) - std::get<0>(start)) * sizeof(char*) );
        std::string start_str, end_str;
        start_str = input_bam.header().ref_ids()[std::get<0>(start)] + ":" +
                    std::to_string(std::get<1>(start));
        region_list.push_back(start_str);
        if (std::get<0>(start) != std::get<0>(end))
        {
            for (size_t j = (std::get<0>(start) + 1); j < std::get<0>(end); j++)
            {
                region_list.push_back(input_bam.header().ref_ids()[j]);
            }
            region_list.push_back(input_bam.header().ref_ids()[std::get<0>(end)] + ":" +
                                  "1-" + std::to_string(std::get<1>(end)));
        }
        else
        {
            start_str += "-";
            start_str += std::to_string(std::get<1>(end));
        }
        region_list[0] = start_str;

        // Copy the region_list from string to region char** array.
        std::vector<char*> region;
        for (const auto & r : region_list)
        {
            region.push_back(strdup(r.data()));
        }
        region.push_back(nullptr);

        seqan3::debug_stream << i << ": " << query;
        _m1 = std::chrono::time_point_cast<std::chrono::microseconds>(std::chrono::system_clock::now());
        bamit::get_overlap_records(input_bam_write, node_list, start, end, false, result_sam_path);
        _m2 = std::chrono::time_point_cast<std::chrono::microseconds>(std::chrono::system_clock::now());
        avg_bit_overlap += std::chrono::duration_cast<std::chrono::microseconds>(_m2 - _m1);

        _m1 = std::chrono::time_point_cast<std::chrono::microseconds>(std::chrono::system_clock::now());
        bamit::get_overlap_file_position(input_bam_offset, node_list, start, end, result);
        _m2 = std::chrono::time_point_cast<std::chrono::microseconds>(std::chrono::system_clock::now());
        avg_bit_offset += std::chrono::duration_cast<std::chrono::microseconds>(_m2 - _m1);

        _m1 = std::chrono::time_point_cast<std::chrono::microseconds>(std::chrono::system_clock::now());
        htslib_overlap_records(index, header, region.data(), (std::get<0>(end) - std::get<0>(start)) + 1, in, out);
        _m2 = std::chrono::time_point_cast<std::chrono::microseconds>(std::chrono::system_clock::now());
        avg_hts_overlap += std::chrono::duration_cast<std::chrono::microseconds>(_m2 - _m1);

        _m1 = std::chrono::time_point_cast<std::chrono::microseconds>(std::chrono::system_clock::now());
        htslib_overlap_file_offset(index, header, region.data(), (std::get<0>(end) - std::get<0>(start)) + 1);
        _m2 = std::chrono::time_point_cast<std::chrono::microseconds>(std::chrono::system_clock::now());
        avg_hts_offset += std::chrono::duration_cast<std::chrono::microseconds>(_m2 - _m1);
        hts_close(out);

        std::filesystem::remove(result_sam_path);
        std::filesystem::remove(result_htslib_path);
    }

    seqan3::debug_stream << "Average for get_overlap_records: " <<
                            std::to_string((avg_bit_overlap.count())/100) << "\n";
    seqan3::debug_stream << "Average for get_overlap_file_position: " <<
                            std::to_string((avg_bit_offset.count())/100) << "\n";
    seqan3::debug_stream << "Average for htslib_overlap_records: " <<
                            std::to_string((avg_hts_overlap.count())/100) << "\n";
    seqan3::debug_stream << "Average for htslib_overlap_file_offset: " <<
                            std::to_string((avg_hts_offset.count())/100) << "\n";
    sam_hdr_destroy(header);
    hts_idx_destroy(index);
    hts_close(in);
}
