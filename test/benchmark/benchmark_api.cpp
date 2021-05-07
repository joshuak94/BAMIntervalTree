#include <gtest/gtest.h>
#include <bamit/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <random>
#include <chrono>
#include <ctime>

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
    std::uniform_int_distribution<> distr_end(rand_chr_start, header.ref_ids().size() - 1); // define the range
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
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path();     // get the temp directory
    std::filesystem::path result_sam_path{tmp_dir/"result.sam"};
    std::filesystem::path large_file{DATADIR"large_file.bam"};
    if (!std::filesystem::exists(large_file))
    {
        return;
    }
    bamit::sam_file_input_type input_bam{large_file};

    // Construct tree.
    std::vector<std::unique_ptr<bamit::IntervalNode>> node_list{};
    RUN((node_list = bamit::index(input_bam)), "Construction");

    // Generate 100 overlaps.
    bamit::Position start, end;
    std::streamoff result{-1};
    bamit::sam_file_input_type input_bam_2{large_file};
    for (int i = 0; i < 100; i++)
    {
        get_random_position(start, end, input_bam.header());
        std::string query{"[" + std::to_string(std::get<0>(start)) + ", " +
                                std::to_string(std::get<1>(start)) + "] - [" +
                                std::to_string(std::get<0>(end)) + ", " +
                                std::to_string(std::get<1>(end)) + "]"};
        RUN(bamit::get_overlap_records(input_bam_2, node_list, start, end, result_sam_path),
            query + " writing");
        RUN(bamit::get_overlap_records(input_bam_2, node_list, start, end),
            query + " offest");
        std::filesystem::remove(result_sam_path);
    }
}
