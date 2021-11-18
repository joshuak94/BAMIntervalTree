#pragma once

#include <algorithm> // For std::sort
#include <random> // For std::random_device
#include <numeric> // For std::reduce
#include <cmath> // For std::sqrt and std::pow

#include <bamit/Record.hpp>
#include <bamit/IntervalNode.hpp>

#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/core/debug_stream.hpp>

namespace bamit
{

/*! Contains statistics which can be calculated over a sample. */
struct EstimationResult
{
    double mean{0}, median{0}, mode{0}, sd{0}, variance{0};

    /*! Print out the statistics stored within. */
    void print()
    {
        seqan3::debug_stream << "Mean: " << mean
                             << "\nMedian: " << median
                             << "\nMode: " << mode
                             << "\nSD: " << sd
                             << "\nVariance: " << variance << "\n";
    }
};

/*!
   \brief A function to sample read depth from an alignment file, given the input file, index, and number of positions.
   \param input_file The input alignment file in sam/bam format.
   \param bamit_index The vector of indices for each chromosome.
   \param sample_value The number of positions to sample.
   \param seed The seed to use for the random generator. Enables reproducibility. Default is 0.

   \details This function samples an alignment file at `sample_value` different random positions according to a uniform
            distribution and gives the read depth of the file as a series of statistics with mean, median, mode, sd and
            variance stored in an EstimationResult struct.

   \return Returns a struct containing statistics over the sampled points.
 */
template <typename traits_type, typename fields_type, typename format_type>
inline EstimationResult sample_read_depth(seqan3::sam_file_input<traits_type, fields_type, format_type> & input_file,
                                          std::vector<std::unique_ptr<IntervalNode>> const & bamit_index,
                                          uint64_t const & sample_value,
                                          uint64_t const & seed = 0)
    {
        if (sample_value <= 1) throw std::invalid_argument("sample_value must be greater than 1.");
        auto & header = input_file.header(); // Obtain the header from the input file.
        // Intialize vector of read depths at # of positions defined by sample_value.
        std::vector<uint64_t> read_depths(sample_value);
        std::map<double, uint64_t> depth_counts{}; // Counts how often each depth occurs.
        std::mt19937 gen(seed); // seed the generator
        std::uniform_int_distribution<> distr_chr(0, header.ref_ids().size() - 1); // define the range

        uint64_t rand_pos, rand_chr;
        std::tuple<uint64_t, uint64_t> pos_tuple;
        std::streamoff position{-1};
        auto it = input_file.begin();

        // For each sample, increment its associated read_depths[i] for each mapped read overlapping the position.
        for (auto & i : read_depths)
        {
            // Obtain a random chromosome from what is listed in the header, then obtain a random
            // position constrained by the chromosome size.
            rand_chr = distr_chr(gen);
            std::uniform_int_distribution<> distr_pos(0, std::get<0>(header.ref_id_info[rand_chr]) - 1);
            rand_pos = distr_pos(gen);
            pos_tuple = std::make_tuple(rand_chr, rand_pos);

            // Get file position of the above location.
            position = get_overlap_records(input_file, bamit_index, pos_tuple, pos_tuple);

            if (position == -1) continue; // The position is not covered by any reads.
            it.seek_to(position);
            // Count number of reads which overlap the position.
            for (; it != input_file.end(); ++it)
            {
                if (unmapped(*it)) continue;
                if (std::make_tuple((*it).reference_id().value(), (*it).reference_position().value()) > pos_tuple) break;
                ++i;
            }
        }
        std::sort(read_depths.begin(), read_depths.end());
        std::for_each(read_depths.begin(), read_depths.end(),
                      [&depth_counts](uint64_t const & value) {
                          ++(depth_counts[value]);
                      });

        EstimationResult result;
        seqan3::debug_stream << "read_depths: " << read_depths << "\ndepth_counts: " << depth_counts << '\n';
        result.mean = std::accumulate(read_depths.begin(), read_depths.end(), 0) / (double) sample_value;
        // If sample_value is even, median is average of two middle numbers. Otherwise, it is just the one number.
        result.median = sample_value % 2 == 0 ? (read_depths[sample_value / 2 - 1] + read_depths[sample_value / 2]) / (double) 2
                                              : read_depths[std::floor(sample_value / 2)];
        result.mode = (std::max_element(depth_counts.begin(), depth_counts.end(),
                                        [](std::pair<double, uint64_t> const & a,
                                           std::pair<double, uint64_t> const & b) {
                                               return a.second < b.second;
                                        }))->first;
        result.variance = std::accumulate(read_depths.begin(), read_depths.end(), 0.0,
                                          [&result](const double & a, const double & b) {
                                              return a + (std::pow((b - result.mean), 2));
                                          }) / (sample_value - 1);
        result.sd = std::sqrt(result.variance);

        return result;
    }
} // namespace bamit
