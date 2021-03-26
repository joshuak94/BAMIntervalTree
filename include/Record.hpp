#pragma once

#include <seqan3/io/sam_file/input.hpp>

namespace bamit
{
struct Record
{
    int32_t start, end, ref_id;

    Record(int32_t start_i, int32_t end_i, int32_t ref_id_i) : start{std::move(start_i)},
                                                               end{std::move(end_i)},
                                                               ref_id{std::move(ref_id_i)} {}
};

struct RecordComparatorStart
{
    bool operator ()(const Record & record1, const Record & record2)
    {
        if(record1.start == record2.start)
            return record1.end < record2.end;
        return record1.start < record2.start;
    }
};

struct RecordComparatorEnd
{
    bool operator ()(const Record & record1, const Record & record2)
    {
        if(record1.end == record2.end)
            return record1.start > record2.start;
        return record1.end > record2.end;
    }
};

int32_t get_length(std::vector<seqan3::cigar> const & cigar)
{
    using seqan3::operator""_cigar_operation;
    using seqan3::get;

    int32_t result{0};
    for (auto const & c : cigar)
    {
        seqan3::cigar::operation op{get<1>(c)};
        uint32_t length{get<0>(c)};
        if (op == 'M'_cigar_operation ||
            op == 'I'_cigar_operation ||
            op == 'D'_cigar_operation ||
            op == '='_cigar_operation ||
            op == 'X'_cigar_operation)
        {
            result += length;
        }
    }
    return result;
}

void parse_file(std::filesystem::path & input_path,
                std::vector<uint32_t> & cumulative_length,
                std::vector<Record> & record_list)
{
    using my_fields = seqan3::fields<seqan3::field::id,
                                     seqan3::field::ref_id,
                                     seqan3::field::ref_offset,
                                     seqan3::field::cigar,
                                     seqan3::field::seq>;
    seqan3::sam_file_input input{input_path, my_fields{}};

    // Calculate total length of genome.
    uint32_t running_sum{0};
    for (auto c : input.header().ref_id_info)
    {
        cumulative_length.push_back(running_sum);
        running_sum += std::get<0>(c);
    }

    for (auto & r : input)
    {
        int32_t ref_id = std::get<1>(r).value();
        int32_t length = get_length(std::get<3>(r));
        int32_t start = std::get<2>(r).value() + cumulative_length[ref_id];
        int32_t end = std::get<2>(r).value() + length + cumulative_length[ref_id];
        Record rec{start, end, ref_id};
        record_list.push_back(rec);
    }
}
}
