#pragma once

#include <seqan3/io/sam_file/input.hpp>

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

void parse_file(std::filesystem::path & input_path,
                std::vector<uint32_t> & cumulative_length,
                std::vector<Record> & record_list);

int32_t get_length(std::vector<seqan3::cigar> const & cigar);
