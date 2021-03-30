#pragma once

#include <seqan3/io/sam_file/input.hpp>

namespace bamit
{
using Position = std::pair<int32_t, int32_t>;
/*! A Record object contains pertinent information about an alignment. */
struct Record
{
    Position start, end;

    Record(Position start_i, Position end_i) :
        start{std::move(start_i)},
        end{std::move(end_i)} {}
};

/*! Used to sort two Record objects in ascending order by start position. */
struct RecordComparatorStart
{
    /*!
       \brief Compare two Record objects and sort them in ascending order by starts.
       \param record1 The first record.
       \param record2 The second record.
       \return Returns `true` if record1 starts before record2, or if they both start at the same position but
               record1 ends before record2.
    */
    bool operator ()(Record const & record1, Record const & record2)
    {
        return (record1.start == record2.start) ? (record1.end < record2.end) : (record1.start < record2.start);
    }
};

/*! Used to sort two Record objects in descending order by end position. */
struct RecordComparatorEnd
{
    /*!
       \brief Compare two Record objects and sort them in descending order by ends.
       \param record1 The first record.
       \param record2 The second record.
       \return Returns `true` if record1 ends after record2, or if they both end at the same position but
               record1 starts after record2.
    */
    bool operator ()(Record const & record1, Record const & record2)
    {
        return (record1.end == record2.end) ? (record1.start > record2.start) : (record1.end > record2.end);
    }
};

/*!
   \brief Get the length of a seqan3::cigar vector based on M/I/D/=/X operations.
   \param cigar The vector of seqan3::cigar characters.
   \return Returns the length of M/I/D/=/X.
*/
int32_t get_length(std::vector<seqan3::cigar> const & cigar)
{
    using seqan3::operator""_cigar_operation;
    using seqan3::get;

    int32_t result{0};
    for (auto const & c : cigar)
    {
        seqan3::cigar::operation const op{get<1>(c)};
        uint32_t const length{get<0>(c)};
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

/*!
   \brief Parse an alignment file and store the records and the cumulative sum across the genome.
   \param input_path The path to the alignment file.
   \param record_list An empty vector which will store all of the records.
*/
void parse_file(std::filesystem::path const & input_path,
                std::vector<Record> & record_list)
{
    using my_fields = seqan3::fields<seqan3::field::ref_id,
                                     seqan3::field::ref_offset,
                                     seqan3::field::cigar>;
    seqan3::sam_file_input input{input_path, my_fields{}};

    // First make sure alingment file is sorted by coordinate.
    if (input.header().sorting != "coordinate")
    {
        throw seqan3::format_error{"ERROR: Input file must be sorted by coordinate (e.g. samtools sort)"};
    }

    for (auto const & r : input)
    {
        int32_t ref_id = r.reference_id().value_or(-1);
        int32_t position = r.reference_position().value_or(-1);
        Position start = std::make_pair(ref_id, position);
        Position end = std::make_pair(ref_id, position + get_length(r.cigar_sequence()));
        Record rec{start, end};
        if (ref_id != -1 && position != -1)
        {
            record_list.push_back(std::move(rec));
        }
    }
}
}
