#pragma once

#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/io/sam_file/output.hpp>
#include <cereal/types/tuple.hpp>

namespace bamit
{

using Position = std::tuple<int32_t, int32_t>;

inline bool unmapped(auto const & rec)
{
    return (rec.reference_id().value_or(-1) == -1) ||
           (rec.reference_position().value_or(-1) == -1) ||
           (static_cast<bool>(rec.flag() & seqan3::sam_flag::unmapped));
}

/*! A Record object contains pertinent information about an alignment. */
struct Record
{
    uint32_t start{}, end{};
    std::streamoff file_position{-1};

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr Record()                        = default; //!< Defaulted.
    Record(Record const &)                    = default; //!< Defaulted.
    Record(Record &&)                         = default; //!< Defaulted.
    Record & operator=(Record const &)        = default; //!< Defaulted.
    Record & operator=(Record &&)             = default; //!< Defaulted.
    ~Record()                                 = default; //!< Defaulted.
     //!\}
    Record(uint32_t start_i, uint32_t end_i, std::streamoff file_position_i) :
        start{std::move(start_i)},
        end{std::move(end_i)},
        file_position{std::move(file_position_i)} {}

    /*!
       \brief Compare two Record objects and return true if they are equal.
       \param rhs The second record.
       \return Returns `true` if the two records are equal, with respect to chromosome and position.
    */
    bool operator==(Record const & rhs) const
    {
        return (start == rhs.start && end == rhs.end);
    }
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
    bool operator()(Record const & record1, Record const & record2)
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
    bool operator()(Record const & record1, Record const & record2)
    {
        return (record1.end == record2.end) ? (record1.start > record2.start) : (record1.end > record2.end);
    }
};

/*!
   \brief Get the length of a seqan3::cigar vector based on M/I/D/=/X operations.
   \param cigar The vector of seqan3::cigar characters.
   \return Returns the length of M/I/D/=/X.
*/
inline int32_t get_length(std::vector<seqan3::cigar> const & cigar)
{
    using seqan3::operator""_cigar_operation;
    using seqan3::get;

    int32_t result{0};
    for (auto const & c : cigar)
    {
        seqan3::cigar::operation const op{get<1>(c)};
        uint32_t const length{get<0>(c)};
        if (op == 'M'_cigar_operation || op == 'I'_cigar_operation || op == 'D'_cigar_operation ||
            op == '='_cigar_operation || op == 'X'_cigar_operation)
            result += length;
    }
    return result;
}
} // namespace bamit
