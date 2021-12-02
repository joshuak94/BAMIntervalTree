#pragma once
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/utility/type_pack/traits.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/vector.hpp>

#include <bamit/Record.hpp>

#include <numeric>

namespace bamit
{
/*! The IntervalNode class stores a single node which is a part of an interval tree. It stores a file position to the
 *  first read which intersects the median, along with pointers to its left and right children.
 *  Additionally, it stores the chromosome it is in, the start of the left-most record and the end of the right-most
 *  record.
 */
class IntervalNode
{
private:
    uint32_t start{}, end{};
    std::streamoff file_position{-1};
    std::unique_ptr<IntervalNode> lNode{nullptr};
    std::unique_ptr<IntervalNode> rNode{nullptr};
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr IntervalNode()                        = default; //!< Defaulted.
    IntervalNode(IntervalNode const &)              = default; //!< Defaulted.
    IntervalNode(IntervalNode &&)                   = default; //!< Defaulted.
    IntervalNode & operator=(IntervalNode const &)  = default; //!< Defaulted.
    IntervalNode & operator=(IntervalNode &&)       = default; //!< Defaulted.
    ~IntervalNode()                                 = default; //!< Defaulted.
     //!\}

     /*!
        \brief Get the left child node of current node.
        \return Returns a std::unique_ptr reference to the left child node.
     */
     std::unique_ptr<IntervalNode> & get_left_node()
     {
         return lNode;
     }

     /*!
        \brief Get the right child node of current node.
        \return Returns a std::unique_ptr reference to the right child node.
     */
     std::unique_ptr<IntervalNode> & get_right_node()
     {
         return rNode;
     }

     /*!
        \brief Get the start for the current node.
        \return Returns the start position of the left-most record which is stored in this node.
     */
     uint32_t const & get_start() const
     {
         return start;
     }

     /*!
        \brief Get the end for the current node.
        \return Returns the end position of the right-most record which is stored in this node.
     */
     uint32_t const & get_end() const
     {
         return end;
     }

     /*!
        \brief Get the file position to the first read stored by this node.
        \return Returns a reference to a std::streamoff which can be used to seek to a file position.
     */
     std::streamoff const & get_file_position() const
     {
         return file_position;
     }

     /*!
        \brief Set the file position to the first read for the current node.
        \param new_file_position The file position based on the records stored by this node.
     */
     void set_file_position(std::streamoff new_file_position)
     {
         this->file_position = std::move(new_file_position);
     }

     /*!
        \brief Set the start value for the current node.
        \param s The start of the left-most record stored by this node.
     */
     void set_start(uint32_t s)
     {
         this->start = std::move(s);
     }

     /*!
        \brief Set the end value for the current node.
        \param e The end of the left-most record stored by this node.
     */
     void set_end(uint32_t e)
     {
         this->end = std::move(e);
     }

     /*!
        \brief Print the Interval Tree starting at this node.
        \param level The level of the current node.
     */
     void print(int32_t level)
     {
         std::string indent(level, '\t');
         seqan3::debug_stream << indent << "Level: " << level << '\n' <<
                                 indent << "Start, end: " << start << ", " << end << '\n' <<
                                 indent << "File position: " << this->get_file_position() << '\n';
         if (lNode)
         {
             seqan3::debug_stream << indent << "left node... \n";
             lNode->print(level + 1);
         }
         if (rNode)
         {
             seqan3::debug_stream << indent << "right node... \n";
             rNode->print(level + 1);
         }
     }

    template <class Archive>
    void serialize(Archive & ar)
    {
        ar(this->start, this->end, this->lNode, this->rNode, this->file_position);
    }

};

/*!
   \brief Calculate the median for a set of records based on the starts and ends of all records.
   \param records_i The list of records from which the median is computed.
   \return Returns a tuple of chromosome of median and the median value.

   The median is calculated by sorting all of the starts and ends from a list of records. Since each record
   has a start and end, the list is an even length and the median is the average of the middle two positions.
*/
inline uint32_t calculate_median(std::vector<Record> const & records_i)
{
    std::vector<uint32_t> values{};
    values.reserve(records_i.size() * 2);
    for (auto const & r : records_i)
    {
        values.push_back(r.start);
        values.push_back(r.end);
    }
    std::sort(values.begin(), values.end());

    return (values[values.size() / 2] + values[(values.size() - 1) / 2]) / 2;
}

/*!
   \brief Construct an interval tree given a set of records.
   \param node The current node to fill.
   \param records_i The list of records to create the tree over.
*/
inline void construct_tree(std::unique_ptr<IntervalNode> & node,
                           std::vector<Record> & records_i)
{
    // If there are no records, exit.
    if (records_i.empty()) return;
    // Set node to an empty IntervalNode pointer.
    node = std::make_unique<IntervalNode>();

    // Calculate and set median.
    auto cur_median = calculate_median(records_i);

    // Get reads which intersect median.
    std::vector<Record> lRecords{};
    std::vector<Record> rRecords{};

    uint32_t start{0}, end{0};
    for (auto & r : records_i)
    {
        // Read ends before the median.
        if (r.end < cur_median) lRecords.push_back(std::move(r));
        // Read starts after the median.
        else if (r.start > cur_median) rRecords.push_back(std::move(r));
        // Read intersects the median. Only store file position and start from the left-most read!
        // End is always updated while the read intersects the median.
        else
        {
            if (node->get_file_position() == -1)
            {
                node->set_file_position(r.file_position);
                start = r.start;
            }
            end = r.end;
        }
    }
    node->set_start(start);
    node->set_end(end);

    // Set left and right subtrees.
    construct_tree(node->get_left_node(), lRecords);
    construct_tree(node->get_right_node(), rRecords);
    return;
}

/*!
   \brief Entry point into the recursive tree construction.
   \param input_file The input file to construct the tree over.
   \param verbose Print verbose output.
   \tparam traits_type The type of the traits for seqan3::sam_file_input
   \tparam fields_type The given fields.
   \tparam format_type The format of the file.
   \return Returns a vector of IntervalNodes, each of which is the root node of an Interval Tree over its respective
           chromosome.
*/
template <typename traits_type, typename fields_type, typename format_type>
inline std::vector<std::unique_ptr<IntervalNode>> index(seqan3::sam_file_input<traits_type,
                                                                               fields_type,
                                                                               format_type> & input_file,
                                                        bool const & verbose = false)
{
    // Very first thing: Check that required fields are non-empty.
    static_assert(fields_type::contains(seqan3::field::ref_id),
                  "Input file must define fields seqan3::field::ref_id, seqan3::field::ref_offset, seqan3::field::cigar, and seqan3::field::flag");
    static_assert(fields_type::contains(seqan3::field::ref_offset),
                  "Input file must define fields seqan3::field::ref_id, seqan3::field::ref_offset, seqan3::field::cigar, and seqan3::field::flag");
    static_assert(fields_type::contains(seqan3::field::cigar),
                  "Input file must define fields seqan3::field::ref_id, seqan3::field::ref_offset, seqan3::field::cigar, and seqan3::field::flag");
    static_assert(fields_type::contains(seqan3::field::flag),
                  "Input file must define fields seqan3::field::ref_id, seqan3::field::ref_offset, seqan3::field::cigar, and seqan3::field::flag");
    // First make sure alingment file is sorted by coordinate.
    if (input_file.header().sorting != "coordinate")
        throw seqan3::format_error{"ERROR: Input file must be sorted by coordinate (e.g. samtools sort)"};

    // Vector containing result, initialized to # of chromosomes in header.
    std::vector<std::unique_ptr<IntervalNode>> result;
    result.reserve(input_file.header().ref_ids().size());
    std::generate_n(std::back_inserter(result), input_file.header().ref_ids().size(),
                    [] { return std::make_unique<IntervalNode>(); });

    // List of records for a single chromosome.
    std::vector<Record> cur_records;

    uint32_t cur_index{0};
    for (auto it = input_file.begin(); it != input_file.end(); ++it)
    {
        if (unmapped(*it)) continue;
        uint32_t ref_id = (*it).reference_id().value();
        uint32_t position = (*it).reference_position().value();
        if (ref_id != cur_index)
        {
            if (verbose) seqan3::debug_stream << "Indexing chr " << input_file.header().ref_ids()[cur_index] << "...";
            construct_tree(result[cur_index], cur_records);
            if (verbose) seqan3::debug_stream << " Done!\n";
            cur_records.clear();
            ++cur_index;
        }
        cur_records.emplace_back(position,
                                 position + get_length((*it).cigar_sequence()),
                                 static_cast<std::streamoff>(it.file_position()));
    }
    if (verbose) seqan3::debug_stream << "Indexing chr " << input_file.header().ref_ids()[cur_index] << "...";
    construct_tree(result[cur_index], cur_records);
    if (verbose) seqan3::debug_stream << " Done!\n";

    return result;
}

/*!
   \brief Find the closest file offset to an overlap query which is stored in the Interval Tree. This may not be the
          record which actually overlaps the query, but it is guaranteed to be to the left of the query.
   \param node The current node to search.
   \param start The start position of the search.
   \param end The end position of the search.
   \param file_position The resulting file position.
   \details This function traverses an interval tree looking for the file position of the closest record to the overlap
            query stored within the tree. Note that this record may not actually overlap the query itself, but is
            guaranteed to be to the left of the overlap query. Thus, this function should be paired with
            bamit::get_correct_position, which will then update the file position to the first record which actually
            overlaps the query. These two functions are used by bamit::get_overlap_file_position to obtain the file
            position of the first record overlapping a query.
 */
inline void get_current_file_position(std::unique_ptr<IntervalNode> const & node,
                                      uint32_t const & start,
                                      uint32_t const & end,
                                      std::streamoff & file_position)
{
    if (!node) return;

    uint32_t cur_start = node->get_start();
    uint32_t cur_end = node->get_end();
    /*
     * There are six possibilities:
     * 1. The query interval is wholly to the left of the current start (end < cur_start).
     * 2. The query interval is partially to the left of the current start (start < cur_start && end >= cur_start).
     * 3. The query interval is contained within the current start and end (start >= cur_start && end <= cur_end).
     * 4. The query interval contains the current start and end (start < cur_start && end > cur_end).
     * 5. The query interval is partially to the right of the current end (start <= cur_end && end > cur_end).
     * 6. The query interval is wholly to the right of the current end (start > cur_end).
     *
     * In case 1, do not store the file position and search the left subtree.
     * In cases 2 through 5, store the file position and search the left subtree.
     * In case 6, do not store the file position and search the right subtree.
    */
    if (end < cur_start)
        get_current_file_position(node->get_left_node(), start, end, file_position);
    else if (start > cur_end)
        get_current_file_position(node->get_right_node(), start, end, file_position);
    else
    {
        file_position = node->get_file_position();
        get_current_file_position(node->get_left_node(), start, end, file_position);
    }
}

/*!
   \brief Move the file position along a BAM/SAM file until it points to the first read overlapping a start position.
   \param input The alignment file to use.
   \param start The start of the given query.
   \param file_position The position to move along.
 */
template <typename traits_type, typename fields_type, typename format_type>
inline void get_correct_position(seqan3::sam_file_input<traits_type, fields_type, format_type> & input,
                                 Position const & start,
                                 std::streamoff & file_position)
{
    // Very first thing: Check that required fields are non-empty.
    static_assert(fields_type::contains(seqan3::field::ref_id),
                  "Input file must define fields seqan3::field::ref_id, seqan3::field::ref_offset, seqan3::field::cigar, and seqan3::field::flag");
    static_assert(fields_type::contains(seqan3::field::ref_offset),
                  "Input file must define fields seqan3::field::ref_id, seqan3::field::ref_offset, seqan3::field::cigar, and seqan3::field::flag");
    static_assert(fields_type::contains(seqan3::field::cigar),
                  "Input file must define fields seqan3::field::ref_id, seqan3::field::ref_offset, seqan3::field::cigar, and seqan3::field::flag");
    static_assert(fields_type::contains(seqan3::field::flag),
                  "Input file must define fields seqan3::field::ref_id, seqan3::field::ref_offset, seqan3::field::cigar, and seqan3::field::flag");

    auto it = input.begin();
    it.seek_to(static_cast<std::streampos>(file_position));
    for (; it != input.end(); ++it)
    {
        if (unmapped(*it)) continue;
        // If the read ends either at or after the start, it is the first read and its position should be returned.
        if (std::make_tuple((*it).reference_id().value(),
                            (*it).reference_position().value() + get_length((*it).cigar_sequence())) >= start)
        {
            file_position = static_cast<std::streamoff>(it.file_position());
            return;
        }
    }
    // Don't think this is possible for the for loop to exit completely, but just in case...
    seqan3::debug_stream << "[ERROR] Improper file position/input file given.\n";
    file_position = -1;
}

/*!
   \brief Obtain the file position of the first record which overlaps a query.
   \param input The sam file input of type bamit::seqan3::sam_file_input.
   \param node_list The list of interval nodes per chromosome
   \param start The start Position of the search.
   \param end The end Position of the search.
   \param file_position The resulting file position.
   \details The main function for obtaining the file position of an overlap query.
 */
template <typename traits_type, typename fields_type, typename format_type>
inline void get_overlap_file_position(seqan3::sam_file_input<traits_type, fields_type, format_type> & input,
                                      std::vector<std::unique_ptr<IntervalNode>> const & node_list,
                                      Position const & start,
                                      Position const & end,
                                      std::streamoff & file_position)
{
    if (std::get<0>(start) == std::get<0>(end)) // Searching in one chromosome.
    {
        get_current_file_position(node_list[std::get<0>(start)], std::get<1>(start), std::get<1>(end), file_position);
        if (file_position != -1) get_correct_position(input, start, file_position);
    }
    else // Searching across multiple chromosomes.
    {
        for (uint32_t i = std::get<0>(start); i < (uint32_t) std::get<0>(end); ++i)
        {
            // Start at given start only for the first chromosome, otherwise start searching from 0.
            // For the first chromosome, we want the actual start position of the query. If no reads
            // in the chromosome are in our query, we want to search the tree of the next chromosome
            // starting from the beginning.
            uint32_t start_position = (uint32_t) std::get<0>(start) == i ? (uint32_t) std::get<1>(start) : 0;
            // End at given end if at the last chromosome, otherwise end at end of chromosome.
            uint32_t end_position = (uint32_t) std::get<0>(end) == i ? (uint32_t) std::get<1>(end) : std::numeric_limits<uint32_t>::max();
            get_current_file_position(node_list[i], start_position, end_position, file_position);

            // If we find the left-most position we can stop. Otherwise we have to check the next tree if
            // the overlap spans more than one chromosome.
            if (file_position != -1)
            {
                get_correct_position(input, start, file_position);
                break;
            }
        }
    }
}

/*!
   \brief Find the records which overlap a given start and end position.
   \param input The sam file input of type bamit::seqan3::sam_file_input.
   \param node_list The list of interval trees.
   \param start The start position of the search.
   \param end The end position of the search.
   \param verbose Print verbose output.
   \param outname The output filename. If not provided the function will only return the file position and
                  not write to any file.

   \return Returns a vector of seqan3::sam_record objects containing records overlapping the query.
   \details The main function for obtaining a vector of records which overlap a query. If just the file position is
            desired, use bamit::get_overlap_file_position instead.
*/
template <typename traits_type, typename fields_type, typename format_type>
inline auto get_overlap_records(seqan3::sam_file_input<traits_type, fields_type, format_type> & input,
                                std::vector<std::unique_ptr<IntervalNode>> const & node_list,
                                Position const & start,
                                Position const & end,
                                bool const & verbose = false,
                                std::filesystem::path const & outname = "")
{
    // Very first thing: Check that required fields are non-empty.
    static_assert(fields_type::contains(seqan3::field::ref_id),
                  "Input file must define fields seqan3::field::ref_id, seqan3::field::ref_offset, seqan3::field::cigar, and seqan3::field::flag");
    static_assert(fields_type::contains(seqan3::field::ref_offset),
                  "Input file must define fields seqan3::field::ref_id, seqan3::field::ref_offset, seqan3::field::cigar, and seqan3::field::flag");
    static_assert(fields_type::contains(seqan3::field::cigar),
                  "Input file must define fields seqan3::field::ref_id, seqan3::field::ref_offset, seqan3::field::cigar, and seqan3::field::flag");
    static_assert(fields_type::contains(seqan3::field::flag),
                  "Input file must define fields seqan3::field::ref_id, seqan3::field::ref_offset, seqan3::field::cigar, and seqan3::field::flag");

    std::streamoff file_position{-1};

    // Get the file position of the first record matching start query.
    get_overlap_file_position(input, node_list, start, end, file_position);

    // Store reads which start before the end of the query, filtering out unmapped reads.
    auto results_list = input | std::views::take_while([file_position](auto & rec) {return file_position != -1;})
                              | std::views::take_while([end](auto & rec) {return std::make_tuple(rec.reference_id().value(), rec.reference_position().value()) < end;})
                              | std::views::filter([](auto & rec) {return !unmapped(rec);})
                              | seqan3::views::to<std::vector>;
    if (results_list.empty() && verbose)
    {
        seqan3::debug_stream << "No overlapping reads found for query "
                             << input.header().ref_ids()[std::get<0>(start)] << ":" << std::get<1>(start) << " through "
                             << input.header().ref_ids()[std::get<0>(end)] << ":" << std::get<1>(end) << "\n";
    }
    if (!outname.empty()) // Outputs an empty file if the list is empty.
    {
        // Need to extract chromosome lengths for the output header file.
        std::vector<int32_t> ref_lengths{};
        std::transform(std::begin(input.header().ref_id_info), std::end(input.header().ref_id_info),
                       std::back_inserter(ref_lengths), [](auto const & pair){ return std::get<0>(pair); });
        seqan3::sam_file_output fout{outname, input.header().ref_ids(), ref_lengths};
        results_list | fout;
    }
    return results_list;
}

template <class Archive>
inline void write(std::vector<std::unique_ptr<IntervalNode>> const & node_list, Archive & archive)
{
    archive(node_list);
}

template <class Archive>
inline void read(std::vector<std::unique_ptr<IntervalNode>> & node_list, Archive & archive)
{
    archive(node_list);
}
} // namespace bamit
