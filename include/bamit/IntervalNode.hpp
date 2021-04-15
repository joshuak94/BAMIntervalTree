#pragma once
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/input.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/vector.hpp>

#include <bamit/Record.hpp>

#include <numeric>

namespace bamit
{
/*! The IntervalNode class stores a single node which is a part of an interval tree. It stores a file offset to the
 *  first read which intersects the median, along with pointers to its left and right children.
 *  Additionally, it stores the chromosome it is in, the start of the left-most record and the end of the right-most
 *  record.
 */
class IntervalNode
{
private:
    uint32_t start{}, end{};
    std::streamoff file_offset{-1};
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
     uint32_t const & get_start()
     {
         return start;
     }

     /*!
        \brief Get the end for the current node.
        \return Returns the end position of the right-most record which is stored in this node.
     */
     uint32_t const & get_end()
     {
         return end;
     }

     /*!
        \brief Get the file offset to the first read stored by this node.
        \return Returns a reference to a std::streamoff which can be used to seek to a file position.
     */
     std::streamoff const & get_file_offset()
     {
         return file_offset;
     }

     /*!
        \brief Set the file offset to the first read for the current node.
        \param new_file_offset The file offset based on the records stored by this node.
     */
     void set_file_offset(std::streamoff new_file_offset)
     {
         this->file_offset = std::move(new_file_offset);
     }

     /*!
        \brief Set the start value for the current node.
        \param s The start of the left-most record stored by this node.
     */
     void set_start(uint32_t s)
     {
         this->start = s;
     }

     /*!
        \brief Set the end value for the current node.
        \param e The end of the left-most record stored by this node.
     */
     void set_end(uint32_t e)
     {
         this->end = e;
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
                                 indent << "File offset: " << this->get_file_offset() << '\n';
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
        ar(this->start, this->end, this->lNode, this->rNode, this->file_offset);
    }

};

/*!
   \brief Calculate the median for a set of records based on the starts and ends of all records.
   \param records_i The list of records to calculate a median off of.
   \return Returns a tuple of chromosome of median and the median value.

   The median is calculated by sorting all of the starts and ends from a list of records. Since each record
   has a start and end, the list is an even length and the median is the average of the middle two positions.
   If the middle two positions are from the same chromsome, then they are averaged. If they're from two different
   chromosomes then averaging is not possible and the right-most is taken.
*/
uint32_t calculate_median(std::vector<Record> const & records_i)
{
    std::vector<uint32_t> values{};
    values.reserve(records_i.size() * 2);
    for (auto it = records_i.cbegin(); it != records_i.cend(); ++it)
    {
        values.push_back((*it).start);
        values.push_back((*it).end);
    }
    std::sort(values.begin(), values.end());

    return (values[values.size() / 2] + values[(values.size() - 1) / 2]) / 2;
}

/*!
   \brief Construct an interval tree given a set of records.
   \param node The current node to fill.
   \param records_i The list of records to create the tree over.
   \param offset How many chromosomes are to the left of the current chromosome?
*/
void construct_tree(std::unique_ptr<IntervalNode> & node,
                    std::vector<Record> & records_i)
{
    // If there are no records, exit.
    if (records_i.empty())
    {
        return;
    }
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
        if (r.end < cur_median)
        {
            // Push it back into the first empty vector.
            lRecords.push_back(std::move(r));
        }
        // Read starts after the median.
        else if (r.start > cur_median)
        {
            rRecords.push_back(std::move(r));
        }
        // Read intersects the median. Only store file offset and start from the left-most read!
        // End is always updated while the read intersects the median.
        else
        {
            if (node->get_file_offset() == -1)
            {
                node->set_file_offset(r.file_offset);
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
   \return Returns a vector of IntervalNodes, each of which is the root node of an Interval Tree over its respective
           chromosome.
*/
std::vector<std::unique_ptr<IntervalNode>> construct_tree(sam_file_input_type & input_file)
{
    // First make sure alingment file is sorted by coordinate.
    if (input_file.header().sorting != "coordinate")
    {
        throw seqan3::format_error{"ERROR: Input file must be sorted by coordinate (e.g. samtools sort)"};
    }

    // Vector containing result.
    std::vector<std::unique_ptr<IntervalNode>> result{};

    // List of records for a single chromosome.
    std::vector<Record> cur_records;

    uint32_t cur_index{0};
    for (auto const & r : input_file | properly_mapped)
    {
        uint32_t ref_id = r.reference_id().value();
        uint32_t position = r.reference_position().value();
        if (ref_id != cur_index)
        {
            std::unique_ptr<IntervalNode> root{nullptr};
            construct_tree(root, cur_records);
            cur_records.clear();
            result.push_back(std::move(root));
            ++cur_index;
        }
        cur_records.emplace_back(position, position + get_length(r.cigar_sequence()), r.file_offset());
    }
    std::unique_ptr<IntervalNode> root{nullptr};
    construct_tree(root, cur_records);
    result.push_back(std::move(root));

    return result;
}

/*!
   \brief Find the file offstream position that is close to the start of the range by traversing the tree. The offstream
          position is guaranteed to be to the left of the start.
   \param node The current node to search.
   \param start The start position of the search.
   \param end The end position of the search.
   \param offset_pos The resulting offset position.
*/
void get_overlap_file_offset(std::unique_ptr<IntervalNode> const & node,
                             uint32_t const & start,
                             uint32_t const & end,
                             std::streamoff & offset_pos)
{
    if (!node)
    {
        return;
    }

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
     * In case 1, do not store the file offset and search the left subtree.
     * In cases 2 through 5, store the file offset and search the left subtree.
     * In case 6, do not store the file offset and search the right subtree.
    */
    if (end < cur_start)
    {
        get_overlap_file_offset(node->get_left_node(), start, end, offset_pos);
    }
    else if (start > cur_end)
    {
        get_overlap_file_offset(node->get_right_node(), start, end, offset_pos);
    }
    else
    {
        offset_pos = node->get_file_offset();
        get_overlap_file_offset(node->get_left_node(), start, end, offset_pos);
    }
}

/*!
   \brief Find the records which overlap a given start and end position.
   \param input The sam file input of type bamit::sam_file_input_type.
   \param node The current node to search.
   \param start The start position of the search.
   \param end The end position of the search.
   \param outname The output filename.
*/
void get_overlap_records(sam_file_input_type & input,
                         std::vector<std::unique_ptr<IntervalNode>> const & node_list,
                         Position const & start,
                         Position const & end,
                         std::filesystem::path const & outname)
{
    std::streamoff offset_pos{-1};

    // Look for the file offset using the tree for the starting chromosome.
    get_overlap_file_offset(node_list[std::get<0>(start)], std::get<1>(start), std::get<1>(end), offset_pos);

    seqan3::sam_file_output fout{outname, bamit::sam_file_output_fields{}};

    if (offset_pos == -1)
    {
        seqan3::debug_stream << "No overlapping reads found.\n";
        return;
    }

    input.seek(offset_pos);

    for (auto & r : input | properly_mapped)
    {
        if (std::make_tuple(r.reference_id(), r.reference_position()) > end)
        {
            break;
        }
        if (std::make_tuple(r.reference_id(), r.reference_position().value() + get_length(r.cigar_sequence())) < start)
        {
            continue;
        }

        fout.push_back(r);
   }
}

template <class Archive>
void write(std::vector<std::unique_ptr<IntervalNode>> const & node_list, Archive & archive)
{
    archive(node_list);
}

template <class Archive>
void read(std::vector<std::unique_ptr<IntervalNode>> & node_list, Archive & archive)
{
    archive(node_list);
}
}
