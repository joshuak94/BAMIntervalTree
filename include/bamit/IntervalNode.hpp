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
    uint32_t start{}, end{}, chr{};
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
        \brief Get the chromosome for the current node.
        \return Returns the chromosom of the current node.
     */
     uint32_t const & get_chr()
     {
         return chr;
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
        \brief Set the chromosome value for the current node.
        \param c The calculated chromosome based on the records stored by this node.
     */
     void set_chr(uint32_t c)
     {
         this->chr = c;
     }

     /*!
        \brief Print the Interval Tree starting at this node.
        \param level The level of the current node.
     */
     void print(int32_t level)
     {
         std::string indent(level, '\t');
         seqan3::debug_stream << indent << "Level: " << level << '\n' <<
                                 indent << "Chromosome, start, end: " << chr << ", " << start << ", " << end << '\n' <<
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
        ar(this->chr, this->start, this->end, this->lNode, this->rNode, this->file_offset);
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
std::tuple<uint32_t, uint32_t> calculate_median(std::vector<std::vector<Record>> const & records_i)
{
    // Find out which chromosome to use for the median.
    auto size_calculator = [](int32_t a, std::vector<Record> b) { return a + b.size(); };
    int32_t middle = std::ceil(std::accumulate(records_i.cbegin(), records_i.cend(), 0, size_calculator) / 2.0);

    uint32_t chr_index{};
    while (middle > 0)
    {
        middle = middle - records_i[chr_index].size();
        ++chr_index;
    }
    --chr_index;

    // How many reads to the left and right of this chromosome?
    uint32_t left = (records_i[chr_index] == records_i.front()) ?
                    0 : std::accumulate(records_i.begin(), records_i.begin() + chr_index, 0, size_calculator);
    uint32_t right = (records_i[chr_index] == records_i.back()) ?
                     0 : std::accumulate(records_i.begin() + chr_index + 1, records_i.end(), 0, size_calculator);
    uint32_t diff_left = left > right ? left - right : 0;
    uint32_t diff_right = right > left ? right - left : 0;

    std::vector<uint32_t> values{};
    // If diff_left > 0, can skip that many reads at end of chromosome.
    // If diff_right > 0, can skip that many reads at start of chromosome.
    // If they are both == 0, skip no reads.
    values.reserve((records_i[chr_index].size() - diff_left - diff_right) * 2);
    for (auto it = records_i[chr_index].cbegin() + diff_right; it != records_i[chr_index].cend() - diff_left; ++it)
    {
        values.push_back((*it).start);
        values.push_back((*it).end);
    }
    std::sort(values.begin(), values.end());

    return std::make_tuple(chr_index, ((values[values.size() / 2] + values[(values.size() - 1) / 2]) / 2));
}

/*!
   \brief Construct an interval tree given a set of records.
   \param node The current node to fill.
   \param records_i The list of records to create the tree over.
   \param offset How many chromosomes are to the left of the current chromosome?
*/
void construct_tree(std::unique_ptr<IntervalNode> & node,
                    std::vector<std::vector<Record>> & records_i,
                    int32_t const offset = 0)
{
    // If there are no records, exit.
    if (records_i.empty() || (records_i.front().empty() && records_i.back().empty()))
    {
        return;
    }
    // Set node to an empty IntervalNode pointer.
    node = std::make_unique<IntervalNode>();

    // Calculate and set chromosome.
    auto [cur_chr, cur_median] = calculate_median(records_i);
    node->set_chr(cur_chr + offset);

    // Get reads which intersect median.
    std::vector<std::vector<Record>> lRecords{};
    std::vector<std::vector<Record>> rRecords{};

    // Add initial empty vector to rRecords.
    rRecords.emplace_back();
    for (size_t i = 0; i < records_i.size(); ++i)
    {
        // Chromosome is less than current chromosome.
        if (i < cur_chr)
        {
            lRecords.push_back(std::move(records_i[i]));
        }
        // Chromosome is greater than current chromosome.
        else if (i > cur_chr)
        {
            rRecords.push_back(std::move(records_i[i]));
        }
        // We are at current chromosome. Reads must be added to current node if they intersect interval.
        // Otherwise, they should be added to the left/right trees.
        else
        {
            std::vector<Record> l_cur_chr{};
            uint32_t start{0}, end{0};
            for (auto & r : records_i[i])
            {
                // Read ends before the median.
                if (r.end < cur_median)
                {
                    // Push it back into the first empty vector.
                    l_cur_chr.push_back(std::move(r));
                }
                // Read starts after the median.
                else if (r.start > cur_median)
                {
                    rRecords[0].push_back(std::move(r));
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
            lRecords.push_back(std::move(l_cur_chr));
        }
    }

    // Set left and right subtrees.
    construct_tree(node->get_left_node(), lRecords, offset);
    construct_tree(node->get_right_node(), rRecords, offset + cur_chr);
    return;
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
                             Position const & start,
                             Position const & end,
                             std::streamoff & offset_pos)
{
    if (!node)
    {
        return;
    }

    Position cur_start = std::make_tuple(node->get_chr(), node->get_start());
    Position cur_end = std::make_tuple(node->get_chr(), node->get_end());
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
                         std::unique_ptr<IntervalNode> const & node,
                         Position const & start,
                         Position const & end,
                         std::filesystem::path const & outname)
{
    std::streamoff offset_pos{-1};
    get_overlap_file_offset(node, start, end, offset_pos);

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
void write(std::unique_ptr<IntervalNode> const & node, Archive & archive)
{
    node->serialize(archive);
}

template <class Archive>
void read(std::unique_ptr<IntervalNode> & node, Archive & archive)
{
    node->serialize(archive);
}
}
