#pragma once
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/input.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/tuple.hpp>
#include <cereal/types/vector.hpp>

#include <bamit/Record.hpp>

namespace bamit
{
/*! The IntervalNode class stores a single node which is a part of an interval tree. It stores a file offset to the
 *  first read which intersects the median, along with pointers to its left and right children.
 */
class IntervalNode
{
private:
    Position median{};
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
        \brief Get the median for the current node.
        \return Returns the median of the current node.
     */
     Position const & get_median()
     {
         return median;
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
        \brief Set the median value for the current node.
        \param m The calculated median based on the records stored by this node.
     */
     void set_median(Position m)
     {
         this->median = std::move(m);
     }

     /*!
        \brief Print the Interval Tree starting at this node.
        \param level The level of the current node.
     */
     void print(int32_t level)
     {
         std::string indent(level, '\t');
         seqan3::debug_stream << indent << "Level: " << level << '\n' <<
                                 indent << "Median(chromosome, position): " << std::get<0>(this->get_median()) <<
                                 ", " << std::get<1>(this->get_median()) << '\n' <<
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
        ar(median, this->lNode, this->rNode, this->file_offset);
    }

};

/*!
   \brief Calculate the median for a set of records based on the starts and ends of all records.
   \param records_i The list of records to calculate a median off of.
   \return Returns the median value and which chromosome it is in.

   The median is calculated by sorting all of the starts and ends from a list of records. Since each record
   has a start and end, the list is an even length and the median is the average of the middle two positions.
   If the middle two positions are from the same chromsome, then they are averaged. If they're from two different
   chromosomes then averaging is not possible and the right-most is taken.
*/
Position calculate_median(std::vector<Record> const & records_i)
{
    std::vector<Position> values{};
    values.reserve(records_i.size() * 2);
    Position median{};
    for (auto const & r : records_i)
    {
        values.push_back(r.start);
        values.push_back(r.end);
    }
    std::sort(values.begin(), values.end());
    size_t size = values.size();

    if (std::get<0>(values[size / 2 - 1]) == std::get<0>(values[size / 2]))
    {
        int32_t median_value = (std::get<1>(values[size / 2 - 1]) + std::get<1>(values[size / 2])) / 2;
        median = std::make_pair(std::get<0>(values[size / 2]), median_value);
    }
    else
    {
        median = values[size / 2];
    }
    return median;
}

/*!
   \brief Construct an interval tree given a set of records.
   \param node The current node to fill.
   \param records_i The list of records to create the tree over.
*/
void construct_tree(std::unique_ptr<IntervalNode> & node, std::vector<Record> const & records_i)
{
    if (records_i.empty())
    {
        return;
    }
    // Set node to an empty IntervalNode pointer.
    node = std::make_unique<IntervalNode>();

    // Calculate and set median.
    node->set_median(calculate_median(records_i));
    Position cur_median = node->get_median();

    // Get reads which intersect median.
    std::vector<Record> lRecords{};
    std::vector<Record> rRecords{};
    for (auto & r : records_i)
    {
        // Median is to the left of the read, so read is in right subtree.
        if (cur_median < r.start)
        {
            rRecords.push_back(std::move(r));
        }
        // Median is to the right of the read, so read is in left subtree.
        else if (cur_median > r.end)
        {
            lRecords.push_back(std::move(r));
        }
        // If median is not to the left or right, it must be within the read.
        // Set the file_offset for this node only if it has not yet been set.
        else if (node->get_file_offset() == -1)
        {
            node->set_file_offset(r.file_offset);
        }
    }

    // Set left and right subtrees.
    construct_tree(node->get_left_node(), lRecords);
    construct_tree(node->get_right_node(), rRecords);
    return;
}

/*!
   \brief Find the file offstream position that is close to the start of the range by traversing the tree. The offstream
          position is guaranteed to be to the left of the start.
   \param input The sam file input of type bamit::sam_file_input_type.
   \param root The root of the tree to search in.
   \param start The start position of the search.
   \param end The end position of the search.
   \param offstream_pos The resulting offstream position.
*/
void overlap(std::unique_ptr<IntervalNode> const & root,
             Position const & start,
             Position const & end,
             std::streamoff & offstream_pos)
{
    if (!root)
    {
        return;
    }

    Position cur_median = root->get_median();
    // If the current median is overlapping the read, take the current offset as result and search the left further.
    if (cur_median >= start && cur_median <= end)
    {
        offstream_pos = root->get_file_offset();
        overlap(root->get_left_node(), start, end, offstream_pos);
    }
    // If current median is to the right of the overlap, go to the left tree.
    else if (end < cur_median)
    {
        offstream_pos = root->get_file_offset();
        overlap(root->get_left_node(), start, end, offstream_pos);
    }
    // If current median is to the left of the overlap, go to the right tree.
    else if (start > cur_median)
    {
        offstream_pos = root->get_file_offset();
        overlap(root->get_right_node(), start, end, offstream_pos);
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
