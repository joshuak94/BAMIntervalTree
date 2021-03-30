#pragma once
#include <seqan3/io/sam_file/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "Record.hpp"

namespace bamit
{
/*! The IntervalNode class stores a single node which is a part of an interval tree. It stores a list of
 *  bamit::Record objects that intersect a given median, along with pointers to its left and right children.
 */
class IntervalNode
{
private:
    Position median{};
    std::vector<Record> records{};
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
        \brief Get the records stored by the current node.
        \return Returns a reference to a vector of bamit::Record objects stored by the current node.
     */
     std::vector<Record> & get_records()
     {
         return records;
     }

     /*!
        \brief Set the median value for the current node.
        \param m The calculated median based on the records stored by this node.
     */
     void set_median(Position const & m)
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
                                 indent << "Reads: ";
         for (auto & r : records)
         {
             seqan3::debug_stream << "[" << r.start << ", " << r.end << "] ";
         }
         seqan3::debug_stream << "\n";

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

};

/*!
   \brief Calculate the median for a set of records based on the starts and ends of all records.
   \param records_i The list of records to calculate a median off of.
   \return Returns the median value and which chromosome it is in.
*/
Position calculate_median(std::vector<Record> const & records_i)
{
    std::vector<Position> values{};
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
     // The size of the vector values will always be an even number, therefore the median is just the leftmost.
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
    for (auto const & r : records_i)
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
        else
        {
            node->get_records().push_back(std::move(r));
        }
    }

    // Set left and right subtrees.
    construct_tree(node->get_left_node(), lRecords);
    construct_tree(node->get_right_node(), rRecords);
    return;
}

/*!
   \brief Find the records which overlap a given start and end position.
   \param root The root of the tree to search in.
   \param start The start position of the search.
   \param end The end position of the search.
   \param results The list of records overlapping the search.
*/
void overlap(std::unique_ptr<IntervalNode> const & root,
             Position const & start,
             Position const & end,
             std::vector<Record> & results)
{
    if (!root)
    {
        return;
    }

    Position cur_median = root->get_median();
    // If the current median is overlapping the read, add all records from this node and search the left and right.
    if (cur_median >= start && cur_median <= end)
    {
        results.insert(std::end(results), std::begin(root->get_records()), std::end(root->get_records()));
        overlap(root->get_left_node(), start, end, results);
        overlap(root->get_right_node(), start, end, results);
    }
    // If current median is to the right of the overlap, sort reads in ascending order and add all reads which
    // start before the overlap ends.
    else if (end < cur_median)
    {
        std::sort(root->get_records().begin(), root->get_records().end(), RecordComparatorStart());
        for (auto const & r : root->get_records())
        {
            if (r.start <= end)
            {
                results.push_back(r);
            }
            else
            {
                break;
            }
        }
        overlap(root->get_left_node(), start, end, results);
    }
    else if (start > cur_median)
    {
        std::sort(root->get_records().begin(), root->get_records().end(), RecordComparatorEnd());
        for (auto const & r : root->get_records())
        {
            if (r.end >= start)
            {
                results.push_back(r);
            }
            else
            {
                break;
            }
        }
        overlap(root->get_right_node(), start, end, results);
    }
}
}
