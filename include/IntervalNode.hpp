#pragma once
#include <seqan3/io/sam_file/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "Record.hpp"

namespace bamit
{
class IntervalNode
{
using my_fields = seqan3::fields<seqan3::field::id,
                                 seqan3::field::ref_id,
                                 seqan3::field::ref_offset,
                                 seqan3::field::cigar>;
using sam_file_input_t = seqan3::sam_file_input<seqan3::sam_file_input_default_traits<>,
                                                my_fields,
                                                seqan3::type_list<seqan3::format_sam, seqan3::format_bam>>;
private:
    bool isLeaf{false};
    int32_t median{};
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

     std::unique_ptr<IntervalNode> & get_left_node()
     {
         return lNode;
     }

     std::unique_ptr<IntervalNode> & get_right_node()
     {
         return rNode;
     }

     int32_t const & get_median()
     {
         return median;
     }

     std::vector<Record> & get_records()
     {
         return records;
     }

     void set_is_leaf()
     {
         this->isLeaf = true;
     }

     void set_median(int32_t const & m)
     {
         this->median = m;
     }

     void print(int32_t level)
     {
         seqan3::debug_stream << "Level: " << level << '\n' << "Median: " << this->get_median() << '\n' << "Reads: ";
         for (auto & r : records)
         {
             seqan3::debug_stream << "[" << r.start << ", " << r.end << "]";
         }
         seqan3::debug_stream << "\n\n";

         if (lNode)
         {
             seqan3::debug_stream << "left node... \n";
             lNode->print(level + 1);
         }
         if (rNode)
         {
             seqan3::debug_stream << "right node... \n";
             rNode->print(level + 1);
         }
     }

};

int32_t calculate_median(std::vector<Record> const & records_i)
{
    std::vector<int32_t> values{};
    int32_t median{};
    for (auto r : records_i)
    {
        values.push_back(r.start);
        values.push_back(r.end);
    }
    std::sort(values.begin(), values.end());
    size_t size = values.size();

     // The size of the vector values will always be an even number, therefore the median is determined by taking the
     // mean of the two values in the middle.
     median = (values[size / 2 - 1] + values[size / 2]) / 2;

    return median;
}


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

    // Get reads which intersect median.
    std::vector<Record> lRecords{};
    std::vector<Record> rRecords{};
    for (const auto & r : records_i)
    {
        if (node->get_median() < r.start) // Median is to the left of the read, so read is in right subtree.
        {
            rRecords.push_back(std::move(r));
        }
        else if (node->get_median() >= r.start && node->get_median() <= r.end) // Median is within the read.
        {
            node->get_records().push_back(std::move(r));
        }
        else if (node->get_median() > r.end) // Median is to the right of the read, so read is in left subtree.
        {
            lRecords.push_back(std::move(r));
        }
    }
    // seqan3::debug_stream << "Median: " << node->get_median() << ", Reads: ";
    // for (auto & r : node->get_records())
    // {
    //     seqan3::debug_stream << "[" << r.start << ", " << r.end << "], ";
    // }
    // seqan3::debug_stream << std::endl;

    // If lRecords and rRecords are empty, it is a leaf.
    if (lRecords.empty() && rRecords.empty())
    {
        node->set_is_leaf();
        // seqan3::debug_stream << "Is leaf! Median: " << median << std::endl;
        return;
    }

    // Set left and right subtrees.
    construct_tree(node->get_left_node(), lRecords);
    construct_tree(node->get_right_node(), rRecords);
    return;
}

void overlap(std::unique_ptr<IntervalNode> & root, int32_t start, int32_t end, std::vector<Record> & results)
{
    if (!root)
    {
        return;
    }

    if (root->get_median() >= start && root->get_median() <= end)
    {
        for (auto & r : root->get_records())
        {
            results.push_back(r);
        }
        overlap(root->get_left_node(), start, end, results);
        overlap(root->get_right_node(), start, end, results);
    }
    else if (end < root->get_median())
    {
        std::sort(root->get_records().begin(), root->get_records().end(), RecordComparatorStart());
        for (auto & r : root->get_records())
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
    else if (start > root->get_median())
    {
        std::sort(root->get_records().begin(), root->get_records().end(), RecordComparatorEnd());
        for (auto & r : root->get_records())
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
