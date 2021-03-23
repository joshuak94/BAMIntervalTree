#include "IntervalNode.hpp"

void IntervalNode::set_left_node(IntervalNode * new_node)
{
    lNode = std::move(new_node);
}

void IntervalNode::set_right_node(IntervalNode * new_node)
{
    rNode = std::move(new_node);
}

IntervalNode * IntervalNode::get_left_node()
{
    return lNode;
}

IntervalNode * IntervalNode::get_right_node()
{
    return rNode;
}

std::vector<Record> & IntervalNode::get_records()
{
    return records;
}

IntervalNode * IntervalNode::construct_tree(std::vector<Record> & records_i)
{
    isLeaf = false;
    // Calculate median.
    int32_t median = this->calculate_median(records_i);
    // Get reads which intersect median.
    std::vector<Record> lRecords{};
    std::vector<Record> rRecords{};
    for (auto r : records_i)
    {
        if (median < r.start) // Median is to the right of the read, so read is in left subtree.
        {
            rRecords.push_back(r);
        }
        else if (median >= r.start && median <= r.end) // Median is within the read.
        {
            get_records().push_back(r);
        }
        else if (median > r.end) // Median is to the left of the read, so read is in right subtree.
        {
            lRecords.push_back(r);
        }
    }
    seqan3::debug_stream << "Median: " << median << ", Reads: ";
    for (auto r : get_records())
    {
        seqan3::debug_stream << "[" << r.start << ", " << r.end << "], ";
    }
    seqan3::debug_stream << std::endl;
    // If lRecords and rRecords are empty, it is a leaf.
    if (lRecords.empty() && rRecords.empty())
    {
        isLeaf = true;
        seqan3::debug_stream << "Is leaf! Median: " << median << std::endl;
        return this;
    }
    // seqan3::debug_stream << lRecords.size() << ", " << rRecords.size() << std::endl;

    // Set left and right subtrees.
    if (!lRecords.empty())
    {
        IntervalNode lInterval{NULL, NULL};
        set_left_node(lInterval.construct_tree(lRecords));
    }
    if (!rRecords.empty())
    {
    }

    // seqan3::debug_stream << lNode->get_records().size() << ", " << rNode->get_records().size() << std::endl;
        IntervalNode rInterval{NULL, NULL};
        set_right_node(rInterval.construct_tree(rRecords));
    return this;
}

// TODO: Fix the printing.
void IntervalNode::print(int32_t level)
{
    seqan3::debug_stream << "Level: " << level << '\n' << "Reads: ";
    for (auto r : get_records())
    {
        seqan3::debug_stream << "[" << r.start << ", " << r.end << "], ";
    }
    seqan3::debug_stream << "\n\n";

    if (get_left_node() != NULL)
    {
        seqan3::debug_stream << "left node... \n";
        (*get_left_node()).print(level + 1);
    }
    if (get_right_node() != NULL)
    {
        seqan3::debug_stream << "right node... \n";
        (*get_right_node()).print(level + 1);
    }
}

int32_t IntervalNode::calculate_median(std::vector<Record> & records_i)
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
    if (size % 2 == 0)
    {
      median = (values[size / 2 - 1] + values[size / 2]) / 2;
    }
    else
    {
      median = values[size / 2];
    }

    return median;
}
