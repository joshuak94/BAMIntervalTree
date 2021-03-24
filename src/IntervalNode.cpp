#include "IntervalNode.hpp"

std::unique_ptr<IntervalNode> & IntervalNode::get_left_node()
{
    return lNode;
}

std::unique_ptr<IntervalNode> & IntervalNode::get_right_node()
{
    return rNode;
}

int32_t & IntervalNode::get_median()
{
    return median;
}

std::vector<Record> & IntervalNode::get_records()
{
    return records;
}

void IntervalNode::set_is_leaf()
{
    this->isLeaf = true;
}

void IntervalNode::set_median(int32_t m)
{
    this->median = m;
}

void IntervalNode::print(int32_t level)
{
    seqan3::debug_stream << "Level: " << level << '\n' << "Reads: ";
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
    for (auto & r : records_i)
    {
        if (node->get_median() < r.start) // Median is to the left of the read, so read is in right subtree.
        {
            rRecords.push_back(r);
        }
        else if (node->get_median() >= r.start && node->get_median() <= r.end) // Median is within the read.
        {
            node->get_records().push_back(r);
        }
        else if (node->get_median() > r.end) // Median is to the right of the read, so read is in left subtree.
        {
            lRecords.push_back(r);
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

std::vector<Record> overlap(std::unique_ptr<IntervalNode> & root, int32_t start, int32_t end)
{
    std::vector<Record> results{};
    return results;
}
