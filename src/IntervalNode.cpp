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

std::vector<std::string> & IntervalNode::get_records()
{
    return records;
}

IntervalNode & IntervalNode::construct_tree(sam_file_input_t & input_file)
{
    size_t length = std::get<0>(input_file.header().ref_id_info[0]);
    seqan3::debug_stream << length << std::endl;
    return *this;
}
