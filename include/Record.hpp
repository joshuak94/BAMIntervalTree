struct Record
{
    int32_t start, end;

    Record(int32_t start_i, int32_t end_i) : start{std::move(start_i)}, end{std::move(end_i)} {}
};
