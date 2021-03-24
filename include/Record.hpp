struct Record
{
    int32_t start, end;

    Record(int32_t start_i, int32_t end_i) : start{std::move(start_i)}, end{std::move(end_i)} {}
};

struct RecordComparatorStart
{
    // Compare 2 Player objects using name
    bool operator ()(const Record & record1, const Record & record2)
    {
        if(record1.start == record2.start)
            return record1.end < record2.end;
        return record1.start < record2.start;
    }
};

struct RecordComparatorEnd
{
    // Compare 2 Player objects using name
    bool operator ()(const Record & record1, const Record & record2)
    {
        if(record1.end== record2.end)
            return record1.start > record2.start;
        return record1.end > record2.end;
    }
};
