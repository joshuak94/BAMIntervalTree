# BAMIntervalTree

The goal of this project is to implement a C++ header library for indexing a BAM file into an interval tree.

**Note: This project is very much a work in progress, and functionality may change between commits.**

There are two ways to use the BAMIntervalTree:
    1. You can make the command line tool and index a BAM file and then perform search queries.
    2. You can use the included header library in your own C++ code.

# Installation and examples

BAMIntervalTree depends on seqan3. It can be downloaded as such:
```sh
git clone --recurse-submodules git@github.com:joshuak94/BAMIntervalTree.git
```
Once downloaded, the command line tool can be built and run as such:
```sh
mkdir BAMIntervalTree_build && cd BAMIntervalTree_build
cmake ../BAMIntervalTree/
make
# Create an index for a given BAM file.
bin/BAMIntervalTree index BAMFILE
# Search for the reads which overlap the region chr1:200-300.
bin/BAMIntervalTree overlap -s chr1,200 -e chr1,300 BAMFILE
# View help pages
bin/BAMIntervalTree -h
bin/BAMIntervalTree index -h
bin/BAMIntervalTree overlap -h
```
Furthermore, the header library is available in `BAMIntervalTree/include/bamit/all.hpp` and once included, can be
used like so:
```cpp
#include <bamit/all.hpp>

#include <filesystem> // For path.
#include <vector>     // For vector.
#include <memory>     // For unique_ptr.

void main()
{
    // Parse the BAM file input.
    std::filesystem::path input{"/PATH/TO/BAMFILE/"};
    std::filesystem::path output{"/PATH/TO/OUTPUT/BAMFILE"};
    bamit::sam_file_input_type input_file{input};

    // Construct the tree over the parsed records.
    std::vector<std::unique_ptr<bamit::IntervalNode>> node_list = bamit::index(input_file);

    // Search for a region. Position is just an alias for std::tuple<int32_t, int32_t> (chromosome, position)
    bamit::Position start{0, 200}, end{0, 300};
    bamit::get_overlap_records(input, root, start, end, output);
}
```
