cmake_minimum_required (VERSION 3.8)

add_api_test (construction_test.cpp)
target_use_datasources (construction_test FILES simulated_chr1_small_golden.bam)
target_use_datasources (construction_test FILES simulated_mult_chr_small_golden.bam)

add_api_test (search_test.cpp)
target_use_datasources (search_test FILES samtools_result.sam)

add_api_test (write_read_test.cpp)

add_api_test (sample_functions_test.cpp)
