cmake_minimum_required (VERSION 3.8)

add_api_test (construction_test.cpp)
target_use_datasources (construction_test FILES simulated_chr1_small_golden.bam)
