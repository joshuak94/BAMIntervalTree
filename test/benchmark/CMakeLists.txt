cmake_minimum_required (VERSION 3.8)

add_benchmark_test (benchmark_api.cpp)
target_use_datasources (construction_test FILES large_file.bam)
