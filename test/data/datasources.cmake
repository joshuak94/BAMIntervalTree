cmake_minimum_required (VERSION 3.8)

include (cmake/app_datasources.cmake)

# copies file to <build>/data/
declare_datasource (FILE simulated_chr1_small_golden.bam
                    URL ${CMAKE_SOURCE_DIR}/test/data/simulated_chr1_small_golden.bam
                    URL_HASH SHA256=c580e7caf6baf0944028a0ba657d7095b8dcb77866bbd9b60aff64b9f43a4f54)
