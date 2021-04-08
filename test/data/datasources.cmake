cmake_minimum_required (VERSION 3.8)

include (cmake/app_datasources.cmake)

# copies file to <build>/data/
declare_datasource (FILE simulated_chr1_small_golden.bam
                    URL ${CMAKE_SOURCE_DIR}/test/data/simulated_chr1_small_golden.bam
                    URL_HASH SHA256=c580e7caf6baf0944028a0ba657d7095b8dcb77866bbd9b60aff64b9f43a4f54)

declare_datasource (FILE simulated_chr1_small_golden.bam.bai
                    URL ${CMAKE_SOURCE_DIR}/test/data/simulated_chr1_small_golden.bam.bai
                    URL_HASH SHA256=63a6a8e1bc3849215058e8af70fd74b6764f73ca11053bd4e65fdde9ef77792e)

declare_datasource (FILE simulated_mult_chr_small_golden.bam
                    URL ${CMAKE_SOURCE_DIR}/test/data/simulated_mult_chr_small_golden.bam
                    URL_HASH SHA256=906cf75a80754d08270780c4722e24097938d4a90bdc1863e5cae78cbc1d457e)

declare_datasource (FILE simulated_mult_chr_small_golden.bam.bai
                    URL ${CMAKE_SOURCE_DIR}/test/data/simulated_mult_chr_small_golden.bam.bai
                    URL_HASH SHA256=77261a82ea19a78518d37039dec473131b76e8f2a692ea2a25e1d9b13265459b)
declare_datasource (FILE samtools_result.sam
                    URL ${CMAKE_SOURCE_DIR}/test/data/samtools_result.sam
                    URL_HASH SHA256=1d5642a8462687e58a65117ed31a9331e0c494b2bb5f97afe4f256b97db3ae71)
