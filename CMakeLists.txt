cmake_minimum_required(VERSION 2.8.0)
project (GRASP2)
set(CMAKE_CXX_FLAGS "-fopenmp -O3 -std=c++0x")
# include directory
include_directories(src)
# handling boost library linking
find_package(Boost 1.46.0 COMPONENTS iostreams filesystem system program_options)
if(NOT Boost_FOUND)
  message(FATAL_ERROR " Cannot find Boost C++ Libraries version 1.46 or higher, abort.")
elseif (Boost_FOUND)
  include_directories(${include_directories} ${Boost_INCLUDE_DIRS})
  link_directories(${link_directories} ${Boost_LIB_DIRS})
  add_executable(bin/grasp2-build ./src/main_build.cc ./src/kmer_filtering.cc ./src/sequence_search.cc ./src/align_batch.cc ./src/string_graph.cc ./src/bwt.cc ./src/bio_alphabet.cc ./src/loader.cc ./src/concatenator.cc ./src/contig_refinement.cc ./src/bwt_search.cc ./src/kmer_unitcoder.cc ./src/minimizer_sort.cc ./src/divsufsort.cc ./src/sssort.cc ./src/trsort.cc ./src/utils.cc)
  target_link_libraries (bin/grasp2-build ${Boost_LIBRARIES})
  add_executable(bin/grasp2-assemble ./src/main_assemble.cc ./src/kmer_filtering.cc ./src/sequence_search.cc ./src/align_batch.cc ./src/string_graph.cc ./src/bwt.cc ./src/bio_alphabet.cc ./src/loader.cc ./src/concatenator.cc ./src/contig_refinement.cc ./src/bwt_search.cc ./src/kmer_unitcoder.cc ./src/minimizer_sort.cc ./src/divsufsort.cc ./src/sssort.cc ./src/trsort.cc ./src/utils.cc)
  target_link_libraries (bin/grasp2-assemble ${Boost_LIBRARIES})
endif (NOT Boost_FOUND)

