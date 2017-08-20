cmake_minimum_required(VERSION 3.8)

project( Matsu )


find_package( OpenMP )

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )

add_subdirectory( nr )
include_directories( nr )

add_library( numfunction numfunctionsplinelog.cc )

add_executable( MatsuRevised theAdlerDavisMatsuOneBetaRevised.cc )
target_link_libraries( MatsuRevised numfunction nr )

add_executable( MatsuVastag theADMatsuVastagnnotRatio.cc )
target_link_libraries( MatsuVastag numfunction nr )

add_executable( Matsu2017_non_extrapolated Matsu2017_non_extrapolated_frequencies.cc )
target_link_libraries( Matsu2017_non_extrapolated numfunction nr )

add_executable( Matsu2017_ratio Matsu2017_ratio.cc )
target_link_libraries( Matsu2017_ratio numfunction nr )

