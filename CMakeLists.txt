cmake_minimum_required (VERSION 3.18)
project (unitig_verify CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
 
# add seqan3 to search path
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../seqan3/build_system")
 
# require seqan3 with a version between >=3.0.0 and <4.0.0
find_package (seqan3 3.0 REQUIRED)
 
# build app with seqan3
add_executable (unitig_verify verify_unitigs.cpp)
add_executable (unitig_canonize canonize_unitigs.cpp)
target_link_libraries (unitig_verify seqan3::seqan3)
target_link_libraries (unitig_canonize seqan3::seqan3)
