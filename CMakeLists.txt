cmake_minimum_required (VERSION 3.14)

project (makegrid C Fortran)

add_executable (mgrid ${mgrid_sources})
add_dependencies (mgrid stell)

target_link_libraries (mgrid stell)
