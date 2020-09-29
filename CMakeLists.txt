cmake_minimum_required (VERSION 3.14)

project (makegrid C Fortran)

find_package (OpenMP)

add_executable (mgrid ${mgrid_sources})
add_dependencies (mgrid stell)

target_link_libraries (mgrid stell)
target_include_directories (mgrid PUBLIC $<TARGET_PROPERTY:stell,BINARY_DIR>)

add_subdirectory (Sources)

if (OpenMP_Found)
    target_link_libraries (mgrid OpenMP::OpenMP)
endif