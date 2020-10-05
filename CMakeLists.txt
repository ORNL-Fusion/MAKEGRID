cmake_minimum_required (VERSION 3.14)

project (makegrid C Fortran)

find_package (OpenMP)

add_executable (mgrid)
add_dependencies (mgrid stell)

target_link_libraries (mgrid stell)
target_include_directories (mgrid PUBLIC $<TARGET_PROPERTY:stell,BINARY_DIR>)

add_subdirectory (Sources)

if (OpenMP_Fortran_FOUND)
    target_link_libraries (mgrid OpenMP::OpenMP_Fortran)
endif ()

################################################################################
#  Testing                                                                     #
################################################################################

#  Build test utilities.
add_executable (xcheck_fields)
add_dependencies (xcheck_fields stell)

target_link_libraries (xcheck_fields stell)
target_include_directories (xcheck_fields PUBLIC $<TARGET_PROPERTY:stell,BINARY_DIR>)

add_subdirectory (Testing)
