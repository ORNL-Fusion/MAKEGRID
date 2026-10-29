cmake_minimum_required (VERSION 3.14)

project (makegrid C Fortran)

find_package (OpenMP)

add_executable (mgrid)
add_dependencies (mgrid pull_makegrid)

target_link_libraries (mgrid
                       
                       PUBLIC
                       stell
                       $<$<BOOL:${OpenMP_Fortran_FOUND}>:OpenMP::OpenMP_Fortran>
)

add_subdirectory (Sources)

install (TARGETS mgrid)

################################################################################
#  Testing                                                                     #
################################################################################

#  Build test utilities.
add_executable (xcheck_fields)
add_dependencies (xcheck_fields pull_makegrid)

target_link_libraries (xcheck_fields PUBLIC stell)

add_subdirectory (Testing)
