#  Copy input file from the source to the build directory.
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/coils.test
                ${CMAKE_CURRENT_BINARY_DIR}/coils.test
                COPYONLY)

add_test (NAME    mgrid_test
          COMMAND $<TARGET_PROPERTY:mgrid,BINARY_DIR>/mgrid coils.test)
add_test (NAME    mgrid_check_fields_test
          COMMAND $<TARGET_PROPERTY:xcheck_fields,BINARY_DIR>/xcheck_fields)
set_tests_properties (mgrid_check_fields_test
                      PROPERTIES
                      DEPENDS mgrid_test)
