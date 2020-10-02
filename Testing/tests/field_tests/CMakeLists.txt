#  Copy input file from the source to the build directory.
add_test (NAME    mgrid_copy_test_file
          COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/coils.test ${CMAKE_CURRENT_BINARY_DIR}/)
add_test (NAME    mgrid_test
          COMMAND $<TARGET_PROPERTY:mgrid,BINARY_DIR>/mgrid coils.test)

