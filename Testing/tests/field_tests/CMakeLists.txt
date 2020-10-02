#  Copy input file from the source to the build directory.
file (COPY        ${CMAKE_CURRENT_SOURCE_DIR}/coils.test
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test (NAME    mgrid_test
          COMMAND $<TARGET_PROPERTY:mgrid,BINARY_DIR>/mgrid coils.test)
add_test (NAME    mgrid_check_fields_test
          COMMAND $<TARGET_PROPERTY:xcheck_fields,BINARY_DIR>/xcheck_fields)
