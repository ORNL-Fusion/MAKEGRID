add_test (NAME mgrid_test
          COMMAND $<TARGET_PROPERTY:mgrid,BINARY_DIR>mgrid coils.test)

