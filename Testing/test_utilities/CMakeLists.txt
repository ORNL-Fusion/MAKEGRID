target_sources(xcheck_fields
               PRIVATE
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/check_fields.f>
)
