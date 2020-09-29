target_sources(mgrid
               PRIVATE
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/makegrid_global.f>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/makegrid.f>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/symcheck.f>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/write_mgrid.f>
)
