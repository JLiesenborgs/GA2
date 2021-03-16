apply_include_paths("${MOGAL2_INTERNAL_INCLUDES}")
apply_include_paths("${MOGAL2_EXTERNAL_INCLUDES}")

add_executable(mpitest mpitest.cpp)
target_link_libraries(mpitest mogal2-static "${MOGAL2_LINK_LIBS}")

add_executable(threadtest threadtest.cpp)
target_link_libraries(threadtest mogal2-static "${MOGAL2_LINK_LIBS}")

add_executable(gatest gatest.cpp )
target_link_libraries(gatest mogal2-static "${MOGAL2_LINK_LIBS}")

add_executable(ndset ndset.cpp)
target_link_libraries(ndset mogal2-static "${MOGAL2_LINK_LIBS}")
