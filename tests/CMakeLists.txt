apply_include_paths("${EATK_INTERNAL_INCLUDES}")
apply_include_paths("${EATK_EXTERNAL_INCLUDES}")

add_executable(mpitest mpitest.cpp)
target_link_libraries(mpitest eatk-static "${EATK_LINK_LIBS}")

add_executable(threadtest threadtest.cpp)
target_link_libraries(threadtest eatk-static "${EATK_LINK_LIBS}")

add_executable(gatest gatest.cpp )
target_link_libraries(gatest eatk-static "${EATK_LINK_LIBS}")

add_executable(ndset ndset.cpp)
target_link_libraries(ndset eatk-static "${EATK_LINK_LIBS}")
