apply_include_paths("${MOGAL2_INTERNAL_INCLUDES}")
apply_include_paths("${MOGAL2_EXTERNAL_INCLUDES}")

add_executable(travelingsalesman travelingsalesman.cpp)
target_link_libraries(travelingsalesman mogal2-static "${MOGAL2_LINK_LIBS}")

