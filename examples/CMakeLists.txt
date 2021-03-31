apply_include_paths("${EATK_INTERNAL_INCLUDES}")
apply_include_paths("${EATK_EXTERNAL_INCLUDES}")

add_executable(travelingsalesman travelingsalesman.cpp)
target_link_libraries(travelingsalesman eatk-static "${EATK_LINK_LIBS}")

add_executable(knapsack knapsack.cpp)
target_link_libraries(knapsack eatk-static "${EATK_LINK_LIBS}")

