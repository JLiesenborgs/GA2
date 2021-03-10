cmake_minimum_required(VERSION 3.0)
project(mogal2)
set(VERSION 0.0.1)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
include(CheckIncludeFile)
include(CheckIncludeFileCXX)
include(${PROJECT_SOURCE_DIR}/cmake/Macros.cmake)
include(CheckCXXCompilerFlag)
include(CheckCXXSourceCompiles)
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
if(COMPILER_SUPPORTS_CXX14)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
else()
	# Test if the compiler already supports C++14
	check_cxx_source_compiles("int main(void)\n{\n\tauto f = [](auto x)\n\t{\n\t\treturn x*x;\n\t};\n\n\treturn 0;\n}\n" C14AVAIL)
	if (NOT C14AVAIL)
		message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++14 support.")
	endif()
endif()

get_install_directory(LIBRARY_INSTALL_DIR)

set(MOGAL2_EXTERNAL_INCLUDES "")
set(MOGAL2_INTERNAL_INCLUDES "")
set(MOGAL2_LINK_LIBS "")

add_additional_stuff(MOGAL2_EXTERNAL_INCLUDES MOGAL2_LINK_LIBS)
get_install_directory(LIBRARY_INSTALL_DIR)

find_package(ErrUt REQUIRED)
list(APPEND MOGAL2_EXTERNAL_INCLUDES "${ERRUT_INCLUDE_DIRS}")
list(APPEND MOGAL2_LINK_LIBS "${ERRUT_LIBRARIES}")

find_package(MPI REQUIRED) # TODO: for now it's required
list(APPEND MOGAL2_LINK_LIBS "${MPI_LIBRARIES}")
list(APPEND MOGAL2_EXTERNAL_INCLUDES "${MPI_INCLUDE_PATH}")

find_package(Threads REQUIRED)
list(APPEND MOGAL2_LINK_LIBS "${CMAKE_THREAD_LIBS_INIT}")

list(APPEND MOGAL2_INTERNAL_INCLUDES "${PROJECT_BINARY_DIR}/src/")
apply_include_paths("${MOGAL_INTERNAL_INCLUDES}")
apply_include_paths("${MOGAL_EXTERNAL_INCLUDES}")

add_executable(mpitest mpitest.cpp vectorgenomefitness.cpp
						mpipopulationfitnesscalculation.cpp
						singlethreadedpopulationfitnesscalculation.cpp
						mpieventdistributor.cpp)
target_link_libraries(mpitest "${MOGAL2_LINK_LIBS}")

add_executable(threadtest threadtest.cpp vectorgenomefitness.cpp
						  multithreadedpopulationfitnesscalculation.cpp)
target_link_libraries(threadtest "${MOGAL2_LINK_LIBS}")

add_executable(gatest gatest.cpp vectorgenomefitness.cpp
			   singlethreadedpopulationfitnesscalculation.cpp
			   multithreadedpopulationfitnesscalculation.cpp
			   mpipopulationfitnesscalculation.cpp
			   mersennerandomnumbergenerator.cpp
			   simplesortedpopulation.cpp
			   rankparentselection.cpp
			   singlethreadedpopulationmutation.cpp
			   singlethreadedpopulationcrossover.cpp
			   mpieventdistributor.cpp
			   singlebestelitism.cpp
			   valuefitness.cpp
			   geneticalgorithm.cpp
			   )
target_link_libraries(gatest "${MOGAL2_LINK_LIBS}")
