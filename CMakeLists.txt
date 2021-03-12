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
apply_include_paths("${MOGAL2_INTERNAL_INCLUDES}")
apply_include_paths("${MOGAL2_EXTERNAL_INCLUDES}")

if (MSVC)
	set(MOGAL2_COMPILE_STATIC ON CACHE BOOL "Flag indicating if a static library should be built, or a dynamic one")
endif (MSVC)
if (NOT MSVC OR MOGAL2_COMPILE_STATIC)
	set(MOGAL2_IMPORT "")
	set(MOGAL2_EXPORT "")
else (NOT MSVC OR MOGAL2_COMPILE_STATIC)
	set(MOGAL2_IMPORT "__declspec(dllimport)")
	set(MOGAL2_EXPORT "__declspec(dllexport)")
endif (NOT MSVC OR MOGAL2_COMPILE_STATIC)

configure_file("${PROJECT_SOURCE_DIR}/src/mogal2config.h.in" "${PROJECT_BINARY_DIR}/src/mogal2config.h")

list(APPEND MOGAL2_INTERNAL_INCLUDES "${PROJECT_SOURCE_DIR}/src/")
list(APPEND MOGAL2_INTERNAL_INCLUDES "${PROJECT_BINARY_DIR}/src/")

add_subdirectory(src)
#add_subdirectory(programs)
add_subdirectory(tests)
add_subdirectory(examples)

if (NOT MSVC)
	set(MOGAL2_LIBS "-L${LIBRARY_INSTALL_DIR}" "-lmogal2")
else (NOT MSVC)
	set(MOGAL2_LIBS optimized "${LIBRARY_INSTALL_DIR}/mogal2.lib" 
		        debug "${LIBRARY_INSTALL_DIR}/mogal2_d.lib")
endif (NOT MSVC)

set(MOGAL2_INCDIRS "${MOGAL2_EXTERNAL_INCLUDES};${CMAKE_INSTALL_PREFIX}/include")
set(MOGAL2_LIBS "${MOGAL2_LIBS};${MOGAL2_LINK_LIBS}")
remove_empty(MOGAL2_INCDIRS)
list(REMOVE_DUPLICATES MOGAL2_INCDIRS)
remove_empty(MOGAL2_LIBS)

foreach(ARG ${MOGAL2_LIBS})
	set(MOGAL2_LIBS_CMAKECONFIG "${MOGAL2_LIBS_CMAKECONFIG} \"${ARG}\"")
endforeach()
foreach(ARG ${MOGAL2_INCDIRS})
	set(MOGAL2_INCDIRS_CMAKECONFIG "${MOGAL2_INCDIRS_CMAKECONFIG} \"${ARG}\"")
endforeach()

configure_file("${PROJECT_SOURCE_DIR}/cmake/MOGAL2Config.cmake.in" "${PROJECT_BINARY_DIR}/cmake/MOGAL2Config.cmake")
install(FILES "${PROJECT_BINARY_DIR}/cmake/MOGAL2Config.cmake" 
	DESTINATION ${LIBRARY_INSTALL_DIR}/cmake/MOGAL2)


