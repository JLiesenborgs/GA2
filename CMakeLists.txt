cmake_minimum_required(VERSION 3.5)
project(eatk)
set(VERSION 0.3.0)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
include(CheckIncludeFile)
include(CheckIncludeFileCXX)
include(${PROJECT_SOURCE_DIR}/cmake/Macros.cmake)
include(CheckCXXCompilerFlag)
include(CheckCXXSourceCompiles)
CHECK_CXX_COMPILER_FLAG("-std=c++17" COMPILER_SUPPORTS_CXX17)
if(COMPILER_SUPPORTS_CXX17)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
else()
	CHECK_CXX_COMPILER_FLAG("/std:c++17" COMPILER_SUPPORTS_CXX17_VS) # For visual studio
	if(COMPILER_SUPPORTS_CXX17_VS)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++17")
	else()
		# Test if the compiler already supports C++17
		check_cxx_source_compiles("#include <utility>\n\nint main(void)\n{\n\tauto f = [](auto x)\n\t{\n\t\treturn x*x;\n\t};\n\tauto [ a, b ] = std::pair<double, int>(1.0,1);\n\treturn 0;\n}\n" C17AVAIL)
		if (NOT C17AVAIL)
			message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++17 support.")
		endif()
	endif()
endif()

get_install_directory(LIBRARY_INSTALL_DIR)

set(EATK_EXTERNAL_INCLUDES "")
set(EATK_INTERNAL_INCLUDES "")
set(EATK_LINK_LIBS "")

add_additional_stuff(EATK_EXTERNAL_INCLUDES EATK_LINK_LIBS)
get_install_directory(LIBRARY_INSTALL_DIR)

find_package(ErrUt REQUIRED)
list(APPEND EATK_EXTERNAL_INCLUDES "${ERRUT_INCLUDE_DIRS}")
list(APPEND EATK_LINK_LIBS "${ERRUT_LIBRARIES}")

find_package(MPI)
if (MPI_FOUND)
	set(TMP TRUE)
else (MPI_FOUND)
	set(TMP FALSE)
endif (MPI_FOUND)
option(EATK_SUPPORT_MPI "Enable MPI support" ${TMP})
if (EATK_SUPPORT_MPI)
	if (MPI_FOUND)
		list(APPEND EATK_LINK_LIBS "${MPI_LIBRARIES}")
		list(APPEND EATK_EXTERNAL_INCLUDES "${MPI_INCLUDE_PATH}")
	endif (MPI_FOUND)
	set(EATKCONFIG_MPISUPPORT "#define EATKCONFIG_MPISUPPORT")
else (EATK_SUPPORT_MPI)
	set(EATKCONFIG_MPISUPPORT "// No MPI support")
endif (EATK_SUPPORT_MPI)

find_package(Threads REQUIRED)
list(APPEND EATK_LINK_LIBS "${CMAKE_THREAD_LIBS_INIT}")

list(APPEND EATK_INTERNAL_INCLUDES "${PROJECT_BINARY_DIR}/src/")
apply_include_paths("${EATK_INTERNAL_INCLUDES}")
apply_include_paths("${EATK_EXTERNAL_INCLUDES}")

if (MSVC)
	set(EATK_COMPILE_STATIC ON CACHE BOOL "Flag indicating if a static library should be built, or a dynamic one")
endif (MSVC)
if (NOT MSVC OR EATK_COMPILE_STATIC)
	set(EATK_IMPORT "")
	set(EATK_EXPORT "")
else (NOT MSVC OR EATK_COMPILE_STATIC)
	set(EATK_IMPORT "__declspec(dllimport)")
	set(EATK_EXPORT "__declspec(dllexport)")
endif (NOT MSVC OR EATK_COMPILE_STATIC)

configure_file("${PROJECT_SOURCE_DIR}/src/eatkconfig.h.in" "${PROJECT_BINARY_DIR}/src/eatkconfig.h")

list(APPEND EATK_INTERNAL_INCLUDES "${PROJECT_SOURCE_DIR}/src/")
list(APPEND EATK_INTERNAL_INCLUDES "${PROJECT_BINARY_DIR}/src/")

add_subdirectory(src)
#add_subdirectory(programs)
add_subdirectory(tests)
add_subdirectory(examples)

if (NOT MSVC)
	set(EATK_LIBS "-L${LIBRARY_INSTALL_DIR}" "-leatk")
else (NOT MSVC)
	set(EATK_LIBS optimized "${LIBRARY_INSTALL_DIR}/eatk.lib" 
		        debug "${LIBRARY_INSTALL_DIR}/eatk_d.lib")
endif (NOT MSVC)

set(EATK_INCDIRS "${EATK_EXTERNAL_INCLUDES};${CMAKE_INSTALL_PREFIX}/include")
set(EATK_LIBS "${EATK_LIBS};${EATK_LINK_LIBS}")
remove_empty(EATK_INCDIRS)
list(REMOVE_DUPLICATES EATK_INCDIRS)
remove_empty(EATK_LIBS)

foreach(ARG ${EATK_LIBS})
	set(EATK_LIBS_CMAKECONFIG "${EATK_LIBS_CMAKECONFIG} \"${ARG}\"")
endforeach()
foreach(ARG ${EATK_INCDIRS})
	set(EATK_INCDIRS_CMAKECONFIG "${EATK_INCDIRS_CMAKECONFIG} \"${ARG}\"")
endforeach()

configure_file("${PROJECT_SOURCE_DIR}/cmake/EATKConfig.cmake.in" "${PROJECT_BINARY_DIR}/cmake/EATKConfig.cmake")
install(FILES "${PROJECT_BINARY_DIR}/cmake/EATKConfig.cmake" 
	DESTINATION ${LIBRARY_INSTALL_DIR}/cmake/EATK)


