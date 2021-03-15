add_definitions(-DMOGAL2_COMPILING)

set(HEADERS
	crossovermutation.h
	geneticalgorithm.h
	genomefitness.h
	mersennerandomnumbergenerator.h
	mpieventdistributor.h
	mpipopulationfitnesscalculation.h
	multithreadedpopulationfitnesscalculation.h
	ndsortedpopulation.h
	permutationordercrossover.h
	permutationswapmutation.h
	population.h
	randomnumbergenerator.h
	rankparentselection.h
	remainingtargetpopulationsizeiteration.h
	simplesortedpopulation.h
	singlebestelitism.h
	singlethreadedpopulationcrossover.h
	singlethreadedpopulationfitnesscalculation.h
	stopcriterion.h
	tournamentparentselection.h
	trackbestonlyselectionpopulation.h
	uniformvectorgenomecrossover.h
	valuefitness.h
	valuevector.h
	vectorgenomefitness.h
	vectorgenomeflipmutation.h
	vectorgenomeuniformmutation.h
	"${PROJECT_BINARY_DIR}/src/mogal2config.h" )

set(SOURCES
	geneticalgorithm.cpp
	mersennerandomnumbergenerator.cpp
	mpieventdistributor.cpp
	mpipopulationfitnesscalculation.cpp
	multithreadedpopulationfitnesscalculation.cpp
	permutationordercrossover.cpp
	permutationswapmutation.cpp
	rankparentselection.cpp
	simplesortedpopulation.cpp
	singlebestelitism.cpp
	singlethreadedpopulationcrossover.cpp
	singlethreadedpopulationfitnesscalculation.cpp
	tournamentparentselection.cpp
	trackbestonlyselectionpopulation.cpp
	valuefitness.cpp
	vectorgenomefitness.cpp)

if (MSVC)
	set(CMAKE_DEBUG_POSTFIX _d)
endif (MSVC)

if (NOT MSVC OR MOGAL2_COMPILE_STATIC)
	add_library(mogal2-static STATIC ${SOURCES} ${HEADERS})
	set_target_properties(mogal2-static PROPERTIES OUTPUT_NAME mogal2)
	set_target_properties(mogal2-static PROPERTIES CLEAN_DIRECT_OUTPUT 1)
	set(MOGAL2_INSTALLTARGETS mogal2-static)
	target_link_libraries(mogal2-static ${MOGAL2_LINK_LIBS})
endif (NOT MSVC OR MOGAL2_COMPILE_STATIC)

if (NOT MSVC OR NOT MOGAL2_COMPILE_STATIC)
	add_library(mogal2-shared SHARED ${SOURCES} ${HEADERS})
	set_target_properties(mogal2-shared PROPERTIES OUTPUT_NAME mogal2)
	set_target_properties(mogal2-shared PROPERTIES VERSION ${VERSION})
	set_target_properties(mogal2-shared PROPERTIES CLEAN_DIRECT_OUTPUT 1)
	set(MOGAL2_INSTALLTARGETS ${MOGAL2_INSTALLTARGETS} mogal2-shared)
	target_link_libraries(mogal2-shared ${MOGAL2_LINK_LIBS})
endif (NOT MSVC OR NOT MOGAL2_COMPILE_STATIC)

apply_include_paths("${MOGAL2_INTERNAL_INCLUDES}")
apply_include_paths("${MOGAL2_EXTERNAL_INCLUDES}")

install(FILES ${HEADERS} DESTINATION include/mogal2)
install(TARGETS ${MOGAL2_INSTALLTARGETS} DESTINATION ${LIBRARY_INSTALL_DIR})

