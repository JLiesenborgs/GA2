add_definitions(-DEATK_COMPILING)

set(HEADERS
	basicnondominatedsetcreator.h
	crossovermutation.h
	duplicateindividualremoval.h
	fitnessbasedduplicateremoval.h
	evolutionaryalgorithm.h
	fasternondominatedsetcreator.h
	genomefitness.h
	mersennerandomnumbergenerator.h
	mpieventdistributor.h
	mpipopulationfitnesscalculation.h
	multithreadedpopulationfitnesscalculation.h
	ndsortedpopulation.h
	ndsortedpopulationinterface.h
	nondominatedsetcreator.h
	permutationordercrossover.h
	permutationswapmutation.h
	population.h
	randomnumbergenerator.h
	rankparentselection.h
	remainingtargetpopulationsizeiteration.h
	simplesortedpopulation.h
	singlebestelitism.h
	singlepopulationcrossover.h
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
	differentialevolutionevolver.h
	vectordifferentialevolution.h
	jadeevolver.h
	individual.h
	calculation.h
	populationevolver.h
	selection.h
	elitism.h
	"${PROJECT_BINARY_DIR}/src/eatkconfig.h" )

set(SOURCES
	basicnondominatedsetcreator.cpp
	fitnessbasedduplicateremoval.cpp
	evolutionaryalgorithm.cpp
	fasternondominatedsetcreator.cpp
	mersennerandomnumbergenerator.cpp
	mpieventdistributor.cpp
	mpipopulationfitnesscalculation.cpp
	multithreadedpopulationfitnesscalculation.cpp
	ndsortedpopulation.cpp
	permutationordercrossover.cpp
	permutationswapmutation.cpp
	rankparentselection.cpp
	simplesortedpopulation.cpp
	singlebestelitism.cpp
	singlepopulationcrossover.cpp
	singlethreadedpopulationfitnesscalculation.cpp
	tournamentparentselection.cpp
	trackbestonlyselectionpopulation.cpp
	valuefitness.cpp
	vectorgenomefitness.cpp
	differentialevolutionevolver.cpp
	jadeevolver.cpp
	)

if (MSVC)
	set(CMAKE_DEBUG_POSTFIX _d)
endif (MSVC)

if (NOT MSVC OR EATK_COMPILE_STATIC)
	add_library(eatk-static STATIC ${SOURCES} ${HEADERS})
	set_target_properties(eatk-static PROPERTIES OUTPUT_NAME eatk)
	set_target_properties(eatk-static PROPERTIES CLEAN_DIRECT_OUTPUT 1)
	set(EATK_INSTALLTARGETS eatk-static)
	target_link_libraries(eatk-static ${EATK_LINK_LIBS})
endif (NOT MSVC OR EATK_COMPILE_STATIC)

if (NOT MSVC OR NOT EATK_COMPILE_STATIC)
	add_library(eatk-shared SHARED ${SOURCES} ${HEADERS})
	set_target_properties(eatk-shared PROPERTIES OUTPUT_NAME eatk)
	set_target_properties(eatk-shared PROPERTIES VERSION ${VERSION})
	set_target_properties(eatk-shared PROPERTIES CLEAN_DIRECT_OUTPUT 1)
	set(EATK_INSTALLTARGETS ${EATK_INSTALLTARGETS} eatk-shared)
	target_link_libraries(eatk-shared ${EATK_LINK_LIBS})
endif (NOT MSVC OR NOT EATK_COMPILE_STATIC)

apply_include_paths("${EATK_INTERNAL_INCLUDES}")
apply_include_paths("${EATK_EXTERNAL_INCLUDES}")

install(FILES ${HEADERS} DESTINATION include/eatk)
install(TARGETS ${EATK_INSTALLTARGETS} DESTINATION ${LIBRARY_INSTALL_DIR})

