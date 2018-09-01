project(invop-1a)

# Compiler
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CPLEX_CXX_FLAGS} ${CPLEX_ILOCPLEX_CXX_FLAGS}")

# Directories, files, includes, output
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${INCLUDE_DIR}
										#${LIB_INCLUDE_DIR}
										#${JSON_INCLUDE_DIR}
										${CPLEX_INCLUDE_DIR}
										${CPLEX_CONCERT_INCLUDE_DIR}
										${CPLEX_ILOCPLEX_INCLUDE_DIR})
link_directories(${CPLEX_LIBRARY}
								 ${CPLEX_CONCERT_LIBRARY}
								 ${CPLEX_ILOCPLEX_LIBRARY})

# Targets
set(SOURCE_FILES main.cpp)
add_executable(1a ${SOURCE_FILES})

# Because CLion can't detect .hpp who includes a .cpp (mimic for
# separated declaration and implementation for template classes),
# I put implementation in .hpp and link it as a library
add_library(1alp INTERFACE)
target_sources(1alp INTERFACE
							 "${INCLUDE_DIR}/problem.hpp"
							 "${CMAKE_CURRENT_SOURCE_DIR}/problem.cpp")
target_link_libraries(1a
											1alp
											${CMAKE_THREAD_LIBS_INIT}
											${CMAKE_DL_LIBS}
											c
											m
											cplex
											concert
											ilocplex)

# Testing
#add_subdirectory(tests)
