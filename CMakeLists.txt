project(invop)
cmake_minimum_required(VERSION 3.6)

# Compiler
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
#set(CMAKE_CXX_COMPILER /usr/bin/g++)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wpedantic -Wall -Werror -fmax-errors=1 -O2")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
#set(CMAKE_C_COMPILER /usr/bin/gcc)


# General directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


# Custom Libraries directories
#set(JSON_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_libraries/json)
#set(JSON_INCLUDE_DIR ${JSON_ROOT_DIR}/src)

#set(LIB_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib/include)


# CPLEX Libraries directories
set(CPLEX_STUDIO_DIR /opt/ibm/ILOG/CPLEX_Studio128)
set(CPLEX_DIR ${CPLEX_STUDIO_DIR}/cplex)
set(CPLEX_CONCERT_DIR ${CPLEX_STUDIO_DIR}/concert)

set(CPLEX_INCLUDE_DIR ${CPLEX_DIR}/include)
find_library(CPLEX_LIBRARY NAMES cplex PATHS ${CPLEX_DIR} PATH_SUFFIXES lib/x86-64_linux/static_pic)
set(CPLEX_LIBRARY_DEBUG ${CPLEX_LIBRARY})
add_library(cplex STATIC IMPORTED GLOBAL)
find_package(Threads)
set_target_properties(cplex PROPERTIES
											IMPORTED_LOCATION "${CPLEX_LIBRARY}"
											IMPORTED_LOCATION_DEBUG "${CPLEX_LIBRARY_DEBUG}"
											INTERFACE_INCLUDE_DIRECTORIES "${CPLEX_INCLUDE_DIR}"
											INTERFACE_INCLUDE_LIBRARIES "${CMAKE_THREAD_LIBS_INIT};m")

set(CPLEX_CONCERT_INCLUDE_DIR ${CPLEX_CONCERT_DIR}/include)
find_library(CPLEX_CONCERT_LIBRARY NAMES concert PATHS ${CPLEX_CONCERT_DIR} PATH_SUFFIXES lib/x86-64_linux/static_pic)
set(CPLEX_CONCERT_LIBRARY_DEBUG ${CPLEX_CONCERT_LIBRARY})
add_library(concert STATIC IMPORTED GLOBAL)
set_target_properties(concert PROPERTIES
											IMPORTED_LOCATION "${CPLEX_CONCERT_LIBRARY}"
											IMPORTED_LOCATION_DEBUG "${CPLEX_CONCERT_LIBRARY_DEBUG}"
											INTERFACE_COMPILE_DEFINITIONS IL_STD # Require standard compliance.
											INTERFACE_INCLUDE_DIRECTORIES "${CPLEX_CONCERT_INCLUDE_DIR}"
											INTERFACE_LINK_LIBRARIES "${CMAKE_THREAD_LIBS_INIT}")

set(CPLEX_ILOCPLEX_INCLUDE_DIR ${CPLEX_INCLUDE_DIR})
find_library(CPLEX_ILOCPLEX_LIBRARY NAMES ilocplex PATHS ${CPLEX_DIR} PATH_SUFFIXES lib/x86-64_linux/static_pic)
set(CPLEX_ILOCPLEX_LIBRARY_DEBUG ${CPLEX_ILOCPLEX_LIBRARY})
add_library(ilocplex STATIC IMPORTED GLOBAL)
set_target_properties(ilocplex PROPERTIES
											IMPORTED_LOCATION "${CPLEX_ILOCPLEX_LIBRARY}"
											IMPORTED_LOCATION_DEBUG "${CPLEX_ILOCPLEX_LIBRARY_DEBUG}"
											INTERFACE_INCLUDE_DIRECTORIES "${CPLEX_ILOCPLEX_INCLUDE_DIR}"
											INTERFACE_LINK_LIBRARIES "concert;cplex")


#set(CPLEX_CXX_FLAGS "-m64 -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD")
#set(CPLEX_CXX_FLAGS "-m64 -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG")
set(CPLEX_CXX_FLAGS "-m64 -fPIC -fno-strict-aliasing -fexceptions")
set(CPLEX_ILOCPLEX_CXX_FLAGS "-Wno-long-long")


# Releases
if (CMAKE_BUILD_TYPE MATCHES Debug)
	add_definitions(-DDEBUG=1)
endif()

if (CMAKE_BUILD_TYPE MATCHES Release)
	set(CPLEX_CXX_FLAGS "${CPLEX_CXX_FLAGS} -DNDEBUG")
endif()

# Projects
#add_subdirectory(external_libraries/json) #nlohmann Json library
#add_subdirectory(lib) #internal libraries
add_subdirectory(CVRP)
add_subdirectory(Food_Manufacture)
