# memcheck.cmake
#
# based on https://github.com/stepcode/stepcode/blob/master/lcov.cmake
# modified by <christoph.roesmann@tu-dortmund.de>
#
# Command:
# `ctest -S memcheck.cmake`

# TODO: chachegrind and callgrind

if(NOT ${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  message(FATAL_ERROR "LCOV is Linux-only")
endif(NOT ${CMAKE_SYSTEM_NAME} STREQUAL "Linux")

set(CTEST_SOURCE_DIRECTORY .)
set(CTEST_BINARY_DIRECTORY ../build)
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
set(CTEST_PROJECT_NAME "CodeCoverage")
set(CTEST_MEMORYCHECK_COMMAND /usr/bin/valgrind)
#set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --leak-check=full" )
#set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE "../src/cmake/valgrind_suppress.txt" ) # relative path does not work here for now...
#set(CTEST_INITIAL_CACHE "
#SITE:STRING=${CTEST_SITE}
#BUILDNAME:STRING=${CTEST_BUILD_NAME}
#SC_ENABLE_TESTING:BOOL=ON
#SC_ENABLE_COVERAGE:BOOL=ON
#SC_BUILD_SCHEMAS:STRING=ALL
#SC_BUILD_TYPE:STRING=Debug
#")

message("val: ${ValgrindCommand}")

ctest_start(memcheck)
#ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})
#message("configuring...")
#ctest_configure(BUILD "${CTEST_BINARY_DIRECTORY}" OPTIONS "${CONFIGURE_OPTIONS}")

#message("building...")
#ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" CONFIGURATION Debug)

#message("running tests...")
#ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}" PARALLEL_LEVEL 1)

message("running memcheck...")
ctest_memcheck(BUILD "${CTEST_BINARY_DIRECTORY}" PARALLEL_LEVEL 1)

message("memory check completed.")

