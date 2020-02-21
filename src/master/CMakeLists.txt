# Library control-box-rst-controllers
project(corbo-master VERSION 0.1 LANGUAGES CXX)

if (NOT RPC_SUPPORT)
   message(STATUS "Skipping corbo-master build, since RPC_SUPPORT is disabled")
   return()
endif (NOT RPC_SUPPORT)

add_library(corbo_master_lib STATIC
    src/master.cpp
)

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build paths.
target_include_directories(corbo_master_lib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/control_box_rst>
    PRIVATE src)

# If we have compiler requirements for this library, list them here
target_compile_features(corbo_master_lib
    PUBLIC cxx_auto_type cxx_range_for cxx_constexpr cxx_lambdas
    PRIVATE cxx_variadic_templates)

# Set compiler definitions
# target_compile_definitions(controllerslib PRIVATE MYDEF=${BLABLA})

# Set compiler optoins/flags
# target_compile_options(corbo_master_lib PUBLIC /MTd)

# Depend on a library that we defined in the top-level file
target_link_libraries(corbo_master_lib
    corbo_communication
    corbo_tasks
)

# Gui application
add_executable(corbo_master
    src/app.cpp
)

target_link_libraries(corbo_master
    corbo_master_lib
)


# 'make install' to the correct location
install(TARGETS corbo_master_lib corbo_master EXPORT corbo_masterConfig
    ARCHIVE  DESTINATION lib/control_box_rst
    LIBRARY  DESTINATION lib/control_box_rst
    RUNTIME  DESTINATION bin/control_box_rst)  # This is for Windows
install(DIRECTORY include/ DESTINATION include/control_box_rst)

# This makes the project importable from the install directory
# Put config file in per-project dir (name MUST match), can also
# just go into <prefix>/cmake.
install(EXPORT corbo_masterConfig DESTINATION share/control_box_rst/corbo_master)

# This makes the project importable from the build directory
export(TARGETS corbo_master_lib corbo_master FILE corbo_masterConfig.cmake)

# Add unit tests
if (BUILD_TESTS)
	add_executable(test_master
	    test/master_test.cpp)

	target_link_libraries(test_master
	    corbo_master_lib
	    gtest
	    #gmock
	)
	add_test(test_master_test test_master)
endif (BUILD_TESTS)

# Add header files as custom target in order to display them in the IDE
# TODO check for a cleaner solution
FILE(GLOB_RECURSE HeaderFiles
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")
add_custom_target(corbo_master_headers SOURCES ${HeaderFiles})
