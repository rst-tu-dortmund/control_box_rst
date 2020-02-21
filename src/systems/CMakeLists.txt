# Library control-box-rst-systems
project(corbo-systems VERSION 0.1 LANGUAGES CXX)

add_library(corbo_systems STATIC
    src/standard_filters.cpp
    src/system_dynamics_interface.cpp
    src/time_value_buffer.cpp
    src/one_step_predictor.cpp
)

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build paths.
target_include_directories(corbo_systems PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/control_box_rst>
    PRIVATE src)

# If we have compiler requirements for this library, list them here
target_compile_features(corbo_systems
    PUBLIC cxx_auto_type cxx_range_for cxx_constexpr cxx_lambdas
    PRIVATE cxx_variadic_templates)

# Set compiler definitions
# target_compile_definitions(controllerslib PRIVATE MYDEF=${BLABLA})
# Set compiler optoins/flags
# target_compile_options(controllerslib PUBLIC -fno-elide-constructors)


# Depend on a library that we defined in the top-level file
target_link_libraries(corbo_systems
    corbo_core
    corbo_communication
    corbo_numerics
)


# 'make install' to the correct location
install(TARGETS corbo_systems EXPORT corbo_systemsConfig
    ARCHIVE  DESTINATION lib/control_box_rst
    LIBRARY  DESTINATION lib/control_box_rst
    RUNTIME  DESTINATION bin/control_box_rst)  # This is for Windows
install(DIRECTORY include/ DESTINATION include/control_box_rst)

# This makes the project importable from the install directory
# Put config file in per-project dir (name MUST match), can also
# just go into <prefix>/cmake.
install(EXPORT corbo_systemsConfig DESTINATION share/control_box_rst/corbo_systems)

# This makes the project importable from the build directory
export(TARGETS corbo_systems FILE corbo_systemsConfig.cmake)

# Add unit tests
if (BUILD_TESTS)
	add_executable(test_systems
	    test/system_test.cpp)

	target_link_libraries(test_systems
	    corbo_systems
	    gtest
	    #gmock
	)
	add_test(test_systems_test test_systems)
endif (BUILD_TESTS)

# Add header files as custom target in order to display them in the IDE
# TODO check for a cleaner solution
FILE(GLOB_RECURSE HeaderFiles
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")
add_custom_target(corbo_systems_headers SOURCES ${HeaderFiles})
