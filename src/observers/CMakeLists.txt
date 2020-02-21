# Library control-box-rst-observers
project(corbo-observers VERSION 0.1 LANGUAGES CXX)

add_library(corbo_observers STATIC
    src/dummy.cpp
)

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build paths.
target_include_directories(corbo_observers PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/control_box_rst>
    PRIVATE src)

# If we have compiler requirements for this library, list them here
target_compile_features(corbo_observers
    PUBLIC cxx_auto_type cxx_range_for cxx_constexpr cxx_lambdas
    PRIVATE cxx_variadic_templates)

# Set compiler definitions
# target_compile_definitions(controllerslib PRIVATE MYDEF=${BLABLA})
# Set compiler optoins/flags
# target_compile_options(controllerslib PUBLIC -fno-elide-constructors)


# Depend on a library that we defined in the top-level file
target_link_libraries(corbo_observers
    corbo_core
    corbo_communication
)


# 'make install' to the correct location
install(TARGETS corbo_observers EXPORT corbo_observersConfig
    ARCHIVE  DESTINATION lib/control_box_rst
    LIBRARY  DESTINATION lib/control_box_rst
    RUNTIME  DESTINATION bin/control_box_rst)  # This is for Windows
install(DIRECTORY include/ DESTINATION include/control_box_rst)

# This makes the project importable from the install directory
# Put config file in per-project dir (name MUST match), can also
# just go into <prefix>/cmake.
install(EXPORT corbo_observersConfig DESTINATION share/control_box_rst/corbo_observers)

# This makes the project importable from the build directory
export(TARGETS corbo_observers FILE corbo_observersConfig.cmake)

# Add unit tests
if (BUILD_TESTS)
	add_executable(test_observers
	    test/observers_test.cpp)

	target_link_libraries(test_observers
	    corbo_observers
	    gtest
	    #gmock
	)
        add_test(test_observers_test test_observers)
endif (BUILD_TESTS)

# Add header files as custom target in order to display them in the IDE
# TODO check for a cleaner solution
FILE(GLOB_RECURSE HeaderFiles
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")
add_custom_target(corbo_observers_headers SOURCES ${HeaderFiles})
