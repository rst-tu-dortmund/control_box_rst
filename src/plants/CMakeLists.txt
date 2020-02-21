# Library control-box-rst-plants
project(corbo-plants VERSION 0.1 LANGUAGES CXX)

add_library(corbo_plants STATIC
    src/plant_interface.cpp
    src/simulated_plant.cpp
    src/simulated_plant_threaded.cpp
    src/disturbances.cpp
)

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build paths.
target_include_directories(corbo_plants PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/control_box_rst>
    PRIVATE src)

# If we have compiler requirements for this library, list them here
target_compile_features(corbo_plants
    PUBLIC cxx_auto_type cxx_range_for cxx_constexpr cxx_lambdas
    PRIVATE cxx_variadic_templates)

# Set compiler definitions
# target_compile_definitions(controllerslib PRIVATE MYDEF=${BLABLA})
# Set compiler optoins/flags
# target_compile_options(controllerslib PUBLIC -fno-elide-constructors)


# Depend on a library that we defined in the top-level file
target_link_libraries(corbo_plants
    corbo_core
    corbo_communication
    corbo_systems
)


# 'make install' to the correct location
install(TARGETS corbo_plants EXPORT corbo_plantsConfig
    ARCHIVE  DESTINATION lib/control_box_rst
    LIBRARY  DESTINATION lib/control_box_rst
    RUNTIME  DESTINATION bin/control_box_rst)  # This is for Windows
install(DIRECTORY include/ DESTINATION include/control_box_rst)

# This makes the project importable from the install directory
# Put config file in per-project dir (name MUST match), can also
# just go into <prefix>/cmake.
install(EXPORT corbo_plantsConfig DESTINATION share/control_box_rst/corbo_plants)

# This makes the project importable from the build directory
export(TARGETS corbo_plants FILE corbo_plantsConfig.cmake)

# Add unit tests
if (BUILD_TESTS)
	add_executable(test_plants
	    test/system_test.cpp)

	target_link_libraries(test_plants
	    corbo_plants
	    gtest
	    #gmock
	)
	add_test(test_plants_test test_plants)
endif (BUILD_TESTS)

# Add header files as custom target in order to display them in the IDE
# TODO check for a cleaner solution
FILE(GLOB_RECURSE HeaderFiles
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")
add_custom_target(corbo_plants_headers SOURCES ${HeaderFiles})
