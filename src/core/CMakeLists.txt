# Library control-box-rst-controllers
project(corbo-core VERSION 0.1 LANGUAGES CXX)

# Include Eigen
find_package(Eigen3 QUIET)
set(EIGEN3_BUILD_INTERFACE_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/../extern/eigen3)
set(EIGEN3_INSTALL_INTERFACE_INCLUDE include/control_box_rst/eigen3)
if (Eigen3_FOUND)
  message(STATUS "Eigen3 found in ${EIGEN3_INCLUDE_DIR}. Using this version instead of the included one.")
  set(EIGEN3_BUILD_INTERFACE_INCLUDE ${EIGEN3_INCLUDE_DIR})
  set(EIGEN3_INSTALL_INTERFACE_INCLUDE ${EIGEN3_BUILD_INTERFACE_INCLUDE})
else (Eigen3_FOUND)
  message(STATUS "Using internal Eigen3 version.")
endif(Eigen3_FOUND)

add_library(corbo_core STATIC
    src/global.cpp
    src/time.cpp
    src/reference_trajectory.cpp
    src/time_series.cpp
    src/signals.cpp
    src/common_signal_target.cpp
    src/data_exporter_interface.cpp
    src/yaml_export.cpp
    src/tsv_export.cpp
)

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build paths.
target_include_directories(corbo_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${EIGEN3_BUILD_INTERFACE_INCLUDE}>
    $<INSTALL_INTERFACE:include/control_box_rst>
    $<INSTALL_INTERFACE:${EIGEN3_INSTALL_INTERFACE_INCLUDE}>
    PRIVATE src)

# If we have compiler requirements for this library, list them here
target_compile_features(corbo_core
    PUBLIC cxx_auto_type cxx_range_for cxx_constexpr cxx_lambdas
    PRIVATE cxx_variadic_templates)

# Set compiler definitions
# e.g., we enable _USE_MATH_DEFINES in order to include M_PI etc. also on visual stuio
target_compile_definitions(corbo_core PUBLIC -D_USE_MATH_DEFINES)
# Set compiler optoins/flags
# target_compile_options(controllerslib PUBLIC -fno-elide-constructors)


# we are currently using the cmake module path in the toplevel cmake directory
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/modules/")

find_package(MATIO)
if(MATIO_FOUND)
    message(STATUS "MatIO found.")
    target_compile_definitions(corbo_core PUBLIC -DMATIO)
    target_include_directories(corbo_core PUBLIC ${MATIO_INCLUDE_DIRS})
    target_link_libraries(corbo_core ${MATIO_LIBRARIES})
endif(MATIO_FOUND)

if (YAML_SUPPORT)
    target_compile_definitions(corbo_core PUBLIC -DYAML_SUPPORT)
    target_link_libraries(corbo_core libyamlcpp)
endif (YAML_SUPPORT)

if(MESSAGE_SUPPORT)
     target_link_libraries(corbo_core corbo_communication)
endif(MESSAGE_SUPPORT)

# 'make install' to the correct location
install(TARGETS corbo_core EXPORT corbo_coreConfig
    ARCHIVE  DESTINATION lib/control_box_rst
    LIBRARY  DESTINATION lib/control_box_rst
    RUNTIME  DESTINATION bin/control_box_rst)  # This is for Windows
install(DIRECTORY include/ DESTINATION include/control_box_rst)

# install eigen if not found
if(NOT Eigen3_FOUND)
install(DIRECTORY ${EIGEN3_BUILD_INTERFACE_INCLUDE} DESTINATION include/control_box_rst/eigen3)
endif(NOT Eigen3_FOUND)

# This makes the project importable from the install directory
# Put config file in per-project dir (name MUST match), can also
# just go into <prefix>/cmake.
install(EXPORT corbo_coreConfig DESTINATION share/control_box_rst/corbo_core)

# This makes the project importable from the build directory
export(TARGETS corbo_core FILE corbo_coreConfig.cmake)

# Add unit tests
if (BUILD_TESTS)
	add_executable(test_core
            test/test_main.cpp
            test/test_yaml_export.cpp)

	target_link_libraries(test_core
	    corbo_core
	    gtest
	    #gmock
	)
	add_test(test_core_test test_core)
endif (BUILD_TESTS)

# Add header files as custom target in order to display them in the IDE
# TODO check for a cleaner solution
FILE(GLOB_RECURSE HeaderFiles
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")
add_custom_target(corbo_core_headers SOURCES ${HeaderFiles})
