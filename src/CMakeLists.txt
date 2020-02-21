# Package control-box-rst toplevel cmake

cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(control-box-rst VERSION 0.1 LANGUAGES CXX)

option(MESSAGE_SUPPORT "Message support (downloads and compiles Protobuf)" ON)
option(RPC_SUPPORT "Network communication (downloads and comiles gRPC, requires MESSAGE_SUPPORT" ON)
option(YAML_SUPPORT "Yaml support" ON)
option(BUILD_GUI "Build gui (requires QT)" ON)
option(ONLY_GUI "If BUILD_GUI is on, build only parts required for running the gui" OFF)
option(BUILD_TESTS "Build unit tests" OFF)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

# Compile static libraries with Visual Studio
if(MSVC)
	include(cmake/msvc_static_runtime.cmake)
	# Add some common library paths to facilitate package finding
	list(APPEND CMAKE_PREFIX_PATH "C:/Qt/5.8/msvc2015_64/;C:/Qt/5.10.0/msvc2015_64/;D:/Qt/5.10.0/msvc2015_64/")
 	# Enable support of large object files and parallel builds
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj /MP /maxcpucount")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /MP /maxcpucount")
elseif (APPLE)
	# Add some common library paths to facilitate package finding
	list(APPEND CMAKE_PREFIX_PATH "$ENV{HOME}/Applications/Qt/5.5/clang_64/")
endif(MSVC)

# Add further compiler specific flags
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	# using GCC
	# enforce error if return type is missing:
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=return-type")
endif()

# Include Boost as an imported target
#find_package(Boost REQUIRED)
#add_library(boost INTERFACE IMPORTED)
#set_property(TARGET boost PROPERTY
#    INTERFACE_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIR})

## External 3rd party libs that we include

# Gtest and Gmock for unit testing
add_subdirectory(extern/gtest EXCLUDE_FROM_ALL)

# Protobuf for serialization, parameters
add_subdirectory(extern/protobuf EXCLUDE_FROM_ALL)

# RPC with protobuf messages
add_subdirectory(extern/grpc EXCLUDE_FROM_ALL)

# yaml-cpp
add_subdirectory(extern/yaml-cpp EXCLUDE_FROM_ALL)

# Configure tests
if (BUILD_TESTS)
	enable_testing()
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/lcov.cmake.in
                       ${CMAKE_CURRENT_BINARY_DIR}/lcov.cmake @ONLY)
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/lcov_skip_build.cmake.in
                       ${CMAKE_CURRENT_BINARY_DIR}/lcov_skip_build.cmake @ONLY)
endif (BUILD_TESTS)


# Targets
add_subdirectory(core)
add_subdirectory(communication)

if (NOT ONLY_GUI)
 add_subdirectory(systems)
 add_subdirectory(plants)
 add_subdirectory(controllers)
 add_subdirectory(observers)
 add_subdirectory(numerics)
 add_subdirectory(optimization)
 add_subdirectory(optimal_control)
 add_subdirectory(tasks)
 add_subdirectory(master)
endif (NOT ONLY_GUI)

add_subdirectory(gui)

# Install top-level cmake config
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/control_box_rstConfig.cmake DESTINATION share/control_box_rst)

# Documentation
find_package(Doxygen)
if (DOXYGEN_FOUND)
	# set input and output files
	set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/../docs/Doxyfile.config)
	set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/doxygen/Doxyfile)

	# request to configure the file
	configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

	# create custom target (add ALL if it should be part of make all)
	add_custom_target(doc
	    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
	    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doxygen
	    COMMENT "Generating API documentation with Doxygen"
	    VERBATIM )
else (DOXYGEN_FOUND)
	message("Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)
