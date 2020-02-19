# <christoph.roesmann@tu-dortmund.de>

if (NOT MESSAGE_SUPPORT)
  return()
endif (NOT MESSAGE_SUPPORT)

project(protobuf_builder VERSION 0.1 LANGUAGES C CXX)

# Enable CMake module for external projects
include(ExternalProject)

# We need thread support
#find_package(Threads REQUIRED)

# Find git since we clone gtest from https://github.com/google/googletest
find_package(Git REQUIRED)

message(STATUS "Loading external project protobuf")

#define library names:
set(NAME_LIBPROTOBUF "${CMAKE_CFG_INTDIR}/libprotobuf${CMAKE_STATIC_LIBRARY_SUFFIX}")
set(NAME_LIBPROTOBUF_LITE "${CMAKE_CFG_INTDIR}/libprotobuf-lite${CMAKE_STATIC_LIBRARY_SUFFIX}")
set(NAME_LIBPROTOC "${CMAKE_CFG_INTDIR}/libprotoc${CMAKE_STATIC_LIBRARY_SUFFIX}")

set(PROTOBUF_CMAKE_ARGS -Dprotobuf_BUILD_TESTS=OFF -Dprotobuf_BUILD_EXAMPLES=OFF -Dprotobuf_DEBUG_POSTFIX= )

# Download and install GoogleTest
ExternalProject_Add(
    protobuf
    # We manage our own fork to have control for new updates
    GIT_REPOSITORY https://github.com/google/protobuf.git
    GIT_TAG "v3.6.1"
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/protobuf
    #SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/protobuf/src/protobuf/cmake
    # Disable download progress bar (better readability in logs)
    DOWNLOAD_NO_PROGRESS 1
    # We want to deactivate building unit tests (it also failes if we do not specify the correct gmock path)
    CMAKE_ARGS ${PROTOBUF_CMAKE_ARGS}
    # We need to switch the directory for cmake to <SOURCE_DIR>/cmake
    # CMAKE 3.7 introduces SOURCE_SUBDIR (https://gitlab.kitware.com/cmake/cmake/merge_requests/63)
    # but up to now the following workaround is used:
    CONFIGURE_COMMAND ${CMAKE_COMMAND} ${PROTOBUF_CMAKE_ARGS} -G${CMAKE_GENERATOR} "<SOURCE_DIR>/cmake"
    # Disable install step
    INSTALL_COMMAND ""
    # Disable update
    UPDATE_COMMAND ""
    BUILD_BYPRODUCTS "<BINARY_DIR>/${NAME_LIBPROTOBUF};<BINARY_DIR>/protoc"
)
# Get source and binary directories
ExternalProject_Get_Property(protobuf source_dir binary_dir)

# Create an interface library for dependent packages
add_library(libprotobuf STATIC IMPORTED GLOBAL)
define_property(TARGET PROPERTY PROTOBUF_PROTOC_EXECUTABLE BRIEF_DOCS "Protobuff Generator Executable" FULL_DOCS "Protobuff Generator Executable")
set_target_properties(libprotobuf PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${source_dir}/src
    IMPORTED_LOCATION ${binary_dir}/${NAME_LIBPROTOBUF}
    PROTOBUF_PROTOC_EXECUTABLE "${binary_dir}/${CMAKE_CFG_INTDIR}/protoc${CMAKE_EXECUTABLE_SUFFIX}"
)

add_dependencies(libprotobuf protobuf)

add_library(libprotoc STATIC IMPORTED GLOBAL)
set_target_properties(libprotoc PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${source_dir}/src
    IMPORTED_LOCATION ${binary_dir}/${NAME_LIBPROTOC}
    PROTOBUF_PROTOC_EXECUTABLE "${binary_dir}/${CMAKE_CFG_INTDIR}/protoc"
)
add_dependencies(libprotoc protobuf)

# Workaround to set the following property: include directories are not created until building
file(MAKE_DIRECTORY ${source_dir}/src)
