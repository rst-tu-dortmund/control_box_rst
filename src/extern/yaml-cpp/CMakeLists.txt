# <christoph.roesmann@tu-dortmund.de>
project(ymalcpp_builder VERSION 0.1 LANGUAGES C CXX)

if (NOT YAML_SUPPORT)
  return()
endif (NOT YAML_SUPPORT)

# Enable CMake module for external projects
include(ExternalProject)


# Find git since we clone yaml-cpp from https://github.com/jbeder/yaml-cpp.git
find_package(Git REQUIRED)

message(STATUS "Loading external project yaml-cpp")

if (MSVC)
    set(YAMLCPPLIBNAME "${CMAKE_CFG_INTDIR}/libyaml-cppmt${CMAKE_STATIC_LIBRARY_SUFFIX}")
else (MSVC)
    set(YAMLCPPLIBNAME "${CMAKE_CFG_INTDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}yaml-cpp${CMAKE_STATIC_LIBRARY_SUFFIX}")
endif (MSVC)

set(YAMLCPP_CMAKE_CUSTOM_ARGS   
    -DYAML_CPP_BUILD_TESTS:BOOL=OFF
    -DYAML_CPP_BUILD_TOOLS:BOOL=OFF
    -DYAML_CPP_BUILD_CONTRIB:BOOL=OFF
    -DMSVC_SHARED_RT:BOOL=OFF
)

# Download and install yaml-cpp
ExternalProject_Add(
    yaml-cpp
    GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
    GIT_TAG ee99c4151c1af794a412b101a75921c086acaac0
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/yaml-cpp
    # Disable download progress bar (better readability in logs)
    DOWNLOAD_NO_PROGRESS 1
	# Apply patches
	PATCH_COMMAND git checkout -- . && git apply "${CMAKE_CURRENT_SOURCE_DIR}/patches/0001-removes-cmake-debug-postfix.patch"
    # Disable install step
    INSTALL_COMMAND ""
    # Disable update
    UPDATE_COMMAND ""
    CMAKE_CACHE_ARGS ${YAMLCPP_CMAKE_CUSTOM_ARGS}
    BUILD_BYPRODUCTS "<BINARY_DIR>/${YAMLCPPLIBNAME}"
)
# Get source and binary directories
ExternalProject_Get_Property(yaml-cpp source_dir binary_dir)

# Get yaml-cpp source and binary directories from CMake project
add_library(libyamlcpp INTERFACE IMPORTED GLOBAL)
set_property(TARGET libyamlcpp PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${source_dir}/include)
set_property(TARGET libyamlcpp PROPERTY INTERFACE_LINK_LIBRARIES ${binary_dir}/${YAMLCPPLIBNAME})
add_dependencies(libyamlcpp yaml-cpp)

# Workaround to set the following property: include directories are not created until building
file(MAKE_DIRECTORY ${source_dir}/include)

