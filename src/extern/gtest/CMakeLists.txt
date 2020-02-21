# <christoph.roesmann@tu-dortmund.de>
project(gtest_builder VERSION 0.1 LANGUAGES C CXX)

if (NOT BUILD_TESTS)
  return()
endif ()

message(STATUS "Loading external projects gtest and gmock")

# Enable CMake module for external projects
include(ExternalProject)

# We need thread support
find_package(Threads REQUIRED)

# Find git since we clone gtest from https://github.com/google/googletest
find_package(Git REQUIRED)

set(GTESTLIBNAME "${CMAKE_CFG_INTDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}gtest${CMAKE_STATIC_LIBRARY_SUFFIX}")
set(GMOCKLIBNAME "${CMAKE_CFG_INTDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}gmock${CMAKE_STATIC_LIBRARY_SUFFIX}")

# CMake arguments forward to the cmake configure command
set(GTEST_CMAKE_CUSTOM_ARGS   
		-DINSTALL_GTEST:BOOL=OFF 
		-DBUILD_SHARED_LIBS:BOOL=OFF
        #-Dgtest_disable_pthreads:BOOL=ON
)

# Download and install GoogleTest
ExternalProject_Add(
    googletest
    #URL https://github.com/google/googletest/archive/master.zip
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG "release-1.8.1"
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gtest
    # Disable download progress bar (better readability in logs)
    DOWNLOAD_NO_PROGRESS 1
	# Apply patches
	PATCH_COMMAND git checkout -- . && git apply "${CMAKE_CURRENT_SOURCE_DIR}/patches/0001-removes-cmake-debug-postfix.patch"
    # Disable install step
    INSTALL_COMMAND ""
    # Disable update
    UPDATE_COMMAND ""
    # We don't need to run update command. Takes time and the version we initially d/l will should be fine
    CMAKE_CACHE_ARGS ${GTEST_CMAKE_CUSTOM_ARGS}
    BUILD_BYPRODUCTS "<BINARY_DIR>/googlemock/gtest/${GTESTLIBNAME};<BINARY_DIR>/googlemock/${GMOCKLIBNAME}"
)
# Get source and binary directories
ExternalProject_Get_Property(googletest source_dir binary_dir)

# Get GTest source and binary directories from CMake project
add_library(gtest INTERFACE IMPORTED GLOBAL)
set_property(TARGET gtest PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${source_dir}/googletest/include)
set_property(TARGET gtest PROPERTY INTERFACE_LINK_LIBRARIES ${binary_dir}/googlemock/gtest/${GTESTLIBNAME} ${CMAKE_THREAD_LIBS_INIT})
add_dependencies(gtest googletest)

# Get GMock source and binary directories from CMake project
# add_library(gmock INTERFACE IMPORTED GLOBAL)
# set_property(TARGET gmock PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${source_dir}/googlemock/include)
# set_property(TARGET gmock PROPERTY INTERFACE_LINK_LIBRARIES ${binary_dir}/googlemock/${GMOCKLIBNAME} ${CMAKE_THREAD_LIBS_INIT})
# add_dependencies(gmock googletest)

# Workaround to set the following property: include directories are not created until building
file(MAKE_DIRECTORY ${source_dir}/googletest/include)
#file(MAKE_DIRECTORY ${source_dir}/googlemock/include)
