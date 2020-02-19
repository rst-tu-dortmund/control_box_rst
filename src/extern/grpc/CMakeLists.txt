# <christoph.roesmann@tu-dortmund.de>

if (NOT RPC_SUPPORT)
  return()
endif (NOT RPC_SUPPORT)
if (NOT MESSAGE_SUPPORT)
  message(STATUS "Skipping gRPC build, since MESSAGE_SUPPORT is disabled")
  return()
endif (NOT MESSAGE_SUPPORT)

project(grpc_builder VERSION 0.1 LANGUAGES C CXX)

# Enable CMake module for external projects
include(ExternalProject)

set(GRPC_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/grpc/src/grpc/include)
set(GRPC_URL https://github.com/grpc/grpc.git)
set(GRPC_BUILD ${CMAKE_CURRENT_BINARY_DIR}/grpc/src/grpc-build)
set(GRPC_TAG b567bb43670f7f2523164fdabdea8786a7595854)

# Find git since we clone gtest from https://github.com/google/googletest
find_package(Git REQUIRED)

# We need thread support
find_package(Threads REQUIRED)

# Get Protobuf Dirs
get_target_property(PROTOBUF_INCLUDE_DIRS libprotobuf INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(protobuf_STATIC_LIBRARIES libprotobuf IMPORTED_LOCATION)
get_target_property(protoc_STATIC_LIBRARIES libprotoc IMPORTED_LOCATION)

message(STATUS "Loading external project grpc")

set(NAME_LIBGRPC_GPR "${CMAKE_CFG_INTDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}gpr${CMAKE_STATIC_LIBRARY_SUFFIX}")
set(NAME_LIBGRPC_CPP_REFLECTION "${CMAKE_CFG_INTDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}grpc++_reflection${CMAKE_STATIC_LIBRARY_SUFFIX}")
set(NAME_LIBGRPC_C_UNSECURE "${CMAKE_CFG_INTDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}grpc_unsecure${CMAKE_STATIC_LIBRARY_SUFFIX}")
set(NAME_LIBGRPC_CPP_UNSECURE "${CMAKE_CFG_INTDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}grpc++_unsecure${CMAKE_STATIC_LIBRARY_SUFFIX}")
set(NAME_LIBGRPC_CARES "third_party/cares/${CMAKE_CFG_INTDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}cares${CMAKE_STATIC_LIBRARY_SUFFIX}")
if (MSVC)
    set(NAME_LIBGRPC_ZLIB "third_party/zlib/${CMAKE_CFG_INTDIR}/zlibstatic${CMAKE_STATIC_LIBRARY_SUFFIX}")
else (MSVC)
    set(NAME_LIBGRPC_ZLIB "third_party/zlib/${CMAKE_CFG_INTDIR}/libz${CMAKE_STATIC_LIBRARY_SUFFIX}")
endif (MSVC)

# CMake arguments forward to the grpc cmake configure command
set(GRPC_CMAKE_CUSTOM_ARGS   
        -DCMAKE_BUILD_TYPE:STRING=Release
        -DCMAKE_VERBOSE_MAKEFILE:BOOL=OFF
		-DgRPC_MSVC_STATIC_RUNTIME:BOOL=ON
        -DPROTOBUF_INCLUDE_DIRS:STRING=${PROTOBUF_INCLUDE_DIRS}
        -DPROTOBUF_LIBRARIES:STRING=${protobuf_STATIC_LIBRARIES}
        -DPROTOC_LIBRARIES:STRING=${protoc_STATIC_LIBRARIES}
        #-DZLIB_ROOT_DIR:STRING=${ZLIB_INSTALL}
)

## the following part loads the visual studio environment variables
## as suggested in the gprc install instructions.
## If enabled: also uncomment the custom configure and build commands in ExternalProject_Add(...)
#if (WIN32 AND MSVC)
#    file(TO_NATIVE_PATH "$ENV{VS140COMNTOOLS}../../VC/vcvarsall.bat" VCVARSALLPATH)#
#    set(PRE_CONFIGURE_COMMAND call ${VCVARSALLPATH} x64 &&)
#endif (WIN32 AND MSVC)

# Download and install GoogleTest
ExternalProject_Add(
    grpc
    DEPENDS libprotobuf libprotoc
    ## We manage our own fork to have control for new updates
    GIT_REPOSITORY ${GRPC_URL}
    GIT_TAG ${GRPC_TAG}
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/grpc
    ## Disable download progress bar (better readability in logs)
    DOWNLOAD_NO_PROGRESS 1
    ## Patch custom CMakeLists.txt
    PATCH_COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/patches/CMakeLists.txt" "<SOURCE_DIR>" && cd <SOURCE_DIR>/third_party/zlib && git checkout -- . && git apply "${CMAKE_CURRENT_SOURCE_DIR}/patches/remove_debug_postfix.patch"
    ## Disable install step
    INSTALL_COMMAND ""
    ## Disable update
    UPDATE_COMMAND ""
	#CONFIGURE_COMMAND ${PRE_CONFIGURE_COMMAND} ${CMAKE_COMMAND} -G${CMAKE_GENERATOR} ${GRPC_CMAKE_CUSTOM_ARGS} "<SOURCE_DIR>"
	#BUILD_COMMAND ${CMAKE_COMMAND} --build .
    CMAKE_CACHE_ARGS ${GRPC_CMAKE_CUSTOM_ARGS}
    BUILD_BYPRODUCTS "<BINARY_DIR>/${NAME_LIBGRPC_GPR};<BINARY_DIR>/${NAME_LIBGRPC_C_UNSECURE};<BINARY_DIR>/${NAME_LIBGRPC_CPP_UNSECURE};<BINARY_DIR>/${NAME_LIBGRPC_ZLIB};<BINARY_DIR>/${NAME_LIBGRPC_CARES}"
)
# Get source and binary directories
ExternalProject_Get_Property(grpc source_dir binary_dir)

set(GRPC_CPP_PLUGIN_EXECUTABLE ${binary_dir}/${CMAKE_CFG_INTDIR}/grpc_cpp_plugin${CMAKE_EXECUTABLE_SUFFIX} PARENT_SCOPE)

# Create interface libraries for dependent packages
add_library(libgrpc_gpr STATIC IMPORTED GLOBAL)
set_target_properties(libgrpc_gpr PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${source_dir}/include
    IMPORTED_LOCATION ${binary_dir}/${NAME_LIBGRPC_GPR}
)
add_dependencies(libgrpc_gpr grpc)

add_library(libgrpc_c STATIC IMPORTED GLOBAL)
set_target_properties(libgrpc_c PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${source_dir}/include
    INTERFACE_LINK_LIBRARIES "libgrpc_gpr;${binary_dir}/${NAME_LIBGRPC_ZLIB};${binary_dir}/${NAME_LIBGRPC_CARES};${CMAKE_THREAD_LIBS_INIT}"
    IMPORTED_LOCATION ${binary_dir}/${NAME_LIBGRPC_C_UNSECURE}
)

if (WIN32 AND MSVC)
	set_property(TARGET libgrpc_c APPEND PROPERTY INTERFACE_LINK_LIBRARIES ";wsock32;ws2_32")
endif (WIN32 AND MSVC)

add_dependencies(libgrpc_c grpc)

add_library(libgrpc STATIC IMPORTED GLOBAL)
set_target_properties(libgrpc PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${source_dir}/include
    INTERFACE_LINK_LIBRARIES libgrpc_c
    IMPORTED_LOCATION ${binary_dir}/${NAME_LIBGRPC_CPP_UNSECURE}
)
add_dependencies(libgrpc libgrpc_c)

## Not yet included in our custom CMakeList.txt for grpc:
#add_library(libgrpc_reflection STATIC IMPORTED GLOBAL)
#set_target_properties(libgrpc_reflection PROPERTIES
#    INTERFACE_INCLUDE_DIRECTORIES ${source_dir}/include
#    INTERFACE_LINK_LIBRARIES libgrpc
#    IMPORTED_LOCATION ${binary_dir}/${NAME_LIBGRPC_CPP_REFLECTION}
#)


#set_property(TARGET libgrpc PROPERTY INTERFACE_LINK_LIBRARIES ${binary_dir}/${NAME_LIBGRPC_CPP_UNSECURE} )

#set_property(TARGET libgrpc PROPERTY GRPC_CPP_PLUGIN_EXECUTABLE "${binary_dir}/${CMAKE_CFG_INTDIR}/grpc_cpp_plugin")
add_dependencies(libgrpc libgrpc_c)

# Workaround to set the following property: include directories are not created until building
file(MAKE_DIRECTORY ${source_dir}/include)

# TODO:
# wrap executable in a target
# see here: https://github.com/IvanSafonov/grpc-cmake-example/blob/master/cmake/FindGRPC.cmake
