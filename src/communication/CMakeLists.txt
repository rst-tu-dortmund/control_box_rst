# Library control-box-rst-communication
project(corbo-communication VERSION 0.1 LANGUAGES CXX)

if (NOT CUSTOM_CORBO_PROTO_MSG_DIR)
      set(CUSTOM_CORBO_PROTO_MSG_DIR ${PROJECT_SOURCE_DIR}/custom_proto)
endif (NOT CUSTOM_CORBO_PROTO_MSG_DIR)

# generate messages and services
if (MESSAGE_SUPPORT)
include(message_generation.cmake)
generate_corbo_proto_msgs(${PROJECT_SOURCE_DIR}/proto
         ${CMAKE_CURRENT_BINARY_DIR}/generated/corbo-communication
         protobuf_msg_autogeneration_target
         protobuf_srv_autogeneration_target
         protobuf_grpc_autogeneration_target
         ${MESSAGE_SUPPORT}
         ${CUSTOM_CORBO_PROTO_MSG_DIR}
)
endif (MESSAGE_SUPPORT)
# check variables

if (RPC_SUPPORT AND NOT MESSAGE_SUPPORT)
	message(FATAL_ERROR "RPC_SUPPORT can only activated if MESSAGE_SUPPORT is on as well")
endif (RPC_SUPPORT AND NOT MESSAGE_SUPPORT)

# Create library

if (MESSAGE_SUPPORT)
# Workaround to set the following property: include directories are not created until building
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/generated/corbo-communication)
endif (MESSAGE_SUPPORT)

add_library(corbo_communication STATIC
    src/main_service_client.cpp
    src/signal_target_rpc.cpp
    src/message_parser.cpp
    src/utilities.cpp
    ${proto_msg_gen_cpp_src}
    ${proto_srv_gen_cpp_src}
    ${proto_grpc_gen_cpp_src}
)
target_link_libraries(corbo_communication
    corbo_core
)

if (MESSAGE_SUPPORT)
    target_compile_definitions(corbo_communication PUBLIC MESSAGE_SUPPORT)
    target_link_libraries(corbo_communication libprotobuf)
    add_dependencies(corbo_communication libprotobuf protobuf_msg_autogeneration_target)
else (MESSAGE_SUPPORT)
    message(STATUS "Message support disabled, since MESSAGE_SUPPORT=False.")
endif (MESSAGE_SUPPORT)


if (RPC_SUPPORT)
   target_compile_definitions(corbo_communication PUBLIC RPC_SUPPORT)
   target_link_libraries(corbo_communication libgrpc)
   add_dependencies(corbo_communication protobuf_srv_autogeneration_target protobuf_grpc_autogeneration_target)
else (RPC_SUPPORT)
    message(STATUS "Rpc support disabled, since RPC_SUPPORT=False.")
endif (RPC_SUPPORT)

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build paths.
target_include_directories(corbo_communication PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/control_box_rst>
    PRIVATE src)

target_compile_features(corbo_communication
    PUBLIC cxx_auto_type cxx_range_for cxx_constexpr cxx_lambdas
    PRIVATE cxx_variadic_templates
)

# Set compiler definitions (required on windows for grpc)
if (MSVC)
    target_compile_definitions(corbo_communication PUBLIC -D_WIN32_WINNT=0x600)
endif (MSVC)


# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build paths.
target_include_directories(corbo_communication PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/generated/corbo-communication>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/generated>
  $<INSTALL_INTERFACE:include/control_box_rst>
)


# 'make install' to the correct location
install(TARGETS corbo_communication EXPORT corbo_communicationConfig
    ARCHIVE  DESTINATION lib/control_box_rst
    LIBRARY  DESTINATION lib/control_box_rst
    RUNTIME  DESTINATION bin/control_box_rst)  # This is for Windows
install(DIRECTORY include/ DESTINATION include/control_box_rst) # TODO install generated files

# This makes the project importable from the install directory
# Put config file in per-project dir (name MUST match), can also
# just go into <prefix>/cmake.
install(EXPORT corbo_communicationConfig DESTINATION share/control_box_rst/corbo_communication)

# This makes the project importable from the build directory
export(TARGETS corbo_communication FILE corbo_communicationConfig.cmake)


# Add unit tests
if (BUILD_TESTS)
        add_executable(test_communication
            test/communication_test.cpp)

        target_link_libraries(test_communication
            corbo_communication
            gtest
            #gmock
        )
        add_test(test_communication_test test_communication)
endif (BUILD_TESTS)

# Add header files as custom target in order to display them in the IDE
# TODO check for a cleaner solution
FILE(GLOB_RECURSE HeaderFiles
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")
add_custom_target(corbo_communication_headers SOURCES ${HeaderFiles})

