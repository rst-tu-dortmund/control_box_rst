cmake_minimum_required(VERSION 2.8.3)
project(control_box_rst)

set(MESSAGE_SUPPORT OFF CACHE BOOL "Message support (downloads and compiles Protobuf)" FORCE)
set(RPC_SUPPORT OFF CACHE BOOL "Network communication (downloads and comiles gRPC, requires MESSAGE_SUPPORT" FORCE)
set(YAML_SUPPORT OFF CACHE BOOL "Yaml support" FORCE)
set(BUILD_GUI OFF CACHE BOOL "Build gui (requires QT)" FORCE)
set(BUILD_TESTS OFF CACHE BOOL "Build unit tests" FORCE)

# Set to Release in order to speed up the program significantly
set(CMAKE_BUILD_TYPE Release) #None, Debug, Release, RelWithDebInfo, MinSizeRel

add_subdirectory(src)

# Install catkin package.xml
install(FILES package.xml DESTINATION share/control_box_rst)
