# Library control-box-rst-controllers
project(corbo-tasks VERSION 0.1 LANGUAGES CXX)

add_library(corbo_tasks STATIC
    src/environment.cpp
    src/task_open_loop_control.cpp
    src/task_closed_loop_control.cpp
    src/benchmark_task_increasing_n_open_loop.cpp
    src/benchmark_task_varying_initial_state.cpp
)

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build paths.
target_include_directories(corbo_tasks PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/control_box_rst>
    PRIVATE src)

# If we have compiler requirements for this library, list them here
target_compile_features(corbo_tasks
    PUBLIC cxx_auto_type cxx_range_for cxx_constexpr cxx_lambdas
    PRIVATE cxx_variadic_templates)

# Set compiler definitions
# target_compile_definitions(controllerslib PRIVATE MYDEF=${BLABLA})
# Set compiler optoins/flags
# target_compile_options(controllerslib PUBLIC -fno-elide-constructors)


# Depend on a library that we defined in the top-level file
target_link_libraries(corbo_tasks
    corbo_core
    corbo_communication
    corbo_controllers
    corbo_observers
    corbo_systems
    corbo_plants
    corbo_numerics
)


# 'make install' to the correct location
install(TARGETS corbo_tasks EXPORT corbo_tasksConfig
    ARCHIVE  DESTINATION lib/control_box_rst
    LIBRARY  DESTINATION lib/control_box_rst
    RUNTIME  DESTINATION bin/control_box_rst)  # This is for Windows
install(DIRECTORY include/ DESTINATION include/control_box_rst)

# This makes the project importable from the install directory
# Put config file in per-project dir (name MUST match), can also
# just go into <prefix>/cmake.
install(EXPORT corbo_tasksConfig DESTINATION share/control_box_rst/corbo_tasks)

# This makes the project importable from the build directory
export(TARGETS corbo_tasks FILE corbo_tasksConfig.cmake)

# Add unit tests
if (BUILD_TESTS)
        add_executable(test_tasks
            test/environment_test.cpp)

        target_link_libraries(test_tasks
            corbo_tasks
	    gtest
	    #gmock
	)
        add_test(test_tasks_test test_tasks)
endif (BUILD_TESTS)

# Add header files as custom target in order to display them in the IDE
# TODO check for a cleaner solution
FILE(GLOB_RECURSE HeaderFiles
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")
add_custom_target(corbo_tasks_headers SOURCES ${HeaderFiles})
