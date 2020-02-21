# Library control-box-rst-optimal-control
project(corbo-optimal-control VERSION 0.1 LANGUAGES CXX)

add_library(corbo_optimal_control STATIC
    src/structured_ocp/structured_optimal_control_problem.cpp
    src/structured_ocp/discretization_grids/full_discretization_grid_base.cpp
    src/structured_ocp/discretization_grids/full_discretization_grid_move_blocking_base.cpp
    src/structured_ocp/discretization_grids/non_uniform_full_discretization_grid_base.cpp
    src/structured_ocp/discretization_grids/shooting_grid_base.cpp
    src/structured_ocp/discretization_grids/non_uniform_shooting_grid_base.cpp
    src/structured_ocp/discretization_grids/finite_differences_grid.cpp
    src/structured_ocp/discretization_grids/finite_differences_grid_move_blocking.cpp
    src/structured_ocp/discretization_grids/finite_differences_variable_grid.cpp
    src/structured_ocp/discretization_grids/non_uniform_finite_differences_variable_grid.cpp
    src/structured_ocp/discretization_grids/multiple_shooting_grid.cpp
    src/structured_ocp/discretization_grids/multiple_shooting_variable_grid.cpp
    src/structured_ocp/discretization_grids/non_uniform_multiple_shooting_variable_grid.cpp
    src/functions/stage_functions.cpp
    src/functions/stage_preprocessor.cpp
    src/functions/nlp_functions.cpp
    src/functions/quadratic_cost.cpp
    src/functions/quadratic_state_cost.cpp
    src/functions/quadratic_control_cost.cpp
    src/functions/final_state_cost.cpp
    src/functions/final_state_constraints.cpp
)

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build paths.
target_include_directories(corbo_optimal_control PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/control_box_rst>
    PRIVATE src)

# If we have compiler requirements for this library, list them here
target_compile_features(corbo_optimal_control
    PUBLIC cxx_auto_type cxx_range_for cxx_constexpr cxx_lambdas
    PRIVATE cxx_variadic_templates)

# Set compiler definitions
# target_compile_definitions(controllerslib PRIVATE MYDEF=${BLABLA})
# Set compiler optoins/flags
# target_compile_options(controllerslib PUBLIC -fno-elide-constructors)

# we are currently using the cmake module path in the toplevel cmake directory
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/modules/")


# Depend on a library that we defined in the top-level file
target_link_libraries(corbo_optimal_control
    corbo_core
    corbo_communication
    corbo_optimization
    corbo_numerics
    corbo_systems
)


# 'make install' to the correct location
install(TARGETS corbo_optimal_control EXPORT corbo_optimal_controlConfig
    ARCHIVE  DESTINATION lib/control_box_rst
    LIBRARY  DESTINATION lib/control_box_rst
    RUNTIME  DESTINATION bin/control_box_rst)  # This is for Windows
install(DIRECTORY include/ DESTINATION include/control_box_rst)

# This makes the project importable from the install directory
# Put config file in per-project dir (name MUST match), can also
# just go into <prefix>/cmake.
install(EXPORT corbo_optimal_controlConfig DESTINATION share/control_box_rst/corbo_optimal_control)

# This makes the project importable from the build directory
export(TARGETS corbo_optimal_control FILE corbo_optimal_controlConfig.cmake)

# Add unit tests
if (BUILD_TESTS)
	add_executable(test_optimal_control
	    test/test_main.cpp
            #test/test_full_discretization_grid.cpp
            #test/test_multiple_shooting_grid.cpp
            #test/test_single_shooting_grid.cpp
	)

	target_link_libraries(test_optimal_control
            corbo_optimal_control
	    gtest
	    #gmock
	)
        add_test(test_optimal_control_test test_optimal_control)
endif (BUILD_TESTS)

# Add header files as custom target in order to display them in the IDE
# TODO check for a cleaner solution
FILE(GLOB_RECURSE HeaderFiles
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")
add_custom_target(corbo_optimal_control_headers SOURCES ${HeaderFiles})
