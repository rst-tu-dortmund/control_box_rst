# Library control-box-rst-optimization
project(corbo-optimization VERSION 0.1 LANGUAGES CXX)

add_library(corbo_optimization STATIC
    src/optimization_problem_interface.cpp
    src/simple_optimization_problem.cpp
    src/hyper_graph/hyper_graph_optimization_problem_base.cpp
    src/hyper_graph/hyper_graph_optimization_problem_edge_based.cpp
    src/hyper_graph/hyper_graph_optimization_problem_vertex_based.cpp
    src/hyper_graph/edge_interface.cpp
    src/hyper_graph/edge_set.cpp
    src/hyper_graph/vertex_set.cpp
    src/hyper_graph/vertex_interface.cpp
    src/hyper_graph/hyper_graph.cpp
    src/solver/levenberg_marquardt_dense.cpp
    src/solver/levenberg_marquardt_sparse.cpp
    src/solver/nlp_solver_ipopt.cpp
    src/solver/nlp_solver_ipopt_wrapper.cpp
    src/solver/qp_solver_osqp.cpp
    src/misc.cpp
)

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build paths.
target_include_directories(corbo_optimization PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/control_box_rst>
    PRIVATE src)

# If we have compiler requirements for this library, list them here
target_compile_features(corbo_optimization
    PUBLIC cxx_auto_type cxx_range_for cxx_constexpr cxx_lambdas
    PRIVATE cxx_variadic_templates)

# Set compiler definitions
# target_compile_definitions(controllerslib PRIVATE MYDEF=${BLABLA})
# Set compiler optoins/flags
# target_compile_options(controllerslib PUBLIC -fno-elide-constructors)

# we are currently using the cmake module path in the toplevel cmake directory
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/modules/")

find_package(IPOPT QUIET)
if(IPOPT_FOUND)
    message(STATUS "IPOPT found.")
    target_compile_definitions(corbo_optimization PUBLIC ${IPOPT_DEFINITIONS} -DIPOPT)
    target_include_directories(corbo_optimization PUBLIC ${IPOPT_INCLUDE_DIRS})
    SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${IPOPT_LINK_FLAGS}" )
    target_link_libraries(corbo_optimization ${IPOPT_LIBRARIES})
else(IPOPT_FOUND)
    message(STATUS "IPOPT not found. Solver IPOPT will not be available.")
endif(IPOPT_FOUND)

find_package(osqp QUIET)
if(osqp_FOUND)
    message(STATUS "OSQP found.")
    target_compile_definitions(corbo_optimization PUBLIC -DOSQP)
    set (CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -ldl")
    target_link_libraries(corbo_optimization osqp::osqpstatic -ldl)
else(osqp_FOUND)
    message(STATUS "OSQP not found. Solver OSQP will not be available.")
endif(osqp_FOUND)

# Depend on a library that we defined in the top-level file
target_link_libraries(corbo_optimization
    corbo_core
    corbo_communication
    corbo_numerics
    corbo_systems
)


# 'make install' to the correct location
install(TARGETS corbo_optimization EXPORT corbo_optimizationConfig
    ARCHIVE  DESTINATION lib/control_box_rst
    LIBRARY  DESTINATION lib/control_box_rst
    RUNTIME  DESTINATION bin/control_box_rst)  # This is for Windows
install(DIRECTORY include/ DESTINATION include/control_box_rst)

# This makes the project importable from the install directory
# Put config file in per-project dir (name MUST match), can also
# just go into <prefix>/cmake.
install(EXPORT corbo_optimizationConfig DESTINATION share/control_box_rst/corbo_optimization)

# This makes the project importable from the build directory
export(TARGETS corbo_optimization FILE corbo_optimizationConfig.cmake)

# Add unit tests
if (BUILD_TESTS)
	add_executable(test_optimization
	    test/test_main.cpp
            test/test_simple_optimization_problem.cpp
            test/test_hyper_graph_optimization_problem_ignore_structure.cpp
            test/test_hyper_graph_optimization_problem_edge_based.cpp
            test/test_hyper_graph_optimization_problem_vertex_based.cpp
            test/test_hyper_graph.cpp
            test/test_solver_ipopt.cpp
            test/test_solver_osqp.cpp
            #test/test_levenberg_marquardt_dense.cpp
            #test/test_levenberg_marquardt_dense_hyper_graph.cpp
            #test/test_levenberg_marquardt_sparse.cpp
	)

	target_link_libraries(test_optimization
            corbo_optimization
	    gtest
	    #gmock
	)
        add_test(test_optimization_test test_optimization)
endif (BUILD_TESTS)

# Add header files as custom target in order to display them in the IDE
# TODO check for a cleaner solution
FILE(GLOB_RECURSE HeaderFiles
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")
add_custom_target(corbo_optimization_headers SOURCES ${HeaderFiles})
