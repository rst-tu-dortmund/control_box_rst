# Library control-box-rst-numerics
project(corbo-numerics VERSION 0.1 LANGUAGES CXX)

add_library(corbo_numerics STATIC
    src/algebraic_riccati_continuous.cpp
    src/algebraic_riccati_discrete.cpp
    src/controllability.cpp
    src/observability.cpp
    src/finite_differences.cpp
    src/lyapunov_continuous.cpp
    src/lyapunov_discrete.cpp
    src/sylvester_continuous.cpp
    src/sylvester_discrete.cpp
    src/schur.cpp
)

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build paths.
target_include_directories(corbo_numerics PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/control_box_rst>
    PRIVATE src)

# If we have compiler requirements for this library, list them here
target_compile_features(corbo_numerics
    PUBLIC cxx_auto_type cxx_range_for cxx_constexpr cxx_lambdas
    PRIVATE cxx_variadic_templates)

# Set compiler definitions
# target_compile_definitions(controllerslib PRIVATE MYDEF=${BLABLA})
# Set compiler optoins/flags
# target_compile_options(controllerslib PUBLIC -fno-elide-constructors)


# Depend on a library that we defined in the top-level file
target_link_libraries(corbo_numerics
    corbo_core
    corbo_systems
)


# 'make install' to the correct location
install(TARGETS corbo_numerics EXPORT corbo_numericsConfig
    ARCHIVE  DESTINATION lib/control_box_rst
    LIBRARY  DESTINATION lib/control_box_rst
    RUNTIME  DESTINATION bin/control_box_rst)  # This is for Windows
install(DIRECTORY include/ DESTINATION include/control_box_rst)

# This makes the project importable from the install directory
# Put config file in per-project dir (name MUST match), can also
# just go into <prefix>/cmake.
install(EXPORT corbo_numericsConfig DESTINATION share/control_box_rst/corbo_numerics)

# This makes the project importable from the build directory
export(TARGETS corbo_numerics FILE corbo_numericsConfig.cmake)

# Add unit tests
if (BUILD_TESTS)
	add_executable(test_numerics
	    test/test_main.cpp
            test/test_finite_differences.cpp
            test/test_algebraic_riccati_continuous.cpp
            test/test_algebraic_riccati_discrete.cpp
            test/test_lyapunov_continuous.cpp
            test/test_lyapunov_discrete.cpp
            test/test_sylvester_continuous.cpp
            test/test_sylvester_discrete.cpp
            test/test_observability.cpp
            test/test_controllability.cpp
            test/test_schur.cpp
            test/test_matrix_utilities.cpp
)

	target_link_libraries(test_numerics
	    corbo_numerics
	    gtest
	    #gmock
	)
        add_test(test_numerics_test test_numerics)
endif (BUILD_TESTS)

# Add header files as custom target in order to display them in the IDE
# TODO check for a cleaner solution
FILE(GLOB_RECURSE HeaderFiles
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")
add_custom_target(corbo_numerics_headers SOURCES ${HeaderFiles})
