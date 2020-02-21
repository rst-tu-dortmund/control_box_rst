# Library control-box-rst-gui
project(corbo-gui VERSION 0.1 LANGUAGES CXX)

## check variables

if (NOT BUILD_GUI)
   message(STATUS "Skipping corbo-gui build, since BUILD_GUI is disabled")
   return()
endif (NOT BUILD_GUI)

if (NOT RPC_SUPPORT)
   message(STATUS "Skipping corbo-gui build, since RPC_SUPPORT is disabled")
   return()
endif (NOT RPC_SUPPORT)

# Add header files as custom target in order to display them in the IDE
# TODO check for a cleaner solution
FILE(GLOB_RECURSE HeaderFiles
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")
add_custom_target(corbo_gui_headers SOURCES ${HeaderFiles})

# Find the QtWidgets library
find_package(Qt5Gui REQUIRED)
find_package(Qt5Widgets REQUIRED)
# Qt5PrintSupport is required for QCustomPlot
find_package(Qt5PrintSupport REQUIRED)

# Include QCustomPlot
set(QCUSTOMPLOT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../extern/qcustomplot" CACHE PATH "")

# Find and Configure QT
# Find includes in corresponding build directories
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
# Fixme: CMake's Automoc features does not find headers automatically yet
#        if located in another location (e.g. include/..) rather than src/).
# Workaround: Parse all header files:
qt5_wrap_cpp(MOC_CPP ${HeaderFiles})

add_library(corbo_gui_lib STATIC
    src/main_window.cpp
    src/toolbox_widget.cpp
    src/signal_widget.cpp
    src/signal_collection_widget.cpp
    src/signal_helper.cpp
    src/scope_widget.cpp
    src/scope_collection_widget.cpp
    src/parameter_widget.cpp
    src/label_edit_widget.cpp
    src/label_slider_widget.cpp
    src/label_combobox_widget.cpp
    src/label_horizontal_line_widget.cpp
    src/horizontal_button_group.cpp
    src/collapsable_groupbox.cpp
    src/one_of_param_widget.cpp
    src/extended_tree_widget.cpp
    src/parameter_cache.cpp
    src/color_manager.cpp
    ${QCUSTOMPLOT_PATH}/qcustomplot.cpp
    ${MOC_CPP}
)

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build paths.
target_include_directories(corbo_gui_lib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${QCUSTOMPLOT_PATH}>
    $<INSTALL_INTERFACE:include/control_box_rst>
    PRIVATE src)

# If we have compiler requirements for this library, list them here
target_compile_features(corbo_gui_lib
    PUBLIC cxx_auto_type cxx_range_for cxx_constexpr cxx_lambdas
    PRIVATE cxx_variadic_templates)


# Depend on a library that we defined in the top-level file
target_link_libraries(corbo_gui_lib
    corbo_communication
    Qt5::Gui
    Qt5::Widgets
    Qt5::PrintSupport
)

# mock only for library target before!
set(CMAKE_AUTOMOC OFF)

# Gui application
add_executable(corbo_gui
    src/app.cpp
)

target_link_libraries(corbo_gui
    corbo_gui_lib
)


# 'make install' to the correct location
install(TARGETS corbo_gui_lib corbo_gui EXPORT corbo_guiConfig
    ARCHIVE  DESTINATION lib/control_box_rst
    LIBRARY  DESTINATION lib/control_box_rst
    RUNTIME  DESTINATION bin/control_box_rst)  # This is for Windows
install(DIRECTORY include/ DESTINATION include/control_box_rst)

# This makes the project importable from the install directory
# Put config file in per-project dir (name MUST match), can also
# just go into <prefix>/cmake.
install(EXPORT corbo_guiConfig DESTINATION share/control_box_rst/corbo_gui)

# This makes the project importable from the build directory
export(TARGETS corbo_gui_lib corbo_gui FILE corbo_guiConfig.cmake)

# Add unit tests
if (BUILD_TESTS)
	add_executable(test_gui
	    test/gui_test.cpp)

	target_link_libraries(test_gui
	    corbo_gui_lib
	    gtest
	    #gmock
	)
	add_test(test_gui_test test_gui)
endif (BUILD_TESTS)


# copy dlls
set(QT5_PATH "${Qt5Core_DIR}/../../../")
#message(STATUS "QT5 Dir: ${QT5_PATH}")
if (WIN32 AND MSVC)
        file(COPY ${QT5_PATH}/bin/Qt5Core.dll DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Release/)
        file(COPY ${QT5_PATH}/bin/Qt5Cored.dll DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Debug/)
        file(COPY ${QT5_PATH}/bin/Qt5Gui.dll DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Release/)
        file(COPY ${QT5_PATH}/bin/Qt5Guid.dll DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Debug/)
        file(COPY ${QT5_PATH}/bin/Qt5Widgets.dll DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Release/)
        file(COPY ${QT5_PATH}/bin/Qt5Widgetsd.dll DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Debug/)
        file(COPY ${QT5_PATH}/bin/Qt5PrintSupport.dll DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Release/)
        file(COPY ${QT5_PATH}/bin/Qt5PrintSupportd.dll DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Debug/)
elseif (WIN32 AND MINGW)
        file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/extern/mingw-dlls/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
endif (WIN32 AND MSVC)

#if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
#    set_target_properties(teb_gui PROPERTIES LINK_FLAGS_DEBUG "/SUBSYSTEM:WINDOWS")
#    set_target_properties(teb_gui PROPERTIES LINK_FLAGS_RELEASE "/SUBSYSTEM:WINDOWS")
#endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
