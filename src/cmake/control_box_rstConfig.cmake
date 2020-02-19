include(CMakeFindDependencyMacro)

foreach(component ${control_box_rst_FIND_COMPONENTS})
  # For requested component, execute its "config" script
   if(EXISTS "${CMAKE_CURRENT_LIST_DIR}/corbo_${component}/corbo_${component}Config.cmake")
     # check install space structure
     include(${CMAKE_CURRENT_LIST_DIR}/corbo_${component}/corbo_${component}Config.cmake)
   else()
     # otherwise check build space structure
     include(${CMAKE_CURRENT_LIST_DIR}/src/${component}/corbo_${component}Config.cmake)
   endif()
endforeach()
