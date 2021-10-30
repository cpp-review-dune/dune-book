# set(ENABLE_BUILD_DYNAMIC ON)
# https://gitlab.onelab.info/gmsh/gmsh/-/blob/master/demos/api/CMakeLists.txt#L84
find_library(GMSH_LIB gmsh)
add_executable("hello-gmsh" hello-gmsh.cc)
target_link_libraries("hello-gmsh" ${GMSH_LIB})
add_executable("l-shape-gmsh" l-shape-gmsh.cc)
target_link_libraries("l-shape-gmsh" ${GMSH_LIB})

#add_custom_command(TARGET "l-shape-gmsh" POST_BUILD
#        COMMAND ${CMAKE_COMMAND} -E copy
#        ${CMAKE_CURRENT_LIST_DIR}/config/config.json
#        ${CMAKE_CURRENT_BINARY_DIR})