# https://aur.archlinux.org/cgit/aur.git/tree/PKGBUILD?h=libmesh
# https://aur.archlinux.org/cgit/aur.git/tree/PKGBUILD?h=mumps
# set(ENV{OMPI_MCA_opal_cuda_support} "0")
# set(ENV{OMPI_MCA_opal_warn_on_missing_libcuda} "0")
# https://stackoverflow.com/a/33630631/9302545
# http://www.idris.fr/eng/simulato/gromacs-eng.html
# https://cgold.readthedocs.io/en/latest/tutorials/variables/environment.html
add_executable("hello-dune" hello-dune.cc)
target_link_dune_default_libraries("hello-dune")

add_executable("recipe-integration" recipe-integration.cc)
target_link_dune_default_libraries("recipe-integration")
add_executable("recipe-iterate-over-grid" recipe-iterate-over-grid.cc)
target_link_dune_default_libraries("recipe-iterate-over-grid")

add_executable("getting-started-poisson-fem" getting-started-poisson-fem.cc) # assembleElementStiffnessMatrix.hh
target_compile_features("getting-started-poisson-fem"
        PUBLIC
            cxx_nullptr
        PRIVATE
            cxx_lambdas)
# cxx_std_17
target_link_dune_default_libraries("getting-started-poisson-fem")
target_compile_options("getting-started-poisson-fem" PRIVATE -Werror)
dune_symlink_to_source_files(FILES l-shape.msh)

add_executable("getting-started-transport-fv" getting-started-transport-fv.cc)
target_compile_features("getting-started-transport-fv"
        PUBLIC
            cxx_nullptr
        PRIVATE
            cxx_lambdas)
target_link_dune_default_libraries("getting-started-transport-fv")
target_compile_options("getting-started-transport-fv" PRIVATE -Werror)

add_executable("grid-adaptivity" grid-adaptivity.cc)
target_compile_features("grid-adaptivity"
        PUBLIC
            cxx_nullptr
        PRIVATE
            cxx_lambdas)
target_link_dune_default_libraries("grid-adaptivity")
target_compile_options("grid-adaptivity" PRIVATE -Werror)

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
