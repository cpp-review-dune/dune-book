# https://aur.archlinux.org/cgit/aur.git/tree/PKGBUILD?h=libmesh
# https://aur.archlinux.org/cgit/aur.git/tree/PKGBUILD?h=mumps
# set(ENV{OMPI_MCA_opal_cuda_support} "0")
# set(ENV{OMPI_MCA_opal_warn_on_missing_libcuda} "0")
# https://stackoverflow.com/a/33630631/9302545
# http://www.idris.fr/eng/simulato/gromacs-eng.html
# https://cgold.readthedocs.io/en/latest/tutorials/variables/environment.html
add_executable("hello-dune" hello-dune.cc)
target_link_dune_default_libraries("hello-dune")

# set(ENABLE_BUILD_DYNAMIC ON)
# https://gitlab.onelab.info/gmsh/gmsh/-/blob/master/demos/api/CMakeLists.txt#L84
find_library(GMSH_LIB gmsh)
add_executable("hello-gmsh" hello-gmsh.cc)
target_link_libraries("hello-gmsh" ${GMSH_LIB})

