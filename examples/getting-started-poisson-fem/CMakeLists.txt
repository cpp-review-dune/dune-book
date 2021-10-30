add_executable("getting-started-poisson-fem"
        getting-started-poisson-fem.cc
        assembleElementStiffnessMatrix.hh
        assembleElementVolumeTerm.hh
        getOccupationPattern.hh
        assemblePoissonProblem.hh)
target_compile_features("getting-started-poisson-fem"
        PUBLIC
        cxx_nullptr
        PRIVATE
        cxx_lambdas)
# cxx_std_17
target_link_dune_default_libraries("getting-started-poisson-fem")
target_compile_options("getting-started-poisson-fem" PRIVATE -Werror)
dune_symlink_to_source_files(FILES l-shape.msh)