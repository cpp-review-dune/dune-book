add_executable("istl-distributed-poisson"
        istl-distributed-poisson.cc
        assemblePoissonProblem.hh)
target_compile_features("istl-distributed-poisson"
        PUBLIC
        cxx_nullptr
        PRIVATE
        cxx_lambdas)
# cxx_std_17
target_link_dune_default_libraries("istl-distributed-poisson")
target_compile_options("istl-distributed-poisson" PRIVATE -Werror)
dune_symlink_to_source_files(FILES l-shape-refined.msh)