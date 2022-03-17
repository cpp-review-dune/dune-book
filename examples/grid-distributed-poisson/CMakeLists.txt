add_executable("grid-distributed-poisson" grid-distributed-poisson.cc)
target_compile_features("grid-distributed-poisson"
        PUBLIC
        cxx_nullptr
        PRIVATE
        cxx_lambdas)
target_link_dune_default_libraries("grid-distributed-poisson")
target_compile_options("grid-distributed-poisson" PRIVATE -Werror)
dune_symlink_to_source_files(FILES l-shape-refined.msh)