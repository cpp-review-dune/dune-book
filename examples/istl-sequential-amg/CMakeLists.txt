add_executable("istl-sequential-amg"
        istl-sequential-amg.cc
        assemblePoissonProblem.hh)
target_compile_features("istl-sequential-amg"
        PUBLIC
        cxx_nullptr
        PRIVATE
        cxx_lambdas)
# cxx_std_17
target_link_dune_default_libraries("istl-sequential-amg")
target_compile_options("istl-sequential-amg" PRIVATE -Werror)
dune_symlink_to_source_files(FILES l-shape-refined.msh)