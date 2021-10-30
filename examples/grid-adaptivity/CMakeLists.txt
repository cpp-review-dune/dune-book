add_executable("grid-adaptivity" grid-adaptivity.cc)
target_compile_features("grid-adaptivity"
        PUBLIC
        cxx_nullptr
        PRIVATE
        cxx_lambdas)
target_link_dune_default_libraries("grid-adaptivity")
target_compile_options("grid-adaptivity" PRIVATE -Werror)
