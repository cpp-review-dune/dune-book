add_executable("getting-started-transport-fv"
        getting-started-transport-fv.cc
        evolve.hh)
target_compile_features("getting-started-transport-fv"
        PUBLIC
        cxx_nullptr
        PRIVATE
        cxx_lambdas)
target_link_dune_default_libraries("getting-started-transport-fv")
target_compile_options("getting-started-transport-fv" PRIVATE -Werror)