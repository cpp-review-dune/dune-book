add_executable("grid-adaptivity-data-transfer"
        grid-adaptivity-data-transfer.cc
        interpolate.hh)
target_compile_features("grid-adaptivity-data-transfer"
        PUBLIC
        cxx_nullptr
        PRIVATE
        cxx_lambdas)
target_link_dune_default_libraries("grid-adaptivity-data-transfer")
target_compile_options("grid-adaptivity-data-transfer" PRIVATE -Werror)
