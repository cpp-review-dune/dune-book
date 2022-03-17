#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <iostream>

int main(int argc, const char **argv)
{
  // Show metada from dune.module
  std::cout << "Module name: " << PACKAGE_NAME << std::endl
            << "Module maintainer: " << PACKAGE_BUGREPORT << std::endl
            << "Module version: " << PACKAGE_VERSION << std::endl;

  // Show dune-module version
  std::cout << "dune-common-version: " << DUNE_COMMON_VERSION << std::endl
            << "dune-geometry-version: " << DUNE_GEOMETRY_VERSION << std::endl
            << "dune-grid-version: " << DUNE_GRID_VERSION << std::endl
            << "dune-uggrid-version: " << DUNE_UGGRID_VERSION << std::endl
            << "dune-istl-version: " << DUNE_ISTL_VERSION << std::endl
            << "dune-typetree-version: " << DUNE_TYPETREE_VERSION << std::endl
            << "dune-localfunctions-version: " << DUNE_LOCALFUNCTIONS_VERSION
            << std::endl
            << "dune-functions-version: " << DUNE_FUNCTIONS_VERSION << std::endl
            << "dune-pdelab-version: " << DUNE_PDELAB_VERSION << std::endl;

  return 0;
}