// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <dune/common/exceptions.hh>         // We use exceptions
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/uggrid.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/matrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/matrixmarket.hh>
#include <dune/istl/preconditioner.hh>
#include <dune/istl/solver.hh>
#include <iostream>

int main(int argc, char **argv) {
  try {
    // Maybe initialize MPI
    Dune::MPIHelper &helper = Dune::MPIHelper::instance(argc, argv);
    std::cout << "Hello World! This is dune-book." << std::endl;

    // https://en.cppreference.com/w/cpp/language/constexpr
    constexpr int dim = 2;

    using Grid = Dune::UGGrid<dim>;
    std::shared_ptr<Grid> grid = Dune::GmshReader<Grid>::read("l-shape.msh");

    grid->globalRefine(2);

    using GridView = Grid::LeafGridView;
    GridView gridView = grid->leafGridView();

    using Matrix = Dune::BCRSMatrix<double>;
    using Vector = Dune::BlockVector<double>;

    Matrix stiffnessMatrix;
    Vector b;

    Dune::Functions::LagrangeBasis<GridView, 1> basis(gridView);

    auto sourceTerm = [](const Dune::FieldVector<double, dim> &x) {
      return -5.0;
    };

    if (Dune::MPIHelper::isFake)
      std::cout << "This is a sequential program." << std::endl;
    else
      std::cout << "I am rank " << helper.rank() << " of " << helper.size()
                << " processes!" << std::endl;
    return 0;
  } catch (Dune::Exception &e) {
    std::cerr << "Dune reported error: " << e << std::endl;
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}