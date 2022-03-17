// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/uggrid.hh>

#include <dune/istl/bccsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/matrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioner.hh>
#include <dune/istl/solver.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>

// { main_begin }
int main(int argc, char *argv[])
{
  // Set up MPI
  const Dune::MPIHelper &mpiHelper = Dune::MPIHelper::instance(argc, argv);

  // Set up the grid
  constexpr int dim = 2;
  using Grid = Dune::UGGrid<dim>;
  using GridView = Grid::LeafGridView;

  std::shared_ptr<Dune::UGGrid<dim>> grid =
      Dune::GmshReader<Dune::UGGrid<dim>>::read("l-shape-refined.msh");

  std::vector<double> dataVector;
  auto gridView = grid->leafGridView();

  if (mpiHelper.rank() == 0) {
    // The initial iterate as a function
    auto initialIterate = [](auto p) { return std::min(p[0], p[1]); };

    // Sample on the grid vertices
    dataVector.resize(gridView.size(dim));
    for (const auto &vertex :
         vertices(gridView, Dune::Partitions::interiorBorder)) {
      auto index = gridView.indexSet().index(vertex);
      dataVector[index] = initialIterate(vertex.geometry().corner(0));
    }
  }

  // Copy vertex data into associative container
  std::map<Grid::LocalIdSet::IdType, double> persistentContainer;
  const auto &idSet = grid->localIdSet();

  for (const auto &vertex : vertices(gridView)) {
    persistentContainer[idSet.id(vertex)] =
        dataVector[gridView.indexSet().index(vertex)];
  }

  // Distribute the grid and the data
  LBVertexDataHandle<Grid, std::map<Grid::LocalIdSet::IdType, double>>
      dataHandle(grid, persistentContainer);

  grid->loadBalance(dataHandle);
}