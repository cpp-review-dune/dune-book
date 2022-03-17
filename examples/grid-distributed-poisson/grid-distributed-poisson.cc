// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/io/file/gmshreader.hh>
// #include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>

#include "LBVertexDataHandle.hh"
#include "assemblePoissonProblem.hh"
#include <dune/grid/uggrid.hh>

// https://stackoverflow.com/a/5192091/9302545

// { main_begin }
int main(int argc, char *argv[])
{

  // Set up MPI
  const Dune::MPIHelper &mpiHelper = Dune::MPIHelper::instance(argc, argv);
  // { mpihelper_end }

  ///////////////////////////
  // Generate the grid
  ///////////////////////////

  // { create_grid_begin }
  constexpr int dim = 2;
  using Grid = Dune::UGGrid<dim>;
  using GridView = Grid::LeafGridView;

  std::shared_ptr<Grid> grid =
      Dune::GmshReader<Grid>::read("l-shape-refined.msh");
  auto gridView = grid->leafGridView();
  // { create_grid_end }

  // { sample_initial_iterate_begin }
  std::vector<double> dataVector;

  if (mpiHelper.rank() == 0) {
    // The initial iterate as a function
    auto initialIterate = [](auto x) { return std::min(x[0], x[1]); };

    // Sample on the grid vertices
    dataVector.resize(gridView.size(dim));
    for (const auto &vertex :
         vertices(gridView, Dune::Partitions::interiorBorder)) {
      auto index = gridView.indexSet().index(vertex);
      dataVector[index] = initialIterate(vertex.geometry().corner(0));
    }
  }
  // { sample_initial_iterate_end }

  // { data_into_map_begin }
  // Copy vertex data into associative container
  using PersistentContainer = std::map<Grid::LocalIdSet::IdType, double>;
  PersistentContainer persistentContainer;
  const auto &idSet = grid->localIdSet();

  for (const auto &vertex : vertices(gridView)) {
    persistentContainer[idSet.id(vertex)] =
        dataVector[gridView.indexSet().index(vertex)];
  }
  // { data_into_map_end }

  // { load_balancing_begin }
  // Distribute the grid and the data
  LBVertexDataHandle<Grid, PersistentContainer> dataHandle(grid,
                                                           persistentContainer);
  grid->loadBalance(dataHandle);
  // { load_balancing_end }

  // { data_from_map_begin }
  // Get gridView again after load-balancing, to make sure it is up-to-date
  gridView = grid->leafGridView();

  // Copy data back into the array
  dataVector.resize(gridView.size(dim));

  for (const auto &vertex : vertices(gridView)) {
    dataVector[gridView.indexSet().index(vertex)] =
        persistentContainer[idSet.id(vertex)];
  }
  // { data_from_map_end }

  /////////////////////////////////////////////////////////
  // Stiffness matrix and right hand side vector
  /////////////////////////////////////////////////////////

  // { create_matrix_vector_begin }
  using Matrix = Dune::BCRSMatrix<double>;
  using Vector = Dune::BlockVector<double>;

  Matrix stiffnessMatrix;
  Vector b;

  auto sourceTerm = [](const Dune::FieldVector<double, dim> &x) {
    return -5.0;
  };

  // Assemble the Poisson system in first-oder lagrange space
  Dune::Functions::LagrangeBasis<GridView, 1> basis(gridView);
  assemblePoissonProblem(basis, stiffnessMatrix, b, sourceTerm);
  // { call_assembler_end }

  // Obtain a consistent representation of the matrix diagonal
  // { make_consistent_diagonal_begin }
  Vector diagonal(basis.size());
  for (std::size_t i = 0; i < basis.size(); ++i) {
    diagonal[i] = stiffnessMatrix[i][i];
  }

  auto consistentDiagonal = diagonal;

  return 0;
}
