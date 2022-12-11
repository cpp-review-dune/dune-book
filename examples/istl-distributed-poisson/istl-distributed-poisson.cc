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
#include <dune/istl/solvers.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>

#include "AdditiveScalarProduct.hh"
#include "JacobiPreconditioner.hh"
#include "LBVertexDataHandle.hh"
#include "assemblePoissonProblem.hh"

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

  // Get gridView again after load-balancing, to make sure it is up-to-date
  gridView = grid->leafGridView();

  // Copy data back into the array
  dataVector.resize(gridView.size(dim));

  for (const auto &vertex : vertices(gridView)) {
    dataVector[gridView.indexSet().index(vertex)] =
        persistentContainer[idSet.id(vertex)];
  }

  /////////////////////////////////////////////////////////
  // Stiffness matrix and right hand side vector
  /////////////////////////////////////////////////////////

  using Vector = Dune::BlockVector<double>;
  using Matrix = Dune::BCRSMatrix<double>;

  Vector rhs;
  Matrix stiffnessMatrix;

  // Assemble the Poisson system in first-order Lagrange space
  Dune::Functions::LagrangeBasis<GridView, 1> basis(gridView);

  auto sourceTerm = [](const Dune::FieldVector<double, dim> &x) {
    return -5.0;
  };

  assemblePoissonProblem(basis, stiffnessMatrix, rhs, sourceTerm);

  // Determine Dirichlet degrees of freedom by marking all those
  // whose Lagrange nodes comply with a given predicate.
  auto dirichletPredicate = [](auto p) {
    return p[0] < 1e-8 || p[1] < 1e-8 || (p[0] > 0.4999 && p[1] > 0.4999);
  };

  // Interpolating the predicate will mark all Dirichlet degrees of freedom
  std::vector<bool> dirichletNodes;
  Dune::Functions::interpolate(basis, dirichletNodes, dirichletPredicate);

  /////////////////////////////////////////////////////////
  // Modify Dirichlet matrix rows
  /////////////////////////////////////////////////////////

  // Loop over the matrix rows
  for (std::size_t i = 0; i < stiffnessMatrix.N(); i++) {
    if (dirichletNodes[i]) {
      auto cIt = stiffnessMatrix[i].begin();
      auto cEndIt = stiffnessMatrix[i].end();
      // Loop over nonzero matrix entries in current row
      for (; cIt != cEndIt; ++cIt)
        *cIt = (i == cIt.index()) ? 1.0 : 0.0;
    }
  }

  // Set Dirichlet values
  for (std::size_t i = 0; i < dirichletNodes.size(); i++)
    if (dirichletNodes[i])
      rhs[i] = dataVector[i];

  Vector x(basis.size());
  std::copy(dataVector.begin(), dataVector.end(), x.begin());

  double reduction = 1e-3;
  int maxIterations = 50;

  Dune::MatrixAdapter<Matrix, Vector, Vector> linearOperator(stiffnessMatrix);
  JacobiPreconditioner<GridView, Matrix, Vector> preconditioner(
      gridView, stiffnessMatrix);
  AdditiveScalarProduct<GridView, Vector> scalarProduct(gridView);

  Dune::CGSolver<Vector> cg(linearOperator, scalarProduct, preconditioner,
                            reduction, maxIterations,
                            (mpiHelper.rank() == 0) ? 2 : 0);

  // Object storing some statistics about the solving process
  Dune::InverseOperatorResult statistics;

  // Solve!
  cg.apply(x, rhs, statistics);

  Dune::VTKWriter<GridView> vtkWriter(gridView);
  vtkWriter.addVertexData(x, "solution");
  vtkWriter.write("istl-distributed-poisson-result");
}