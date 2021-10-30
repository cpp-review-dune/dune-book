// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <dune/common/parallel/mpihelper.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/uggrid.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/matrixmarket.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#include "assemblePoissonProblem.hh"

int main(int argc, char **argv) {
  Dune::MPIHelper::instance(argc, argv);

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

  assemblePoissonProblem(basis, stiffnessMatrix, b, sourceTerm);

  auto predicate = [](auto x) {
    return x[0] < 1e-8 || x[1] < 1e-8 || (x[0] > 0.4999 && x[1] > 0.4999);
  };

  std::vector<bool> dirichletNodes;
  Dune::Functions::interpolate(basis, dirichletNodes, predicate);

  for (std::size_t i = 0; i < stiffnessMatrix.N(); i++) {
    if (dirichletNodes[i]) {
      auto cIt = stiffnessMatrix[i].begin();
      auto cEndIt = stiffnessMatrix[i].end();

      for (; cIt != cEndIt; ++cIt)
        *cIt = (cIt.index() == i) ? 1.0 : 0.0;
    }
  }

  auto dirichletValues = [](auto x) {
    return (x[0] < 1e-8 || x[1] < 1e-8) ? 0 : 0.5;
  };

  Dune::Functions::interpolate(basis, b, dirichletValues, dirichletNodes);

  std::string baseName = "getting-started-poisson-fem-" +
                         std::to_string(grid->maxLevel()) + "-refinements";
  Dune::storeMatrixMarket(stiffnessMatrix, baseName + "-matrix.mtx");
  Dune::storeMatrixMarket(b, baseName + "-rhs.mtx");

  Vector x(basis.size());
  x = b;

  Dune::MatrixAdapter<Matrix, Vector, Vector> linearOperator(stiffnessMatrix);

  Dune::SeqILU<Matrix, Vector, Vector> preconditioner(stiffnessMatrix, 1.0);

  Dune::CGSolver<Vector> cg(linearOperator, preconditioner, 1e-5, 50, 2);

  Dune::InverseOperatorResult statistics;

  cg.apply(x, b, statistics);

  Dune::VTKWriter<GridView> vtkWriter(gridView);
  vtkWriter.addVertexData(x, "solution");
  vtkWriter.write("getting-started-poisson-fem-result");

  return 0;
}