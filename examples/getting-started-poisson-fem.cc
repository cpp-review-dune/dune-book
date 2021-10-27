// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <dune/common/parallel/mpihelper.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/uggrid.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/matrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/matrixmarket.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
//#include <iostream>
//#include "assembleElementStiffnessMatrix.hh"

template <class LocalView, class Matrix>
void assembleElementStiffnessMatrix(const LocalView &localView,
                                    Matrix &elementMatrix) {
  using Element = typename LocalView::Element;
  constexpr int dim = Element::dimension;
  auto element = localView.element();
  auto geometry = element.geometry();
  const auto &localFiniteElement = localView.tree().finiteElement();
  elementMatrix.setSize(localView.size(), localView.size());
  elementMatrix = 0;

  int order = 2 * (localFiniteElement.localBasis().order() - 1);
  const auto &quadRule =
      Dune::QuadratureRules<double, dim>::rule(element.type(), order);
  for (const auto &quadPoint : quadRule) {
    const auto quadPos = quadPoint.position();
    const auto jacobian = geometry.jacobianInverseTransposed(quadPos);
    const auto integrationElement = geometry.integrationElement(quadPos);
    std::vector<Dune::FieldMatrix<double, 1, dim>> referenceGradients;
    localFiniteElement.localBasis().evaluateJacobian(quadPos,
                                                     referenceGradients);
    std::vector<Dune::FieldVector<double, dim>> gradients(
        referenceGradients.size());
    for (std::size_t i = 0; i < gradients.size(); i++)
      jacobian.mv(referenceGradients[i][0], gradients[i]);

    for (std::size_t p = 0; p < elementMatrix.M(); p++) {
      auto localRow = localView.tree().localIndex(p);
      for (std::size_t q = 0; q < elementMatrix.M(); q++) {
        auto localCol = localView.tree().localIndex(q);
        elementMatrix[localRow][localCol] += (gradients[p] * gradients[q]) *
                                             quadPoint.weight() *
                                             integrationElement;
      }
    }
  }
};

template <class LocalView>
void assembleElementVolumeTerm(
    const LocalView &localView, Dune::BlockVector<double> &localB,
    const std::function<
        double(Dune::FieldVector<double, LocalView::Element::dimension>)>
        volumeTerm) {
  using Element = typename LocalView::Element;
  auto element = localView.element();
  constexpr int dim = Element::dimension;

  const auto &localFiniteElement = localView.tree().finiteElement();

  localB.resize(localFiniteElement.size());
  localB = 0;

  int order = dim;
  const auto &quadRule =
      Dune::QuadratureRules<double, dim>::rule(element.type(), order);

  for (const auto &quadPoint : quadRule) {
    const Dune::FieldVector<double, dim> &quadPos = quadPoint.position();
    const double integrationElement =
        element.geometry().integrationElement(quadPos);

    double functionValue = volumeTerm(element.geometry().global(quadPos));

    std::vector<Dune::FieldVector<double, 1>> shapeFunctionValues;
    localFiniteElement.localBasis().evaluateFunction(quadPos,
                                                     shapeFunctionValues);

    for (std::size_t p = 0; p < localB.size(); p++) {
      auto localIndex = localView.tree().localIndex(p);
      localB[localIndex] += shapeFunctionValues[p] * functionValue *
                            quadPoint.weight() * integrationElement;
    }
  }
}

template <class Basis>
void getOccupationPattern(const Basis &basis, Dune::MatrixIndexSet &nb) {
  nb.resize(basis.size(), basis.size());
  auto gridView = basis.gridView();
  auto localView = basis.localView();
  for (const auto &element : elements(gridView)) {
    localView.bind(element);
    for (std::size_t i = 0; i < localView.size(); i++) {
      auto row = localView.index(i);
      for (std::size_t j = 0; j < localView.size(); j++) {
        auto col = localView.index(j);
        nb.add(row, col);
      }
    }
  }
}

template <class Basis>
void assemblePoissonProblem(
    const Basis &basis, Dune::BCRSMatrix<double> &matrix,
    Dune::BlockVector<double> &b,
    const std::function<
        double(Dune::FieldVector<double, Basis::GridView::dimension>)>
        volumeTerm) {
  auto gridView = basis.gridView();
  Dune::MatrixIndexSet occupationPattern;
  getOccupationPattern(basis, occupationPattern);
  occupationPattern.exportIdx(matrix);
  matrix = 0;
  b.resize(basis.dimension());
  b = 0;
  auto localView = basis.localView();

  for (const auto &element : elements(gridView)) {

    localView.bind(element);
    Dune::Matrix<double> elementMatrix;
    assembleElementStiffnessMatrix(localView, elementMatrix);
    for (std::size_t p = 0; p < elementMatrix.N(); p++) {
      auto row = localView.index(p);
      for (std::size_t q = 0; q < elementMatrix.M(); q++) {
        auto col = localView.index(q);
        matrix[row][col] += elementMatrix[p][q];
      }
    }
    Dune::BlockVector<double> localB;
    assembleElementVolumeTerm(localView, localB, volumeTerm);

    for (std::size_t p = 0; p < localB.size(); p++) {
      auto row = localView.index(p);
      b[row] += localB[p];
    }
  }
}

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
