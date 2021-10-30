// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
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
#include <map>
#include <vector>

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

  int order = 2 * (dim * localFiniteElement.localBasis().order() - 1);
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

    for (std::size_t p = 0; p < elementMatrix.N(); p++) {
      auto localRow = localView.tree().localIndex(p);
      for (std::size_t q = 0; q < elementMatrix.M(); q++) {
        auto localCol = localView.tree().localIndex(q);
        elementMatrix[localRow][localCol] += (gradients[p] * gradients[q]) *
                                             quadPoint.weight() *
                                             integrationElement;
      }
    }
  }
}

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

  for (const auto &element : elements(gridView, Dune::Partitions::interior)) {
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

template <class Grid, class AssociativeContainer>
struct LBVertexDataHandle : public Dune::CommDataHandleIF<
                                LBVertexDataHandle<Grid, AssociativeContainer>,
                                typename AssociativeContainer::mapped_type> {
  LBVertexDataHandle(const std::shared_ptr<Grid> &grid,
                     AssociativeContainer &dataContainer)
      : idSet_(grid->localIdSet()), dataContainer_(dataContainer) {}

  bool contains(int dim, int codim) const {
    assert(dim == Grid::dimension);
    return (codim == dim);
  }

  bool fixedSize(int dim, int codim) const { return true; }

  template <class Entity> std::size_t size(const Entity &entity) const {
    return 1;
  }

  template <class MessageBuffer, class Entity>
  void gather(MessageBuffer &buffer, const Entity &entity, std::size_t n) {
    assert(n == 1);

    auto id = idSet_.id(entity);
    buffer.read(dataContainer_[id]);
  }

  template <class MessageBuffer, class Entity>
  void scatter(MessageBuffer &buffer, const Entity &entity, std::size_t n) {
    assert(n == 1);

    auto id = idSet_.id(entity);
    buffer.read(dataContainer_[id]);
  }

private:
  const typename Grid::LocalIdSet &idSet_;
  AssociativeContainer &dataContainer_;
};

template <class GridView, class Vector>
struct VertexDataUpdate
    : public Dune::CommDataHandleIF<VertexDataUpdate<GridView, Vector>,
                                    typename Vector::value_type> {
  using DataType = typename Vector::value_type;

  VertexDataUpdate(const GridView &gridView, const Vector &userDataSend,
                   Vector &userDataReceive)
      : gridView_(gridView), userDataSend_(userDataSend),
        userDataReceive_(userDataReceive) {}

  bool contains(int dim, int codim) const { return (codim == dim); }

  bool fixedSize(int dim, int codim) const { return true; }

  template <class Entity> std::size_t size(const Entity &e) const { return 1; }

  template <class MessageBuffer, class Entity>
  void gather(MessageBuffer &buffer, const Entity &entity) const {
    auto index = gridView_.indexSet().index(entity);
    buffer.write(userDataSend_[index]);
  }

  template <class MessageBuffer, class Entity>
  void scatter(MessageBuffer &buffer, const Entity &entity, std::size_t n) {
    assert(n == 1);
    DataType x;
    buffer.read(x);

    userDataReceive_[gridView_.indexSet().index(entity)] += x;
  }

private:
  const GridView gridView_;
  const Vector &userDataSend_;
  Vector &userDataReceive_;
};

// https://stackoverflow.com/a/5192091/9302545
int main(int argc, char *argv[]) {

  const Dune::MPIHelper &mpiHelper = Dune::MPIHelper::instance(argc, argv);

  constexpr int dim = 2;
  using Grid = Dune::UGGrid<dim>;
  using GridView = Grid::LeafGridView;

  std::shared_ptr<Grid> grid =
      Dune::GmshReader<Grid>::read("l-shape-refined.msh");
  auto gridView = grid->leafGridView();

  std::vector<double> dataVector;

  if (mpiHelper.rank() == 0) {
    auto initialIterate = [](auto x) { return std::min(x[0], x[1]); };
    dataVector.resize(gridView.size(dim));
    for (const auto &vertex :
         vertices(gridView, Dune::Partitions::interiorBorder)) {
      auto index = gridView.indexSet().index(vertex);
      dataVector[index] = initialIterate(vertex.geometry().corner(0));
    }
  }

  using PersistentContainer = std::map<Grid::LocalIdSet::IdType, double>;
  PersistentContainer persistentContainer;
  const auto &idSet = grid->localIdSet();

  for (const auto &vertex : vertices(gridView)) {
    persistentContainer[idSet.id(vertex)] =
        dataVector[gridView.indexSet().index(vertex)];
  }

  return 0;
}
