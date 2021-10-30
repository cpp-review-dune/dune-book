//
// Created by carlosal1015 on 10/26/21.
//

#ifndef DUNE_BOOK_ASSEMBLEELEMENTSTIFFNESSMATRIX_HH
#define DUNE_BOOK_ASSEMBLEELEMENTSTIFFNESSMATRIX_HH

#include <dune/geometry/quadraturerules.hh>

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
}

#endif // DUNE_BOOK_ASSEMBLEELEMENTSTIFFNESSMATRIX_HH
