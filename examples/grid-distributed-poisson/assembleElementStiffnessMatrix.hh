//
// Created by carlosal1015 on 03/17/22.
//

#pragma once

#include <dune/geometry/quadraturerules.hh>

// Compute the stiffness matrix for a single element
template <class LocalView, class Matrix>
void assembleElementStiffnessMatrix(const LocalView &localView,
                                    Matrix &elementMatrix)
{
  using Element = typename LocalView::Element;
  constexpr int dim = Element::dimension;
  auto element = localView.element();
  auto geometry = element.geometry();

  // Get set of shape functions for this element
  const auto &localFiniteElement = localView.tree().finiteElement();

  // Set all matrix entries to zero
  elementMatrix.setSize(localView.size(), localView.size());
  elementMatrix = 0; // Fill the entire matrix with zeros

  // Get a quadrature rule
  int order = 2 * (dim * localFiniteElement.localBasis().order() - 1);
  const auto &quadRule =
      Dune::QuadratureRules<double, dim>::rule(element.type(), order);

  // Loop over all quadrature points
  for (const auto &quadPoint : quadRule) {
    // Position of the current quadrature point in the reference element
    const auto quadPos = quadPoint.position();

    // The transposed inverse Jacobian of the map from the reference element
    // to the grid element
    const auto jacobian = geometry.jacobianInverseTransposed(quadPos);

    // The multiplicative factor in the integral transformation formula
    const auto integrationElement = geometry.integrationElement(quadPos);

    // The gradients of the shape functions on the reference element
    std::vector<Dune::FieldMatrix<double, 1, dim>> referenceGradients;
    localFiniteElement.localBasis().evaluateJacobian(quadPos,
                                                     referenceGradients);

    // Compute the shape function gradients on the grid element
    std::vector<Dune::FieldVector<double, dim>> gradients(
        referenceGradients.size());
    for (std::size_t i = 0; i < gradients.size(); i++)
      jacobian.mv(referenceGradients[i][0], gradients[i]);

    // Compute the actual matrix entries
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