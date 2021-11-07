//
// Created by carlosal1015 on 10/26/21.
//

#ifndef DUNE_BOOK_ASSEMBLEELEMENTSTIFFNESSMATRIX_HH
#define DUNE_BOOK_ASSEMBLEELEMENTSTIFFNESSMATRIX_HH

#include <dune/geometry/quadraturerules.hh>

// Compute the stiffness matrix for a single element
// { local_assembler_signature_begin }
template <class LocalView, class Matrix>
void assembleElementStiffnessMatrix(const LocalView &localView,
                                    Matrix &elementMatrix)
// { local_assembler_signature_end }
{
  // { local_assembler_get_geometry_begin }
  using Element = typename LocalView::Element;
  constexpr int dim = Element::dimension;
  auto element = localView.element();
  auto geometry = element.geometry();
  // { local_assembler_get_geometry_end }

  // Get set of shape functions for this element
  // { get_shapefunctions_begin }
  const auto &localFiniteElement = localView.tree().finiteElement();
  // { get_shapefunctions_end }

  // Set all matrix entries to zero
  // { init_element_matrix_begin }
  elementMatrix.setSize(localView.size(), localView.size());
  elementMatrix = 0; // Fill the entire matrix with zeros
  // { init_element_matrix_end }

  // Get a quadrature rule
  // { get_quadrature_begin }
  int order = 2 * (localFiniteElement.localBasis().order() - 1);
  const auto &quadRule =
      Dune::QuadratureRules<double, dim>::rule(element.type(), order);
  // { get_quadrature_end }

  // Loop over all quadrature points
  // { loop_over_quad_points_begin }
  for (const auto &quadPoint : quadRule) {
    // { loop_over_quad_points_end }

    // { get_quad_points_info_begin }
    // Position of the current quadrature point in the reference element
    const auto quadPos = quadPoint.position();

    // The transposed inverse Jacobian of the map from the reference element
    // to the grid elements
    const auto jacobian = geometry.jacobianInverseTransposed(quadPos);

    // The determinant term in the integral transformation formula
    // { get_quad_points_info_end }

    // { compute_gradients_begin }
    const auto integrationElement = geometry.integrationElement(quadPos);
    std::vector<Dune::FieldMatrix<double, 1, dim>> referenceGradients;
    localFiniteElement.localBasis().evaluateJacobian(quadPos,
                                                     referenceGradients);

    // Compute the shape function gradients on the grid element
    std::vector<Dune::FieldVector<double, dim>> gradients(
        referenceGradients.size());
    for (std::size_t i = 0; i < gradients.size(); i++)
      jacobian.mv(referenceGradients[i][0], gradients[i]);
    // { compute_gradients_end }

    // Compute the actual matrix entries
    // { compute_matrix_entries_begin }
    for (std::size_t p = 0; p < elementMatrix.M(); p++) {
      auto localRow = localView.tree().localIndex(p);
      for (std::size_t q = 0; q < elementMatrix.M(); q++) {
        auto localCol = localView.tree().localIndex(q);
        elementMatrix[localRow][localCol] += (gradients[p] * gradients[q]) *
                                             quadPoint.weight() *
                                             integrationElement;
      }
    }
    // { compute_matrix_entries_end }
  }
}

#endif // DUNE_BOOK_ASSEMBLEELEMENTSTIFFNESSMATRIX_HH
