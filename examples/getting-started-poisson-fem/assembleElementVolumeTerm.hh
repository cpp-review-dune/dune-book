//
// Created by carlosal1015 on 10/29/21.
//

#ifndef DUNE_BOOK_ASSEMBLEELEMENTVOLUMETERM_HH
#define DUNE_BOOK_ASSEMBLEELEMENTVOLUMETERM_HH

#include <dune/geometry/quadraturerules.hh>
#include <dune/istl/bvector.hh>

// Compute the source term for a single element
template <class LocalView>
void assembleElementVolumeTerm(
    const LocalView &localView, Dune::BlockVector<double> &localB,
    const std::function<
        double(Dune::FieldVector<double, LocalView::Element::dimension>)>
        volumeTerm)
{
  using Element = typename LocalView::Element;
  auto element = localView.element();
  constexpr int dim = Element::dimension;

  // Set a shape functions for a single element
  const auto &localFiniteElement = localView.tree().finiteElement();

  // Set all entries to zero
  localB.resize(localFiniteElement.size());
  localB = 0;

  // A quadrature rule
  int order = dim;
  const auto &quadRule =
      Dune::QuadratureRules<double, dim>::rule(element.type(), order);

  // Loop over all quadrature points
  for (const auto &quadPoint : quadRule) {
    // Position of the current quadrature point in the reference element
    const Dune::FieldVector<double, dim> &quadPos = quadPoint.position();

    // The multiplicative factor in the integral transformation formula
    const double integrationElement =
        element.geometry().integrationElement(quadPos);

    double functionValue = volumeTerm(element.geometry().global(quadPos));

    // Evaluate all shape function values at this point
    std::vector<Dune::FieldVector<double, 1>> shapeFunctionValues;
    localFiniteElement.localBasis().evaluateFunction(quadPos,
                                                     shapeFunctionValues);

    // Actually compute the vector entries
    for (std::size_t p = 0; p < localB.size(); p++) {
      auto localIndex = localView.tree().localIndex(p);
      localB[localIndex] += shapeFunctionValues[p] * functionValue *
                            quadPoint.weight() * integrationElement;
    }
  }
}

#endif // DUNE_BOOK_ASSEMBLEELEMENTVOLUMETERM_HH
