//
// Created by carlosal1015 on 03/17/22.
//

#pragma once

#include <dune/istl/matrixindexset.hh>

// Get the occupation pattern of the stiffness matrix
template <class Basis>
void getOccupationPattern(const Basis &basis, Dune::MatrixIndexSet &nb)
{
  nb.resize(basis.size(), basis.size());

  auto gridView = basis.gridView();

  // A loop over all elements of the grid
  auto localView = basis.localView();

  for (const auto &element : elements(gridView)) {
    localView.bind(element);

    for (std::size_t i = 0; i < localView.size(); i++) {
      // The global index of the i-th vertex of the element
      auto row = localView.index(i);

      for (std::size_t j = 0; j < localView.size(); j++) {
        // The global index of the j-th vertex of the element
        auto col = localView.index(j);
        nb.add(row, col);
      }
    }
  }
}