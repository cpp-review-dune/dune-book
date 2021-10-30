//
// Created by carlosal1015 on 10/29/21.
//

#ifndef DUNE_BOOK_GETOCCUPATIONPATTERN_HH
#define DUNE_BOOK_GETOCCUPATIONPATTERN_HH

#include <dune/istl/matrixindexset.hh>

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

#endif // DUNE_BOOK_GETOCCUPATIONPATTERN_HH
