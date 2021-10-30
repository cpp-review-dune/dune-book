//
// Created by carlosal1015 on 10/29/21.
//

#ifndef DUNE_BOOK_ASSEMBLEPOISSONPROBLEM_HH
#define DUNE_BOOK_ASSEMBLEPOISSONPROBLEM_HH

#include "assembleElementStiffnessMatrix.hh"
#include "assembleElementVolumeTerm.hh"
#include "getOccupationPattern.hh"
#include <dune/istl/bcrsmatrix.hh>

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

#endif // DUNE_BOOK_ASSEMBLEPOISSONPROBLEM_HH
