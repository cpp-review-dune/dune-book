//
// Created by carlosal1015 on 10/29/21.
//

#ifndef DUNE_BOOK_ASSEMBLEPOISSONPROBLEM_HH
#define DUNE_BOOK_ASSEMBLEPOISSONPROBLEM_HH

#include "assembleElementStiffnessMatrix.hh"
#include "assembleElementVolumeTerm.hh"
#include "getOccupationPattern.hh"

// { include_matrix_begin }
#include <dune/istl/bcrsmatrix.hh>
// Included by bcrsmatrix.hh
// #include <dune/istl/bvector.hh>
// { include_matrix_end }
/** \brief Assemble the Laplace stiffness matrix on the given grid view */
// { global_assembler_signature_begin }
template <class Basis>
void assemblePoissonProblem(
    const Basis &basis, Dune::BCRSMatrix<double> &matrix,
    Dune::BlockVector<double> &b,
    const std::function<
        double(Dune::FieldVector<double, Basis::GridView::dimension>)>
        volumeTerm)
// { global_assembler_signature_end }
{
  // { assembler_get_grid_info_begin }
  auto gridView = basis.gridView();
  // { assembler_get_grid_info_end }

  // MatrixIndexSets store the occupation pattern of a sparse matrix.
  // They are not particularly efficient, but simple to use.
  // { assembler_matrix_pattern_begin }
  Dune::MatrixIndexSet occupationPattern;
  getOccupationPattern(basis, occupationPattern);
  occupationPattern.exportIdx(matrix);
  // { assembler_matrix_pattern_end }

  // Set all entries to zero
  // { assembler_zero_matrix_begin }
  matrix = 0;
  // { assembler_zero_matrix_end }

  // { assembler_zero_vector_begin }
  // Set b to correct length
  b.resize(basis.dimension());

  // Set all entries to zero
  b = 0;
  // { assembler_zero_vector_end }

  // A loop over all elements of the grid
  // { assembler_element_loop_begin }
  auto localView = basis.localView();

  for (const auto &element : elements(gridView)) {
    // { assembler_element_loop_end }

    // Now let's get the element stiffness matrix
    // A dense matrix is used for the element stiffness matrix
    // { assembler_assemble_element_matrix_begin }
    localView.bind(element);

    Dune::Matrix<double> elementMatrix;
    assembleElementStiffnessMatrix(localView, elementMatrix);
    // { assembler_assemble_element_matrix_end }

    // { assembler_add_element_matrix_begin }
    for (std::size_t p = 0; p < elementMatrix.N(); p++) {
      // The global index of the p-th degree of freedom of the element
      auto row = localView.index(p);

      for (std::size_t q = 0; q < elementMatrix.M(); q++) {
        // The global index of the q-th degree of freedom of the element
        auto col = localView.index(q);
        matrix[row][col] += elementMatrix[p][q];
      }
    }
    // { assembler_add_element_matrix_end }

    // Now get the local contribution to the right-hand side vector
    Dune::BlockVector<double> localB;
    assembleElementVolumeTerm(localView, localB, volumeTerm);

    for (std::size_t p = 0; p < localB.size(); p++) {
      // The global index of the p-th vertex of the element
      auto row = localView.index(p);
      b[row] += localB[p];
    }
  }
}

#endif // DUNE_BOOK_ASSEMBLEPOISSONPROBLEM_HH
