// Assemble the Laplace stiffness matrix on the given grid view
template <class Basis>
void assemblePoissonProblem(
    const Basis &basis, BCRSMatrix<double> &matrix,
    Dune::BlockVector<double> &b,
    const std::function<
        double(Dune::FieldVector<double, Basis::GridView::dimension>)>
        volumeTerm)
{
  auto gridView = basis.gridView();

  // MatrixIndexSets store the occupation pattern of a sparse matrix.
  // They are not particularly efficient, but simple to use.
  Dune::MatrixIndexSet occupationPattern;
  getOccupationPattern(basis, occupationPattern);

  // ... and give it the occupation pattern we want.
  occupationPattern.exportIdx(matrix);

  // Set all entries to zero
  matrix = 0;

  // Set b to correct length
  b.resize(basis.dimension());

  // Set all entries to zero
  b = 0;

  // A loop over all elements of the grid
  auto localView = basis.localView();

  for (const auto &element : elements(gridView, Dune::Partitions::interior)) {
    // Now let's get the element stiffness matrix
    // A dense matrix is used for the element stiffness matrix
    localView.bind(element);

    Dune::Matrix<double> elementMatrix;
    assembleElementStiffnessMatrix(localView, elementMatrix);

    for (std::size_t p = 0; p < elementMatrix.N(); p++) {
      // The global index of the p-th degree of freedom of the element
      auto row = localView.index(p);

      for (std::size_t q = 0; q < elementMatrix.M(); q++) {
        // The global index of the q-th degree of freedom of the element
        auto col = localView.index(q);
        matrix[row][col] += elementMatrix[p][q];
      }
    }

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