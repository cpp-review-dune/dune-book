// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

/**
 * @file getting-started-poisson-fem.cc
 * @brief
 * Equation: Poisson equation
 * Discretization: First-order Lagrange finite elements
 * Grid: Unstructured triangle grid, implemented with UGGrid
 * Solver: CG Solver with ILU preconditioner
 * Distributed: no
 */

// Included by uggrid.hh
// #include <dune/common/parallel/mpihelper.hh>

// { include_uggrid_begin }
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/uggrid.hh>
// { include_uggrid_end }
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/matrixmarket.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>

// { using_namespace_dune_begin }
// { using_namespace_dune_end }

#include "assemblePoissonProblem.hh"

int main(int argc, char **argv)
{
  // { mpi_setup_begin }
  // Set up MPI, if available
  Dune::MPIHelper::instance(argc, argv);
  // { mpi_setup_end }

  ////////////////////////
  // Generate the grid
  ////////////////////////

  // { create_grid_begin }
  // https://en.cppreference.com/w/cpp/language/constexpr
  constexpr int dim = 2;
  using Grid = Dune::UGGrid<dim>;
  std::shared_ptr<Grid> grid = Dune::GmshReader<Grid>::read("l-shape.msh");

  grid->globalRefine(2);

  using GridView = Grid::LeafGridView;
  GridView gridView = grid->leafGridView();
  // { create_grid_end }

  ////////////////////////
  // Stiffness matrix and right-hand side vector
  ////////////////////////

  // { create_matrix_vector_begin }
  using Matrix = Dune::BCRSMatrix<double>;
  using Vector = Dune::BlockVector<double>;

  Matrix stiffnessMatrix;
  Vector b;
  // { create_matrix_vector_end }

  ////////////////////////
  // Assemble the system
  ////////////////////////

  // { setup_basis_begin }

  Dune::Functions::LagrangeBasis<GridView, 1> basis(gridView);

  auto sourceTerm = [](const Dune::FieldVector<double, dim> &x) {
    return -5.0;
  };
  // { setup_basis_end }
  // { call_assembler_begin }
  assemblePoissonProblem(basis, stiffnessMatrix, b, sourceTerm);
  // { call_assembler_end }

  // Determine Dirichlet dofs by marking all degrees of freedom whose Lagrange
  // nodes comply with a given predicate
  // { dirichlet_marking_begin }
  auto predicate = [](auto x) {
    return x[0] < 1e-8 || x[1] < 1e-8 || (x[0] > 0.4999 && x[1] > 0.4999);
  };

  // Evaluating the predicate will mark all Dirichlet degrees of freedom
  std::vector<bool> dirichletNodes;
  Dune::Functions::interpolate(basis, dirichletNodes, predicate);
  // { dirichlet_marking_end }

  ////////////////////////
  // Modify Dirichlet rows
  ////////////////////////
  // { dirichlet_matrix_modification_begin }
  // Loop over the matrix rows
  for (std::size_t i = 0; i < stiffnessMatrix.N(); i++) {
    if (dirichletNodes[i]) {
      auto cIt = stiffnessMatrix[i].begin();
      auto cEndIt = stiffnessMatrix[i].end();
      // Loop over nonzero matrix entries in current row
      for (; cIt != cEndIt; ++cIt)
        *cIt = (cIt.index() == i) ? 1.0 : 0.0;
    }
  }
  // { dirichlet_matrix_modification_end }

  // Set Dirichlet rows
  // { dirichlet_rhs_modification_begin }
  auto dirichletValues = [](auto x) {
    return (x[0] < 1e-8 || x[1] < 1e-8) ? 0 : 0.5;
  };

  Dune::Functions::interpolate(basis, b, dirichletValues, dirichletNodes);
  // { dirichlet_rhs_modification_end }

  ////////////////////////
  // Write matrix and load vector to files, to be used in later examples
  ////////////////////////
  // { matrix_rhs_writing_begin }

  std::string baseName = "getting-started-poisson-fem-" +
                         std::to_string(grid->maxLevel()) + "-refinements";
  Dune::storeMatrixMarket(stiffnessMatrix, baseName + "-matrix.mtx");
  Dune::storeMatrixMarket(b, baseName + "-rhs.mtx");
  // { matrix_rhs_writing_end }

  ////////////////////////
  // Compute solution
  ////////////////////////
  // { algebraic_solving_begin }
  // Choose an initial iterate that fulfills the Dirichlet conditions
  Vector x(basis.size());
  x = b;

  // Turn the matrix into a linear operator
  Dune::MatrixAdapter<Matrix, Vector, Vector> linearOperator(stiffnessMatrix);

  // Sequential incomplete LU decomposition as the preconditioner
  Dune::SeqILU<Matrix, Vector, Vector> preconditioner(stiffnessMatrix,
                                                      1.0); // Relaxation factor

  // Preconditioned conjugate gradient solver
  Dune::CGSolver<Vector> cg(linearOperator, preconditioner,
                            1e-5, // Desired residual reduction factor
                            50,   // Maximum number of iterations
                            2);   // Verbosity of the solver

  // Object storing some statistics about the solving process
  Dune::InverseOperatorResult statistics;

  // Solve!
  cg.apply(x, b, statistics);
  // { algebraic_solving_end }

  // Output result
  // { vtk_output_begin }
  Dune::VTKWriter<GridView> vtkWriter(gridView);
  vtkWriter.addVertexData(x, "solution");
  vtkWriter.write("getting-started-poisson-fem-result");
  // { vtk_output_end }

  return 0;
}