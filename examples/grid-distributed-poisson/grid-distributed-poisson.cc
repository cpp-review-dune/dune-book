// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>

#include "LBVertexDataHandle.hh"
#include "VertexDataUpdate.hh"
#include "assemblePoissonProblem.hh"
#include <dune/grid/uggrid.hh>
#include <iomanip>

// https://stackoverflow.com/a/5192091/9302545

// { main_begin }
int main(int argc, char *argv[])
{

  // Set up MPI
  const Dune::MPIHelper &mpiHelper = Dune::MPIHelper::instance(argc, argv);
  // { mpihelper_end }

  ///////////////////////////
  // Generate the grid
  ///////////////////////////

  // { create_grid_begin }
  constexpr int dim = 2;
  using Grid = Dune::UGGrid<dim>;
  using GridView = Grid::LeafGridView;

  std::shared_ptr<Grid> grid =
      Dune::GmshReader<Grid>::read("l-shape-refined.msh");
  auto gridView = grid->leafGridView();
  // { create_grid_end }

  // { sample_initial_iterate_begin }
  std::vector<double> dataVector;

  if (mpiHelper.rank() == 0) {
    // The initial iterate as a function
    auto initialIterate = [](auto x) { return std::min(x[0], x[1]); };

    // Sample on the grid vertices
    dataVector.resize(gridView.size(dim));
    for (const auto &vertex :
         vertices(gridView, Dune::Partitions::interiorBorder)) {
      auto index = gridView.indexSet().index(vertex);
      dataVector[index] = initialIterate(vertex.geometry().corner(0));
    }
  }
  // { sample_initial_iterate_end }

  // { data_into_map_begin }
  // Copy vertex data into associative container
  using PersistentContainer = std::map<Grid::LocalIdSet::IdType, double>;
  PersistentContainer persistentContainer;
  const auto &idSet = grid->localIdSet();

  for (const auto &vertex : vertices(gridView)) {
    persistentContainer[idSet.id(vertex)] =
        dataVector[gridView.indexSet().index(vertex)];
  }
  // { data_into_map_end }

  // { load_balancing_begin }
  // Distribute the grid and the data
  LBVertexDataHandle<Grid, PersistentContainer> dataHandle(grid,
                                                           persistentContainer);
  grid->loadBalance(dataHandle);
  // { load_balancing_end }

  // { data_from_map_begin }
  // Get gridView again after load-balancing, to make sure it is up-to-date
  gridView = grid->leafGridView();

  // Copy data back into the array
  dataVector.resize(gridView.size(dim));

  for (const auto &vertex : vertices(gridView)) {
    dataVector[gridView.indexSet().index(vertex)] =
        persistentContainer[idSet.id(vertex)];
  }
  // { data_from_map_end }

  /////////////////////////////////////////////////////////
  // Stiffness matrix and right hand side vector
  /////////////////////////////////////////////////////////

  // { create_matrix_vector_begin }
  using Matrix = Dune::BCRSMatrix<double>;
  using Vector = Dune::BlockVector<double>;

  Matrix stiffnessMatrix;
  Vector b;

  auto sourceTerm = [](const Dune::FieldVector<double, dim> &x) {
    return -5.0;
  };

  // Assemble the Poisson system in first-order lagrange space
  Dune::Functions::LagrangeBasis<GridView, 1> basis(gridView);
  assemblePoissonProblem(basis, stiffnessMatrix, b, sourceTerm);
  // { call_assembler_end }

  // Obtain a consistent representation of the matrix diagonal
  // { make_consistent_diagonal_begin }
  Vector diagonal(basis.size());
  for (std::size_t i = 0; i < basis.size(); ++i) {
    diagonal[i] = stiffnessMatrix[i][i];
  }

  auto consistentDiagonal = diagonal;
  VertexDataUpdate<GridView, Vector> matrixDataHandle(gridView, diagonal,
                                                      consistentDiagonal);

  gridView.communicate(matrixDataHandle,
                       Dune::InteriorBorder_InteriorBorder_Interface,
                       Dune::ForwardCommunication);
  // { make_consistent_diagonal_end }

  // Determine Dirichlet degrees of freedom by marking all degrees
  // of freedom whose Lagrange nodes comply with a given predicate.
  // { dirichlet_marking_begin }
  auto dirichletPredicate = [](auto p) {
    return p[0] < 1e-8 || p[1] < 1e-8 || (p[0] > 0.4999 && p[1] > 0.4999);
  };

  // Interpolating the predicate with mark
  // all desired Dirichlet degrees of freedom
  std::vector<bool> dirichletNodes;
  Dune::Functions::interpolate(basis, dirichletNodes, dirichletPredicate);
  // { dirichlet_marking_end }

  /////////////////////////////////////////////////////////
  // Modify Dirichlet matrix rows and load vector entries
  /////////////////////////////////////////////////////////

  // { dirichlet_modification_begin }
  // Loop over the matrix rows
  for (std::size_t i = 0; i < stiffnessMatrix.N(); i++) {
    if (dirichletNodes[i]) {
      auto cIt = stiffnessMatrix[i].begin();
      auto cEndIt = stiffnessMatrix[i].end();
      // Loop over nonzero matrix entries in current row
      for (; cIt != cEndIt; ++cIt) {
        *cIt = (i == cIt.index()) ? 1.0 : 0.0;
      }
      // Modify corresponding load vector entry
      b[i] = dataVector[i];
    }
  }
  // { dirichlet_modification_end }

  /////////////////////////////////////////////////////////
  // Compute solution
  /////////////////////////////////////////////////////////

  // { algebraic_solving_preprocess_begin }
  // Set the initial iterate
  Vector x(basis.size());
  std::copy(dataVector.begin(), dataVector.end(), x.begin());

  // Solver parameters
  double reduction = 1e-3; // Desired residual reduction factor
  int maxIterations = 50;  // Maximum number of iterations
  // { algebraic_solving_preprocess_end }

  // Solve!
  // { algebraic_solving_begin }
  auto r = b;
  stiffnessMatrix.mmv(x, r); // r -= Ax
  // { initial_residual_end }

  // { make_residual_consistent_begin }
  // Construct consistent representation of the data in r
  auto rConsistent = r;

  VertexDataUpdate<GridView, Vector> vertexUpdateHandle(gridView, r,
                                                        rConsistent);

  gridView.communicate(vertexUpdateHandle,
                       Dune::InteriorBorder_InteriorBorder_Interface,
                       Dune::ForwardCommunication);
  // { make_residual_consistent_end }

  // { global_initial_residual_begin }
  double defect0 = r.dot(rConsistent); // Norm on the local process
  defect0 = grid->comm().sum(defect0);
  defect0 = sqrt(defect0);
  // { global_initial_residual_end}

  // { output_header_begin }
  if (mpiHelper.rank() == 0) {
    std::cout << "Iteration Defect Rate" << std::endl;
    std::cout << "0" << std::setw(16) << defect0 << std::endl;
  }
  // { output_header_end }

  // { initial_direction_begin }
  // Construct initial search direction in variable d by applying
  // the Jacobi preconditioner to the residual in r.
  Vector d(r.size());
  for (std::size_t i = 0; i < stiffnessMatrix.N(); ++i) {
    d[i] = 0;
    if (std::abs(consistentDiagonal[i]) > 1e-5) // Degree of freedom
      d[i] = rConsistent[i] / consistentDiagonal[i];
  }
  // { initial_direction_end }

  // { orthogonalization_begin }
  double rho = d.dot(r);
  rho = grid->comm().sum(rho);
  // { orthogonalization_end }

  // { loop_and_alpha_begin }
  // Current residual norm
  double defect = defect0;

  for (int k = 0; k < maxIterations; ++k) {
    // Step length in search direction d
    Vector tmp(d.size());
    stiffnessMatrix.mv(d, tmp);     // tmp = Ad^k
    double alphaDenom = d.dot(tmp); // Scalar product
    alphaDenom = grid->comm().sum(alphaDenom);
    double alpha = rho / alphaDenom;
    // { loop_and_alpha_end }

    // { update_iterate_begin }
    x.axpy(alpha, d); // Update iterate
    // { update_iterate_end }
    // { update_residual_begin }
    r.axpy(-alpha, tmp); // Update residual
    // { update_residual_end }

    // Convergence test
    // { check_convergence_begin }
    // Compute residual norm again
    rConsistent = r;
    gridView.communicate(vertexUpdateHandle,
                         Dune::InteriorBorder_InteriorBorder_Interface,
                         Dune::ForwardCommunication);

    auto residualNorm = r.dot(rConsistent);
    residualNorm = grid->comm().sum(residualNorm);
    residualNorm = sqrt(residualNorm);

    if (mpiHelper.rank() == 0) {
      std::cout << std::setw(5) << k + 1 << " ";
      std::cout << std::setw(16) << residualNorm << " ";
      // Estimated convergence rate
      std::cout << std::setw(16) << residualNorm / defect << std::endl;
    }

    defect = residualNorm; // Update norm

    if (defect < defect0 * reduction) // Convergence check
      break;
    // { check_convergence_end }

    // Determine new search direction
    // { compute_new_direction_begin }
    // Precondition the residual
    Vector preconditionedResidual(d.size());
    for (std::size_t i = 0; i < stiffnessMatrix.N(); i++) {
      preconditionedResidual[i] = 0;
      if (std::abs(consistentDiagonal[i]) > 1e-5) // Degree of freedom
        // not on ghost vertex
        preconditionedResidual[i] = rConsistent[i] / consistentDiagonal[i];
    }

    double rhoNext = preconditionedResidual.dot(r);
    rhoNext = grid->comm().sum(rhoNext);
    double beta = rhoNext / rho;

    // Compute new search direction
    d *= beta;
    d += preconditionedResidual;
    rho = rhoNext; // Remember rho for the next iterate
  }
  // { algebraic_solving_end }

  // Output result
  // { vtk_output_begin }
  // For visualization: Write the rank number for each element
  Dune::MultipleCodimMultipleGeomTypeMapper<GridView> elementMapper(
      gridView, Dune::mcmgElementLayout());
  std::vector<int> ranks(elementMapper.size());
  for (const auto &element : elements(gridView))
    ranks[elementMapper.index(element)] = mpiHelper.rank();

  Dune::VTKWriter<GridView> vtkWriter(gridView);
  vtkWriter.addVertexData(x, "solution");
  vtkWriter.addCellData(ranks, "ranks");
  vtkWriter.write("grid-distributed-poisson-result");
  // { vtk_output_end }

  return 0;
}
