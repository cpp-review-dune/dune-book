// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "evolve.hh"
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/yaspgrid.hh>

// { using_namespace_dune_begin }
// { using_namespace_dune_end }

// { main_begin }
int main(int argc, char **argv)
{

  // Set up MPI, if available
  Dune::MPIHelper::instance(argc, argv);
  // { main_signature_end }

  // { create_grid_begin }
  constexpr int dim = 2;
  using Grid = Dune::YaspGrid<dim>;
  Grid grid({1.0, 1.0}, // Upper right corner, the lower left one is (0,0)
            {80, 80});  // Number of elements per direction

  using GridView = Grid::LeafGridView;
  GridView gridView = grid.leafGridView();
  // { create_grid_end }

  // Assigns a unique number to each element
  // { create_concentration_begin }
  Dune::MultipleCodimMultipleGeomTypeMapper<GridView> mapper(
      gridView, Dune::mcmgElementLayout());

  // Allocate a vector for the concentration
  std::vector<double> c(mapper.size());
  // { create_concentration_end }

  // Initial concentration
  // { lambda_initial_concentration_begin }
  auto c0 = [](const Dune::FieldVector<double, dim> &x) {
    return (x.two_norm() > 0.125 && x.two_norm() < 0.5) ? 1.0 : 0.0;
  };
  // { lambda_initial_concentration_end }

  // { sample_initial_concentration_begin }
  // Iterate over grid elements and evaluate c0 at element centers
  for (const auto &element : elements(gridView)) {
    // Get element geometry
    auto geometry = element.geometry();

    // Get global coordinate of element center
    auto global = geometry.center();

    // Sample initial concentration \c_0 at the element center
    c[mapper.index(element)] = c0(global);
  }
  // { sample_initial_concentration_end }

  // Construct VTK writer
  // { construct_vtk_writer_begin }
  auto vtkWriter = std::make_shared<Dune::VTKWriter<GridView>>(gridView);
  Dune::VTKSequenceWriter<GridView> vtkSequenceWriter(
      vtkWriter, "getting-started-transport-fv-result"); // File name

  // Write the initial values
  vtkWriter->addCellData(c, "concentration");
  vtkSequenceWriter.write(0.0); // 0.0 is the current time
  // { construct_vtk_writer_end }

  // Now do the time steps
  // { time_loop_begin }
  double t = 0;            // Initial time
  const double tend = 0.6; // Final time
  const double dt = 0.006; // Time step size
  int k = 0;               // Time step counter

  // Inflow boundary values
  auto inflow = [](const Dune::FieldVector<double, dim> &x) { return 0.0; };

  // Velocity field
  auto v = [](const Dune::FieldVector<double, dim> &x) {
    return Dune::FieldVector<double, dim>(1.0);
  };

  while (t < tend) {
    // Apply finite volume scheme
    evolve(gridView, mapper, dt, c, v, inflow);

    // Augment time and time step counter
    t += dt;
    ++k;

    // Write data. We do not have call addCellData again!
    vtkSequenceWriter.write(t);

    // Print iteration number, time, and time step size
    std::cout << "k = " << k << "t = " << t << std::endl;
  }
  // { time_loop_end }
  return 0;
}
// { main_end }
