// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include "Sphere.hh"

int main(int argc, char **argv)
{

  Dune::MPIHelper::instance(argc, argv);

  constexpr int dim = 2;
  using Grid = Dune::UGGrid<dim>;

  const std::array<unsigned, dim> n = {8, 25};

  const Dune::FieldVector<double, dim> lower = {0, 0};
  const Dune::FieldVector<double, dim> upper = {6, 15};

  std::shared_ptr<Grid> grid =
      Dune::StructuredGridFactory<Grid>::createSimplexGrid(lower, upper, n);

  using GridView = Grid::LeafGridView;
  const GridView gridView = grid->leafGridView();

  Sphere<dim> sphere({3.0, 2.5}, 1.0);

  const int steps = 30;
  const Dune::FieldVector<double, dim> stepDisplacement = {0, 0.5};

  const double epsilon = 0.4;

  const int levels = 3;

  for (int i = 0; i < steps; ++i) {
    std::cout << "Step " << i << std::endl;

    for (int k = 0; k < levels - 1; ++k) {
      for (const auto &element : elements(gridView))
        grid->mark(-1, element);
      grid->preAdapt();
      grid->adapt();
      grid->postAdapt();
    }

    for (int k = 0; k < levels - 1; ++k) {
      for (const auto &element : elements(gridView))
        if (sphere.distanceTo(element.geometry().center()) < epsilon)
          grid->mark(1, element);
      grid->preAdapt();
      grid->adapt();
      grid->postAdapt();
    }
    Dune::VTKWriter<GridView> vtkWriter(gridView);
    vtkWriter.write("refined_grid_" + std::to_string(i));

    sphere.displace(stepDisplacement);
  }

  return 0;
}
