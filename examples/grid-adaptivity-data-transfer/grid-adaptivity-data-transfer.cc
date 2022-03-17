// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <map>

#include "Sphere.hh"
#include "interpolate.hh"

// https://stackoverflow.com/a/5192091/9302545
// { grid_setup_begin }
int main(int argc, char *argv[])
{
  // Set up MPI if available
  Dune::MPIHelper::instance(argc, argv);

  constexpr int dim = 2; // Grid dimension
  using Grid = Dune::UGGrid<dim>;

  // Create UGGrid from structured triangle grid
  const std::array<unsigned, dim> n = {8, 25};

  const Dune::FieldVector<double, dim> lower = {0, 0};
  const Dune::FieldVector<double, dim> upper = {6, 15};

  std::shared_ptr<Grid> grid =
      Dune::StructuredGridFactory<Grid>::createSimplexGrid(lower, upper, n);

  using GridView = Grid::LeafGridView;
  const GridView gridView = grid->leafGridView();
  const auto &indexSet = gridView.indexSet();
  const auto &idSet = grid->localIdSet();
  // { grid_setup_end }

  // Create sphere
  // { sphere_setup_begin }
  Sphere<dim> sphere({3.0, 2.5}, 1.0);

  // Set parameters
  const int steps = 30;
  const Dune::FieldVector<double, dim> stepDisplacement = {0, 0.5};

  const double epsilon = 0.4; // Thickness of the refined region
                              // around the sphere

  const int levels = 2;
  // { sphere_setup_end }

  // Construct a piecewise linear function on the grid
  // { function_setup_begin }
  auto dataFunction = [](const Dune::FieldVector<double, dim> &x) {
    return std::sin(x[1]);
  };

  std::vector<double> data(gridView.size(dim));
  for (auto &&v : vertices(gridView)) {
    data[indexSet.index(v)] = dataFunction(v.geometry().corner(0));
  }
  // { function_setup_end }

  // { coarsening_refinement_begin }
  for (int i = 0; i < steps; ++i) {
    std::map<Grid::LocalIdSet::IdType, double> persistentContainer;

    // Coarsen everything
    for (int j = 0; j < levels; ++j) {
      for (const auto &element : elements(gridView)) {
        grid->mark(-1, element); // Mark element for coarsening
      }

      grid->preAdapt();

      for (const auto &vertex : vertices(gridView)) {
        persistentContainer[idSet.id(vertex)] = data[indexSet.index(vertex)];
      }

      grid->adapt();

      data.resize(gridView.size(dim));

      for (const auto &v : vertices(gridView)) {
        data[indexSet.index(v)] = persistentContainer[idSet.id(v)];
      }

      grid->postAdapt();
    }
    // { coarsening_end }

    // Refine near the sphere
    // { refinement_begin }
    for (int j = 0; j < levels - 1; ++j) {
      // Select elements that are close to the sphere for grid refinement
      for (const auto &element : elements(gridView)) {
        if (sphere.distanceTo(element.geometry().center()) < epsilon) {
          grid->mark(1, element);
        }
      }

      grid->preAdapt();

      for (const auto &vertex : vertices(gridView)) {
        persistentContainer[idSet.id(vertex)] = data[indexSet.index(vertex)];
      }

      grid->adapt();
      // { refinement_after_modification }

      // { refinement_data_interpolation_begin }
      data.resize(gridView.size(dim));

      for (const auto &element : elements(gridView)) {
        if (element.isNew()) {
          for (std::size_t k = 0; k < element.subEntities(dim); k++) {
            auto father = element;
            auto positionInFather =
                Dune::ReferenceElements<double, dim>::general(element.type())
                    .position(k, dim);

            do {
              positionInFather =
                  father.geometryInFather().global(positionInFather);
              father = father.father();
            } while (father.isNew());

            // Extract corner values
            std::vector<double> values(father.subEntities(dim));
            for (std::size_t l = 0; l < father.subEntities(dim); l++)
              values[l] = persistentContainer[idSet.subId(father, l, dim)];

            // Interpolate linearly on the ancestor complex
            data[indexSet.subIndex(element, k, dim)] =
                interpolate(values, positionInFather);
          }
        }
        else
          for (std::size_t k = 0; k < element.subEntities(dim); k++)
            data[indexSet.subIndex(element, k, dim)] =
                persistentContainer[idSet.subId(element, k, dim)];
      }

      grid->postAdapt();
      // { refinement_data_interpolation_end }
    }
    // { coarsening_refinement_end }

    // Write grid to file
    // { writing_moving_begin }
    Dune::VTKWriter<GridView> vtkWriter(gridView);
    vtkWriter.addVertexData(data, "data");
    vtkWriter.write("refined_grid_" + std::to_string(i));

    // Move sphere
    sphere.displace(stepDisplacement);
  }
  // { writing_moving_end }
  return 0;
}
