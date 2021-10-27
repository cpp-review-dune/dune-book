// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/yaspgrid.hh>
#include <iostream>
#include <vector>

template <class GridView, class Mapper>
void evolve(
    const GridView &gridView, const Mapper &mapper, double dt,
    std::vector<double> &c,
    const std::function<Dune::FieldVector<double, GridView::dimension>(
        Dune::FieldVector<double, GridView::dimension>)>
        v,
    const std::function<double(Dune::FieldVector<double, GridView::dimension>)>
        inflow) {
  constexpr int dim = GridView::dimension;

  std::vector<double> update(c.size());
  std::fill(update.begin(), update.end(), 0.0);

  for (const auto &element : elements(gridView)) {
    auto geometry = element.geometry();

    double elementVolume = geometry.volume();

    typename Mapper::Index i = mapper.index(element);

    for (const auto &intersection : intersections(gridView, element)) {
      auto intersectionGeometry = intersection.geometry();

      Dune::FieldVector<double, dim> intersectionCenter =
          intersectionGeometry.center();

      Dune::FieldVector<double, dim> velocity = v(intersectionCenter);

      const auto &intersectionReferenceElement =
          Dune::ReferenceElements<double, dim - 1>::general(
              intersection.type());

      Dune::FieldVector<double, dim - 1> intersectionLocalCenter =
          intersectionReferenceElement.position(0, 0);

      Dune::FieldVector<double, dim> integrationOuterNormal =
          intersection.integrationOuterNormal(intersectionLocalCenter);

      double intersectionFlow = velocity * integrationOuterNormal;

      update[i] -= c[i] * std::max(0.0, intersectionFlow) / elementVolume;

      if (intersectionFlow <= 0) {
        if (intersection.neighbor()) {
          auto j = mapper.index(intersection.outside());
          update[i] -= c[j] * intersectionFlow / elementVolume;
        }
        if (intersection.boundary())
          update[i] -=
              inflow(intersectionCenter) * intersectionFlow / elementVolume;
      }
    }
  }
  for (std::size_t i = 0; i < c.size(); ++i)
    c[i] += dt * update[i];
}

int main(int argc, char **argv) {

  Dune::MPIHelper::instance(argc, argv);

  constexpr int dim = 2;
  using Grid = Dune::YaspGrid<dim>;
  Grid grid({1.0, 1.0}, {80, 80});
  using GridView = Grid::LeafGridView;
  GridView gridView = grid.leafGridView();

  Dune::MultipleCodimMultipleGeomTypeMapper<GridView> mapper(
      gridView, Dune::mcmgElementLayout());

  std::vector<double> c(mapper.size());

  auto c0 = [](const Dune::FieldVector<double, dim> &x) {
    return (x.two_norm() > 0.125 && x.two_norm() < 0.5) ? 1.0 : 0.0;
  };

  for (const auto &element : elements(gridView)) {
    auto geometry = element.geometry();
    auto global = geometry.center();
    c[mapper.index(element)] = c0(global);
  }

  auto vtkWriter = std::make_shared<Dune::VTKWriter<GridView>>(gridView);
  Dune::VTKSequenceWriter<GridView> vtkSequenceWriter(
      vtkWriter, "getting-started-transport-fv-result");

  vtkWriter->addCellData(c, "concentration");
  vtkSequenceWriter.write(0.0);

  double t = 0;
  const double tend = 0.6;
  const double dt = 0.006;
  int k = 0;

  auto inflow = [](const Dune::FieldVector<double, dim> &x) { return 0.0; };

  auto v = [](const Dune::FieldVector<double, dim> &x) {
    return Dune::FieldVector<double, dim>(1.0);
  };

  while (t < tend) {
    evolve(gridView, mapper, dt, c, v, inflow);
    t += dt;
    ++k;
    vtkSequenceWriter.write(t);
    std::cout << "k = " << k << "t = " << t << std::endl;
  }
  return 0;
}
