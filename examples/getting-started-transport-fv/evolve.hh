//
// Created by carlosal1015 on 10/29/21.
//

#ifndef DUNE_BOOK_EVOLVE_HH
#define DUNE_BOOK_EVOLVE_HH

#include <dune/common/fvector.hh>
#include <dune/geometry/referenceelements.hh>
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

#endif // DUNE_BOOK_EVOLVE_HH
