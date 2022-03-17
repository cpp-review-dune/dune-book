//
// Created by carlosal1015 on 10/29/21.
//

#pragma once

#include <dune/common/fvector.hh>
#include <dune/geometry/referenceelements.hh>

// { evolve_signature_begin }
template <class GridView, class Mapper>
void evolve(
    const GridView &gridView, const Mapper &mapper,
    double dt, // Time step size
    std::vector<double> &c,
    const std::function<Dune::FieldVector<double, GridView::dimension>(
        Dune::FieldVector<double, GridView::dimension>)>
        v,
    const std::function<double(Dune::FieldVector<double, GridView::dimension>)>
        inflow)
// { evolve_signature_end }
{
  // { evolve_init_begin }
  // Grid dimension
  constexpr int dim = GridView::dimension;

  // Allocate a temporary vector for the update
  std::vector<double> update(c.size());
  std::fill(update.begin(), update.end(), 0.0);
  // { evolve_init_end }

  // Compute update vector
  // { element_loop_begin }
  for (const auto &element : elements(gridView)) {
    // Element geometry
    auto geometry = element.geometry();

    // Element volume
    double elementVolume = geometry.volume();

    // Unique element number
    typename Mapper::Index i = mapper.index(element);
    // { element_loop_end }

    // Loop over all intersections \gamma_{ij} with neighbors and boundary
    // { intersection_loop_begin }
    for (const auto &intersection : intersections(gridView, element)) {
      // Geometry of the intersection
      auto intersectionGeometry = intersection.geometry();

      // Center of intersection in global coordinates
      Dune::FieldVector<double, dim> intersectionCenter =
          intersectionGeometry.center();

      // Velocity at intersection center \vector{v}_{ij}
      Dune::FieldVector<double, dim> velocity = v(intersectionCenter);

      // Center of the intersection in local coordinates
      const auto &intersectionReferenceElement =
          Dune::ReferenceElements<double, dim - 1>::general(
              intersection.type());

      Dune::FieldVector<double, dim - 1> intersectionLocalCenter =
          intersectionReferenceElement.position(0, 0);

      // Normal vector scaled with intersection area:
      // \vector{n}_{ij}|\gamma_{ij}|
      Dune::FieldVector<double, dim> integrationOuterNormal =
          intersection.integrationOuterNormal(intersectionLocalCenter);

      // Compute factor occurring in flux formula:
      // \langle\vector{v}_{ij},\vector{n}_{ij}\rangle|\gamma_{ij}|
      double intersectionFlow = velocity * integrationOuterNormal;
      // { intersection_loop_initend }

      // { intersection_loop_mainbegin }
      // Outflow contributions
      update[i] -= c[i] * std::max(0.0, intersectionFlow) / elementVolume;

      // Inflow contributions
      if (intersectionFlow <= 0) {
        // Handle interior intersection
        if (intersection.neighbor()) {
          // Access neighbor
          auto j = mapper.index(intersection.outside());
          update[i] -= c[j] * intersectionFlow / elementVolume;
        }

        // Handle boundary intersection
        if (intersection.boundary())
          update[i] -=
              inflow(intersectionCenter) * intersectionFlow / elementVolume;
      }
      // { intersection_loop_end }
    } // End loop over all intersections
  }   // End loop over the grid elments
  // { element_loop_end }

  // { evolve_laststeps }
  // Update the concentration vector
  for (std::size_t i = 0; i < c.size(); ++i)
    c[i] += dt * update[i];
}
// { evolve_end }