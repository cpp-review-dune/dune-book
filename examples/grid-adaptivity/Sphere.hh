//
// Created by carlosal1015 on 10/29/21.
//

#pragma once

#include <dune/common/fvector.hh>

// { sphere_begin }
template <int dim> class Sphere {
  double radius_;
  Dune::FieldVector<double, dim> center_;

public:
  Sphere(const Dune::FieldVector<double, dim> &center, const double &radius)
      : radius_(radius), center_(center)
  {
  }

  double distanceTo(const Dune::FieldVector<double, dim> &point) const
  {
    return std::abs((center_ - point).two_norm() - radius_);
  }

  void displace(const Dune::FieldVector<double, dim> &increment)
  {
    center_ += increment;
  }
};
// { sphere_end }