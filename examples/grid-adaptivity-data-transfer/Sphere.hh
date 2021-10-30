//
// Created by carlosal1015 on 10/29/21.
//

#ifndef DUNE_BOOK_SPHERE_HH
#define DUNE_BOOK_SPHERE_HH

#include <dune/common/fvector.hh>

template <int dim> class Sphere {
  double radius_;
  Dune::FieldVector<double, dim> center_;

public:
  Sphere(const Dune::FieldVector<double, dim> &c, const double &r)
      : radius_(r), center_(c) {}

  double distanceTo(const Dune::FieldVector<double, dim> &point) const {
    return std::abs((center_ - point).two_norm() - radius_);
  }

  void displace(const Dune::FieldVector<double, dim> &increment) {
    center_ += increment;
  }
};

#endif // DUNE_BOOK_SPHERE_HH
