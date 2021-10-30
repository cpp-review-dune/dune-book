// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/yaspgrid.hh>

int main(int argc, char **argv) {
  Dune::MPIHelper::instance(argc, argv);

  constexpr int dim = 4;
  using Grid = Dune::YaspGrid<dim>;
  Dune::FieldVector<double, dim> len;
  for (auto &l : len)
    l = 1.0;
  std::array<int, dim> cells{};
  for (auto &c : cells)
    c = 5;
  Grid grid(len, cells);

  Dune::FieldVector<double, 4> x({1, 2, 3, 4});
  auto y(x);
  y *= 1.0 / 3.0;
  auto s = x * y;
  auto norm = x.two_norm();
  Dune::FieldMatrix<double, 4, 4> A(
      {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}});
  A.mv(x, y);
  A.usmv(0.5, x, y);

  auto u = [](const auto &x) { return std::exp(x.two_norm()); };

  double integral = 0.0;
  auto gv = grid.leafGridView();
  for (const auto &e : Dune::elements(gv))
    integral += u(e.geometry().center() * e.geometry().volume());
  std::cout << "integral = " << integral << std::endl;

  double integral2 = 0.0;
  using QR = Dune::QuadratureRules<Grid::ctype, dim>;
  for (const auto &e : elements(gv)) {
    auto geo = e.geometry();
    auto quadrature = QR::rule(geo.type(), 5);
    for (const auto &qp : quadrature)
      integral2 += u(geo.global(qp.position())) *
                   geo.integrationElement(qp.position()) * qp.weight();
  }
  std::cout << "integral2 = " << integral2 << std::endl;

  auto f = [](const auto &x) { return x; };
  double divergence = 0.0;
  for (const auto &i : elements(gv)) {
    for (const auto &I : intersections(gv, i))
      if (!I.neighbor()) {
        auto geoI = I.geometry();
        divergence +=
            f(geoI.center()) * I.centerUnitOuterNormal() * geoI.volume();
      }
  }
  std::cout << "divergence  = " << divergence << std::endl;
}
