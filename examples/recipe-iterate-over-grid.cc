// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <dune/grid/yaspgrid.hh>
int main(int argc, char **argv) {
  Dune::MPIHelper &helper = Dune::MPIHelper::instance(argc, argv);

  constexpr int dim = 4;
  using Grid = Dune::YaspGrid<dim>;
  Dune::FieldVector<double, dim> len;
  for (auto &l : len)
    l = 1.0;
  std::array<int, dim> cells{};
  for (auto &c : cells)
    c = 4;
  Grid grid(len, cells);
  auto gv = grid.leafGridView();

  constexpr int codim = 2;
  for (const auto &e : entities(gv, Dune::Codim<codim>{}))
    if (!e.type().isCube())
      std::cout << "not a cube" << std::endl;

  for (const auto &e : elements(gv))
    ;
  for (const auto &e : vertices(gv))
    ;
  for (const auto &e : edges(gv))
    ;
  for (const auto &e : facets(gv))
    ;

  constexpr int mycodim = 2;
  for (const auto &e : elements(gv))
    for (unsigned int i = 0; i < e.subEntities(mycodim); ++i)
      auto v = e.template subEntity<codim>(i);
}
