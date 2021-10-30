//
// Created by carlosal1015 on 10/29/21.
//

#ifndef DUNE_BOOK_INTERPOLATE_HH
#define DUNE_BOOK_INTERPOLATE_HH

#include <dune/common/fvector.hh>

template <int dim>
double interpolate(const std::vector<double> values,
                   Dune::FieldVector<double, dim> p) {
  assert(values.size() == dim + 1);
  double result = values[0];
  for (std::size_t i = 0; i < p.size(); i++)
    result += p[i] * (values[i + 1] - values[0]);
  return result;
}

#endif // DUNE_BOOK_INTERPOLATE_HH
