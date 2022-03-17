//
// Created by carlosal1015 on 10/29/21.
//

#pragma once

#include <dune/common/fvector.hh>

// Linear interpolation on a simplex
// { linear_interpolation_begin }
template <int dim>
double interpolate(const std::vector<double> values,
                   Dune::FieldVector<double, dim> p)
{
  assert(values.size() == dim + 1);
  double result = values[0];
  for (std::size_t i = 0; i < p.size(); i++)
    result += p[i] * (values[i + 1] - values[0]);
  return result;
}
// { linear_interpolation_end }