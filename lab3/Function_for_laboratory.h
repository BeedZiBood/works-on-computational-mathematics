#ifndef LAB3_FUNCTION_FOR_LABORATORY_H
#define LAB3_FUNCTION_FOR_LABORATORY_H
#include <cmath>

namespace mashkin
{
template< class X, class DX >
  void func(double t, X x, DX dx)
  {
    dx[0] = -130 * x[0] + 900 * x[1] + std::exp(-10 * t);
    dx[1] = 30 * x[0] - 300 * x[1] + std::log(1 + 100 * t * t);
  }
}
#endif
