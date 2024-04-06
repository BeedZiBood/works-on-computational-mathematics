#include "Function_for_laboratory.h"
#include <cmath>

namespace mashkin
{
  void func(double t, double* x, double* dx)
  {
    dx[0] = -130 * x[0] + 900 * x[1] + std::exp(-10 * t);
    dx[1] = 30 * x[0] - 300 * x[1] + std::log(1 + 100 * t * t);
  }
}
