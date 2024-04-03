#include "Lagrange.h"

REAL Lagrange(REAL* f, REAL* x, REAL z, int size)
{
  REAL var = 1.0;
  REAL result = 0.0;
  for (auto k = 0; k <= size; k++)
  {
    for (auto i = 0; i <= size; i++)
    {
      if (k == i)
      {
        continue;
      }
      var *= (z - x[i]) / (x[k] - x[i]);
    }
    result += var * f[k];
    var = 1.0;
  }
  return result;
}