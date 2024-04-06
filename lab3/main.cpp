#include <iostream>
#include "RKF45.h"
#include "Function_for_laboratory.h"

int main()
{
  double x[2]{3.0, -1.0};
  double t(0);
  double tOut(0);
  double relerr(0.0001);
  double abserr(0.0001);
  double work[3+6*2];
  int flag = 1;
  int line = 1;
  for (double i = 0.0075; i <= 0.1501; i += 0.0075)
  {
    tOut = i;
    std::cout << line++ << " ";
    RKF45(mashkin::func, 2, x, t, tOut, relerr, abserr, work, flag);
    std::cout << t << " " << tOut << " " << x[0] << " " << x[1] << " flag = " << flag << "\n";
  }
  return 0;
}
