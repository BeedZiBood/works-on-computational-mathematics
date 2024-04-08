#include <iostream>
#include <iomanip>
#include "RKF45.h"
#include "Function_for_laboratory.h"
#include "RungeKutta.h"

int main()
{
  double x[2]{3.0, -1.0};
  double t(0);
  double tOut(0);
  double relerr(0.0001);
  double abserr(0.0001);
  double work[3+6*2];
  int flag = 1;
  std::vector< double > x_vect_first{3.0, -1.0};
  std::vector< double > x_vect_second{3.0, -1.0};
  std::cout << "RKF45:                                    | RungeKutta h = 0.0075    | h = 0.037\n";
  for (double i = 0.0075; i <= 0.1501; i += 0.0075)
  {
    tOut = i;
    RKF45(mashkin::func, 2, x, t, tOut, relerr, abserr, work, flag);
    std::cout << std::setw(6) << tOut << " ";
    std::cout << std::setprecision(5) << std::setw(11) << x[0] << " ";
    std::cout << std::setw(12) << x[1] << "  flag = " << flag << " |" << std::setw(12);
    std::cout << mashkin::rungeKutta3degree(tOut, x_vect_first, 0.0075) << " |" << std::setw(10);
    std::cout << mashkin::rungeKutta3degree(tOut, x_vect_second, 0.0037) << "\n";
  }
  return 0;
}
