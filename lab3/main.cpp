#include <iostream>
#include <iomanip>
#include <fstream>
#include "RKF45.h"
#include "Function_for_laboratory.h"
#include "RungeKutta.h"

int main(int argc, char** argv)
{
  if (argc != 2)
  {
    std::cout << "Not enough argument\n";
    return 1;
  }
  std::ofstream outFile;
  outFile.exceptions(std::ofstream::badbit | std::ofstream::failbit);
  try
  {
    outFile.open(argv[1]);
  }
  catch (const std::exception& ex)
  {
    std::cerr << ex.what() << "\n";
    return 1;
  }
  double x[2]{3.0, -1.0};
  double x_to_csv[2]{3.0, -1.0};
  double t(0);
  double tOutFirst(0);
  double tOutSecond(0);
  double relerr(0.0001);
  double abserr(0.0001);
  double work[3+6*2];
  int flag = 1;
  std::vector< double > x_vect_first{3.0, -1.0};
  std::vector< double > x_vect_second{3.0, -1.0};
  outFile << "t,rkf1,rkf2,rk11,rk12,rk21,rk22\n";
  std::cout << "RKF45:                                    | RungeKutta h = 0.0075    |\n";
  for (double i = 0.0075; i <= 0.1501; i += 0.0075)
  {
    tOutFirst = i;
    RKF45(mashkin::func< double*, double* >, 2, x, t, tOutFirst, relerr, abserr, work, flag);
    std::cout << std::setw(6) << tOutFirst << " ";
    outFile << tOutFirst << "," << x[0] << "," << x[1] << ",";
    std::cout << std::setprecision(5) << std::setw(11) << x[0] << " ";
    std::cout << std::setw(12) << x[1] << "  flag = " << flag << " |" << std::setw(12);
    std::cout << mashkin::rungeKutta3degree(tOutFirst, x_vect_first, 0.0075) << " |" << std::setw(10) << "\n";
    outFile << x_vect_first[0] << "," << x_vect_first[1] << "\n";
  }
  std::cout << "\nRKF45 with h = 0.00375:                       | RungeKutta h = 0.00375   |\n";
  t = 0.0;
  for (double i = 0.00375; i <= 0.1501; i += 0.00375)
  {
    tOutSecond = i;
    RKF45(mashkin::func, 2, x, t, tOutSecond, relerr, abserr, work, flag);
    std::cout << std::setw(8) << tOutSecond << " ";
    outFile << tOutSecond << "," << x[0] << "," << x[1] << ",";
    std::cout << std::setprecision(5) << std::setw(13) << x[0] << " ";
    std::cout << std::setw(13) << x[1] << "  flag = " << flag << "|" << std::setw(12);
    std::cout << mashkin::rungeKutta3degree(tOutSecond, x_vect_second, 0.00375) << " |\n";
    outFile << x_vect_second[0] << "," << x_vect_second[1] << "\n";
  }
  return 0;
}
