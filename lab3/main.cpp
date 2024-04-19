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
  double x_first[2]{3.0, -1.0};
  double x_to_csv[2]{3.0, -1.0};
  double t(0);
  double tOutFirst(0);
  double relerr(0.0001);
  double abserr(0.0001);
  double work_first[3+6*2];
  int flag = 1;
  std::vector< double > x_vect_first{3.0, -1.0};
  outFile << "t,rkf1,rkf2,rk1,rk2\n";
  std::cout << "RKF45:                                    | RungeKutta h = 0.0075    |\n";
  for (double i = 0.0075; i <= 0.1501; i += 0.0075)
  {
    tOutFirst = i;
    RKF45(mashkin::func< double*, double* >, 2, x_first, t, tOutFirst, relerr, abserr, work_first, flag);
    std::cout << std::setw(6) << tOutFirst << " ";
    outFile << tOutFirst << "," << x_first[0] << "," << x_first[1] << ",";
    std::cout << std::setprecision(5) << std::setw(11) << x_first[0] << " ";
    std::cout << std::setw(12) << x_first[1] << "  flag = " << flag << " |" << std::setw(12);
    std::cout << mashkin::rungeKutta3degree(tOutFirst, x_vect_first, 0.0075) << " |" << std::setw(10) << "\n";
    outFile << x_vect_first[0] << "," << x_vect_first[1] << "\n";
  }

  double x_second[2]{3.0, -1.0};
  double tOutSecond(0);
  relerr = 0.0001;
  abserr = 0.0001;
  double work_second[3+6*2];
  flag = 1;
  std::vector< double > x_vect_second{3.0, -1.0};
  std::cout << "\nRKF45 with h = 0.00375:                       | RungeKutta h = 0.00375   |\n";
  t = 0.0;
  for (double i = 0.00375; i <= 0.1501; i += 0.00375)
  {
    tOutSecond = i;
    RKF45(mashkin::func, 2, x_second, t, tOutSecond, relerr, abserr, work_second, flag);
    std::cout << std::setw(8) << tOutSecond << " ";
    outFile << tOutSecond << "," << x_second[0] << "," << x_second[1] << ",";
    std::cout << std::setprecision(5) << std::setw(13) << x_second[0] << " ";
    std::cout << std::setw(13) << x_second[1] << "  flag = " << flag << "|" << std::setw(12);
    std::cout << mashkin::rungeKutta3degree(tOutSecond, x_vect_second, 0.00375) << " |\n";
    outFile << x_vect_second[0] << "," << x_vect_second[1] << "\n";
  }
  return 0;
}
