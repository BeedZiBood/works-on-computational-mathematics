#include <iostream>
#include <functional>
#include <function.h>
#include <fstream>
#include <iomanip>
#include "Lagrange.h"
#include "FORSYTHE.H"

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
  double a, b, epsabs, epsrel, result, errest;
  int nfe = 0;
  double flag = 0.0;
  REAL x[] = {-1.000, -0.960, -0.860, -0.790, 0.220, 0.500, 0.930};
  REAL y[] = {-1.000, -0.151, 0.894, 0.986, 0.895, 0.500, -0.306};
  auto size = sizeof(x) / sizeof(double);
  SPLINE spline(7, x, y);
  std::cout << "      x      |    Lagrange   |       Spline\n";
  for (auto i = -1.0; i < -1.0 + 0.1 * 19; i += 0.1)
  {
    outFile << i << "," << Lagrange(y, x, i, size) << "," << spline.Eval(i) << "\n";
    std::cout << std::setw(12) << i << " | ";
    std::cout << std::setw(13) << Lagrange(y, x, i, size) << " | ";
    std::cout << std::setw(12) << spline.Eval(i) << "\n";
  }
  outFile.close();

  a = 1.0;
  b = 2.0;
  epsrel = 0.0001;
  epsabs = 0.0;

  QUANC8(f, a, b, epsabs, epsrel, result, errest, nfe, flag);
  std::cout << "\n\n";
  std::cout << "Result: " << result << "\n";
  std::cout << "Error: " << errest << "\n";
  std::cout << "NoFun: " << nfe << "\n";
  std::cout << "Flag: " << flag << "\n";
}