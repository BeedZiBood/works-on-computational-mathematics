#include <iostream>
#include "functions_for_course_work.h"
#include "FORSTYTHE.h"
#include <math.h>

int main()
{
  //find min x for E
  double x;
  int flag;
  x = zeroin(3.75, 5.0, mashkin::function_for_search_min_x, 0.0000001, &flag);
  std::cout << "Min X is " << x << "\n";
  double l = 0.2231271 * x;
  std::cout << "l = " << l << "\n";
  double a = 0.0;
  double b = M_PI/2;
  double epsrel = 0.0000001;
  double epsabs = 0.0;
  double errest = 0.0;
  double flag_quanc8 = 0.0;
  int nfe = 0;
  double E;
  QUANC8(mashkin::function_for_integrate, a, b, epsabs, epsrel, E, errest, nfe, flag_quanc8);
  E *= 0.6436369;
  std::cout << "E = " << E << "\n";
}