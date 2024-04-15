#include <iostream>
#include "functions_for_course_work.h"
#include "FORSTYTHE.h"
#include <cstddef>
#include "D:/Labs/CompMath/works-on-computational-mathematics/lab2/decomp.h"
#include "D:/Labs/CompMath/works-on-computational-mathematics/lab2/solve.h"
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
  std::cout << "E = " << E << "\n\n";
  long double A[19 * 19] = {0};
  long double B[19] = {0};
  int pivot[19], flag_for_matrix;
  long double cond = 0;
  for (int i = 1; i <= 2; i++)
  {
    int n = 10 * i;
    double h = l / n;
    B[n - 2] = - E + h * (n - 1) * h / 2;
    A[0] = -2 * E - h * h * i * h;
    A[1] = E + h * i * h / 2;
    for (size_t j = 1; j < n - 2; j++)
    {
      //std::cout << j * (n - 1) + j - 1 << " " << j * (n - 1) + j << " " << j * (n - 1) + j + 1;
      A[j * (n - 1) + j - 1] = E - h * ((j + 1) * h) / 2;
      A[j * (n - 1) + j] = -2 * E - h * h * (j + 1) * h;
      A[j * (n - 1) + j + 1] = E + h * (j + 1) * h / 2;
      //std::cout << "\n";
    }
    A[(n - 1) * (n - 1) - 2] = E - h * (9 * h) / 2;
    A[(n - 1) * (n - 1) - 1] = -2 * E - h * h * 9 * h;

    decomp(n - 1, n - 1, A, &cond, pivot, &flag_for_matrix);
    solve(n - 1, n - 1, A, B, pivot);

    std::cout << "result for n = " << n << ":\n";
    for (size_t k = 0; k < n - 1; k++)
    {
      std::cout << B[k] << " ";
      B[k] = 0.0;
      for (size_t j = 0; j < n - 1; j++) {
        //std::cout << A[k * (n - 1) + j] << " ";
        //std::cout << k * (n - 1) + j << " ";
        A[k * (n - 1) + j] = 0.0;
      }
      //std::cout << "\n";
    }
    std::cout << "\n";
  }
}