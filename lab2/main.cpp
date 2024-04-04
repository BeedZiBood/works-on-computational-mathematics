#include <iostream>
#include <iomanip>
#include <algorithm>
#include "decomp.h"
#include "solve.h"
#include "matrix_function.h"

int main() {
  const auto default_precision{std::cout.precision()};
  static long double cond = 0;
  double det = 0;
  const int ndim = 12;
  int n = 2;
  int p = 4;
  int pivot[ndim], i, flag;
  long double B[ndim * ndim] = {0};
  long double E[ndim * ndim] = {0};
  long double result_for_multiplication[ndim * ndim] = {0};
  long double result_for_subtraction[ndim * ndim] = {0};
  long double copy_of_b[ndim * ndim] = {0};
  long double copy_of_e[ndim * ndim] = {0};
  for (int k = 0; k < 5; k++)
  {
    n += 2;
    for (int i = 0; i < ndim; i++)
    {
      for (int j = 0; j < ndim; j++)
      {
        B[i * ndim + j] = 1.0 / (p + i + j - 1);
        if (i == j)
        {
          E[i * ndim + j] = 1;
        }
        else
        {
          E[i * ndim + j] = 0;
        }
      }
    }
    if (k == 0)
    {
      std::copy(B, B + ndim * ndim, copy_of_b);
      std::copy(E, E + ndim * ndim, copy_of_e);
    }
    decomp(n, ndim, B, &cond, pivot, &flag);
    std::cout << "----------------------------------------------------------------------\n";
    std::cout << "Result for order of matrix = " << n << "\n";
    if (flag == 0)
    {
      for (int i = 0; i < n; i++)
      {
        solve(n, ndim, B, (E + i * ndim), pivot);
      }

      //            std::cout << "Solution (inverse matrix) = [\n";
      //            for (int i = 0; i < n; i++) {
      //                for (int j = 0; j < n; j++)
      //                {
      //                    std::cout << std::setw(13) << std::fixed << std::setprecision(0) << E[i * ndim + j] << " ";
      //                }
      //                std::cout << "\n";
      //            }
      //            std::cout << "]\n";
      std::cout << "Condition number = " << cond << "\n" << std::scientific << std::setprecision(1);

      std::cout << "R = [\n";
      mashkin::multiply(n, ndim, copy_of_b, E, result_for_multiplication);
      mashkin::subtract(n, ndim, result_for_multiplication, copy_of_e, result_for_subtraction);
      mashkin::print_matrix(n, ndim, result_for_subtraction);
      std::cout << "]\n";

      det = pivot[n-1];
      for (i = 0; i < n; i++)
        det = det * B[i * ndim + i];
    }
    else
    {
      std::cout << "Ehd with error flag = " << flag << "\n";
    }
    std::cout << "----------------------------------------------------------------------\n\n";
  }
  return 0;
}
