#include "matrix_function.h"
#include <iomanip>
#include <iostream>

namespace mashkin {
    void multiply(int n, int ndim, long double *a, long double *b, long double *result) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                long double temp = 0.0;
                for (int k = 0; k < n; k++) {
                    temp += a[i * ndim + k] * b[k * ndim + j];
                }
                result[i * ndim + j] = temp;
            }
        }
    }

    void subtract(int n, int ndim, long double *a, long double *b, long double *result)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                result[i * ndim + j] = a[i * ndim + j] - b[i * ndim + j];
            }
        }
    }

    void print_matrix(int n, int ndim, long double *matrix)
    {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++)
            {
                std::cout << std::setw(8) << matrix[i * ndim + j] << " ";
            }
            std::cout << "\n";
        }
    }
}
