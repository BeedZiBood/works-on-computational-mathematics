#ifndef LAB2_MATRIX_FUNCTION_H
#define LAB2_MATRIX_FUNCTION_H

namespace mashkin
{
    void multiply(int n, int ndim, long double *a, long double *b, long double *result);
    void subtract(int n, int ndim, long double *a, long double *b, long double *result);
    void print_matrix(int n, int ndim, long double *matrix);
    double calculate_the_norm_of_a_matrix(long double *matrix, int ndim);
}
#endif
