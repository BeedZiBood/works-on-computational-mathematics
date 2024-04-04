#include "solve.h"

#define AINDEX(i,j) (i * ndim + j)

int solve (int n, int ndim, long double *a, long double b[], int pivot[])
{   /* --- begin function solve() --- */

    int    i, j, k, m;
    long double t;

    if (n == 1)
    {
        /* trivial */
        b[0] /= a[0];
    }
    else
    {
        /* Forward elimination: apply multipliers. */
        for (k = 0; k < n-1; k ++)
        {
            m = pivot[k];
            t = b[m]; b[m] = b[k]; b[k] = t;
            for (i = k+1; i < n; ++i) b[i] += a[AINDEX(i,k)] * t;
        }

        /* Back substitution. */
        for (k = n-1; k >= 0; --k)
        {
            t = b[k];
            for (j = k+1; j < n; ++j) t -= a[AINDEX(k,j)] * b[j];
            b[k] = t / a[AINDEX(k,k)];
        }
    }

    return(0);
}  /* --- end function solve() --- */
