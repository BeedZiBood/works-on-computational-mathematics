#include "decomp.h"
#include <stdexcept>
#include <cmath>
#include "solve.h"

#define AINDEX(i,j) (i * ndim + j)

int decomp(int n, int ndim, long double *a, long double *cond, int pivot[], int *flag)
{   /* --- function decomp() --- */
  long double EPSILON = 2.2e-16;
  long double ek, t, pvt, anorm, ynorm, znorm;
  int    i, j, k, m;
  long double *pa, *pb;      /* temporary pointers */
  long double *work;

  *flag = 0;
  work = (long double *) nullptr;

  if (a == nullptr || pivot == nullptr || n < 1 || ndim < n)
  {
    *flag = 2;
    return (0);
  }

  pivot[n-1] = 1;
  if (n == 1)
  {
    /* One element only */
    *cond = 1.0;
    if (*a == 0.0)
    {
      *cond = 1.0e+32;  /* singular */
      *flag = 3;
      return (0);
    }
    return (0);
  }

  work = (long double *) malloc(n * sizeof(long double));
  if (work == nullptr)
  {
    *flag = 1;
    return (0);
  }

  /* --- compute 1-norm of a --- */

  anorm = 0.0;
  for (j = 0; j < n; ++j)
  {
    t = 0.0;
    for (i = 0; i < n; ++i) t += fabs(a[AINDEX(i,j)]);
    if (t > anorm) anorm = t;
  }

  /* Apply Gaussian elimination with partial pivoting. */

  for (k = 0; k < n-1; ++k)
  {
    /* Find pivot and label as row m.
  This will be the element with largest magnitude in
  the lower part of the kth column. */
    m = k;
    pvt = fabs(a[AINDEX(m,k)]);
    for (i = k+1; i < n; ++i)
    {
      t = fabs(a[AINDEX(i,k)]);
      if ( t > pvt )  { m = i; pvt = t; }
    }
    pivot[k] = m;
    pvt = a[AINDEX(m,k)];

    if (m != k)
    {
      pivot[n-1] = -pivot[n-1];
      /* Interchange rows m and k for the lower partition. */
      for (j = k; j < n; ++j)
      {
        pa = a+AINDEX(m,j); pb = a+AINDEX(k,j);
        t = *pa; *pa = *pb; *pb = t;
      }
    }
    /* row k is now the pivot row */

    /* Bail out if pivot is too small */
    if (fabs(pvt) < anorm * EPSILON)
    {
      /* Singular or nearly singular */
      *cond = 1.0e+32;
      *flag = 3;
      goto DecompExit;
    }

    /* eliminate the lower matrix partition by rows
  and store the multipliers in the k sub-column */
    for (i = k+1; i < n; ++i)
    {
      pa = a+AINDEX(i,k);          /* element to eliminate */
      t = -( *pa / pvt );          /* compute multiplier   */
      *pa = t;                     /* store multiplier     */
      for (j = k+1; j < n; ++j)    /* eliminate i th row */
      {
        if (fabs(t) > anorm * EPSILON)
          a[AINDEX(i,j)] += a[AINDEX(k,j)] * t;
      }
    }

  }  /* End of Gaussian elimination. */

  /* cond = (1-norm of a)*(an estimate of 1-norm of a-inverse)
 estimate obtained by one step of inverse iteration for the
 small singular vector. This involves solving two systems
 of equations, (a-transpose)*y = e and a*z = y where e
 is a vector of +1 or -1 chosen to cause growth in y.
 estimate = (1-norm of z)/(1-norm of y)

 Solve (a-transpose)*y = e   */

  for (k = 0; k < n; ++k)
  {
    t = 0.0;
    if (k != 0)
    {
      for (i = 0; i < k; ++i)  t += a[AINDEX(i,k)] * work[i];
    }
    if (t < 0.0) ek = -1.0; else  ek = 1.0;
    pa = a+AINDEX(k,k);
    if (fabs(*pa) < anorm * EPSILON)
    {
      /* Singular */
      *cond = 1.0e+32;
      *flag = 3;
      goto DecompExit;
    }

    work[k] = -(ek + t) / *pa;
  }

  for (k = n-2; k >= 0; --k)
  {
    t = 0.0;
    for (i = k+1; i < n; i++)
      t += a[AINDEX(i,k)] * work[i];
    /* we have used work[i] here, however the use of work[k]
     makes some difference to cond */
    work[k] = t;
    m = pivot[k];
    if (m != k) { t = work[m]; work[m] = work[k]; work[k] = t; }
  }

  ynorm = 0.0;
  for (i = 0; i < n; ++i) ynorm += fabs(work[i]);

  /* --- solve a * z = y */
  solve(n, ndim, a, work, pivot);

  znorm = 0.0;
  for (i = 0; i < n; ++i) znorm += fabs(work[i]);

  /* --- estimate condition --- */
  *cond = anorm * znorm / ynorm;
  if (*cond < 1.0) *cond = 1.0;
  if (*cond + 1.0 == *cond) *flag = 3;

DecompExit:
  if (work != NULL) { free (work); work = (long double *) NULL; }
  return (0);
}   /* --- end of function decomp() --- */
