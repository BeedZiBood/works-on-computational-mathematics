#ifndef LAB3_RKF45_H
#define LAB3_RKF45_H

#include <math.h>
#include <stdlib.h>
#include <functional>
#include <vector>
#define REAL double
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define EXP(x) exp(x)
#define SQRT(x) sqrt(x)
#define SIN(x) sin(x)
#define COS(x) cos(x)
#define TAN(x) tan(x)
#define ATAN(x) atan(x)


inline double absval(double x)
{
  return fabs(x);
}
inline int absval(int x)
{
  return abs(x);
}
#define ABS(x) absval(x)
#define SIGN(a,b) (b<0) ? (-absval(a)) : absval(a)

void RKF45(void(F)(REAL T,REAL*Y,REAL*YP),int NEQN,REAL *Y,REAL &T,REAL TOUT,REAL &RELERR,REAL &ABSERR,REAL *WORK,int &IFLAG);
#endif
