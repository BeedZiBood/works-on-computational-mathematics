#include "RungeKutta.h"
#include <iostream>
#include <iomanip>
#include <cstddef>
#include <utility>
#include <cmath>
#include "outputStructs.h"

namespace mashkin
{
  MyArr::MyArr(std::vector< double >& vect):
    x(vect)
  {
  }

  MyArr MyArr::operator=(MyArr&& rhs)
  {
    x = rhs.x;
    return *this;
  }

  MyArr MyArr::operator=(MyArr& rhs)
  {
    *this = std::move(rhs);
    return *this;
  }

  double& MyArr::operator[](size_t ind)
  {
    return this->x[ind];
  }

  MyArr MyArr::operator/(double num)
  {
    MyArr result = *this;
    for (size_t i = 0; i < result.x.size(); i++)
    {
      result[i] /= num;
    }
    return result;
  }

  MyArr MyArr::operator+(MyArr&& rhs)
  {
    MyArr result = *this;
    for (size_t i = 0; i < result.x.size(); i++)
    {
      result[i] += rhs[i];
    }
    return result;
  }

  MyArr MyArr::operator+(double num)
  {
    MyArr result = *this;
    for (size_t i = 0; i < result.x.size(); i++)
    {
      result[i] += num;
    }
    return result;
  }

  MyArr MyArr::operator*(double num)
  {
    MyArr result = *this;
    for (size_t i = 0; i < result.x.size(); i++)
    {
      result[i] *= num;
    }
    return result;
  }

  double func1(double t, MyArr&& x)
  {
    return -130 * x[0] + 900 * x[1] + std::exp(-10 * t);
  }

  double func2(double t, MyArr&& x)
  {
    double var = 30 * x[0] - 300 * x[1] + std::log(1 + 100 * t * t);
    return var;
  }

  size_t MyArr::size()
  {
    return x.size();
  }

  MyArr rungeKutta3degree(double t, std::vector< double >& x, double h)
  {
    std::vector< double > zero_vect{0, 0};
    MyArr k1(zero_vect);
    MyArr k2(zero_vect);
    MyArr k3(zero_vect);
    MyArr zn(x);
    k1[0] = h * func1(t, std::move(zn));
    k1[1] = h * func2(t, std::move(zn));
    k2[0] = h * func1(t + h / 2, zn + k1 / 2);
    k2[1] = h * func2(t + h / 2, zn + k1 / 2);
    k3[0] = h * func1(t + 3 * h / 4, zn + (k2 * 3) / 4);
    k3[1] = h * func2(t + 3 * h / 4, zn + (k2 * 3) / 4);
    zn = zn + (k1 * 2 + k2 * 3 + k3 * 4) / 9;
    for (size_t i = 0; i < x.size(); i++)
    {
      x[i] = zn[i];
    }
    return zn;
  }

  std::ostream& operator<<(std::ostream& out, MyArr&& dest)
  {
    std::ostream::sentry sentry(out);
    if (!sentry)
    {
      return out;
    }
    iofmtguard fmtguard(out);
    for (int i = 0; i < dest.size(); i++)
    {
      out << dest[i];
      if (i + 1 != dest.size())
      {
        out << std::setw(13);
      }
    }
    return out;
  }
}
