#ifndef LAB3_RUNGEKUTTA_H
#define LAB3_RUNGEKUTTA_H
#include "Function_for_laboratory.h"
#include <cstddef>
#include <iostream>
#include <vector>

namespace mashkin
{
  class MyArr
  {
  public:
    MyArr() = default;
    ~MyArr() = default;
    MyArr(double* arr);
    MyArr(MyArr& rhs) = default;
    MyArr(MyArr&& rhs) = default;

    MyArr operator/(double num);
    MyArr operator=(MyArr& rhs);
    MyArr operator=(MyArr&& rhs);
    MyArr operator*(double num);
    MyArr operator+(MyArr&& rhs);
    MyArr operator+(double num);

    double& operator[](size_t ind);
  private:
    std::vector< double > x;
  };

  std::ostream& operator<<(std::ostream& out, MyArr& dest);
  MyArr rungeKutta3degree(double t, double* x, double h);
}
#endif
