#include "functions_for_course_work.h"
#include <cmath>

namespace mashkin
{
  double function_for_integrate(double x)
  {
    return std::sqrt(3 -x) * std::cos(x);
  }

  double function_for_search_min_x(double x)
  {
    return  std::cos(x)/std::sin(x) + (x/(1 - std::pow(x, 2)));
  }
}
