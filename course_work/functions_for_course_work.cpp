#include "functions_for_course_work.h"
#include <cmath>

namespace mashkin
{
  double function_for_search_min_x(double x)
  {
    return  1.0/std::tan(x) + x/(1 - x * x);
  }
}
