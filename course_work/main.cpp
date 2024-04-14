#include <iostream>
#include "functions_for_course_work.h"
#include "FORSTYTHE.h"

int main()
{
  //find min x for E
  double ans;
  int flag;
  ans = zeroin(3.75, 5.0, mashkin::function_for_search_min_x, 0.0000001, &flag);
  std::cout << "Min ans is " << ans << "\n";
}