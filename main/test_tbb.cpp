#include "counter.h"
#include <timer.h>
#include <oneapi/dpl/algorithm>
#include <oneapi/dpl/execution>
#include <oneapi/tbb/global_control.h>
#include <vector>

double 
kernel (const double& x, const double& y) {
  return x + y;
};

int main () {

  constexpr size_t N = 100000l * 5000l;
  std::vector<double> a(N);
  std::vector<double> b(N);
  std::vector<double> c(N);
  
  cdf::timer::timer_t timer;

  std::fill (a.begin (), a.end (), 1.0);
  std::fill (b.begin (), b.end (), 2.0);
  std::fill (c.begin (), c.end (), 0.0);
  
  timer.tic ("seriale");
  for (size_t ii = 0; ii < N; ++ii) {
    c[ii] += a[ii] + b[ii];
  }
  timer.toc ("seriale");
  std::cout << a[N-1] << " " << b[N-1] << " " << c[N-1] << std::endl;
  
  std::fill (a.begin (), a.end (), 1.0);
  std::fill (b.begin (), b.end (), 2.0);
  std::fill (c.begin (), c.end (), 0.0);
  
  timer.tic ("par index based");
  range<size_t> r (0, N);
  dpl::for_each (std::execution::par, r.begin (), r.end (), [&] (const size_t& ii) { c[ii] += kernel (a[ii], b[ii]); });
  timer.toc ("par index based");
  std::cout << a[N-1] << " " << b[N-1] << " " << c[N-1] << std::endl;
  
  std::fill (a.begin (), a.end (), 1.0);
  std::fill (b.begin (), b.end (), 2.0);
  std::fill (c.begin (), c.end (), 0.0);
  
  timer.tic ("par transform");
  dpl::transform (std::execution::par, a.begin (), a.end (), b.begin (), c.begin (), kernel);
  timer.toc ("par transform");
  std::cout << a[N-1] << " " << b[N-1] << " " << c[N-1] << std::endl;

  std::fill (a.begin (), a.end (), 1.0);
  std::fill (b.begin (), b.end (), 2.0);
  std::fill (c.begin (), c.end (), 0.0);
  
  timer.tic ("par_unseq index based");
  dpl::for_each (std::execution::par_unseq, r.begin (), r.end (), [&] (const size_t& ii) { c[ii] += kernel (a[ii], b[ii]); });
  timer.toc ("par_unseq index based");
  std::cout << a[N-1] << " " << b[N-1] << " " << c[N-1] << std::endl;

  std::fill (a.begin (), a.end (), 1.0);
  std::fill (b.begin (), b.end (), 2.0);
  std::fill (c.begin (), c.end (), 0.0);
  
  timer.tic ("par_unseq transform");
  dpl::transform (std::execution::par, a.begin (), a.end (), b.begin (), c.begin (), kernel);
  timer.toc ("par_unseq transform");
  std::cout << a[N-1] << " " << b[N-1] << " " << c[N-1] << std::endl;

  std::fill (a.begin (), a.end (), 1.0);
  std::fill (b.begin (), b.end (), 2.0);
  std::fill (c.begin (), c.end (), 0.0);
  
  timer.tic ("unseq transform");
  dpl::transform (std::execution::unseq, a.begin (), a.end (), b.begin (), c.begin (), kernel);
  timer.toc ("unseq transform");
  std::cout << a[N-1] << " " << b[N-1] << " " << c[N-1] << std::endl;

  
  return 0;
};
