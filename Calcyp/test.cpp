#include <test.h>


Uniform::Uniform(double min_, double max_)
{
  num =  calculate(5);
}
int Uniform::calculate(int num)
{
  return(num*2);
}

RCPP_MODULE(unif_module) {
  class_<Uniform>( "Uniform" )
  .constructor<double,double>()
  .field( "min", &Uniform::min )
  .field( "max", &Uniform::max )
  .field( "num", &Uniform::num )
  ;
}