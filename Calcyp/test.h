#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp14)]]
class Uniform {
public:
  Uniform(double , double );  
  int calculate(int);
  
  NumericVector draw(int n) const {
    RNGScope scope;
    return runif( n, min, max );
  }
  double min, max;
  int num;
};

