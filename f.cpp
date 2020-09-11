#include "f.h"


void f(int R, Matrix<mpreal, Dynamic, 1>& x, Matrix<mpreal, Dynamic, 1>& f) {

  //std::cout << "Entered f4." << std::endl;
  
  mpreal mu = x(R-1);

  //printf("x = \n");
  //for (int i=0; i<R; i++) {
  //  mpfr_printf ("%10Rf \n", x(i));
  //}
  //mpfr_printf("mu = %10Rf \n", mu);
  
  for (int i=0; i<(R-2); i++) {
    f(i) = mu*x(i)*(1.0-x(i)) - x(i+1);
  }
  f(R-2) = mu*x(R-2)*(1.0-x(R-2)) - x(0);
    
  mpreal mp = 1.0;
  for (int i=0; i<(R-1); i++) {
    mp = mp*mu*(1-2*x(i));
  }

  f(R-1) = mp+1;
  
  return;

}
