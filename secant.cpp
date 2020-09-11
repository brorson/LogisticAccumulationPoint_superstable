#include "globals.h"
#include <cmath> 


mpreal secant(mpreal lam_0, mpreal lam_1, unsigned long long R) {

  mpreal::set_default_prec(PREC);
  mpreal tol = TOL;
  
  mpreal s, lam_2;
  mpreal x_0, x_1, x_2;

  // Now initialize computation by computing x_l and x_r
  printf("Entered secant ....\n");
  mpfr_printf("lam_0 = %10Re \n", lam_0);
  mpfr_printf("lam_1 = %10Re \n", lam_1);  
  x_0 = f(lam_0, R);
  x_1 = f(lam_1, R);  
  mpfr_printf("x_0 = %10Re \n", x_0);
  mpfr_printf("x_1 = %10Re \n", x_1);  
  
  // Must verify x_0 and x_1 have different signs
  if (sgn(x_0) == sgn(x_1)) {
    printf("Left and right walls have the same sign!  Throwing exception.\n");
    throw;
  }
  
  // Now start secant iteration
  for (int cnt=0; cnt<250; cnt++) {
    
    // Check if difference between lam_l and lam_r is small enough.
    if (abs(lam_1 - lam_0) < tol) { 
      // Converged!
      printf("Converged after %d iterations.\n", cnt);
      return (lam_1 + lam_0)/2.0;
    }
    
    // Else iterate again.
    s = (x_1-x_0)/(lam_1-lam_0);
    lam_2 = lam_1 - (x_1/s);
    x_2 = f(lam_2, R);
    
    // Update values
    lam_0 = lam_1;
    lam_1 = lam_2;
    x_0 = x_1;
    x_1 = x_2;
    
  } // for (int cnt=0; cnt<25; cnt++)
  
  printf("Convergence failed!  Throwing exception.\n");
  throw;
  
}
