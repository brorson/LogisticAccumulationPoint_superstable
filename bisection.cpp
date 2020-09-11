#include "globals.h"
#include <cmath> 


mpreal bisection(mpreal lam_l, mpreal lam_r, unsigned long long R) {

  mpreal tol = TOL;
  
  mpreal lam_m;
  mpreal x_l, x_r, x_m;

  // Now initialize computation by computing x_l and x_r
  x_l = f(lam_l, R);
  x_r = f(lam_r, R);  
  mpfr_printf("x_l = %10Re \n", x_l);
  mpfr_printf("x_r = %10Re \n", x_r);  
  
  // Must verify x_l and x_r have different signs
  if (sgn(x_l) == sgn(x_r)) {
    printf("Left and right walls have the same sign!  Throwing exception.\n");
    throw;
  }
  
  
  // Now start midpoint iteration
  for (int cnt=0; cnt<1000; cnt++) {
    //mpfr_printf("laml_l- lam_r = %10Re \n", lam_l-lam_r);    
    
    lam_m = (lam_r + lam_l)/2.0;
    x_m = f(lam_m, R);
    
    //mpfr_printf("lam_m = %10Rf, x_m = %10Rf \n", lam_m, x_m);      
    
    if (sgn(x_l) == sgn(x_m)) {
      // Move in left wall
      lam_l = lam_m;
      x_l = x_m;
    } else if (sgn(x_m) == sgn(x_r)) {
      // Move in right wall
      lam_r = lam_m;
      x_r = x_m;
    } else {
      // Error!
      printf("Can't figure out which wall to move!  Throwing exception.\n");
      throw;
    }
    
    // Check for convergence
    if ((lam_r - lam_l) < tol) {
      // Converged!
      printf("Converged after %d iterations.\n", cnt);
      return lam_m;
    }
    
    
  } // for (int cnt=0; cnt<25; cnt++)
  
  printf("Convergence failed!  Throwing exception.\n");
  throw;

}
