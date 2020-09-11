#include "globals.h"
#include <cmath> 


mpreal secant_accelerated(mpreal xm1, mpreal x0, unsigned long long R) {
  // This implements the 2nd order accelerated Secant method presented in the
  // paper XXXXXXXXXXXXXXXXXXXXXXXXXXXX
  // x is the location on the x axis, and f is the value of the fcn there.
  
  mpreal::set_default_prec(PREC);
  mpreal tol = TOL;
  
  mpreal x1, x2, x3;
  mpreal fm1, f0, f1, f2;
  mpreal fm1inv, f0inv, f1inv, f2inv;
  mpreal am1, a0, a1, a2;
  mpreal tm1, t0, t1, t2;  
  mpreal num, den;

  // Now initialize computation by computing x_l and x_r
  printf("Entered secant ....\n");
  mpfr_printf("xm1 = %10Re \n", xm1);
  mpfr_printf("x0 = %10Re \n", x0);  
  fm1 = f(xm1, R);
  f0 = f(x0, R);  
  mpfr_printf("fm1 = %10Re \n", fm1);
  mpfr_printf("f0 = %10Re \n", f0);  
  
  // Must verify fm1 and f0 have different signs
  if (sgn(fm1) == sgn(f0)) {
    printf("Function f at left and right walls have the same sign!  Throwing exception.\n");
    throw;
  }
  
  // First pass, do a normal secant step and get first estimated
  // root x1.
  fm1inv = mpreal(1.0)/fm1;
  f0inv = mpreal(1.0)/f0;
  tm1 = fm1inv;
  t0 = f0inv;
  num = -xm1*tm1 + x0*t0;
  den = tm1 + t0;
  x1 = num/den;

  // Second pass, do 2nd order step, get next estimated root x2.
  f1 = f(x1, R);
  f1inv = mpreal(1.0)/f1;
  am1 = -(x0-x1);
  a0 = (xm1-x1);
  a1 = -(xm1-x0);
  tm1 = am1*fm1inv;
  t0 = a0*f0inv;
  t1 = a1*f1inv;

  num = xm1*tm1 + x0*t0 + x1*t1;
  den = tm1 + t0 - t1;
  
  x2 = num/den;
  
  // Now start 3nd order secant iteration
  for (int cnt=0; cnt<20; cnt++) {
    f2 = f(x2,R);
    f2inv = mpreal(1.0)/f2;    

    am1 = -(x0-x1)*(x0-x2)*(x1-x2);
    a0 = (xm1-x1)*(xm1-x2)*(x1-x2);
    a1 = -(xm1-x0)*(xm1-x2)*(x0-x2);
    a2 = (xm1-x0)*(xm1-x1)*(x0-x1);

    tm1 = am1*fm1inv;
    t0 = a0*f0inv;
    t1 = a1*f1inv;
    t2 = a2*f2inv;    

    num = xm1*tm1 + x0*t0 + x1*t1 + x2*t2;
    den = tm1 + t0 + t1 + t2;;
  
    x3 = num/den;


    // Check if difference between x3 and x2 is small enough.
    if (abs(x3 - x2) < tol) { 
      // Converged!
      printf("Converged after %d iterations.\n", cnt);
      return (x3 + x2)/mpreal(2.0);
    }

    // Else update values and run another loop.
    xm1 = x0;
    x0 = x1;
    x1 = x2;
    x2 = x3;
    fm1inv = f0inv;
    f0inv = f1inv;
    f1inv = f2inv;
    
  } // for (int cnt=0; cnt<20; cnt++)
  
  printf("Convergence failed!  Throwing exception.\n");
  throw;
  
}
