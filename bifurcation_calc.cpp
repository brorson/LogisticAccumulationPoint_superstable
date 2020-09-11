#include <algorithm> 
#include "globals.h"
#include "f.h"
#include "J.h"


#define PREC 200

mpreal tol = 1e-50;

void bifurcation_calc(int, mpreal&, mpreal&);



//-------------------------------------------------
int main(int argc, char* argv[])
{

  mpreal::set_default_prec(PREC);
  mpreal x0;
  mpreal lam;
  int R;
  int ord;

  // Initialize calc
  ord = 4;
  x0 = 8.0e-01;
  lam = 3.4;

  FILE *outfile;
  outfile = fopen("results.txt", "w"); 

  fprintf(outfile, "PREC = %d, tol = %e\n", PREC, tol);
  // outfile << "PREC = " << PREC << ", tol = " << tol << "\n";
  
  // This just loops through each bifurcation point.
  for (int i=0; i<21; i++) {
    printf("=============================================\n");
    mpfr_printf("Calling with order = %d, lam = %10Rf, x0 = %10Rf \n", ord, lam, x0);
    R = ord+1;
    bifurcation_calc(R, lam, x0);
    mpfr_printf ("Computed order = %d, lam = %10Re, x0 = %10Rf \n", ord, lam, x0);

    mpfr_fprintf(outfile, "order = %d, lambda = %10Re \n", ord, lam);
    fflush(outfile);
    
    // Go to next bifurcation order
    ord = 2*ord;
  }

  // All done.
  fclose(outfile);
  
}


//-------------------------------------------------
void bifurcation_calc(int R, mpreal& lam, mpreal& x0) {
  mpreal::set_default_prec(PREC);  

  int i, j;
  Matrix<mpreal, Dynamic, 1> fn(R,1);
  DynamicSparseMatrix<mpreal> Jn(R, R);
  Matrix<mpreal, Dynamic, 1> delta(R,1);
  Matrix<mpreal, Dynamic, 1> xn(R, 1);
  SparseLU<SparseMatrix<mpreal>, COLAMDOrdering<int> >  solver;

  // Initialize J
  initJ(R, Jn);
  
  // Iterate create starting vector for Newton's method.
  mpreal ex = 3.3*log10(mpreal(R-1));   // log base 2 of R-1
  mpreal base = 0.21;
  mpreal dlam = pow(base, ex);
  lam = lam + dlam;
  printf("dlam = \n");
  mpfr_printf("%10Re \n", dlam);
  printf("starting lam = \n");
  mpfr_printf("%10Re \n", lam);
  for (i=0; i<10*R; i++) {
    x0 = lam*x0*(1-x0);
  }

  // Now iterate R-1 times more and collect x values
  for (i=0; i<R-1; i++) {
    x0 = lam*x0*(1-x0);
    xn(i) = x0;
  }
  xn(R-1) = lam;
  
  // Now compute the difference between elements in this vector
  // and print out the smallest difference
  // Must do this operation on a copy
  Matrix<mpreal, Dynamic, 1> x1(R,1);
  std::copy(xn.begin(), xn.end()-1, x1.begin());
  std::sort(x1.begin(), x1.end());
  // Just do diff in place to save memory
  for (i=0; i<R-1; i++) {
    x1[i] = abs(x1[i]-x1[i+1]);
  }
  mpreal diff = 1.0;
  for (i=0; i<R-1; i++) {
    if (x1[i]<diff) {
      diff = x1[i];
    }
  }
  mpfr_printf("Minimum difference between zeros = %10Re \n", diff);  
  
  
  
  //--------------------------------------
  // Run a Newton's method loop
  for (int cnt=0; cnt<75; cnt++) {

    //printf("------------------------------------\n");
    
    // Compute f
    f(R, xn, fn);

    // Compute J
    J(R, xn, Jn);
    
    // Now compute delta = Newton step
    if (cnt == 0) {
      // configure solver on first iteration.
      solver.analyzePattern(Jn); 
    }
    solver.factorize(Jn);
    delta = solver.solve(fn); 

    // Take step
    for (i=0; i<R; i++) {    
      xn(i) = xn(i) - delta(i);
    }

    // Check if we're close enough to quit
    //mpfr_printf("norm(delta) = %10Re\n", delta.norm());
    if (delta.norm() < tol) {
      printf("  ---  Terminating after %d iterations because norm(delta) < tol.\n", cnt);
      lam = xn(R-1);
      x0 = xn(0);
      return;
    } else {
      //printf("Must take another loop iteration.\n");
    }
    
  } // for (int cnt=0; cnt<25; cnt++)

  printf("Convergence failed!  Throwing exception.\n");
  throw;

}
