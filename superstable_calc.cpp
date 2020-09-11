#include "globals.h"

using namespace mpfr;
using namespace std;
using namespace Eigen;


//-----------------------------------------------
// Sign function stolen from
// https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
// Used in bisection calc.
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


//-----------------------------------------------
mpreal f(mpreal lam, unsigned long long R) {
  
  mpreal::set_default_prec(PREC);
  mpreal x = mpreal("0.5");
  mpreal x0 = mpreal("0.5");
  mpreal one = mpreal("1.0");
  
  for (unsigned long i=1; i<R+1; i++) {
    // std::cout << i << std::endl;
    x = lam*x*(one-x);
  }

  mpreal diff = x-x0;
  return diff;
}
  

//-------------------------------------------------
int main(int argc, char* argv[])
{
  mpreal tol = TOL;
  mpreal::set_default_prec(PREC);
  mpreal lam;   // lam is the control parameter
  mpreal newlam;
  int N;
  int R;        // R = 2^N -- order of fixed point

  // Initialize calc
  N = 2;
  R = (unsigned long long) std::pow(2.0, N);
  lam = mpreal("3.25");  // Approx location of superstable pt for R=2.

  // Open output file.  Give it name related to requested PREC.
  string filenamebase = "results_superstable_prec";
  int iprec = PREC;
  string sprec= to_string(iprec);
  string filename = filenamebase + sprec + ".txt";
  cout << "Opening filename = " << filename << endl;
  std::ofstream outfile;
  outfile.open (filename.c_str(), std::ofstream::out);

  // Set precision for mpfr output.  Divide by 3 since PREC is in bits
  // but I want decimal digits.
  outfile.precision(PREC/3);
  
  outfile << "PREC = " << PREC << ", tol = " << tol << "\n";
  
  // This just loops through each supertable point.
  for (int i=0; i<38; i++) {

    printf("=============================================\n");
    mpfr_printf("Calling with order = %d, lam = %10Rf \n", N, lam);
    auto t0 = std::chrono::high_resolution_clock::now();
    newlam = superstable_calc(N, lam);
    mpfr_printf ("Computed order = %d, lam = %10Re \n", R, newlam);
    auto t1 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(t1-t0);
    cout << "Calculation time = " << duration.count() << std::endl;
    outfile << "order = " << N << ", lambda = " << newlam << endl;
    outfile.flush();
    
    // Go to next period
    lam = newlam;
    N = N+1;
    R = (unsigned long) std::pow(2.0, N);
  }

  // All done.
  outfile.close();
  
}


//-------------------------------------------------
mpreal superstable_calc(int N, mpreal& lam0) {
  // The goal of this program is to find the lambda corresponding
  // to the superstable point using root finding.
  // lam_l is left hand lambda.

  mpreal::set_default_prec(PREC);  

  int i, j;
  mpreal lam_l, lam_r, lam_m;
  mpreal x_l, x_r, x_m;
  unsigned long long R;
  mpreal four = mpreal("4.000");    
  
  R = (unsigned long long) std::pow(2.0, N);
  
  // Compute left and right lambdas.
  // First estimate width of new interval
  mpreal dlam = mpfr::pow(mpreal("0.216"), static_cast<mpreal>(N-1), MPFR_RNDN);
  // Next estimate pos of superstable point in new interval.
  lam_m = lam0+dlam;
  mpfr_printf("Starting laml_m = %10Re\n dlam = %10Re \n", lam_m, dlam);
  
  // Now set left and right search walls based on these estimations.
  lam_l = lam_m - dlam/four;
  lam_r = lam_m + dlam/four;
  mpfr_printf("Starting laml_l = %10Rf, lam_r = %10Rf \n", lam_l, lam_r);

  // Now call preferred solver method.
  //lam_m = bisection(lam_l, lam_r, R);
  //lam_m = secant(lam_l, lam_r, R);
  lam_m = secant_accelerated(lam_l, lam_r, R);    

  return lam_m;
  
}
