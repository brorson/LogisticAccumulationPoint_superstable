#include <stdio.h>
#include <iostream>
#include <cstdint>
#include <Eigen/Core>
#include <Eigen/LU>
#include <unsupported/Eigen/MPRealSupport>
#include <unsupported/Eigen/SparseExtra>
#include <mpfr.h>
#include <gmp.h>
#include "mpreal.h"
#include <Eigen/Sparse>

#include <chrono>

using namespace mpfr;
using namespace std;
using namespace Eigen;

// Set number of bits used by MPFR
#define PREC 600

// This is my stopping tolerance for iteration
#define TOL mpreal(1e-150)

// ----------------------------------------------
// Function declarations
mpreal superstable_calc(int, mpreal&);
mpreal f(mpreal lam, unsigned long long R);
//mpreal bisection(mpreal lam_l, mpreal lam_r, unsigned long long R);
//mpreal secant(mpreal lam_l, mpreal lam_r, unsigned long long R);
mpreal secant_accelerated(mpreal lam_l, mpreal lam_r, unsigned long long R);
template <typename T> int sgn(T val);

