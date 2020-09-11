=== What is this?

This is a high-precision computation of the accumulation point,
lambda_inf, of the period 2^N cycles of the logistic map.  Basic info
about the accumulation point is presented at
https://mathworld.wolfram.com/LogisticMap.html so I won't repeat it
here. 

This computation finds lambda_inf by computing the bifurcation point
for increasing values of n.  It finds the lambda_n by solving the
nonlinear system, 

x2 = lam*x1*(1-x1)
x3 = lam*x2*(1-x2)
...
x1 = lam*xn*(1-xn)
lam^n*(1-2*x1)*(1-2*x2)*....*(1-2*xn) = -1

where n=2, 4, 8, 16, etc is the number of fixed points.  The first n
equations are found by simply unrolling the logistic map iteration,
and requiring that after n iterations I get the original fixed point
back again.  The last equation expresses the fact that at the
bifurcation point, the slope of the overall iterated map f^(n)(x) = -1
for any fixed point x.  Note that this system has n+1 equations and
n+1 unknowns.

I solve this system using Newton's method.  I first solve the n=2
cycle, then use the resulting x and lam as starting guesses for the
n=4 Newton iteration.  I repeat this procedure for n = 4, 8, 16,
... etc.  The computation takes longer for increasing n.  I can get up
to around period 2^17 but after that Newton's method fails to
converge and/or I run out of memory depending upon my precision and
tolerance settings.

The computation uses MPFR for high-precision numerics, and Eigen to
provide solvers and other linear algebra functions.  The Jacobian
matrix used in Newton's method is sparse to help save memory.  A
dynamically-sized sparse solver is used when solving for the Newton
step. 

Results of the computation are dumped into a file called
"results.txt".  After a run I move the results file into the "results"
directory for post-processing.  Specific post-processing steps involve
checking how many digits match between different runs using different
precision and tolerance settings.  I also post-process using Aitken
sequence acceleration to squeeze more digits out of the estimated
value of lambda_inf.  More information about that procedure is given
in the results directory.

=== Instructions 

The code is written in easy C++.  If you want to repeat these
computations, you need to download and install the following libraries
(and their dependencies): 

MPFR++: http://www.holoborodko.com/pavel/mpfr/
Eigen:  http://eigen.tuxfamily.org/index.php?title=Main_Page

Place both libraries into the main directory as sub-directories.  Then
do "make" and if all goes well you should have an executable called
bifurcation_calc which you can play with.  If you manage to compute
lambda_n values for n > 2^17, please drop me a line and tell me how
you did it.

Stuart Brorson
6.17.2020

