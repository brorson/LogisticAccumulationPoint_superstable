#include "J.h"

//-------------------------------------------------------------
//void initJ(int R, SparseMatrix<mpreal, Dynamic, Dynamic> J) {
void initJ(int R, DynamicSparseMatrix<mpreal>& J) {  
  // The purpose of this fcn is to initialize J the first
  // time, so from now on out I only need to update its
  // elements when doing Newton's method.
  
  //std::cout << "Entered initJ" << std::endl;

  // Init main part
  for (int i=0; i<(R-2); i++) {
    J.insert(i,i) = 1.0;
    J.insert(i,i+1) = -1.0;
    J.insert(i,R-1) = 1.0;
  }
  J.insert(R-2,R-2) = 1.0;
  J.insert(R-2,0) = -1;
  J.insert(R-2,R-1) = 1.0;

  // Init last row
  for (int k=0; k<(R-1); k++) {
    J.insert(R-1,k) = 1.0;
  }

  // Init last element in last row
  J.insert(R-1,R-1) = 1.0;

  return;
}

//-------------------------------------------------------------
//void J(int R, Matrix<mpreal, Dynamic, 1>& x, SparseMatrix<mpreal, Dynamic, Dynamic> J) {
void J(int R, Matrix<mpreal, Dynamic, 1>& x, DynamicSparseMatrix<mpreal>& J) {  

  //std::cout << "Entered J4" << std::endl;
  
  mpreal mu = x(R-1);
  mpreal z;  // Temp variable
  Matrix<mpreal, Dynamic, 1> xx(R,1);

  // Compute main part
  for (int i=0; i<(R-2); i++) {
    xx(i) = 1-2*x(i);
    J.coeffRef(i,i) = mu*xx(i);
    J.coeffRef(i,i+1) = -1;
    J.coeffRef(i,R-1) = x(i)*(1-x(i));
  }
  xx(R-2) = 1-2*x(R-2);
  J.coeffRef(R-2,R-2) = mu*xx(R-2);
  J.coeffRef(R-2,0) = -1;
  J.coeffRef(R-2,R-1) = x(R-2)*(1-x(R-2));

  // Last row
  // This is O(N^2).  Lots of redundancy.
  // Can I find a better way to do this?
  for (int k=0; k<(R-1); k++) {
    z = 1.0;
    for (int i=0; i<(R-1); i++) {
      if (k != i) {
        z = z*mu*xx(i);
      }
    }
    J.coeffRef(R-1,k) = -2.0*z*mu;
  }

  // Last element in last row
  z = 1.0;
  for (int k=0; k<(R-1); k++) {
    z = z*mu*xx(k);
  }
  J.coeffRef(R-1,R-1) = (R-1)*z/mu;

  return;
}
