#include "globals.h"

//void initJ(int R, SparseMatrix<mpreal, Dynamic, Dynamic> J);
//void J(int R, Matrix<mpreal, Dynamic, 1>& x, SparseMatrix<mpreal, Dynamic, Dynamic> J);

void initJ(int R, DynamicSparseMatrix<mpreal>& J);
void J(int R, Matrix<mpreal, Dynamic, 1>& x, DynamicSparseMatrix<mpreal>& J);
