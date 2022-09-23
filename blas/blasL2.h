/*
* BLAS Level 2 funcition
* level2: matrix-vector, ��������������
*
* x, y, z: ����
* a, b, c: ����
* A, B, C: ����
*/
#ifndef BLASL2_H
#define BLASL2_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

class CSRMatrix;

// y := A*x
void dgemvCSR(const CSRMatrix* A, const double* x, double* y);
// y := A*x + b*y
void dAxPbyCSR(const CSRMatrix* A, const double* x, const double* b, double* y);
//


#endif