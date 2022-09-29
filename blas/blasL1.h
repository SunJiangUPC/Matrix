/*
* BLAS Level 1 funcition
* level1: vector-vector, ��������������
* 
* x, y, z: ����
* a, b, c: ����
* A, B, C: ����
*/
#ifndef BLASL1_H
#define BLASL1_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C"
{
#endif
	// y = x
	void dcopy(const int* N, const double* x, double* y);
	// y = a * x
	void dscal(const int* N, const double* a, const double* x, double* y);
	// a = x dot y
	double ddot(const int* N, const double* x, const double* y);
	// a = sum( x(i) * x(i))
	double dnorm(const int* N, const double* x);
	// y := a*x + b*y
	void daxpby(const int* N, const double* a, const double* x, const double* b, double* y);
	//
	

#ifdef __cplusplus
}
#endif



#endif