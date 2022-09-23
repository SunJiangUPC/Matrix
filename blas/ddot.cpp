#include "blasL1.h"

// a = x dot y
double ddot(const int* N, const double* x, const double* y)
{
	double dDot = 0.0;
	for (int i = 0; i < *N; i++)
	{
		dDot += x[i] * y[i];
	}
	return dDot;
}