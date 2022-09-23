#include "blasL1.h"

// y := a*x + b*y
void daxpby(const int* N, const double* a, const double* x, const double* b, double* y)
{
	double tempa = *a, tempb = *b;
	for (int i = 0; i < *N; i++)
	{
		y[i] = tempa * x[i] + tempb * y[i];
	}
}
