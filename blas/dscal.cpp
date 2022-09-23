#include "blasL1.h"

// y = a * x
void dscal(const int* N, const double* a, const double* x, double* y)
{
	double tempa = *a;
	for (int i = 0; i < *N; i++)
	{
		y[i] = tempa * x[i];
	}
}