#include "blasL1.h"

// a = sum( x(i) * x(i))
double dnorm(const int* N, const double* x)
{
	double dNorm = 0.0;
	for (int i = 0; i < *N; i++)
	{
		dNorm += x[i] * x[i];
	}
	return dNorm;
}