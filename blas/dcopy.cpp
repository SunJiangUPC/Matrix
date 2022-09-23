#include "blasL1.h"

// y = x
void dcopy(const int* N, const double* x, double* y)
{
	if (!y)
	{
		printf("����δ�����ڴ�!(File: %s, Line: %d)\n", __FILE__, __LINE__);
		return;
	}
	for (int i = 0; i < *N; i++)
	{
		y[i] = x[i];
	}
}