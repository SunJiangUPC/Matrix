#include "vector.h"

void dzeros(const int* N, double* x)
{
	double zero = 0.0;
	if (x)
	{
		for (int i = 0; i < *N; i++)
		{
			x[i] = zero;
		}
	}
	else
	{
		printf("ERROR: 数组为NULL,请先分配内存!(File:%s, Line:%d)\n", __FILE__, __LINE__);
	}
}
void fzeros(const int* N, float* x)
{
	float zero = 0.0;
	if (x)
	{
		for (int i = 0; i < *N; i++)
		{
			x[i] = zero;
		}
	}
	else
	{
		printf("ERROR: 数组为NULL,请先分配内存!(File:%s, Line:%d)\n", __FILE__, __LINE__);
	}
}
void izeros(const int* N, int* x)
{
	int zero = 0;
	if (x)
	{
		for (int i = 0; i < *N; i++)
		{
			x[i] = zero;
		}
	}
	else
	{
		printf("ERROR: 数组为NULL,请先分配内存!(File:%s, Line:%d)\n", __FILE__, __LINE__);
	}
}

void dones(const int* N, double* x)
{
	double one = 1.0;
	if (x)
	{
		for (int i = 0; i < *N; i++)
		{
			x[i] = one;
		}
	}
	else
	{
		printf("ERROR: 数组为NULL,请先分配内存!(File:%s, Line:%d)\n", __FILE__, __LINE__);
	}
}
void fones(const int* N, float* x)
{
	float one = 1.0;
	if (x)
	{
		for (int i = 0; i < *N; i++)
		{
			x[i] = one;
		}
	}
	else
	{
		printf("ERROR: 数组为NULL,请先分配内存!(File:%s, Line:%d)\n", __FILE__, __LINE__);
	}
}
void iones(const int* N, int* x)
{
	int one = 1;
	if (x)
	{
		for (int i = 0; i < *N; i++)
		{
			x[i] = one;
		}
	}
	else
	{
		printf("ERROR: 数组为NULL,请先分配内存!(File:%s, Line:%d)\n", __FILE__, __LINE__);
	}
}

void resize(const int* N, double* x, const int* size)
{
	if (*size == *N)
	{
		return;
	}
	else
	{
		x = (double*)realloc(x, *size * sizeof(double));
		if (!x)
		{
			printf("ERROR: 内存分配失败!(File:%s, Line:%d)\n", __FILE__, __LINE__);
		}
	}
}

void cut(const int* N, const double* x, const int* size, double* y)
{
	if (!x || !y)
	{
		printf("ERROR: 未分配内存!(File:%s, Line:%d)\n", __FILE__, __LINE__);
		return;
	}
	if (*size <= *N)
	{
		for (int i = 0; i < *size; i++)
		{
			y[i] = x[i];
		}
	}
	else
	{
		for (int i = 0; i < *N; i++)
		{
			y[i] = x[i];
		}
		for (int i = *N; i < *size; i++)
		{
			y[i] = 0.0;
		}
	}
}


bool isfiniteall(const int* N, const double* x)
{
	for (int i = 0; i < *N; i++)
	{
		if (!isfinite(x[i]))
		{
			return false;
		}
	}
	return true;
}