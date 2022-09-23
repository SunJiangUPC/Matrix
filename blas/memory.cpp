#include "memory.h"

// double数组分配内存
void dvecmal(const int* N, double** x)
{
	*x = (double*)malloc(*N * sizeof(double));
	if (*x == NULL)
	{
		printf("内存分配失败.\n");
		return;
	}
}
// float数组分配内存
void fvecmal(const int* N, float** x)
{
	*x = (float*)malloc(*N * sizeof(float));
	if (*x == NULL)
	{
		printf("内存分配失败.\n");
		return;
	}
}
// int数组分配内存
void ivecmal(const int* N, int** x)
{
	*x = (int*)malloc(*N * sizeof(int));
	if (*x == NULL)
	{
		printf("内存分配失败.\n");
		return;
	}
}
// double/float/int数组内存释放
void dvecfree(double** x)
{
	if (*x)
	{
		free(*x);
		*x = NULL;
	}
}
void fvecfree(float** x)
{
	if (*x)
	{
		free(*x);
		*x = NULL;
	}
}
void ivecfree(int** x)
{
	if (*x)
	{
		free(*x);
		*x = NULL;
	}
}