#include "memory.h"

// double��������ڴ�
void dvecmal(const int* N, double** x)
{
	*x = (double*)malloc(*N * sizeof(double));
	if (*x == NULL)
	{
		printf("�ڴ����ʧ��.\n");
		return;
	}
}
// float��������ڴ�
void fvecmal(const int* N, float** x)
{
	*x = (float*)malloc(*N * sizeof(float));
	if (*x == NULL)
	{
		printf("�ڴ����ʧ��.\n");
		return;
	}
}
// int��������ڴ�
void ivecmal(const int* N, int** x)
{
	*x = (int*)malloc(*N * sizeof(int));
	if (*x == NULL)
	{
		printf("�ڴ����ʧ��.\n");
		return;
	}
}
// double/float/int�����ڴ��ͷ�
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