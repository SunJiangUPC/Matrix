#include "information.h"
#include "../blas/CSRmatrix.h"
#include "../blas/vector.h"

void error(const char* str)
{
	printf("%s\n", str);
}

void warning(const char* str)
{
	printf("%s\n", str);
}


void printMatrix(CSRMatrix* A)
{
	if (A->N > 20)
	{
		printf("矩阵维度过大,不输出显示!\n");
		return;
	}
	double* tempcol;
	tempcol = (double*)malloc(A->N * sizeof(double));
	if (!tempcol)
		printf("内存分配失败!\n");
	else
	{
		
		for (int i = 0; i < A->N; i++)
		{
			dzeros(&A->N, tempcol);
			for (int j = A->I[i]; j < A->I[i + 1]; j++)
			{
				tempcol[A->J[j]] = A->A[j];
			}
			for (int j = 0; j < A->N; j++)
			{
				cout.width(10);
				cout << tempcol[j];
			}
			cout << endl;
		}
	}

	free(tempcol);
}