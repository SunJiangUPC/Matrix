/*
* ���ܾ�����㷽��
*/
#include <stdio.h>
#include <stdlib.h>

class DenseMatrix
{
public:
	DenseMatrix(int M, int N);
	DenseMatrix(int M, int N, double constval);
	DenseMatrix(int M, int N, double* vec);
	DenseMatrix(DenseMatrix& B);
	~DenseMatrix();
	
	int M = 0;
	int N = 0;
	double* A = NULL;

	DenseMatrix operator=(DenseMatrix& B);
};

// �����ܾ��󷽳�: A*x = b
void solveDenseMatrix(const int* N, const DenseMatrix* A, const double* b, double* x);
