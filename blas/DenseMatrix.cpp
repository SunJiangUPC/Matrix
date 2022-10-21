#include "DenseMatrix.h"

DenseMatrix::DenseMatrix(int M, int N)
{
}

DenseMatrix::DenseMatrix(int M, int N, double constval)
{


}

DenseMatrix::DenseMatrix(int M, int N, double* vec)
{

}

DenseMatrix::DenseMatrix(DenseMatrix& B)
{

}

DenseMatrix::~DenseMatrix()
{
	M = 0;
	N = 0;
	free(A);
}

DenseMatrix DenseMatrix::operator=(DenseMatrix& B)
{
	if (this != &B) { // 是否为本身赋值
		this->M = B.M;
		this->N = B.N;
		this->A = B.A;
		return *this;
	}

}


// 求解稠密矩阵方程: A*x = b
void solveDenseMatrix(const int* N, const DenseMatrix* A, const double* b, double* x);