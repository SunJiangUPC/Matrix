#ifndef SOLVE_H
#define SOLVE_H
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>

class CSRMatrix;
// ���������Ϣ
class IterInfo
{
public:
    int flag = 0;// ���������
    double relres = 0.0;// ��Բв�
    float iter = 0;// ��������: ���ڰ벽
    double* resvec = NULL;// ÿ���в� resvec = norm(b - A * x)
    void initial(int nLength)
    {
        resvec = (double*)malloc(nLength * sizeof(double));
        if (resvec == NULL) printf("�ڴ����ʧ��!");
    }
    void clear()
    {
        free(resvec);
        resvec = NULL;
    }
    //
    int outeriter = 0;// ������gmres
    int inneriter = 0;
};

// Solve lower matrix: L*y = b
void solveLCSR(const CSRMatrix* L, const double* b, double* y);
// Solve upper matrix: U*x = y
void solveUCSR(const CSRMatrix* U, const double* y, double* x);
// BicgStab algorithm
void bicgstab(IterInfo* info, double* x, const CSRMatrix* A, const double* b, const double* Tol, const int* MaxIter,
	const CSRMatrix* L = NULL, const CSRMatrix* U = NULL, const double* x0 = NULL);
// gmres algorithm
void gmres(IterInfo* info, double* x, const CSRMatrix* A, const double* b, const double* Tol,const int* MaxIter,
    const int* Restart = NULL, const CSRMatrix* L = NULL, const CSRMatrix* U = NULL, const double* x0 = NULL);

// Ԥ������
void ilu0(const CSRMatrix* A, CSRMatrix* L, CSRMatrix* U);

#endif // !SOLVE_H
