#include "solve.h"
#include "../blas/CSRmatrix.h"
#include "../config/information.h"
#include "../blas/blasL1.h"

// Solve upper matrix: U*x = y
void solveUCSR(const CSRMatrix* U, const double* y, double* x)
{
    // ��������Ǿ���CSR��ʽ
    // U * x = y
    // ����Ƿ�Ϊ�����Ǿ���
    int N = U->N;
    for (int i = 0; i < N; i++)
    {
        for (int j = U->I[i]; j < U->I[i + 1]; j++)
        {
            if (U->J[j] < i)
                error("No Upper Matrix.");
        }
    }
    //
    double sum = 0.0;
    dcopy(&N, y, x);// x = y
    x[N - 1] = y[N - 1] / U->A[U->I[N] - 1];
    for (int i = (N - 2); i >= 0; i--)
    {
        sum = 0.0;
        for (int j = U->I[i] + 1; j < U->I[i + 1]; j++)
        {
            sum += U->A[j] * x[U->J[j]];
        }
        x[i] = (y[i] - sum) / U->A[U->I[i]];
    }
}
