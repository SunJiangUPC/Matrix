#include "solve.h"
#include "../blas/CSRmatrix.h"
#include "../config/information.h"
#include "../blas/blasL1.h"

// Solve lower matrix: L*y = b
void solveLCSR(const CSRMatrix* L, const double* b, double* y)
{
    // 求解下三角矩阵：CSR格式
    // L* y = b
    // 
    int N = L->N;
    // 检查是否为下三角矩阵
    for (int i = 0; i < N; i++)
    {
        for (int j = L->I[i]; j < L->I[i + 1]; j++)
        {
            if (L->J[j] > i)
                error("No Lower Matrix.");
        }
    }
    //
    double sum = 0.0;
    dcopy(&N, b, y);// y = b
    y[0] = b[0] / L->A[0];
    for (int i = 1; i < N; i++)
    {
        sum = 0.0;
        for (int j = L->I[i]; j < L->I[i + 1] - 1; j++)
        {
            sum += L->A[j] * y[L->J[j]];
        }
        y[i] = (b[i] - sum) / L->A[L->I[i + 1] - 1];
    }
}
