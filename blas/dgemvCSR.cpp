#include "blasL1.h"
#include "blasL2.h"
#include "CSRmatrix.h"

// y := A*x
/*
* Example:
* A = [1  2  0  0
    * 3  4  5  0
    * 0  6  7  8
    * 0  0  9  10]
* AI = [0  2        5        8     10];
* AJ = [1  2  1  2  3  2  3  4  3  4];
* AA = [1  2  3  4  5  6  7  8  9  10];
*/
void dgemvCSR(const CSRMatrix* A, const double* x, double* y)
{
    int N = A->N;
    int* I = A->I;
    int* J = A->J;
    double* TempA = A->A;

    double sum = 0.0;
    // dcopy(&N, x, y);// y = x

    for (int i = 0; i < N; i++)
    {
        sum = 0.0;
        for (int j = I[i]; j < I[i + 1]; j++)
        {
            sum += TempA[j] * x[J[j]];
        }
        y[i] = sum;
    }
}