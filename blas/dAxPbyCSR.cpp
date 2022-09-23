#include "blasL2.h"
#include "CSRmatrix.h"

// y := A*x + b*y
void dAxPbyCSR(const CSRMatrix* A, const double* x, const double* b, double* y)
{
    int N = A->N;
    int* I = A->I;
    int* J = A->J;
    double* TempA = A->A;
    double tempb = *b;

    double sum = 0.0;

    for (int i = 0; i < N; i++)
    {
        sum = 0.0;
        for (int j = I[i]; j < I[i + 1]; j++)
        {
            sum += TempA[j] * x[J[j]];
        }
        y[i] = sum + tempb * y[i];
    }
}