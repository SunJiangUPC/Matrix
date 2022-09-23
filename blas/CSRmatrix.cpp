#include "CSRmatrix.h"

void CSRMatrix::clear()
{
	free(I);
	free(J);
	free(A);
	Type = '\0';
	N = 0;
	NNZ = 0;
	I = NULL;
	J = NULL;
	A = NULL;
}

void CSRMatrix::initial(char cT, int nN, int nNNZ, int* nI, int* nJ, double* dA)
{
	Type = cT;
//#pragma warning(suppress : 4996)
//	strcpy(&Type, &cT);
	N = nN;
	NNZ = nNNZ;
	nN += 1;
	ivecmal(&nN, &I);
	ivecmal(&nNNZ, &J);
	dvecmal(&nNNZ, &A);
	for (int i = 0; i < nN; i++)
	{
		I[i] = nI[i];
	}
	for (int i = 0; i < nNNZ; i++)
	{
		J[i] = nJ[i];
		A[i] = dA[i];
	}
	
}

void CSRMatrix::initial(char cT, int nN, int nNNZ)
{
	Type = cT;
	N = nN;
	NNZ = nNNZ;
	nN += 1;
	ivecmal(&nN, &I);
	ivecmal(&nNNZ, &J);
	dvecmal(&nNNZ, &A);
	for (int i = 0; i < nN; i++)
	{
		I[i] = 0;
	}
	for (int i = 0; i < nNNZ; i++)
	{
		J[i] = 0;
		A[i] = 0.0;
	}

}