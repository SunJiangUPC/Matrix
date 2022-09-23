/*
* ����ϡ�����: CSRѹ����ʽ
* Example:
* A = [1  2  0  0
*      3  4  5  0
*      0  6  7  8
*      0  0  9  10]
* N = 4;
* NNZ = 10;
* AI = [0  2        5        8     10];
* AJ = [1  2  1  2  3  2  3  4  3  4];
* AA = [1  2  3  4  5  6  7  8  9  10];
*/
#ifndef CSRMATRIX_H
#define CSRMATRIX_H
#include <stdlib.h>
#include <string.h>

#include "memory.h"
#include "blas.h"

class CSRMatrix
{
public:
	void initial(char cT, int nN, int nNNZ, int* nI, int* nJ, double* dA);
	void initial(char cT, int nN, int nNNZ);
	void clear();
	char Type;//S-����; L-��������; U-��������
	int N;//�������
	int NNZ;//����Ԫ�ظ���
	int* I;
	int* J;
	double* A;
};

/* 
* ���ϡ�����
* type:
*    1: �Ƿ�Ϊ�����Ǿ���
*    2: �Ƿ�Ϊ�����Ǿ���
*/
bool check_matrix(const CSRMatrix* A, int type);

#endif