/*
* 定义稀疏矩阵: CSR压缩格式
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
	char Type;//S-方阵; L-下三角阵; U-上三角阵
	int N;//矩阵阶数
	int NNZ;//非零元素个数
	int* I;
	int* J;
	double* A;
};

/* 
* 检查稀疏矩阵
* type:
*    1: 是否为上三角矩阵
*    2: 是否为下三角矩阵
*/
bool check_matrix(const CSRMatrix* A, int type);

#endif