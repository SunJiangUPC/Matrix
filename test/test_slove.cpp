/*
* x, y, z: 向量
* a, b, c : 常数
* A, B, C : 矩阵
*/
#pragma warning(disable : 4996)

#include <stdio.h>
#include "../blas/blas.h"
#include "../blas/CSRmatrix.h"
#include "../blas/memory.h"
#include "../blas/vector.h"
#include "../config/information.h"
#include "../solve/solve.h"
#include <iostream>
#include <Windows.h>
#include <string>
#include <fstream>
#include <vector>

using namespace std;

void readCSRMatrixVector(string filename, CSRMatrix* A, vector<double>& b);// 读取压缩矩阵与右边项
void readCSRMatrix(string filename, CSRMatrix* A);// 读取压缩矩阵
void readVector(string filename, vector<double>& b);// 读取右边项

int main(int argnum, char* argvar[])
{
	/*
	* A =
	*[4 3 0 0
	* 2 6 5 0
	* 0 3 7 2
	* 0 0 4 5]
	* 
	* L =
	*[4 0 0 0
	* 2 6 0 0
	* 0 3 7 0
	* 0 0 4 5]
	* 
	* U = 
	*[4 3 0 0
	* 0 6 5 0
	* 0 0 7 2
	* 0 0 0 5]
	* 
	*/
	// ---- 定义矩阵 ---- //
	int N = 4;
	int NNZ = 10;
	int AI[] = { 0,2,5,8,10 };
	int AJ[] = { 0,1,0,1,2,1,2,3,2,3 };
	double AA[] = { 4,3,2,6,5,3,7,2,4,5 };
	double Rhs[] = { 4,3,5,1 };
	//
	int LNNZ = 7;
	int LI[] = { 0,1,3,5,7 };
	int LJ[] =    { 0,0,1,1,2,2,3 };
	double LA[] = { 4,2,6,3,7,4,5 };
	//
	int UNNZ = 7;
	int UI[] = { 0,2,4,6,7 };
	int UJ[] =    { 0,1,1,2,2,3,3 };
	double UA[] = { 4,3,6,5,7,2,5 };
	//
	double x[] = { 3,6,2,1 };
	double y[] = { 3,2,1,5 };
	double z[] = { 0,0,0,0 };
	double a = 4, b = 3;
	// ---- 稀疏矩阵 ---- //
	CSRMatrix A;
	A.initial('S', N, NNZ, AI, AJ, AA);
	CSRMatrix L;
	L.initial('L', N, LNNZ, LI, LJ, LA);
	CSRMatrix U;
	U.initial('U', N, UNNZ, UI, UJ, UA);
	printf("A:\n");
	printMatrix(&A);
	printf("L:\n");
	printMatrix(&L);
	printf("U:\n");
	printMatrix(&U);
	printf("RHS:\n");
	printvecRow(N, Rhs);

	// ---- 测试上下三角矩阵 ---- //
	solveLCSR(&L, Rhs, z);
	printf("x = L/Rhs:\n");
	printvecRow(N, z);

	solveUCSR(&U, Rhs, z);
	printf("x = U/Rhs:\n");
	printvecRow(N, z);

	// ---- 测试bicgstab ---- //
	// 读取矩阵
	
	CSRMatrix A1;
	CSRMatrix L1;
	CSRMatrix U1;
	vector<double> b1;
	//string filename1 = "test/CSR_Data2.dat";
	//string filename_A_b = "test/Ex1_CSR_A_b.dat";
	string filename_A = "test/Ex1_CSR_A.dat";
	string filename_b = "test/Ex1_CSR_b.dat";
	string filename_L = "test/Ex1_CSR_L.dat";
	string filename_U = "test/Ex1_CSR_U.dat";

	//readCSRMatrixVector(filename_A_b, &A1, b1);
	readCSRMatrix(filename_A, &A1);
	readCSRMatrix(filename_L, &L1);
	readCSRMatrix(filename_U, &U1);
	readVector(filename_b, b1);

	double start = 0.0, stop = 0.0;

	IterInfo info;
	double tol = 1.0e-8;
	int maxit = 100;
	//bicgstab(&info, z, &A, Rhs, &tol, &maxit);
	//printf("A:\n");	printMatrix(&A);
	//printf("RHS:\n"); printvecRow(N, Rhs);
	//printf("Solution:\n"); printvecRow(N, z);

	//
	double* sol = NULL;
	sol = (double*)malloc(A1.N * sizeof(double));
	start = clock();
	bicgstab(&info, sol, &A1, b1.data(), &tol, &maxit);
	stop = clock();
	printf("BicgStab Time: %f ms\n", stop - start);
	printf("Flag: %d; Iter: %.2f; RelRes: %.4e; Solution:\n",info.flag, info.iter, info.relres);
	//printvecCol(A1.N, sol);

	//
	printf("\n\n预处理bicgstab:\n");
	start = clock();
	bicgstab(&info, sol, &A1, b1.data(), &tol, &maxit, &L1, &U1);
	stop = clock();
	printf("BicgStab Time: %f ms\n", stop - start);
	printf("Flag: %d; Iter: %.2f; RelRes: %.4e; Solution:\n", info.flag, info.iter, info.relres);
	//printvecCol(A1.N, sol);

	printvecCol(2 * info.iter + 1, info.resvec);

	A.clear();
	L.clear();
	U.clear();
	free(sol);
	system("pause");
	return 0;
}



// 读取压缩矩阵与右边项
// 文本文件格式
// N: 矩阵维度
// NNZ : 非零个数
// I : 数据行, 共N + 1行
// J : 数据行, 共NNZ行
// A : 数据行, 共NNZ行
// b : 数据行, 共N行
void readCSRMatrixVector(string filename, CSRMatrix* A, vector<double>& b)
{
	int N = 0, NNZ = 0;
	int nTemp = 0;
	double dTemp = 0.0;
	ifstream infile(filename);
	string line;
	if (infile) // 有该文件
	{
		//while (getline(in, line)) // line中不包括每行的换行符
		//{
		//	cout << line << endl;
		//}
		// N
		getline(infile, line);
		N = atoi(line.c_str());
		// NNZ
		getline(infile, line);
		NNZ = atoi(line.c_str());
		A->initial('S', N, NNZ);
		b.resize(N);
		// I
		for (int i = 0; i < N + 1; i++)
		{
			getline(infile, line);
			A->I[i] = atoi(line.c_str());
		}
		// J
		for (int i = 0; i < NNZ; i++)
		{
			getline(infile, line);
			A->J[i] = atoi(line.c_str());
		}
		// A
		for (int i = 0; i < NNZ; i++)
		{
			getline(infile, line);
			A->A[i] = atof(line.c_str());
		}
		// b
		for (int i = 0; i < N; i++)
		{
			getline(infile, line);
			b[i] = atof(line.c_str());
		}
		
	}
	else // 没有该文件
	{
		cout << "no such file" << endl;//温馨小提示，木有此文件
	}
}

void readCSRMatrix(string filename, CSRMatrix* A)// 读取压缩矩阵
{
	// 文本文件格式
	// N: 矩阵维度
	// NNZ : 非零个数
	// I : 数据行, 共N + 1行
	// J : 数据行, 共NNZ行
	// A : 数据行, 共NNZ行
	int N = 0, NNZ = 0;
	int nTemp = 0;
	double dTemp = 0.0;
	ifstream infile(filename);
	string line;
	if (infile) // 有该文件
	{
		// N
		getline(infile, line);
		N = atoi(line.c_str());
		// NNZ
		getline(infile, line);
		NNZ = atoi(line.c_str());
		A->initial('S', N, NNZ);
		// I
		for (int i = 0; i < N + 1; i++)
		{
			getline(infile, line);
			A->I[i] = atoi(line.c_str());
		}
		// J
		for (int i = 0; i < NNZ; i++)
		{
			getline(infile, line);
			A->J[i] = atoi(line.c_str());
		}
		// A
		for (int i = 0; i < NNZ; i++)
		{
			getline(infile, line);
			A->A[i] = atof(line.c_str());
		}
	}
	else // 没有该文件
	{
		cout << "no such file" << endl;//温馨小提示，木有此文件
	}
}
void readVector(string filename, vector<double>& b)// 读取右边项
{
	// 文本文件格式
	// N: 数组维度
	// b: 数据行, 共N行
	int N = 0;
	ifstream infile(filename);
	string line;
	if (infile) // 有该文件
	{
		// N
		getline(infile, line);
		N = atoi(line.c_str());
		b.resize(N);	
		// b
		for (int i = 0; i < N; i++)
		{
			getline(infile, line);
			b[i] = atof(line.c_str());
		}
	}
	else // 没有该文件
	{
		cout << "no such file" << endl;//温馨小提示，木有此文件
	}
}
