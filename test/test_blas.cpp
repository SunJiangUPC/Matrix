/*
* x, y, z: 向量
* a, b, c : 常数
* A, B, C : 矩阵
*/
#include <stdio.h>
#include "../blas/blas.h"
#include "../blas/CSRmatrix.h"
#include "../blas/memory.h"
#include "../blas/vector.h"
#include "../config/information.h"
#include <iostream>
#include <Windows.h>

using namespace std;



int main()
{
	/*
	* A = 
	*[4 3 0 0
	* 2 6 5 0
	* 0 3 7 2
	* 0 0 4 5]
	*/
	int N = 4;
	int NNZ = 10;
	int AI[] = { 0,2,5,8,10 };
	int AJ[] = { 0,1,0,1,2,1,2,3,2,3 };
	double AA[] = { 4,3,2,6,5,3,7,2,4,5 };
	double Rhs[] = { 4,3,5,1 };
	double x[] = { 3,6,2,1 };
	double y[] = { 3,2,1,5 };
	double z[] = { 0,0,0,0 };
	double a = 4, b = 3;

	// -------- 测试memory.h -------- //
	int M = 10;
	double* dx = NULL;
	float* fx = NULL;
	int* ix = NULL;
	dvecmal(&M, &dx);// double数组分配内存
	fvecmal(&M, &fx);// float数组分配内存
	ivecmal(&M, &ix);// int数组分配内存
	cout << "分配内存:" << endl;
	cout << "dx = "; printvecCol(M, dx);
	cout << "fx = "; printvecCol(M, fx);
	cout << "ix = "; printvecCol(M, ix);

	/*
	double* d1 = NULL;
	double* d2 = NULL;
	double* d3 = NULL;
	double* d4 = NULL;
	double* d5 = NULL;
	d1 = (double*)malloc(M * sizeof(double));
	d2 = (double*)malloc(M * sizeof(double));
	d3 = (double*)malloc(M * sizeof(double));
	d4 = (double*)malloc(M * sizeof(double));
	free(d1);
	free(d2);
	free(d3);
	free(d4);
	*/

	/*
	double* p = NULL;
	for (int i = 0; i < 100; i++)
	{
		p = (double*)malloc(100);
		free(p);
		cout << "i=" << i << ":" << p << endl;

		//Sleep(100);
	}
	*/

	// -------- 测试vector.h -------- //
	
	dzeros(&M, dx);
	fzeros(&M, fx);
	izeros(&M, ix);
	cout << "初始化0:" << endl;
	cout << "dx = "; printvecRow(M, dx);
	cout << "fx = "; printvecRow(M, fx);
	cout << "ix = "; printvecRow(M, ix);
	dones(&M, dx);
	fones(&M, fx);
	iones(&M, ix);
	cout << "初始化1:" << endl;
	cout << "dx = "; printvecRow(M, dx);
	cout << "fx = "; printvecRow(M, fx);
	cout << "ix = "; printvecRow(M, ix);
	resize(&M, dx, &N);
	cout << "重定义数组长度:" << endl;
	cout << "dx = "; printvecRow(N, dx);
	double* dy2 = NULL;
	dvecmal(&M, &dy2);
	cut(&N, dx, &M, dy2);
	cout << "裁减数组:" << endl;
	cout << "dx = "; printvecRow(N, dx);
	cout << "y2 = "; printvecRow(M, dy2);


	dvecfree(&dx);// double/float/int数组内存释放
	fvecfree(&fx);
	ivecfree(&ix);
	dvecfree(&dy2);

	// -------- 测试matrixCSR.h -------- //
	CSRMatrix A;
	//A.Type = 'S';
	//A.N = N;
	//A.NNZ = NNZ;
	//A.I = AI;
	//A.J = AJ;
	//A.A = AA;
	A.initial('S', N, NNZ, AI, AJ, AA);
	//A.clear();

	// -------- 测试blas.h -------- //
	// blas level 1
	printf("x: ");
	printvecRow(N, x);
	dcopy(&N, x, y);// y = x
	printf("y=x: ");
	printvecRow(N, y);
	dscal(&N, &a, x, y);// y = a * x
	printf("y=a*x: ");
	printvecRow(N, y);
	double dotval = ddot(&N, x, y);// a = x dot y
	printf("x.*y: %f\n", dotval);
	double normval = dnorm(&N, x);// a = sum( x(i) * x(i))
	printf("sum(x(i)*x(i)): %f\n", normval);
	daxpby(&N, &a, x, &b, y);// y := a*x + b*y
	printf("y := a*x + b*y: ");
	printvecRow(N, y);

	// blas level 2
	// y := A*x
	dgemvCSR(&A, x, z);
	printf("y := A*x: ");
	printvecRow(N, z);
	// y := A*x + b*y
	dAxPbyCSR(&A, x, &b, y);
	printf("y := A*x + b*y: ");
	printvecRow(N, y);


	A.clear();

	cout << "Please input any character to stop:";
	char cmd;
	//cmd = getchar();
	system("pause");
	return 0;
}