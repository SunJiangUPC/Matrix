#ifndef INFORMATION_H
#define INFORMATION_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

using namespace std;

class CSRMatrix;

void* message(char* info);
void error(const char* str);
void warning(const char* str);


// ���д�ӡvector
template<typename T>
void printvecCol(const int N, T* x)
{
	//cout << "Print Vector:" << endl;
	for (int i = 0; i < N; i++)
	{
		std::cout << x[i] << std::endl;
	}
}

// ���д�ӡvector
template<typename T>
void printvecRow(const int N, T* x)
{
	//cout << "Print Vector:" << endl;
	for (int i = 0; i < N - 1; i++)
	{
		std::cout << x[i] << ",";
	}
	std::cout << x[N - 1] << std::endl;
}

// ��ӡ����
void printMatrix(CSRMatrix* A);



#endif