/*
* 数组与矩阵内存管理
*/
#ifndef MEMORY_H
#include <stdlib.h>
#include <stdio.h>

// double数组分配内存
void dvecmal(const int* N, double** x);
// float数组分配内存
void fvecmal(const int* N, float** x);
// int数组分配内存
void ivecmal(const int* N, int** x);
// double/float/int数组内存释放
void dvecfree(double** x);
void fvecfree(float** x);
void ivecfree(int** x);


#endif // !MEMORY_H

