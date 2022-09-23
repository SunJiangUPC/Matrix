/*
* 数组简易运算: 初始化
*/
#ifndef VECTOR_H
#define VECTOR_H

#include <stdlib.h>
#include <stdio.h>
#include <algorithm>

// double数组元素全为0
void dzeros(const int* N, double* x);
// float数组元素全为0
void fzeros(const int* N, float* x);
// int数组元素全为0
void izeros(const int* N, int* x);
// double数组元素全为1
void dones(const int* N, double* x);
// float数组元素全为1
void fones(const int* N, float* x);
// int数组元素全为1
void iones(const int* N, int* x);

// 重新设置数组长度
void resize(const int* N, double* x, const int* size);

// 从数组x中裁减出数组y,数组x长度不足时补0
void cut(const int* N, const double* x, const int* size, double* y);

// 所有元素都是有限值
bool isfiniteall(const int* N, const double* x);


#endif // !VECTOR_H
