/*
* �����������: ��ʼ��
*/
#ifndef VECTOR_H
#define VECTOR_H

#include <stdlib.h>
#include <stdio.h>
#include <algorithm>

// double����Ԫ��ȫΪ0
void dzeros(const int* N, double* x);
// float����Ԫ��ȫΪ0
void fzeros(const int* N, float* x);
// int����Ԫ��ȫΪ0
void izeros(const int* N, int* x);
// double����Ԫ��ȫΪ1
void dones(const int* N, double* x);
// float����Ԫ��ȫΪ1
void fones(const int* N, float* x);
// int����Ԫ��ȫΪ1
void iones(const int* N, int* x);

// �����������鳤��
void resize(const int* N, double* x, const int* size);

// ������x�вü�������y,����x���Ȳ���ʱ��0
void cut(const int* N, const double* x, const int* size, double* y);

// ����Ԫ�ض�������ֵ
bool isfiniteall(const int* N, const double* x);


#endif // !VECTOR_H
