/*
* ����������ڴ����
*/
#ifndef MEMORY_H
#include <stdlib.h>
#include <stdio.h>

// double��������ڴ�
void dvecmal(const int* N, double** x);
// float��������ڴ�
void fvecmal(const int* N, float** x);
// int��������ڴ�
void ivecmal(const int* N, int** x);
// double/float/int�����ڴ��ͷ�
void dvecfree(double** x);
void fvecfree(float** x);
void ivecfree(int** x);


#endif // !MEMORY_H

