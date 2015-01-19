#include<iostream>
#include"algebra.h"

#ifndef MULTI_MATRIC_1
#define MULTI_MATRIC_1

using namespace std;

int cal_app_number(int number){
	int mid = sqrt(number);
	if (mid*mid < number){
		++mid;
	}
	for (int i = mid; i <= number; i++){
		if (number%i == 0){
			return i;
		}
	}
	return -1;
}

matric2D * multi(matric2D * pA, matric2D * pB){
	if (pA->x_max != pB->x_max || pA->y_max != pB->y_max){
		cout << "矩阵行列数不一!无法相减！\tA:" << pA->x_max << "x" << pA->y_max << "\tB:" << pB->x_max << "x" << pB->y_max << endl;
		return NULL;
	}

	int m = cal_app_number(pA->x_max);
	int k = cal_app_number(pA->y_max);
	int n = cal_app_number(pB->y_max);

	matric2D * * ppA;
	ppA = new matric2D *[m];
	for (int i = 0; i < m; ++i){
		ppA[i] = new matric2D(pA->x_max/m,pA->y_max/k);
	}
}
#endif