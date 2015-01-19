#include"matric1D.h"
#include"matric2D.h"
#include<omp.h>

#ifndef ALGEBRA_1
#define ALGEBRA_1
#define ZERO 10e-8

using namespace std;



class algebra{
public:
	static matric1D * get_cloumn_vector(matric2D * pA,const int column);
	static matric2D * add(matric2D * pA,matric2D * pB);
	static matric2D * add_e(matric2D * pA,const double value);
	static matric2D * minus(matric2D * pA,matric2D * pB);
	static matric2D * multi(matric2D * pA,matric2D * pB);
	static matric2D * multi(matric2D * pA,const double times);
	static matric2D * transpose(matric2D * pA);
	static matric2D * non_nagetive_map(matric2D * pA);
	static double norm_1(matric2D * pA);
	static double norm_f(matric2D * pA);
	//static double norm_max(matric2D * pA);
public:
	static matric1D * add(matric1D * pA,matric1D * pB);
	static matric1D * add_e(matric1D * pA,const double value);
	static matric1D * minus(matric1D * pA,matric1D * pB);
	static double multi(matric1D * pA,matric1D * pB);
	static matric1D * non_nagetive_map(matric1D * pA);
	static double norm_1(matric1D * pA);
	static double norm_2(matric1D * pA);
	//static double norm_max(matric1D * pA);
};

#endif