#include<iostream>
#include<cmath>
#include<ctime>
#include<map>
#include<iomanip> 
#include<omp.h>

#ifndef MATRIC_1D
#define MATRIC_1D

#define ZERO 10e-8
#define END this->matric.end()

using namespace std;

class matric1D{
public:
	//map<int,double> matric;
	double * matric;
	int x_max;
public:
	matric1D(const int x);
	~matric1D();
	void set(const int x,const double value);
	double get(const int x);
	void random_init(const bool normalization);
	void normalization();
	void print();
	matric1D * get_copy();
private:
	bool validate(const double value);
};

#endif