#include<iostream>
#include<cmath>
#include<ctime>
#include<map>
#include<iomanip> 
#include<omp.h>

#ifndef MATRIC_2D
#define MATRIC_2D

#define ZERO 10e-8
#define END this->matric.end()
#define ROW_END this->matric[x].end()

using namespace std;

class matric2D{
public:
	//map<int,map<int,double>> matric;
	double ** matric;
	int x_max;
	int y_max;
public:
	matric2D(const int x,const int y);
	~matric2D();
	void set(const int x,const int y,const double value);
	double get(const int x,const int y);
	void random_init(bool normalization);
	void normalization();
	void print();
	matric2D * get_copy();
private:
	bool validate(const double value);
};

#endif