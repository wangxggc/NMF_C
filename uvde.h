#include"algebra.h"
#include<fstream>
#ifndef UVDE_H
#define UVDE_H

using namespace std;

struct uint{
	int x, y;
	double value;
};

class uv_decompose{
public:
	matric2D * pD;
	matric2D * pU;
	matric2D * pV;

	double u_eps;
	double v_eps;

	int u_loop_max;
	int v_loop_max;
	int d_loop_max;

	int m;
	int n;
	int k;

	double a;
	double b;
	double c;

	string log_file_name;
public:
	uv_decompose(int m,int n,int k);
	uv_decompose(double ** data,int m,int n,int k);
	uv_decompose(struct uint * data, int length, int m, int n, int k);
	void update_u();
	void update_v();
	void decompose();
	matric2D * uutu_u(double times);
	matric2D * uvvt_dvt(matric2D * pVVt,matric2D * pDVt,double times);
	int update_v_column(int column,matric2D * pS,matric2D * pR);
	void save(string filename, matric2D * matric);
};


#endif