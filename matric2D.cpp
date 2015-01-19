#include"matric2D.h"

matric2D::matric2D(const int x,const int y){
	this->x_max = x;
	this->y_max = y;
	matric = new double*[x];
	for(int i=0;i<x_max;++i){
		matric[i] = new double[y];
	}
	//矩阵初始化
	for(int i=0;i<this->x_max;++i){
		for(int j=0;j<this->y_max;++j){
			this->matric[i][j] = 0.0;
		}
	}
}
matric2D::~matric2D(){
	//int core = omp_get_num_procs();
	//#pragma omp parallel for num_threads(core)
	//for(int i=0;i<this->x_max;++i){
	//	if(this->matric.find(i)!=END){
	//		this->matric[i].clear();
	//	}
	//}
	//this->matric.clear();
	for(int i=0;i<this->x_max;i++){
		delete[] matric[i];
	}
	delete[] matric;
}
void matric2D::set(const int x,const int y,const double value){
	if(x<0||x>=this->x_max
		||y<0||y>=this->y_max){
		cout<<"数组越界"<<endl;
		return ;
	}
	if(this->validate(value))
		this->matric[x][y] = value;
}
double matric2D::get(const int x,const int y){
	/*if(this->matric.find(x)!=END){
		if(this->matric[x].find(y)!=ROW_END){
			return this->matric[x][y];
		}
	}*/
	if(x<0||x>=this->x_max
		||y<0||y>=this->y_max){
		cout<<"数组越界"<<endl;
		return -1;
	}
	return this->matric[x][y];
}
void matric2D::random_init(bool normalization){
	srand(time(0));
	
	//int core = omp_get_num_procs();
	//#pragma omp parallel for num_threads(core)
	for(int i=0;i<this->x_max;++i){
		for(int j=0;j<this->y_max;++j){
			this->set(i,j,rand()%10+1);
			//this->matric[i][j] = rand()%10+1;
		}
	}
	if(normalization)
		this->normalization();
}
void matric2D::normalization(){
	int core = omp_get_num_procs();
	#pragma omp parallel for num_threads(core)
	for(int j=0;j<this->y_max;++j){
		double value = 0.0;
		for(int i=0;i<this->x_max;++i){
			value += pow(this->get(i,j),2);
		}
		value = sqrt(value);
		if(this->validate(value)){
			for(int i=0;i<this->x_max;++i){
				this->set(i,j,this->get(i,j)/value);
			}
		}
	}
}
//打印不能并行
void matric2D::print(){
	cout<<"x-dim:"<<this->x_max<<"\ty-dim:"<<this->y_max<<endl;
	for(int i=0;i<this->x_max;++i){
		for(int j=0;j<this->y_max;++j){
			cout<<setiosflags(ios::fixed)<<setprecision(10)<<this->get(i,j)<<" ";
		}
		cout<<endl;
	}
}
matric2D * matric2D::get_copy(){
	matric2D * pm = new matric2D(this->x_max,this->y_max);
	
	int core = omp_get_num_procs();
	#pragma omp parallel for num_threads(core)
	for(int i=0;i<pm->x_max;++i){
		for(int j=0;j<pm->y_max;++j){
			pm->set(i,j,this->get(i,j));
		}
	}
	return pm;
}

bool matric2D::validate(const double value){
	if(abs(value)>ZERO){
		return true;
	}else{
		return false;
	}
}