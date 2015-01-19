#include"matric1D.h"

matric1D::matric1D(const int x){
	this->x_max = x;
	this->matric = new double[x];
	for(int i=0;i<this->x_max;++i){
		this->matric[i] = 0.0;
	}
}
matric1D::~matric1D(){
	//this->matric.clear();
	delete matric;
}
void matric1D::set(const int x,const double value){
	if(this->validate(value))
		this->matric[x] = value;
}
double matric1D::get(const int x){
	//if(this->matric.find(x)!=END){
	//	return this->matric[x];
	//}
	return this->matric[x];
}
void matric1D::random_init(const bool normalization){
	srand(time(0));
	
	for(int i=0;i<this->x_max;i++){
		this->set(i,rand()%10+1);
	}
	if(normalization){
		this->normalization();
	}
}
void matric1D::normalization(){
	double value = 0.0;

	int core = omp_get_num_procs();
	#pragma omp parallel for num_threads(core) reduction(+:value)
	for(int i=0;i<this->x_max;i++){
		value+=pow(this->get(i),2);
	}
	value = sqrt(value);
	
	if(this->validate(value)){
		#pragma omp parallel for num_threads(core)
		for(int i=0;i<this->x_max;i++){
			this->set(i,this->get(i)/value);
		}
	}
}
//打印不能并行
void matric1D::print(){
	cout<<"x-dim:"<<this->x_max<<endl;
	for(int i=0;i<this->x_max;i++){
		cout<<setiosflags(ios::fixed)<<setprecision(10)<<this->get(i)<<" ";
	}
	cout<<endl;
}
matric1D * matric1D::get_copy(){
	matric1D * pm = new matric1D(this->x_max);
	
	int core = omp_get_num_procs();
	#pragma omp parallel for num_threads(core)
	for(int i=0;i<this->x_max;i++){
		pm->set(i,this->get(i));
	}
	return pm;
}
bool matric1D::validate(const double value){
	if(abs(value)>ZERO){
		return true;
	}else{
		return false;
	}
}