#include"algebra.h"

matric1D * algebra::get_cloumn_vector(matric2D * pA,const int column){
	if(pA->y_max<=column){
		cout<<"行数应在0~"<<pA->y_max-1<<"之间"<<endl;
		return NULL;
	}
	matric1D * pm = new matric1D(pA->x_max);

	//int core = omp_get_num_procs();
	//#pragma omp parallel for num_threads(core)
	for(int i=0;i<pA->x_max;++i){
		pm->set(i,pA->get(i,column));
	}

	return pm;
}

matric2D * algebra::add(matric2D * pA,matric2D * pB){
	if(pA->x_max!=pB->x_max||pA->y_max!=pB->y_max){
		cout<<"矩阵行列数不一!无法相加！\tA:"<<pA->x_max<<"x"<<pA->y_max<<"\tB:"<<pB->x_max<<"x"<<pB->y_max<<endl;
		return NULL;
	}

	matric2D * pm = new matric2D(pA->x_max,pA->y_max);

	int core = omp_get_num_procs();
	#pragma omp parallel for num_threads(core)
	for(int i=0;i<pm->x_max;++i){
		for(int j=0;j<pm->y_max;++j){
			pm->set(i,j,pA->get(i,j)+pB->get(i,j));
		}
	}

	return pm;
}

matric2D * algebra::add_e(matric2D * pA,double value){
	matric2D * pm = new matric2D(pA->x_max,pA->y_max);

	int core = omp_get_num_procs();
	#pragma omp parallel for num_threads(core)
	for(int i=0;i<pm->x_max;++i){
		for(int j=0;j<pm->y_max;++j){
			if(i==j)
				pm->set(i,j,pA->get(i,j)+value);
			else
				pm->set(i,j,pA->get(i,j));
		}
	}

	return pm;
}

matric2D * algebra::minus(matric2D * pA,matric2D * pB){
	if(pA->x_max!=pB->x_max||pA->y_max!=pB->y_max){
		cout<<"矩阵行列数不一!无法相减！\tA:"<<pA->x_max<<"x"<<pA->y_max<<"\tB:"<<pB->x_max<<"x"<<pB->y_max<<endl;
		return NULL;
	}

	matric2D * pm = new matric2D(pA->x_max,pA->y_max);

	int core = omp_get_num_procs();
	#pragma omp parallel for num_threads(core)
	for(int i=0;i<pm->x_max;++i){
		for(int j=0;j<pm->y_max;++j){
			pm->set(i,j,pA->get(i,j)-pB->get(i,j));
		}
	}

	return pm;
}

matric2D * algebra::multi(matric2D * pA,matric2D * pB){
	if(pA->y_max!=pB->x_max){
		cout<<"矩阵行列数不一!无法相乘！\tA:"<<pA->x_max<<"x"<<pA->y_max<<"\tB:"<<pB->x_max<<"x"<<pB->y_max<<endl;
		return NULL;
	}
	matric2D * pBt = algebra::transpose(pB);
	matric2D * pm = new matric2D(pA->x_max,pB->y_max);
	
	int core = omp_get_num_procs();
	#pragma omp parallel for num_threads(core)
	for(int i=0;i<pm->x_max;++i){
		for(int j=0;j<pm->y_max;++j){
			double value = 0.0;
			for(int k=0;k<pA->y_max;k++){
				//value += pA->get(i,k)*pB->get(k,j);
				/*if (pA->get(i, k) == 0 || pB->get(j, k) == 0)
					continue;
				value += pA->get(i, k)*pBt->get(j, k);*/
				value += pA->matric[i][k] * pBt->matric[j][k];
			}
			pm->set(i,j,value);
		}
	}

	delete pBt;
	return pm;
}

matric2D * algebra::multi(matric2D * pA,const double times){
	matric2D * pm = new matric2D(pA->x_max,pA->y_max);	
	
	int core = omp_get_num_procs();
	#pragma omp parallel for num_threads(core)
	for(int i=0;i<pm->x_max;++i){
		for(int j=0;j<pm->y_max;++j){
			pm->set(i,j,pA->get(i,j)*times);
		}
	}
	return pm;
}

matric2D * algebra::transpose(matric2D * pA){
	matric2D * pm = new matric2D(pA->y_max,pA->x_max);	
	
	int core = omp_get_num_procs();
	#pragma omp parallel for num_threads(core)
	for(int i=0;i<pm->x_max;++i){
		for(int j=0;j<pm->y_max;++j){
			pm->set(i,j,pA->get(j,i));
		}
	}
	return pm;
}

matric2D * algebra::non_nagetive_map(matric2D * pA){
	matric2D * pm = new matric2D(pA->x_max,pA->y_max);	
	
	int core = omp_get_num_procs();
	#pragma omp parallel for num_threads(core)
	for(int i=0;i<pm->x_max;++i){
		for(int j=0;j<pm->y_max;++j){
			double value = pA->get(i,j);
			if(value>ZERO){
				pm->set(i,j,value);
			}else{
				pm->set(i,j,0);
			}
		}
	}
	return pm;
}

double algebra::norm_1(matric2D * pA){
	double value = 0.0;

	int core = omp_get_num_procs();
	#pragma omp parallel for num_threads(core) reduction(+:value)
	for(int i=0;i<pA->x_max;++i){
		double max = 0.0;
		for(int j=0;j<pA->y_max;++j){
			if(max<abs(pA->get(i,j))){
				max = abs(pA->get(i,j));
			}
		}
		value+=max;
	}
	return value;
}

double algebra::norm_f(matric2D * pA){
	double value = 0.0;

	int core = omp_get_num_procs();
	#pragma omp parallel for num_threads(core) reduction(+:value)
	for(int i=0;i<pA->x_max;++i){
		for(int j=0;j<pA->y_max;++j){
			value+=pow(pA->get(i,j),2);
		}
	}
	return sqrt(value);
}

matric1D * algebra::add(matric1D * pA,matric1D * pB){
	if(pA->x_max!=pB->x_max){
		cout<<"向量规模不一，无法相加!\tA:"<<pA->x_max<<"\tB:"<<pB->x_max<<endl;
		return NULL;
	}
	matric1D * pm = new matric1D(pA->x_max);

	//int core = omp_get_num_procs();
	//#pragma omp parallel for num_threads(core)
	for(int i=0;i<pm->x_max;++i){
		pm->set(i,pA->get(i)+pB->get(i));
	}
	return pm;
}

matric1D * algebra::add_e(matric1D * pA,const double value){
	matric1D * pm = new matric1D(pA->x_max);

	//int core = omp_get_num_procs();
	//#pragma omp parallel for num_threads(core)
	for(int i=0;i<pm->x_max;++i){
		pm->set(i,pA->get(i)+value);
	}
	return pm;
}

matric1D * algebra::minus(matric1D * pA,matric1D * pB){
	if(pA->x_max!=pB->x_max){
		cout<<"向量规模不一，无法相加!\tA:"<<pA->x_max<<"\tB:"<<pB->x_max<<endl;
		return NULL;
	}
	matric1D * pm = new matric1D(pA->x_max);

	//int core = omp_get_num_procs();
	//#pragma omp parallel for num_threads(core)
	for(int i=0;i<pm->x_max;++i){
		pm->set(i,pA->get(i)-pB->get(i));
	}
	return pm;
}

double algebra::multi(matric1D * pA,matric1D * pB){
	if(pA->x_max!=pB->x_max){
		cout<<"向量规模不一，无法相乘!\tA:"<<pA->x_max<<"\tB:"<<pB->x_max<<endl;
		return NULL;
	}
	double value = 0.0;

	//int core = omp_get_num_procs();
	//#pragma omp parallel for num_threads(core) reduction(+:value)
	for(int i=0;i<pA->x_max;++i){
		value+=pA->get(i)*pB->get(i);
	}

	return value;
}

matric1D * algebra::non_nagetive_map(matric1D * pA){
	matric1D * pm = new matric1D(pA->x_max);

	//int core = omp_get_num_procs();
	//#pragma omp parallel for num_threads(core)
	for(int i=0;i<pm->x_max;++i){
		double value = pA->get(i);
		if(value>ZERO){
			pm->set(i,value);
		}else{
			pm->set(i,0.0);
		}
	}
	return pm;
}

double algebra::norm_1(matric1D * pA){
	double max = 0.0;
	for(int i=0;i<pA->x_max;++i){
		if(max<abs(pA->get(i))){
			max = abs(pA->get(i));
		}
	}
	return max;
}

double algebra::norm_2(matric1D * pA){
	double value = 0.0;

	//int core = omp_get_num_procs();
	//#pragma omp parallel for num_threads(core) reduction(+:value)
	for(int i=0;i<pA->x_max;i++){
		value+=pow(pA->get(i),2);
	}
	return sqrt(value);
}