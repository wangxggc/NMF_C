#include"uvde.h"

uv_decompose::uv_decompose(int m,int n,int k){
	this->m = m;
	this->n = n;
	this->k = k;

	//this->pD = new matric2D(this->m,this->n);
	this->pU = new matric2D(this->m,this->k);
	this->pV = new matric2D(this->k,this->n);
	this->pD = new matric2D(this->m,this->n);
	//this->pU->random_init(true);
	//this->pV->random_init(true);
	//this->pD = algebra::multi(this->pU,this->pV);

	this->u_eps = 1e-8;
	this->v_eps = 1e-8;
	
	this->u_loop_max = 1000;
	this->v_loop_max = 1000;
	this->d_loop_max = 500;

	this->c = 0.01;
	this->log_file_name = "result";
}
uv_decompose::uv_decompose(struct uint * data, int length, int m, int n, int k){
	this->m = m;
	this->n = n;
	this->k = k;

	//this->pD = new matric2D(this->m,this->n);
	this->pU = new matric2D(this->m, this->k);
	this->pV = new matric2D(this->k, this->n);
	this->pD = new matric2D(this->m, this->n);

	for (int i = 0; i < this->pD->x_max; i++){
		for (int j = 0; j < this->pD->y_max; j++){
			this->pD->set(i, j, 0);
		}
	}
	for (int i = 0; i < length; i++){
		this->pD->set(data[i].x, data[i].y, data[i].value);
	}

	this->pD->normalization();
	this->u_eps = 1e-8;
	this->v_eps = 1e-8;

	this->u_loop_max = 1000;
	this->v_loop_max = 1000;
	this->d_loop_max = 500;

	this->c = 0.01;
	this->log_file_name = "result";
}
uv_decompose::uv_decompose(double ** data,int m,int n,int k){
	this->m = m;
	this->n = n;
	this->k = k;

	//this->pD = new matric2D(this->m,this->n);
	this->pU = new matric2D(this->m,this->k);
	this->pV = new matric2D(this->k,this->n);
	this->pD = new matric2D(this->m,this->n);

	for(int i=0;i<this->m;++i){
		for(int j=0;j<this->n;++j){
			this->pD->set(i,j,data[i][j]);
		}
	}
	this->pD->normalization();
	this->u_eps = 1e-8;
	this->v_eps = 1e-8;
	
	this->u_loop_max = 1000;
	this->v_loop_max = 1000;
	this->d_loop_max = 500;

	this->c = 0.01;
	this->log_file_name = "result";
}
void uv_decompose::update_u(){
	cout<<"update_u();";
	matric2D * pVt = algebra::transpose(this->pV);
	matric2D * pVVt = algebra::multi(this->pV,pVt);
	matric2D * pDVt = algebra::multi(this->pD,pVt);
	delete pVt;

	matric2D * pUU = NULL;
	double t = 1.0;
	do{
		if (pUU != NULL){
			delete pUU;
		}
		double tt = this->c / sqrt(t);

		matric2D * pUVVt = algebra::multi(this->pU, pVVt);
		matric2D * pUVVt_DVt = algebra::minus(pUVVt, pDVt);
		matric2D * pR = algebra::multi(pUVVt_DVt, 2.0*tt);
		delete pUVVt;
		delete pUVVt_DVt;

		matric2D * pgradientU = NULL;

		pgradientU = pR->get_copy();
		
		delete pR;

		matric2D * pU1 = algebra::minus(this->pU,pgradientU);
		delete pgradientU;

		matric2D * pU_new = algebra::non_nagetive_map(pU1);
		delete pU1;

		pUU = algebra::minus(pU_new,pU);
		delete pU;
		
		pU = pU_new;
		t++;

	}while(algebra::norm_f(pUU)>this->u_eps&&t<=this->u_loop_max);
	if (this->log_file_name.size()>0){
		fstream out(log_file_name, fstream::app);
		out << "update_u():" << t-1;
		out.close();
	}
	cout<<t-1<<"\t";

	if(pUU!=NULL){
		delete pUU;
	}
	delete pVVt;
	delete pDVt;
}
void uv_decompose::update_v(){
	cout<<"update_v();";
	matric2D * pUt = algebra::transpose(this->pU);
	matric2D * pS = algebra::multi(pUt,pU);
	matric2D * pR = algebra::multi(pUt,pD);
	delete pUt;

	int t = 1;
	//int core = omp_get_num_procs();
	//#pragma omp parallel for num_threads(core)
	double value = 0.0;
	for(int i=0;i<this->pV->y_max;++i){
		value+=update_v_column(i,pS,pR);
	}
	cout<<setiosflags(ios::fixed)<<setprecision(5)<<1+value/this->pV->y_max;

	if (this->log_file_name.size()>0){
		fstream out(log_file_name,fstream::app);
		out<<"\tupdate_v():"<<setiosflags(ios::fixed)<<setprecision(5)<<(1+value/this->pV->y_max);
		out.close();
	}
	

	delete pS;
	delete pR;
}

int uv_decompose::update_v_column(int column,matric2D * pS,matric2D * pR){
	matric1D * ppv = algebra::get_cloumn_vector(this->pV,column);
	matric1D * pv = ppv->get_copy();

	int t=0;
	for(t=0;t<this->v_loop_max;t++){
		//int core = omp_get_num_procs();
		//#pragma omp parallel for num_threads(core)
		
		for(int i=0;i<pv->x_max;++i){
			double value = pR->get(i,column);
			for(int l=0;l<pv->x_max;++l){
				if(l!=i)
					value-=pS->get(i,l)*pv->get(l);
			}
			value/=pS->get(i,i);
			if(value>ZERO)
				pv->set(i,value);
			else
				pv->set(i,0);
		}
		matric1D * pdeltv = algebra::minus(ppv,pv);
		delete ppv;
		ppv = pv->get_copy();

		if(algebra::norm_1(pdeltv)<this->v_eps){
			/*for(int i=0;i<pv->x_max;++i){
				this->pV->set(i,column,pv->get(i));
			}*/
			delete pdeltv;
			break;
		}
		delete pdeltv;
	}
	//迭代次数到了也返回
	for(int i=0;i<pv->x_max;++i){
		this->pV->set(i,column,pv->get(i));
	}
	delete ppv;
	delete pv;
	return t;
}

void uv_decompose::decompose(){
	this->pU->random_init(true);
	this->pV->random_init(true);

	for(int i=0;i<this->d_loop_max;++i){
		clock_t begin = clock();
		this->update_u();
		this->update_v();
		clock_t end = clock();
		
		if (this->log_file_name.size()>0){
			fstream out(log_file_name, fstream::app);
			out << "\titeration:" << i + 1 <<"\ttime:"<< end-begin <<"ms"<<endl;
			out.close();
		}
		cout << "\titeration:" << i + 1 << "\ttime:" << end - begin << "ms" << endl;
		/*if (i % 100 == 99){
			string s = "iteration-" + (i + 1);
			save(this->log_file_name + "-D-"+s, this->pD);
			save(this->log_file_name + "-U-"+s, this->pU);
			save(this->log_file_name + "-V-"+s, this->pV);
		}*/
	}
}

void uv_decompose::save(string filename, matric2D * pm){
	ofstream out(filename);
	out << "-------------------------------x_max:" << pm->x_max;
	out << "\ty_max:" << pm->y_max << "-------------------------------\n";
	for (int i = 0; i<pm->x_max; ++i){
		for (int j = 0; j<pm->y_max; ++j){
			out << setiosflags(ios::fixed) << setprecision(5) << pm->get(i, j) << "  ";
		}
		out << endl;
	}
	out.close();
}
matric2D * uv_decompose::uutu_u(double times){
	matric2D * pUt = algebra::transpose(this->pU);
	matric2D * pUtU = algebra::multi(pUt,this->pU);
	matric2D * pUtU_E = algebra::add_e(pUtU,-1.0);
	matric2D * pUUtU_U = algebra::multi(this->pU,pUtU_E);
	matric2D * pUUtU_U_4b = algebra::multi(pUUtU_U,4*b*times);

	delete pUt;
	delete pUtU;
	delete pUtU_E;
	delete pUUtU_U;

	return pUUtU_U_4b;
}
matric2D * uv_decompose::uvvt_dvt(matric2D * pVVt,matric2D * pDVt,double times){
	matric2D * pUVVt = algebra::multi(this->pU,pVVt);
	matric2D * pUVVt_DVt = algebra::minus(pUVVt,pDVt);
	matric2D * pUVVt_DVt_2 = algebra::multi(pUVVt_DVt,2.0*times);

	delete pUVVt;
	delete pUVVt_DVt;

	return pUVVt_DVt_2;
}