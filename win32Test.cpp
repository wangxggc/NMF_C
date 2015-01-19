#include"uvde.h"
#include<ctime>
#include<cstdio>
#include<fstream>
#include<sstream>

using namespace std;

void test_update_u(){
	for(int i=0;i<10;i++){
		system("cls");
		cout<<"test:"<<i<<endl;
		uv_decompose * de = new uv_decompose(50,40000,5);
		de->pU->random_init(true);
		de->pV->random_init(true);
		de->pD = algebra::multi(de->pU,de->pV);
		de->c = 0.0005;
		de->pU->print();
		getchar();
		de->pU->random_init(true);
		de->update_u();
		de->pU->print();
		getchar();
	}
}

void test_update_v(){
	for(int i=0;i<10;i++){
		system("cls");
		cout<<"test:"<<i<<endl;
		uv_decompose * de = new uv_decompose(40000,100,5);
		de->pU->random_init(true);
		de->pV->random_init(true);
		//de->pV->print();
		de->pD = algebra::multi(de->pU,de->pV);
		de->b = 1;
		de->c = 0.0001;
		de->pV->print();
		getchar();
		de->pV->random_init(true);
		//de->pV->print();
		de->update_v();
		de->pV->print();
		getchar();
	}
}

void test(){
	system("cls");
	uv_decompose * de = new uv_decompose(20,8,3);
	de->pU->random_init(true);
	de->pV->random_init(true);
	de->pD = algebra::multi(de->pU,de->pV);
	getchar();
	de->pU->random_init(true);
	de->pV->random_init(true);

	de->a = 0;
	de->b = 1;
	de->c = 0.1;
	de->u_loop_max = 1000;
	de->v_loop_max = 1000;
	de->d_loop_max = 300;
	de->u_eps = 1e-6;
	de->v_eps = 1e-6;

	de->decompose();

	de->pD->print();
	algebra::multi(de->pU, de->pV)->print();
	cout << "UtU==S" << endl;
	algebra::multi(algebra::transpose(de->pU), de->pU)->print();
	cout << "UtD==R" << endl; 
	algebra::multi(algebra::transpose(de->pU), de->pD)->print();
	cout << "V" << endl;
	de->pV->print();
		
}

void save(string filename,matric2D * pm){
	ofstream out(filename);
	out<<"-------------------------------x_max:"<<pm->x_max;
	out<<"\ty_max:"<<pm->y_max<<"-------------------------------\n";
	for(int i=0;i<pm->x_max;++i){
		for(int j=0;j<pm->y_max;++j){
			out<<setiosflags(ios::fixed)<<setprecision(8)<<pm->get(i,j)<<"  ";
		}
		out<<endl;
	}
	out.close();
}

void main(){
	/*读入文件*/
	//test();
	//getchar();
	ifstream input;
	/*
	saved in sparse like
	0 0 0.35 
	1 1 0.23
	1 2 1.02
	where the 1th number means x, 2ed number means y, 3rd number means value 
	*/
	input.open("tfidf");
	//input.open("ORLMatrix.txt");
	//input.open("senti/p-t20-book-V");
	if (!input){
		cout << "文件打开失败" << endl;
		return;
	}

	int row_num = 0;
	int column_num = 0;
	int length = 0;
	/*double ** filedata = new double *[6000];
	for (int i = 0; i < 6000; i++){
		filedata[i] = new double[20000];
		memset(filedata[i], 0, 20000 * sizeof(double));
	}*/
	struct uint * data = new uint[500000];
	cout << "init data" << endl;
	while (!input.eof()){
		string str;
		getline(input, str);
		istringstream strin(str);
		double value = 0.0;
		int row = 0;
		int column = 0;
		
		while (!strin.eof()){
			strin >> row;
			strin >> column;
			double value;
			strin >> value;
			data[length].x = row;
			data[length].y = column;
			data[length].value = value;
			length++;
		}
		if (column_num < column){
			column_num = column;
		}
		if (row_num < row){
			row_num = row;
		}
		str.clear();
	}

	input.close();

	uv_decompose * de = new uv_decompose(data, length-1,row_num+1, column_num+1, 160);

	cout << de->pD->get(0, 0) << endl;;

	//for (int i = 0; i < 10500; ++i){
	//	delete[] filedata[i];
	//}
	//delete[] filedata;
	delete[] data;

	char c_filename[1000];
	sprintf(c_filename, "NMF");

	string file_name(c_filename);
	//string file_name = "senti/p-t20-book-V";
	/*参数设置*/
	de->c = 0.0003;
	de->u_eps = 5e-6;
	de->v_eps = 5e-6;
	de->u_loop_max = 1000;
	de->v_loop_max = 1000;
	de->d_loop_max = 500;
	de->log_file_name = file_name;

	/*分解*/
	cout << "decompose() " << de->m << " " << de->k << " " << de->n << endl;
	de->decompose();

	/*存储*/
	//save(file_name + "-D", de->pD);
	save(file_name + "U", de->pU);
	save(file_name + "V", de->pV);
	cout << file_name << endl;
}