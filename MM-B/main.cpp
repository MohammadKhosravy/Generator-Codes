#include "sel.h"

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <sstream>

using namespace std;

int main (int argc, char* argv[])
{
	srand(atoi(argv[8]));

	int n = atoi(argv[1]);
	int p = atoi(argv[2]);
	int gamma = atoi(argv[3]);
	int param = atoi(argv[4]);
	int R = atoi(argv[5]);          //to choose what random method
	int t = atoi(argv[6]);
    float b = atof(argv[7]);



	Select sel;

	//double start = clock();
	if (param == 0)
		sel.generate_rand(n,p,gamma,R);
	else if (param == 1)
		sel.generate_hard_c_1(n,p,gamma,R,b,t);
    else if (param == 2)
		sel.generate_hard_d_1(n,p,gamma,R,b,t);
    else if (param == 3)
		sel.generate_hard_c_d_1(n,p,gamma,R,b,t);
	else if (param == 4)
		sel.generate_hard_c_2(n,p,gamma,R,b,t);
    else if (param == 5)
		sel.generate_hard_d_2(n,p,gamma,R,b,t);
    else if (param == 6)
		sel.generate_hard_c_d_2(n,p,gamma,R,b,t);


	vector<double> c = sel.get_c();
	vector<double> d = sel.get_d();

	stringstream ss;
	ss<<"instance";
	for (int i=0; i<8; ++i)
		ss<<"-"<<argv[i+1];
	ss<<".dat";

	ofstream out(ss.str().c_str());
	out<<n<<";"<<p<<";"<<gamma<<"\n";
	for (int i=0; i<n; ++i)
        out<<c[i]<<";";
    out<<"\n";
    for (int i=0; i<n; ++i)
        out<<d[i]<<";";
	out.close();


//	stringstream ss;
//	ss<<"modifiedc-"<<n<<"-"<<p<<"-"<<nomgamma<<"-"<<num<<"-"<<budget<<".dat";
//	ofstream out (ss.str().c_str());
//	out<< n << ";" << p << ";" << nomgamma <<"\n";
//    for (int i=0; i<n; ++i)
//        out<< cplex.getValue(c[i]) <<";";
//	out<<"\n";
//	for (int i=0; i<n; ++i)
//	{
//        out<< nomd[i] <<";";
//	}
//
//	out.close();
//	env.end();

	return 0;
}
