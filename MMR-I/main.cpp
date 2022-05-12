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
	srand(atoi(argv[7]));

	int n = atoi(argv[1]);
	int p = atoi(argv[2]);
	int param = atoi(argv[3]);
	int R = atoi(argv[4]);          //to choose what random method
	int t = atoi(argv[5]);
	double b = atof(argv[6]);



	Select sel;

	double start = clock();
	if (param == 0)
		sel.generate_rand(n,p,R);
	else if (param == 1)
		sel.generate_hard_c_d(n,p,R,b,t);
	else
	{
		cout << "WARNING: THE THIRD ARGUMENT COULD ONLY BE 0 OR 1" << endl;
		return 0;
	}

	vector<double> c = sel.get_c();
	vector<double> d = sel.get_d();

	stringstream ss;
	ss<<"instance";
	for (int i=0; i<7; ++i)
		ss<<"-"<<argv[i+1];
	ss<<".dat";

	ofstream out(ss.str().c_str());
	out<<n<<";"<<p<<"\n";
	for (int i=0; i<n; ++i)
        out<<c[i]<<";";
    out<<"\n";
    for (int i=0; i<n; ++i)
        out<<d[i]<<";";
	out.close();

	return 0;
}
