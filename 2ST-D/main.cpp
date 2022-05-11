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
	int N = atoi(argv[3]);
	int param = atoi(argv[4]);
	int R = atoi(argv[5]);

	Select sel;
	sel.set_timelimit(atoi(argv[6]));
	sel.set_budget(atoi(argv[7]));

	double start = clock();
	if (param == 0)
		sel.generate_rand(n,p,N,R);
	else if (param == 1)
		sel.generate_hard_C(n,p,N,R);
    else if (param == 2)
		sel.generate_hard_C_and_c(n,p,N,R);
	else
	{
		cout << "WARNING: THE FOURTH ARGUMENT COULD ONLY BE 0,1 OR 2" << endl;
		return 0;
	}

    //sel.print();

	cout<<"ALG-TI;"<<(clock()-start)/CLOCKS_PER_SEC<<"\n";

    //instance evaluation  has time limit 60
	start = clock();
	Solution sol = sel.solve_ip(60);
	double soltime = (clock()-start)/CLOCKS_PER_SEC;
	cout<<"SOLUTION-TIME;"<<soltime<<"\n";
	cout<<"NODES;"<<sol.nodes<<"\n";
	cout<<"OBJ;"<<sol.ub<<"\n";

	vector<vector<double> > c = sel.get_c();
	vector<double> C = sel.get_C();
	stringstream ss;
	ss<<"instance";
	for (int i=0; i<8; ++i)
		ss<<"-"<<argv[i+1];
	ss<<".dat";

	ofstream out(ss.str().c_str());
	out<<n<<";"<<p<<";"<<N<<"\n";
	for (int i=0; i<n; ++i)
        out<<C[i]<<";";
    out<<"\n";
	for (int i=0; i<N; ++i)
	{
        for (int j=0; j<n; ++j)
			out<<c[i][j]<<";";
		out<<"\n";
	}
	out.close();


	return 0;
}
