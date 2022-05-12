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
		sel.generate_hard_c(n,p,N,R);
	else
	{
		cout << "WARNING: THE FOURTH ARGUMENT COULD ONLY BE 0 OR 1" << endl;
		return 0;
	}

	cout<<"ALGTIME;"<<(clock()-start)/CLOCKS_PER_SEC<<"\n";

    //instance evaluation  has time limit 60
	start = clock();
	Solution sol = sel.solve_ip(60);
	double soltime = (clock()-start)/CLOCKS_PER_SEC;
	cout<<"./TIME;"<<soltime<<"\n";
	cout<<"NODES;"<<sol.nodes<<"\n";
	cout<<"OBJ;"<<sol.ub<<"\n";

	vector<vector<double> > c = sel.get_c();

	stringstream ss;
	ss<<"instance";
	for (int i=0; i<8; ++i)
		ss<<"-"<<argv[i+1];
	ss<<".dat";

	ofstream out(ss.str().c_str());
	out<<n<<";"<<p<<";"<<N<<"\n";
	for (int i=0; i<N; ++i)
	{
        for (int j=0; j<n; ++j)
			out<<c[i][j]<<";";
		out<<"\n";
	}
	out.close();


	return 0;
}
