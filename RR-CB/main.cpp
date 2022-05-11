#include "sel.h"
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>

using namespace std;

int main(int argc, char* argv[])
{
    Sel sel;

    srand(atoi(argv[5]));

    int n = atoi(argv[1]);
    int p = atoi(argv[2]);
    int gamma = atoi(argv[3]);
    int param = atoi(argv[4]);


    if (param == 0)
		sel.generate_RR_CB_U(n,p,gamma);
    else if (param == 1)
		sel.generate_RR_CB_1(n,p,gamma);
  	else if (param == 2)
    sel.generate_RR_CB_2(n,p,gamma);
    else
    cout << "Warning: the fourth argument should be 0, 1 or 2" << endl;

    vector<int> C = sel.get_C();
    vector<int> c = sel.get_c();
    vector<int> d = sel.get_d();

    stringstream ss;
	  ss<<"instance";
	  for (int i=0; i<5; ++i)
		ss<<"-"<<argv[i+1];
	  ss<<".dat";

  	ofstream out (ss.str().c_str());
	  out<<n<<";"<<p<<";"<<gamma<<"\n";
	  for (int i=0; i<n; ++i)
        out<<C[i]<<";";
    out<<"\n";
	  for (int i=0; i<n; ++i)
        out<<c[i]<<";";
    out<<"\n";
    for (int i=0; i<n; ++i)
        out<<d[i]<<";";
    out<<"\n";
	  out.close();

    return 0;
}
