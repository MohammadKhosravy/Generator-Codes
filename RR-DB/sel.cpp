
#include "sel.h"
#include <set>
#include <iostream>
#include "ilcplex/ilocplex.h"

ILOSTLBEGIN

using namespace std;

void Sel::generate_RR_DB_U(int _n, int _p, int _gamma)
{
    n = _n;
    p = _p;
    gamma = _gamma;

    C.resize(n);
    c.resize(n);
    d.resize(n);

	for(int i=0; i<n; ++i)
    {
        C[i] = rand()%100+1;
        c[i] = rand()%100+1;
        d[i] = rand()%100+1;
    }

}

void Sel::generate_RR_DB_1(int _n, int _p, int _gamma)
{
    n = _n;
    p = _p;
    gamma = _gamma;

    C.resize(n);
    c.resize(n);
    d.resize(n);

    for(int i=0; i<n; ++i)
    {
        C[i] = rand()%100+1;
        c[i] = rand()%10+1;
    }

	for(int i=0; i<n; ++i)
        d[i] = rand()%(c[i]+1)+(100-c[i]);

}

void Sel::generate_RR_DB_2(int _n, int _p, int _gamma)
{
    n = _n;
    p = _p;
    gamma = _gamma;

    C.resize(n);
    c.resize(n);
    d.resize(n);

    for(int i=0; i<n; ++i)
        C[i] = rand()%100+1;

    for(int i=0; i<n; ++i)
        c[i] = 100 - C[i];

	for(int i=0; i<n; ++i)
        d[i] = rand()%(101-c[i])+(c[i]);

}



vector<int> Sel::get_C()
{
	return C;
}



vector<int> Sel::get_c()
{
	return c;
}



vector<int> Sel::get_d()
{
	return d;
}

