#include <ctime>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include "ilcplex/ilocplex.h"

ILOSTLBEGIN

using namespace std;

int main (int argc, char* argv[])
{

	//read file
	ifstream in(argv[1]);
	string line;

	getline(in,line);
	int n = atoi(line.c_str());
	size_t pos = line.find(";");
	line=line.substr(pos+1);

	int p = atoi(line.c_str());
	pos = line.find(";");
	line=line.substr(pos+1);

	int gamma = atoi(line.c_str());
	pos = line.find(";");
	line=line.substr(pos+1);

	int delta = atoi(line.c_str());


	vector<double> C(n);          //first stage scenario
    getline(in,line);

    for (int j=0; j<n; ++j)
    {
        C[j] = atof(line.c_str());
        pos = line.find(";");
        line=line.substr(pos+1);
    }

    vector<double> c(n);
    getline(in,line);

    for (int j=0; j<n; ++j)
    {
        c[j] = atof(line.c_str());
        pos = line.find(";");
        line=line.substr(pos+1);
    }

    vector<double> d(n);
    getline(in,line);

    for (int j=0; j<n; ++j)
    {
        d[j] = atof(line.c_str());
        pos = line.find(";");
        line=line.substr(pos+1);
    }

    //solve problem

	IloEnv env;
	IloModel model(env);

	//definition of the variables

	vector<IloNumVar> cplexx(n);
	for (int i=0; i<n; ++i)
		cplexx[i] = IloNumVar(env, 0, 1, ILOBOOL);

    vector<IloNumVar> cplexy(n);
	for (int i=0; i<n; ++i)
		cplexy[i] = IloNumVar(env, 0, 1, ILOFLOAT);

    vector<IloNumVar> cplexz(n);
	for (int i=0; i<n; ++i)
		cplexz[i] = IloNumVar(env, 0, 1, ILOBOOL);

    vector<IloNumVar> rho(n);
	for (int i=0; i<n; ++i)
		rho[i] = IloNumVar(env,0,IloInfinity,ILOFLOAT);

	IloNumVar pi(env,0,IloInfinity,ILOFLOAT);

	//definition of constraints

    {
        IloExpr obj(env);
        for (int i=0; i<n; ++i)
            obj += (C[i] * cplexx[i]) + (c[i] * cplexy[i]) + (d[i] * rho[i]);
        obj += gamma * pi;

        model.add(IloMinimize(env, obj));
    }

	{
		IloExpr con(env);
		for (int i=0; i<n; ++i)
			con += (cplexx[i]);
		model.add(con == p);
	}

	{
		IloExpr con(env);
		for (int i=0; i<n; ++i)
			con += (cplexy[i]);
		model.add(con == p);
	}

	{
		IloExpr con(env);
		for (int i=0; i<n; ++i)
			con += (cplexz[i]);
		model.add((p-delta) <= con);
	}

	for (int i=0; i<n; ++i)
	{
        model.add(cplexz[i] <= cplexx[i]);
        model.add(cplexz[i] <= cplexy[i]);
	}

	for (int i=0; i<n; ++i)
	{
        IloExpr con(env);
            con = (pi + rho[i]);

        model.add((cplexy[i]) <= con);
	}

	IloCplex cplex(model);

	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::Threads, 1);
	cplex.setParam(IloCplex::PreInd,0);
	cplex.setParam(IloCplex::TiLim, 600);

	double dstart = cplex.getDetTime();
	double tstart = clock();

	bool result = cplex.solve();

	double dend = cplex.getDetTime();
	double tend = clock();

	cout<<"ALGTIME;"<<(dend-dstart)/CLOCKS_PER_SEC<<"\n";
	cout<<"./TIME;"<<(tend-tstart)/CLOCKS_PER_SEC<<"\n";
	cout<<"NODES;"<<cplex.getNnodes()<<"\n";

	env.end();

}
