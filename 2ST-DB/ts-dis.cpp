#include <ctime>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include "ilcplex/ilocplex.h"
#include <bits/stdc++.h>
#include <cmath>
#include <utility>

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

	int k = atoi(line.c_str());


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

    ///////////////////////////////////////////////////////////

    set<double> temp;
    temp.insert(0);
    for (int j=0; j<n; ++j)
    {
        temp.insert(c[j]);
        temp.insert(c[j]+d[j]);
    }

    vector<double> alpha;
    for (std::set<double>::iterator it=temp.begin(); it!=temp.end(); ++it)
    {
        alpha.push_back(*it);
    }

    int S = alpha.size();

    ///////////////////////////////////////////////////////////

    //solve problem

	IloEnv env;
	IloModel model(env);

	//definition of the variables

	vector<IloNumVar> cplexx(n);
	for (int i=0; i<n; ++i)
		cplexx[i] = IloNumVar(env, 0, 1, ILOBOOL);

    vector<vector<IloNumVar> > rho(n);
	for (int i=0; i<n; ++i)
	{
        rho[i].resize(S);
        for (int j=0; j<S; ++j)
            rho[i][j] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
    }

    vector<IloNumVar> pi(S);
    for (int j=0; j<S; ++j)
        pi[j] = IloNumVar(env,0,IloInfinity,ILOFLOAT);

	IloNumVar lambda(env,0,IloInfinity,ILOFLOAT);

	//definition of constraints

    for (int j=0; j<S; ++j)
    {
        IloExpr con(env);
        for (int i=0; i<n; ++i)
            con += (C[i] * cplexx[i]) + rho[i][j] - (cplexx[i] * alpha[j]) - ((1 - cplexx[i]) * max(0.0,(alpha[j] - c[i])));
        con += (gamma * pi[j]) + (p * alpha[j]);

        model.add(con <= lambda);
    }

    {
		IloExpr con(env);
		for (int i=0; i<n; ++i)
			con += (cplexx[i]);
		model.add(con <= p);
	}

    for (int j=0; j<S; ++j)
    {
        for (int i=0; i<n; ++i)
        {
            IloExpr con(env);
            con =  (1 - cplexx[i]) * (max(0.0,(alpha[j] - c[i])) - max(0.0,(alpha[j] - c[i] - d[i])));

            IloExpr lhs(env);
            lhs = (pi[j] + rho[i][j]);

            model.add(con <= lhs);
        }
    }

	model.add(IloMinimize(env, lambda));

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
	//cout<<"objective:"<<cplex.getValue(lambda)<<"\n";

	env.end();
}
