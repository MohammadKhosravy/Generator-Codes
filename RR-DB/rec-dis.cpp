#include <ctime>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
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

    ////////////////////////////////////////////////////////////////////////////////

    set<double> setc;
    //setc.insert(0);
    for (int j=0; j<n; ++j)
    {
        setc.insert(c[j]);
        setc.insert(c[j]+d[j]);
    }

    vector<double> sizec;
    for (std::set<double>::iterator it=setc.begin(); it!=setc.end(); ++it)
    {
        sizec.push_back(*it);
    }

    int T = sizec.size();


    set<pair<double,double> > sets;
    for (int i=0; i<T; ++i)
        for (int j=0; j<T; ++j)
        {
            sets.insert(make_pair(sizec[i],max(0.0, sizec[j]-sizec[i])));
        }

    vector<pair<double,double> > sizes;
    for (std::set<pair<double,double> >::iterator it=sets.begin(); it!=sets.end(); ++it)
    {
        sizes.push_back(*it);
    }

    int S = sizes.size();

    vector<double> alpha(S);
    for (int i=0; i<S; ++i)
        alpha[i] = (sizes[i].first);

    vector<double> beta(S);
    for (int i=0; i<S; ++i)
        beta[i] = (sizes[i].second);


/*


    vector<double> pi;
    for (std::set<double>::iterator it=tempi.begin(); it!=tempi.end(); ++it)
    {
        pi.push_back(*it);
    }

    int S = pi.size();




    vector<double> cu(n);
    for (int i=0; i<n; ++i)
        cu[i] = (cl[i] + d[i]);


    int _n = 0;
    _n = 2 * n;

    vector<double> setc(_n);
    for(int i=0; i<(_n); ++i)
    {
        if (i<n)
            setc[i] = cl[i];
        else
            setc[i] = cu[i-n];
    }

    vector<vector<double> > sets1(_n);
    for(int i=0; i<(_n); ++i)
    {
        sets1[i].resize(_n);
        for(int j=0; j<(_n); ++j)
        {
            if(setc[i] < setc[j])
                sets1[i][j] = (setc[j]-setc[i]);
            else
                sets1[i][j] = 0;
        }
    }

    int m = 0;
    m = (_n) * (_n);

    vector<pair<double,double> > sets(m);
    for(int i=0; i<m; ++i)
    {
        int i1 = 0;
        int i2 = 0;
        i1 = int(i/(_n));
        i2 = i%(_n);

        sets[i] = make_pair(setc[i1],sets1[i1][i2]);
    }

    vector<double> alpha(m);
    for (int i=0; i<m; ++i)
        alpha[i] = (sets[i].first);

    vector<double> beta(m);
    for (int i=0; i<m; ++i)
        beta[i] = (sets[i].second);
*/
    ////////////////////////////////////////////////////////////////////////////////

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
            con += (C[i] * cplexx[i]) + rho[i][j];
        for (int i=0; i<n; ++i)
            con += -(((max(0.0,(alpha[j]+beta[j]-c[i]))-max(0.0,(alpha[j]-c[i])))*cplexx[i]) + (max(0.0,(alpha[j]-c[i]))));
        con += (gamma * (pi[j])) + (p * alpha[j]) + ((p-delta) * beta[j]);

        model.add(con <= lambda);
    }

    for (int j=0; j<S; ++j)
    {
        for (int i=0; i<n; ++i)
        {
            IloExpr con1(env);
            con1 = (((max(0.0,(alpha[j]+beta[j]-c[i]))-max(0.0,(alpha[j]-c[i])))*cplexx[i])+(max(0.0,(alpha[j]-c[i]))));

            IloExpr con2(env);
            con2 = (((max(0.0,(alpha[j]+beta[j]-c[i])-d[i])-max(0.0,(alpha[j]-c[i]-d[i])))*cplexx[i])+(max(0.0,(alpha[j]-c[i]-d[i]))));

            IloExpr lhs(env);
            lhs = (pi[j] + rho[i][j]);

            model.add((con1-con2) <= lhs);
        }
    }

	{
		IloExpr con(env);
		for (int i=0; i<n; ++i)
			con += (cplexx[i]);
		model.add(con == p);
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

	env.end();
}
