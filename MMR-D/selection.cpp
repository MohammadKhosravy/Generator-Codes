#include "sel.h"
#include <algorithm>
#include <random>
#include <set>
#include "ilcplex/ilocplex.h"

ILOSTLBEGIN

using namespace std;


double timelimit = 30;
double globalstart;


void Select::generate_MMR_D_U()
{
    //costs
	c.resize(N);
	nomc.resize(N);
	for (int i=0; i<N; ++i)
	{
		c[i].resize(n,-1);
		nomc[i].resize(n,-1);
	}

	// find random lower and upper bounds
	cl.resize(N);
	cu.resize(N);

	for (int i=0; i<N; ++i)
	{
		cl[i].resize(n);
		cu[i].resize(n);

		for (int k=0; k<n; ++k)
		{
			c[i][k] = rand()%101;
			nomc[i][k] = c[i][k];

			cl[i][k] = c[i][k]-scenbudget;
			if (cl[i][k] < 0)
				cl[i][k] = 0;

			cu[i][k] = c[i][k]+scenbudget;
			if (cu[i][k] > 100)
				cu[i][k] = 100;
		}
	}

}

void Select::generate_MMR_D_1()
{
    //costs
	c.resize(N);
	nomc.resize(N);
	for (int i=0; i<N; ++i)
	{
		c[i].resize(n,-1);
		nomc[i].resize(n,-1);
	}

	// find random lower and upper bounds
	cl.resize(N);
	cu.resize(N);

	for (int i=0; i<N; ++i)
	{
		cl[i].resize(n);
		cu[i].resize(n);

		for (int k=0; k<n; ++k)
		{
			if (rand()%2==0)
                c[i][k] = rand()%10+1;
            else
                c[i][k] = rand()%10+91;

			nomc[i][k] = c[i][k];

			cl[i][k] = c[i][k]-scenbudget;
			if (cl[i][k] < 0)
				cl[i][k] = 0;

			cu[i][k] = c[i][k]+scenbudget;
			if (cu[i][k] > 100)
				cu[i][k] = 100;
		}
	}

}

void Select::generate_MMR_D_2()
{
    //costs
	c.resize(N);
	nomc.resize(N);
	for (int i=0; i<N; ++i)
	{
		c[i].resize(n,-1);
		nomc[i].resize(n,-1);
	}

	// find random lower and upper bounds
	cl.resize(N);
	cu.resize(N);

	for (int i=0; i<N; ++i)
	{
		cl[i].resize(n);
		cu[i].resize(n);

		for (int k=0; k<n; ++k)
		{
            if (k<(n/2))
                c[i][k] = rand()%100+1;
            else
                c[i][k] = 100 - c[i][k-(n/2)];

			nomc[i][k] = c[i][k];

			cl[i][k] = c[i][k]-scenbudget;
			if (cl[i][k] < 0)
				cl[i][k] = 0;

			cu[i][k] = c[i][k]+scenbudget;
			if (cu[i][k] > 100)
				cu[i][k] = 100;
		}
	}

}


void  Select::generate_hard_c(int _n, int _p, int _N, int _R)
{
    globalstart = clock();

	n = _n;
	p = _p;
	N = _N;
	R = _R;

    if (R==0)
        generate_MMR_D_U();
    else if (R==1)
        generate_MMR_D_1();
    else if (R==2)
        generate_MMR_D_2();

	Solution startsol = solve_ip();
    if (!status)
    {
        cout<<"Timelimit.\n";
        return;
    }

	x.push_back(startsol.x);

	double bestobj = startsol.ub;
	vector<vector<double> > bestc = c;
    //vector<double> bestC = C;

	double cstart = clock();

	//solve master
	double obj = solve_master();
    if (!status)
    {
        cout<<"Timelimit.\n";
        return;
    }

    double lowerbound = obj;

	//solve sub
	double start = clock();
	Solution sol = solve_ip();
	double end = clock();
    if (!status)
    {
        cout<<"Timelimit.\n";
        return;
    }

	if (sol.ub > bestobj - 0.01)
	{
		bestobj = sol.ub;
		bestc = c;
        //bestC = C;
	}


	int itnr = 1;
	cout<<itnr<<";"<<obj<<";"<<bestobj<<";"<<(end-start)/CLOCKS_PER_SEC<<"\n"<<flush;

	while(obj > bestobj + 0.01 && (clock()-cstart)/CLOCKS_PER_SEC < timelimit)
	{
		x.push_back(sol.x);

      double masterstart = clock();
      obj = solve_master();
      double masterend = clock();

    if (!status)
    {
        cout<<"Timelimit.\n";
        return;
    }

		start = clock();
		sol = solve_ip();
		end = clock();

    if (!status)
    {
        cout<<"Timelimit.\n";
        return;
    }

		if (sol.ub > bestobj - 0.01)
		{
			cout<<"*";
			bestobj = sol.ub;
			bestc = c;
            //bestC = C;
		}

		++itnr;
		cout<<itnr<<";"<<obj<<";"<<bestobj<<";"<<sol.ub<<";";
      cout<<(end-start)/CLOCKS_PER_SEC<<";"<<(masterend-masterstart)/CLOCKS_PER_SEC<<"\n"<<flush;
	}

//	c = bestc;

}

double Select::solve_master()
{
    int K = x.size();

    IloEnv env;
    IloModel model(env);

    vector<vector<IloNumVar> > cplexc(N);
    for (int i=0; i<N; ++i)
    {
        cplexc[i].resize(n);
        for (int k=0; k<n; ++k)
            cplexc[i][k] = IloNumVar(env,cl[i][k],cu[i][k],ILOFLOAT);
    }

    vector<vector<IloNumVar> > cplexy(N);
    for (int i=0; i<N; ++i)
    {
        cplexy[i].resize(n);
        for (int k=0; k<n; ++k)
            cplexy[i][k] = IloNumVar(env,0,1,ILOBOOL);
    }

    vector<vector<IloNumVar> > cplexlambda(K);
    for (int j=0; j<K; ++j)
    {
        cplexlambda[j].resize(N);
        for (int i=0; i<N; ++i)
            cplexlambda[j][i] = IloNumVar(env,0,1,ILOBOOL);
    }

    vector<vector<vector<IloNumVar> > > cplexd(N);
    for (int i=0; i<N; ++i)
    {
        cplexd[i].resize(K);
        for (int j=0; j<K; ++j)
        {
            cplexd[i][j].resize(n);
            for (int k=0; k<n; ++k)
                cplexd[i][j][k] = IloNumVar(env,0,cu[i][k],ILOFLOAT);
        }
    }

    vector<vector<vector<IloNumVar> > > cplexalpha(N);
    for (int i=0; i<N; ++i)
    {
        cplexalpha[i].resize(K);
        for (int j=0; j<K; ++j)
        {
            cplexalpha[i][j].resize(n);
            for (int k=0; k<n; ++k)
                cplexalpha[i][j][k] = IloNumVar(env,0,cu[i][k],ILOFLOAT);
        }
    }

    IloNumVar cplext(env,0,IloInfinity,ILOFLOAT);

    for (int j=0; j<K; ++j)
    {
        IloExpr con(env);
        for (int i=0; i<N; ++i)
            for (int k=0; k<n; ++k)
                con += ((x[j][k] * cplexd[i][j][k]) - cplexalpha[i][j][k]);
        model.add(cplext <= con);
    }

    for (int j=0; j<K; ++j)
    {
        IloExpr con(env);
        for (int i=0; i<N; ++i)
            con += cplexlambda[j][i];
        model.add(con == 1);
    }

    for (int i=0; i<N; ++i)
    {
        IloExpr con(env);
        for (int k=0; k<n; ++k)
            con += cplexy[i][k];
        model.add(con == p);
    }

    for (int i=0; i<N; ++i)
        for (int j=0; j<K; ++j)
            for (int k=0; k<n; ++k)
            {
                model.add(cplexd[i][j][k] <= cplexc[i][k]);
                model.add(cplexd[i][j][k] <= cu[i][k]*cplexlambda[j][i]);
            }

    for (int i=0; i<N; ++i)
        for (int j=0; j<K; ++j)
            for (int k=0; k<n; ++k)
            {
                model.add(cplexalpha[i][j][k] >= (cplexc[i][k] - (cu[i][k]*(2-cplexlambda[j][i]-cplexy[i][k]))));
            }

    for (int i=0; i<N; ++i)
    {
        double csum = 0;
        IloExpr con(env);
        for (int k=0; k<n; ++k)
        {
            csum += nomc[i][k];
            con += cplexc[i][k];
        }
        model.add(con <= csum);
    }

    model.add(IloMaximize(env, cplext));

    IloCplex cplex(model);

	//warmstart
	if (c[0][0] > -0.5)
	{
		IloNumVarArray startVar(env);
		IloNumArray startVal(env);
        for (int i=0; i<N; ++i)
            for (int k=0; k<n; ++k)
            {
               startVar.add(cplexc[i][k]);
               startVal.add(c[i][k]);
            }

		cplex.addMIPStart(startVar, startVal);
		startVal.end();
		startVar.end();
	}

	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::Threads, 1);
    cplex.setParam(IloCplex::TiLim, timelimit - (clock()-globalstart)/CLOCKS_PER_SEC);

	bool result = cplex.solve();

    double obj = 0;
    //cout<<cplex.getCplexStatus()<<"\n";
    if (cplex.getCplexStatus() == IloCplex::Optimal || cplex.getCplexStatus() == IloCplex::OptimalTol )
    {
        status = true;
        for (int i=0; i<N; ++i)
            for (int k=0; k<n; ++k)
                c[i][k] = cplex.getValue(cplexc[i][k]);
        obj = cplex.getObjValue();
    }
    else
        status = false;

	env.end();

	return obj;
}

void Select::generate_rand(int _n, int _p, int _N, int _R)
{
	n = _n;
	p = _p;
	N = _N;
    R = _R;

    if (R==0)
        generate_MMR_D_U();
    else if (R==1)
        generate_MMR_D_1();
    else if (R==2)
        generate_MMR_D_2();

};



Solution Select::solve_ip(double tlim)
{
	Solution sol;

	//int n = numelements;
    //int N = numscenarios;
    //vector< vector<double> > c = scenariocosts;

    vector<vector<pair<double,int> > > cindex(N);
    for (int i=0; i<N; ++i)
    {
        cindex[i].resize(n);
        for (int j=0; j<n; ++j)
            cindex[i][j] = make_pair(c[i][j],j);

        sort(cindex[i].begin(), cindex[i].end());
    }

    vector<vector<double> > sortc(N);
	for (int i=0; i<N; ++i)
    {
        sortc[i].resize(n);
        for (int j=0; j<n; ++j)
            sortc[i][j] = c[i][cindex[i][j].second];
    }

    vector<double> optc(N,0);
    for (int i=0; i<N; ++i)
    {
        for (int j=0; j<p; ++j)
            optc[i] += sortc[i][j];
    }

    IloEnv env;
    IloModel model(env);

    vector<IloNumVar> cplexx(n);
    for (int i=0; i<n; ++i)
        cplexx[i] = IloNumVar(env, 0, 1, ILOBOOL);

    IloNumVar cplexz(env,0,IloInfinity,ILOFLOAT);

    for (int i=0; i<N; ++i)
    {
        IloExpr con(env);
        for (int j=0; j<n; ++j)
            con += c[i][j] * cplexx[j];
        model.add((cplexz+optc[i]) >= con);
    }

    {
        IloExpr con(env);
        for (int j=0; j<n; ++j)
            con += cplexx[j];
        model.add(con == p);
    }

    model.add(IloMinimize(env, cplexz));

    IloCplex cplex(model);

	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::Threads, 1);
	cplex.setParam(IloCplex::PreInd, 0);
    if (tlim < 0)
        cplex.setParam(IloCplex::TiLim, timelimit - (clock()-globalstart)/CLOCKS_PER_SEC);
    else
        cplex.setParam(IloCplex::TiLim, tlim);

	bool result = cplex.solve();

    if (cplex.getCplexStatus() == IloCplex::Optimal || cplex.getCplexStatus() == IloCplex::OptimalTol )
    {
        status = true;
        sol.nodes = cplex.getNnodes();
        sol.ub = cplex.getObjValue();
        sol.x.resize(n);
        for (int j=0; j<n; ++j)
            sol.x[j] = cplex.getValue(cplexx[j]);
    }
    else
        status = false;

	env.end();

	return sol;
}

vector<vector<double> > Select::get_c()
{
	return c;
}

void Select::set_timelimit(double _timelimit)
{
	timelimit = _timelimit;
}

void Select::set_budget(int _scenbudget)
{
	scenbudget = _scenbudget;
}

void Select::print()
{
	for (int i=0; i<N; ++i)
	{
		for (int k=0; k<n; ++k)
			cout<<c[i][k]<<";";
		cout<<"\n";
	}
}
