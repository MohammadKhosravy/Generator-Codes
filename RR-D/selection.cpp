#include "sel.h"
#include <algorithm>
#include <random>
#include <set>

ILOSTLBEGIN

using namespace std;


double timelimit = 30;
double globalstart;


void Select::gen_RU()
{
    //first-stage costs
    C.resize(n);
    nomC.resize(n);
    Cl.resize(n);
	Cu.resize(n);
    for (int k=0; k<n; ++k)
    {
        C[k] = rand()%101;
        nomC[k] = C[k];

        Cl[k] = max(0.0, C[k] - scenbudget);
        Cu[k] = min(100.0, C[k] + scenbudget);
    }


    //second-stage costs
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

void Select::gen_R17()
{
    //first-stage costs
    C.resize(n);
    nomC.resize(n);
    Cl.resize(n);
	Cu.resize(n);
    for (int k=0; k<n; ++k)
    {
        if (rand()%2==0)
            C[k] = rand()%11+45;
        else
            C[k] = rand()%51+25;

        nomC[k] = C[k];

        Cl[k] = max(0.0, C[k] - scenbudget);
        Cu[k] = min(100.0, C[k] + scenbudget);
    }


    //second-stage costs
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
                c[i][k] = rand()%11+C[k]-5;
            else
            {
                if (rand()%2==0)
                    c[i][k] = rand()%10+1;
                else
                    c[i][k] = rand()%10+91;
            }

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

void Select::gen_R21()
{
    //first-stage costs
    C.resize(n);
    nomC.resize(n);
    Cl.resize(n);
	Cu.resize(n);
    for (int k=0; k<n; ++k)
    {
        if (rand()%2==0)
            C[k] = rand()%100+1;
        else
            C[k] = 50;

        nomC[k] = C[k];

        Cl[k] = max(0.0, C[k] - scenbudget);
        Cu[k] = min(100.0, C[k] + scenbudget);
    }


    //second-stage costs
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
			if (C[k]==50)
            {
                if (rand()%2==0)
                    c[i][k]= rand()%10+1;
                else
                    c[i][k]= rand()%10+91;
            }

            else
                c[i][k]= rand()%11+C[k]-5;

			nomc[i][k] = max(0.0 , c[i][k]);

			cl[i][k] = c[i][k]-scenbudget;
			if (cl[i][k] < 0)
				cl[i][k] = 0;

			cu[i][k] = c[i][k]+scenbudget;
			if (cu[i][k] > 100)
				cu[i][k] = 100;
		}
	}

}

void  Select::generate_hard_C(int _n, int _p, int _N, int _delta, int _R)
{
    globalstart = clock();

	n = _n;
	p = _p;
	N = _N;
	delta = _delta;
	R = _R;

	if (R == 0)
        gen_RU();
    else if (R == 17)
        gen_R17();
    else if (R == 21)
        gen_R21();

	Solution startsol = solve_ip();
    if (!status)
    {
        cout<<"Timelimit.\n";
        return;
    }

	x.push_back(startsol.x);

	double bestobj = startsol.ub;
	vector<vector<double> > bestc = c;
    vector<double> bestC = C;

	double cstart = clock();

	//solve master
	double obj = solve_master_C();
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
        bestC = C;
	}


	int itnr = 1;
	cout<<itnr<<";"<<obj<<";"<<bestobj<<";"<<(end-start)/CLOCKS_PER_SEC<<"\n"<<flush;

	while(obj > bestobj + 0.01 && (clock()-cstart)/CLOCKS_PER_SEC < timelimit)
	{
		x.push_back(sol.x);

      double masterstart = clock();
      obj = solve_master_C();
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
            bestC = C;
		}

		++itnr;
		cout<<itnr<<";"<<obj<<";"<<bestobj<<";"<<sol.ub<<";";
      cout<<(end-start)/CLOCKS_PER_SEC<<";"<<(masterend-masterstart)/CLOCKS_PER_SEC<<"\n"<<flush;
	}

//	c = bestc;

}


void  Select::generate_hard_C_and_c(int _n, int _p, int _N, int _delta, int _R)
{
	n = _n;
	p = _p;
	N = _N;
	delta = _delta;
    R = _R;

	if (R == 0)
        gen_RU();
    else if (R == 17)
        gen_R17();
    else if (R == 21)
        gen_R21();


	Solution startsol = solve_ip();
	x.push_back(startsol.x);

    if (!status)
    {
        cout<<"Timelimit.\n";
        return;
    }

	double bestobj = startsol.ub;
	vector<vector<double> > bestc = c;
    vector<double> bestC = C;

	double cstart = clock();

	//solve master
	double obj = solve_master_C_and_c();
    double lowerbound = obj;

    if (!status)
    {
        cout<<"Timelimit.\n";
        return;
    }

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
        bestC = C;
	}


	int itnr = 1;
	cout<<itnr<<";"<<obj<<";"<<bestobj<<";"<<(end-start)/CLOCKS_PER_SEC<<"\n"<<flush;

	while(obj > bestobj + 0.01 && (clock()-cstart)/CLOCKS_PER_SEC < timelimit)
	{
		x.push_back(sol.x);

      double masterstart = clock();
      obj = solve_master_C_and_c();
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
            bestC = C;
		}

		++itnr;
		cout<<itnr<<";"<<obj<<";"<<bestobj<<";"<<sol.ub<<";";
      cout<<(end-start)/CLOCKS_PER_SEC<<";"<<(masterend-masterstart)/CLOCKS_PER_SEC<<"\n"<<flush;
	}

//	c = bestc;

}

double Select::solve_master_C()
{
	int K = x.size();

	IloEnv env;
	IloModel model(env);

    //variables

    vector<IloNumVar> cplexC(n);
    for (int k=0; k<n; ++k)
        cplexC[k] = IloNumVar(env,Cl[k],Cu[k],ILOFLOAT);

	vector<vector<IloNumVar> > cplexlambda(K);
	for (int j=0; j<K; ++j)
	{
		cplexlambda[j].resize(N);
        for (int i=0; i<N; ++i)
            cplexlambda[j][i] = IloNumVar(env,0,1,ILOBOOL);
	}

	IloNumVar cplext(env,0,IloInfinity,ILOFLOAT);

    vector<IloNumVar> cplexbeta(K);
    for (int j=0; j<K; ++j)
        cplexbeta[j] = IloNumVar(env,0,IloInfinity,ILOFLOAT);

    vector<IloNumVar> cplexeta(K);
    for (int j=0; j<K; ++j)
        cplexeta[j] = IloNumVar(env,0,IloInfinity,ILOFLOAT);

    vector<vector<IloNumVar> > cplexgamma(K);
    for (int j=0; j<K; ++j)
    {
        cplexgamma[j].resize(n);
        for (int k=0; k<n; ++k)
            cplexgamma[j][k] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
    }


    //constraints

	for (int j=0; j<K; ++j)
	{
		IloExpr con(env);
        for (int k=0; k<n; ++k)
            con += cplexC[k]*x[j][k];

        //double betaval=p;
        //for (int k=0; k<n; ++k)
            //betaval -= x[j][k];
        con += (p*cplexbeta[j]) + ((p-delta)*cplexeta[j]);

        for (int k=0; k<n; ++k)
            con -= cplexgamma[j][k];

        model.add(cplext <= con);
	}

	for (int j=0; j<K; ++j)
	{
		IloExpr con(env);
		for (int i=0; i<N; ++i)
			con += cplexlambda[j][i];
		model.add(con == 1);
	}

    for (int j=0; j<K; ++j)
        for (int k=0; k<n; ++k)
        {
            IloExpr lhs(env);
            lhs = cplexbeta[j] + (cplexeta[j] * x[j][k]);

            IloExpr con(env);
            for (int i=0; i<N; ++i)
                con += nomc[i][k]*cplexlambda[j][i];
            con += cplexgamma[j][k];

            model.add(lhs <= con);
        }


    {
        double csum = 0;
		IloExpr con(env);
		for (int k=0; k<n; ++k)
		{
			csum += nomC[k];
			con += cplexC[k];
		}
		model.add(con == csum);
    }

	model.add(IloMaximize(env, cplext));

	IloCplex cplex(model);


	//warmstart
	if (C[0] > -0.5)
	{
		IloNumVarArray startVar(env);
		IloNumArray startVal(env);
        for (int k=0; k<n; ++k)
        {
            startVar.add(cplexC[k]);
            startVal.add(C[k]);
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
        for (int k=0; k<n; ++k)
            C[k] = cplex.getValue(cplexC[k]);
        obj = cplex.getObjValue();
    }
    else
        status = false;

	env.end();

	return obj;
}

double Select::solve_master_C_and_c()
{
	int K = x.size();

	IloEnv env;
	IloModel model(env);

    //variables

    vector<IloNumVar> cplexC(n);
    for (int k=0; k<n; ++k)
        cplexC[k] = IloNumVar(env,Cl[k],Cu[k],ILOFLOAT);

	vector<vector<IloNumVar> > cplexc(N);
	for (int i=0; i<N; ++i)
	{
		cplexc[i].resize(n);
		for (int k=0; k<n; ++k)
			cplexc[i][k] = IloNumVar(env,cl[i][k],cu[i][k],ILOFLOAT);
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

	IloNumVar cplext(env,0,IloInfinity,ILOFLOAT);

    vector<IloNumVar> cplexbeta(K);
    for (int j=0; j<K; ++j)
        cplexbeta[j] = IloNumVar(env,0,IloInfinity,ILOFLOAT);

    vector<IloNumVar> cplexeta(K);
    for (int j=0; j<K; ++j)
        cplexeta[j] = IloNumVar(env,0,IloInfinity,ILOFLOAT);

    vector<vector<IloNumVar> > cplexgamma(K);
    for (int j=0; j<K; ++j)
    {
        cplexgamma[j].resize(n);
        for (int k=0; k<n; ++k)
            cplexgamma[j][k] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
    }


    //constraints

	for (int j=0; j<K; ++j)
	{
		IloExpr con(env);
        for (int k=0; k<n; ++k)
            con += cplexC[k]*x[j][k];

        con += (p*cplexbeta[j]) + ((p-delta)*cplexeta[j]);

        for (int k=0; k<n; ++k)
            con -= cplexgamma[j][k];

        model.add(cplext <= con);
	}

	for (int j=0; j<K; ++j)
	{
		IloExpr con(env);
		for (int i=0; i<N; ++i)
			con += cplexlambda[j][i];
		model.add(con == 1);
	}

	for (int j=0; j<K; ++j)
        for (int k=0; k<n; ++k)
        {
            IloExpr lhs(env);
            lhs = cplexbeta[j] + (cplexeta[j] * x[j][k]);

            IloExpr con(env);
            for (int i=0; i<N; ++i)
                con += cplexd[i][j][k];
            con += cplexgamma[j][k];

            model.add(lhs <= con);
        }

	for (int i=0; i<N; ++i)
		for (int j=0; j<K; ++j)
			for (int k=0; k<n; ++k)
			{
				model.add(cplexd[i][j][k] <= cplexc[i][k]);
				model.add(cplexd[i][j][k] <= cu[i][k]*cplexlambda[j][i]);
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
		model.add(con == csum);
	}

    {
        double csum = 0;
		IloExpr con(env);
		for (int k=0; k<n; ++k)
		{
			csum += nomC[k];
			con += cplexC[k];
		}
		model.add(con == csum);
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
        for (int k=0; k<n; ++k)
        {
            startVar.add(cplexC[k]);
            startVal.add(C[k]);
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

    if (cplex.getCplexStatus() == IloCplex::Optimal || cplex.getCplexStatus() == IloCplex::OptimalTol )
    {
        status = true;

        for (int k=0; k<n; ++k)
            C[k] = cplex.getValue(cplexC[k]);

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



void Select::generate_rand(int _n, int _p, int _N, int _delta, int _R)
{
	n = _n;
	p = _p;
	N = _N;
	delta = _delta;
    R = _R;

	if (R == 0)
        gen_RU();
    else if (R == 17)
        gen_R17();
    else if (R == 21)
        gen_R21();

}



Solution Select::solve_ip(double tlim)
{
	Solution sol;

	IloEnv env;
	IloModel model(env);

	vector<IloNumVar> cplexx(n);
	for (int i=0; i<n; ++i)
		cplexx[i] = IloNumVar(env, 0, 1, ILOBOOL);

    vector<vector<IloNumVar> > cplexy(N);
    for (int i=0; i<N; ++i)
    {
        cplexy[i].resize(n);
        for (int j=0; j<n; ++j)
            cplexy[i][j] = IloNumVar(env, 0, 1, ILOBOOL);
    }

    vector<vector<IloNumVar> > cplexz(N);
    for (int i=0; i<N; ++i)
    {
        cplexz[i].resize(n);
        for (int j=0; j<n; ++j)
            cplexz[i][j] = IloNumVar(env,0,1,ILOBOOL);
    }

	IloNumVar cplexs(env,0,IloInfinity,ILOFLOAT);

	for (int i=0; i<N; ++i)
	{
		IloExpr con(env);
		for (int j=0; j<n; ++j)
			con += c[i][j] * cplexy[i][j];

		model.add(cplexs >= con);
	}

	{
        IloExpr con(env);
        for(int i=0; i<n; i++)
            con += cplexx[i];

        model.add(con == p);
    }

    for(int i=0; i<N; i++)
    {
        IloExpr con(env);
        for(int j=0; j<n; j++)
            con += cplexy[i][j];

        model.add(con == p);
    }

    for(int i=0; i<N; i++)
    {
        IloExpr con(env);
        for(int j=0; j<n; j++)
            con += cplexz[i][j];

        model.add(p-delta <= con);
    }

    for(int i=0; i<N; i++)
    {
        for(int j=0; j<n; j++)
            model.add(cplexz[i][j] <= cplexx[j]);
    }

    for(int i=0; i<N; i++)
    {
        for(int j=0; j<n; j++)
            model.add(cplexz[i][j] <= cplexy[i][j]);
    }

    IloExpr obj(env);
    for (int j=0; j<n; ++j)
        obj += C[j]*cplexx[j];

	model.add(IloMinimize(env, obj + cplexs));

	IloCplex cplex(model);

	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::Threads, 1);
	cplex.setParam(IloCplex::PreInd, 0);
    if (tlim < 0)
        cplex.setParam(IloCplex::TiLim, timelimit - (clock()-globalstart)/CLOCKS_PER_SEC);
    else
        cplex.setParam(IloCplex::TiLim, tlim);

	bool result = cplex.solve();

	//cout<<cplex.getCplexStatus()<<"\n";
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



vector<double> Select::get_C()
{
	return C;
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
    for (int k=0; k<n; ++k)
        cout<<C[k]<<";";
    cout<<"\n";

	for (int i=0; i<N; ++i)
	{
		for (int k=0; k<n; ++k)
			cout<<c[i][k]<<";";
		cout<<"\n";
	}
}
