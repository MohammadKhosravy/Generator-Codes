#include "sel.h"
#include <algorithm>
#include <random>
#include <set>
#include "ilcplex/ilocplex.h"

ILOSTLBEGIN

using namespace std;


void Select::generate_MMR_I_U()
{
    //costs
	c.resize(n);
	d.resize(n);
	nomc.resize(n);
	nomd.resize(n);
	for (int i=0; i<n; ++i)
	{
		c[i] = rand()%101;
		d[i] = rand()%101;

		nomc[i] = c[i];
		nomd[i] = d[i];
	}

}

void Select::generate_MMR_I_1()
{
    //costs
	c.resize(n);
	d.resize(n);
	nomc.resize(n);
	nomd.resize(n);
	for (int i=0; i<n; ++i)
	{
		if (rand()%2==0)
            c[i] = rand()%10+1;
        else
            c[i] = rand()%10+91;

        if (c[i]<=10)
            d[i] = rand()%10+91;
        else
            d[i] = rand()%10+1;

		nomc[i] = c[i];
		nomd[i] = d[i];
	}

}

void Select::generate_MMR_I_2()
{
    //costs
	c.resize(n);
	d.resize(n);
	nomc.resize(n);
	nomd.resize(n);
	for (int i=0; i<n; ++i)
	{
		if (rand()%2==0)
            c[i] = rand()%10+1;
        else
            c[i] = rand()%10+91;

        if (c[i]<=10)
            d[i] = rand()%10+1;
        else
            d[i] = rand()%10+91;

		nomc[i] = c[i];
		nomd[i] = d[i];
	}

}

void Select::generate_rand(int _n, int _p, int _R)
{
	n = _n;
	p = _p;
	R = _R;

	if (R == 0)
        generate_MMR_I_U();
    else if (R == 1)
        generate_MMR_I_1();
    else if (R == 2)
        generate_MMR_I_2();
}

void  Select::generate_hard_c_d(int _n, int _p, int _R, double _b, int _t)
{
    n = _n;
    p = _p;
    R = _R;
    b = _b;
    t = _t;

	if (R == 0)
        generate_MMR_I_U();
    else if (R == 1)
        generate_MMR_I_1();
    else if (R == 2)
        generate_MMR_I_2();


    set<double> tempi;
    tempi.insert(0);
    for (int j=0; j<n; ++j)
    {
        tempi.insert(nomc[j]);
        tempi.insert(nomc[j]+nomd[j]);
    }

    vector<double> pi;
    for (std::set<double>::iterator it=tempi.begin(); it!=tempi.end(); ++it)
    {
        pi.push_back(*it);
    }

    int S = pi.size();

    IloEnv env;
    IloModel model(env);

    vector<IloNumVar> cplexc(n);
    for (int i=0; i<n; ++i)
        cplexc[i] = IloNumVar(env,max(nomc[i] - b,0.0),nomc[i] + b,ILOFLOAT);

    vector<IloNumVar> cplexd(n);
    for (int i=0; i<n; ++i)
        cplexd[i] = IloNumVar(env,max(nomd[i]-b,0.0),nomd[i]+b,ILOFLOAT);

    vector<IloNumVar> cplexalpha(S);
    for (int k=0; k<S; ++k)
        cplexalpha[k] = IloNumVar(env,0,IloInfinity,ILOFLOAT);

    vector<vector<IloNumVar> > cplexbeta(S);
    for (int k=0; k<S; ++k)
    {
        cplexbeta[k].resize(n);
        for (int i=0; i<n; ++i)
            cplexbeta[k][i] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
    }

    vector<vector<IloNumVar> > cplexs(S);
    for (int k=0; k<S; ++k)
    {
        cplexs[k].resize(n);
        for (int i=0; i<n; ++i)
            cplexs[k][i] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
    }

    vector<vector<IloNumVar> > cplexz(S);
    for (int k=0; k<S; ++k)
    {
        cplexz[k].resize(n);
        for (int i=0; i<n; ++i)
            cplexz[k][i] = IloNumVar(env,0,1,ILOBOOL);
    }

    vector<vector<IloNumVar> > cplexzhat(S);
    for (int k=0; k<S; ++k)
    {
        cplexzhat[k].resize(n);
        for (int i=0; i<n; ++i)
            cplexzhat[k][i] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
    }

    vector<vector<IloNumVar> > cplexq(S);
    for (int k=0; k<S; ++k)
    {
        cplexq[k].resize(n);
        for (int i=0; i<n; ++i)
            cplexq[k][i] = IloNumVar(env,0,1,ILOBOOL);
    }

    vector<vector<IloNumVar> > cplexqhat(S);
    for (int k=0; k<S; ++k)
    {
        cplexqhat[k].resize(n);
        for (int i=0; i<n; ++i)
            cplexqhat[k][i] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
    }

    vector<vector<IloNumVar> > cplexqtilde(S);
    for (int k=0; k<S; ++k)
    {
        cplexqtilde[k].resize(n);
        for (int i=0; i<n; ++i)
            cplexqtilde[k][i] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
    }


    IloNumVar cplext(env,0,IloInfinity,ILOFLOAT);

    //constraints
    for (int k=0; k<S; ++k)
    {
        IloExpr con(env);
        for (int i=0; i<n; ++i)
            con += ((cplexz[k][i] * pi[k]) - (cplexzhat[k][i]));
        for (int i=0; i<n; ++i)
            con -= (cplexbeta[k][i]);
        con += (p * (cplexalpha[k] - pi[k]));

        model.add(cplext <= con);
    }

    for (int k=0; k<S; ++k)
        for (int i=0; i<n; ++i)
        {
            IloExpr lhs(env);
            lhs = cplexc[i] + cplexd[i] + pi[k]*cplexq[k][i] - cplexqhat[k][i] - cplexqtilde[k][i] - cplexs[k][i];

            IloExpr rhs(env);
            rhs = cplexalpha[k] - cplexbeta[k][i];

            model.add(rhs <= lhs);
        }

    for (int k=0; k<S; ++k)
        for (int i=0; i<n; ++i)
        {
            IloExpr con(env);
            //con = cplexs[k][i] - pi[k] + cplexc[i];

            model.add(cplexs[k][i] >= pi[k] - cplexc[i]);
        }

    for (int k=0; k<S; ++k)
        for (int i=0; i<n; ++i)
        {
            IloExpr con(env);
            con = cplexzhat[k][i] - cplexc[i] + (nomc[i]+b) * (1 - cplexz[k][i]);

            model.add(0 <= con);
        }


    for (int k=0; k<S; ++k)
        for (int i=0; i<n; ++i)
        {
            IloExpr con(env);
            con = cplexqhat[k][i] - cplexc[i] + (nomc[i]+b) * (1 - cplexq[k][i]);

            model.add(0 <= con);
        }


    for (int k=0; k<S; ++k)
        for (int i=0; i<n; ++i)
        {
            IloExpr con(env);
            con = cplexqtilde[k][i] - cplexd[i] + (nomd[i]+b) * (1 - cplexq[k][i]);

            model.add(0 <= con);
        }

/*    for (int i=0; i<n; ++i)
    {
        model.add(cplexc[i] <= ((1 + b) * nomc[i]));
        model.add(((1 - b) * nomc[i]) <= cplexc[i]);

        model.add(cplexd[i] <= ((1 + b) * nomd[i]));
        model.add(((1 - b) * nomd[i]) <= cplexd[i]);
    }*/
/*
    for (int i=0; i<n; ++i)
    {
        model.add(cplexc[i] <= (nomc[i] + b));
        model.add((nomc[i] - b) <= cplexc[i]);

        model.add(cplexd[i] <= (nomd[i] + b));
        model.add((nomd[i] - b) <= cplexd[i]);
    }
*/
    {
        double csum = 0;
		IloExpr con(env);
		for (int i=0; i<n; ++i)
		{
			csum += nomc[i];
			con += cplexc[i];
		}
		model.add(con <= csum);
    }

    {
        double dsum = 0;
		IloExpr con(env);
		for (int i=0; i<n; ++i)
		{
			dsum += nomd[i];
			con += cplexd[i];
		}
		model.add(con <= dsum);
    }

	model.add(IloMaximize(env, cplext));

    IloCplex cplex(model);

	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::Threads,1);
	cplex.setParam(IloCplex::PreInd,0);
    cplex.setParam(IloCplex::TiLim,t);

	double dstart = cplex.getDetTime();
	double tstart = clock();

	bool result = cplex.solve();

	double dend = cplex.getDetTime();
	double tend = clock();


    cout<<"ALG-TI;"<<(dend-dstart)/CLOCKS_PER_SEC<<"\n";
	cout<<"TIME;"<<(tend-tstart)/CLOCKS_PER_SEC<<"\n";
	cout<<"NODES;"<<cplex.getNnodes()<<"\n";
	cout<<"OBJ;"<<cplex.getValue(cplext)<<"\n";

	for (int i=0; i<n; ++i)
    {
        c[i] = cplex.getValue(cplexc[i]);
        d[i] = cplex.getValue(cplexd[i]);
    }

	env.end();

}

vector<double> Select::get_c()
{
	return c;
}

vector<double> Select::get_d()
{
	return d;
}
