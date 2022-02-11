#include "sel.h"
#include <algorithm>
#include <random>
#include <set>
#include "ilcplex/ilocplex.h"

ILOSTLBEGIN

using namespace std;

//void Select::generate_RU_dis1()
//{
//    c.resize(n);
//    d.resize(n);
//    nomc.resize(n);
//    nomd.resize(n);
//
//    for (int i=0; i<n; ++i)
//    {
//        c[i] = rand()%101;
//        d[i] = rand()%101;
//
//        nomc[i] = c[i];
//        nomd[i] = d[i];
//    }
//
//}
//
//void Select::generate_RU_con1()
//{
//
//    c.resize(n);
//    d.resize(n);
//    nomc.resize(n);
//    nomd.resize(n);
//
//    for (int i=0; i<n; ++i)
//    {
//        c[i] = fRand(0.0 , 100.0);
//        d[i] = fRand(0.0 , 100.0);
//
//        nomc[i] = c[i];
//        nomd[i] = d[i];
//    }
//
//}

void Select::generate_RU_dis2()
{

    c.resize(n);
    d.resize(n);
    nomc.resize(n);
    nomd.resize(n);

    for (int i=0; i<n; ++i)
    {
        c[i] = rand()%100+1;
        d[i] = rand()%100+1;

        nomc[i] = c[i];
        nomd[i] = d[i];
    }

//    c.resize(n);
//    d.resize(n);
//    nomc.resize(n);
//    nomd.resize(n);
//
//    for (int i=0; i<n; ++i)
//    {
//        c[i] = 30;
//        d[i] = 3;
//
//        nomc[i] = c[i];
//        nomd[i] = d[i];
//    }

}

//void Select::generate_RU_con2()
//{
//
//    c.resize(n);
//    d.resize(n);
//    nomc.resize(n);
//    nomd.resize(n);
//
//    for (int i=0; i<n; ++i)
//    {
//        c[i] = fRand(1.0 , 100);
//        d[i] = fRand(1.0 , 100);
//
//        nomc[i] = c[i];
//        nomd[i] = d[i];
//    }
//
//}
//
//void Select::generate_R1_dis1()
//{
//
//    c.resize(n);
//    d.resize(n);
//    nomc.resize(n);
//    nomd.resize(n);
//
//    for (int i=0; i<n; ++i)
//    {
//        c[i] = rand()%101;
//        d[i] = 100 - c[i];
//
//        nomc[i] = c[i];
//        nomd[i] = d[i];
//    }
//
//}
//
//void Select::generate_R1_con1()
//{
//
//    c.resize(n);
//    d.resize(n);
//    nomc.resize(n);
//    nomd.resize(n);
//
//    for (int i=0; i<n; ++i)
//    {
//        c[i] = fRand(0.0 , 100.0);
//        d[i] = 100.0 - c[i];
//
//        nomc[i] = c[i];
//        nomd[i] = d[i];
//    }
//
//}

void Select::generate_R1_dis2()
{

    c.resize(n);
    d.resize(n);
    nomc.resize(n);
    nomd.resize(n);

    for (int i=0; i<n; ++i)
    {
        c[i] = rand()%100+1;
        d[i] = 100 - c[i];

        nomc[i] = c[i];
        nomd[i] = d[i];
    }

}

//void Select::generate_R1_con2()
//{
//
//    c.resize(n);
//    d.resize(n);
//    nomc.resize(n);
//    nomd.resize(n);
//
//    for (int i=0; i<n; ++i)
//    {
//        c[i] = fRand(1.0 , 100.0);
//        d[i] = 100.0 - c[i];
//
//        nomc[i] = c[i];
//        nomd[i] = d[i];
//    }
//
//}
//
//void Select::generate_R5_dis1()
//{
//
//    c.resize(n);
//    d.resize(n);
//    nomc.resize(n);
//    nomd.resize(n);
//
//    for (int i=0; i<n; ++i)
//    {
//        c[i] = rand()%11;
//        d[i] = rand()%(int(c[i])+1)+(100-c[i]);
//
//        nomc[i] = c[i];
//        nomd[i] = d[i];
//    }
//
//}
//
//void Select::generate_R5_con1()
//{
//
//    c.resize(n);
//    d.resize(n);
//    nomc.resize(n);
//    nomd.resize(n);
//
//    for (int i=0; i<n; ++i)
//    {
//        c[i] = fRand(0.0 , 10.0);
//        d[i] = fRand(100-c[i] , 100.0);
//
//        nomc[i] = c[i];
//        nomd[i] = d[i];
//    }
//
//}

void Select::generate_R5_dis2()
{

    c.resize(n);
    d.resize(n);
    nomc.resize(n);
    nomd.resize(n);

    for (int i=0; i<n; ++i)
    {
        c[i] = rand()%10+1;
        d[i] = 100 - rand()%(int (c[i])+1);

        nomc[i] = c[i];
        nomd[i] = d[i];
    }

}

//void Select::generate_R5_con2()
//{
//
//    c.resize(n);
//    d.resize(n);
//    nomc.resize(n);
//    nomd.resize(n);
//
//    for (int i=0; i<n; ++i)
//    {
//        c[i] = fRand(1.0 , 10.0);
//        d[i] = fRand(100-c[i] , 100.0);
//
//        nomc[i] = c[i];
//        nomd[i] = d[i];
//    }
//
//}


void Select::generate_rand(int _n, int _p, int _gamma, int _R)
{
	n = _n;
	p = _p;
	gamma = _gamma;
	R = _R;

	if (R == 0)
        generate_RU_dis2();
    else if (R == 1)
        generate_R1_dis2();
    else if (R == 5)
        generate_R5_dis2();


//	if (R == 1)
//        generate_RU_dis1();
//    else if (R == 2)
//        generate_RU_con1();
//    else if (R == 3)
//        generate_RU_dis2();
//    else if (R == 4)
//        generate_RU_con2();
//    else if (R == 5)
//        generate_R1_dis1();
//    else if (R == 6)
//        generate_R1_con1();
//    else if (R == 7)
//        generate_R1_dis2();
//    else if (R == 8)
//        generate_R1_con2();
//    else if (R == 9)
//        generate_R5_dis1();
//    else if (R == 10)
//        generate_R5_con1();
//    else if (R == 11)
//        generate_R5_dis2();
//    else if (R == 12)
//        generate_R5_con2();
}

void  Select::generate_hard_c_1(int _n, int _p, int _gamma, int _R, float _b, int _t)
{
    n = _n;
    p = _p;
    gamma = _gamma;
    R = _R;
    b = _b;
    t = _t;

   	if (R == 0)
        generate_RU_dis2();
    else if (R == 1)
        generate_R1_dis2();
    else if (R == 5)
        generate_R5_dis2();


//	if (R == 1)
//        generate_RU_dis1();
//    else if (R == 2)
//        generate_RU_con1();
//    else if (R == 3)
//        generate_RU_dis2();
//    else if (R == 4)
//        generate_RU_con2();
//    else if (R == 5)
//        generate_R1_dis1();
//    else if (R == 6)
//        generate_R1_con1();
//    else if (R == 7)
//        generate_R1_dis2();
//    else if (R == 8)
//        generate_R1_con2();
//    else if (R == 9)
//        generate_R5_dis1();
//    else if (R == 10)
//        generate_R5_con1();
//    else if (R == 11)
//        generate_R5_dis2();
//    else if (R == 12)
//        generate_R5_con2();

    //////////////////////////////////////

    int K = n + 1;

    vector<double> pi(K,0);
	for (int k=0; k<(K-1); ++k)
		pi[k] = nomd[k];

//    pi[K-1] = 0;

    vector<double> minc(n);
	for (int i=0; i<n; ++i)
		minc[i] = max (0.0 , nomc[i] - b);

	vector<double> maxc(n);
	for (int i=0; i<n; ++i)
		maxc[i] = min (100.0 , nomc[i] + b);

    //////////////////////////////////////

//    set<double> tempi;
//    tempi.insert(0);
//    for (int j=0; j<n; ++j)
//        tempi.insert(nomd[j]);
//
//
//    vector<double> pi;
//    for (std::set<double>::iterator it=tempi.begin(); it!=tempi.end(); ++it)
//    {
//        pi.push_back(*it);
//    }
//
//    int K = pi.size();

    //////////////////////////////////////

	IloEnv env;
	IloModel model(env);


    vector<IloNumVar> cplexc(n);
	for (int i=0; i<n; ++i)
		cplexc[i] = IloNumVar(env, minc[i], maxc[i], ILOFLOAT);


	vector<IloNumVar> alpha(K);
	for (int k=0; k<K; ++k)
		alpha[k] = IloNumVar(env, 0, IloInfinity, ILOFLOAT);


    vector<vector<IloNumVar> > beta(K);
	for (int k=0; k<K; ++k)
    {
        beta[k].resize(n);
        for (int i=0; i<n; ++i)
            beta[k][i] = IloNumVar(env, 0, IloInfinity, ILOFLOAT);
    }

	IloNumVar cplext(env,0,IloInfinity,ILOFLOAT);

	//constraints

    for (int k=0; k<K; ++k)
    {
        IloExpr con(env);
            con = p * alpha[k];
        for (int i=0; i<n; ++i)
            con += (-beta[k][i]);

        model.add(cplext <= (con + (gamma * pi[k])));
	}


	for (int k=0; k<K; ++k)
    {
        for (int i=0; i<n; ++i)
        {
            IloExpr con(env);
                con = cplexc[i] + max (0.0 , nomd[i] - pi[k]);
            IloExpr lhs(env);
                lhs = alpha[k] - beta[k][i] ;

            model.add(lhs <= con);

        }
    }

//    for (int i=0; i<n; ++i)
//        model.add((nomc[i] - b) <= cplexc[i]);
//
//	for (int i=0; i<n; ++i)
//        model.add(cplexc[i] <= (nomc[i] + b));
//
//    for (int i=0; i<n; ++i)
//        model.add(((1 - b) * nomc[i]) <= cplexc[i]);
//
//	for (int i=0; i<n; ++i)
//        model.add(cplexc[i] <= ((1 + b) * nomc[i]));

	{
        IloExpr con(env);
        for (int i=0; i<n; ++i)
            con += nomc[i];
        IloExpr lhs(env);
        for (int i=0; i<n; ++i)
            lhs += cplexc[i];

        model.add(lhs <= con);
	}

	model.add(IloMaximize(env, cplext));

	IloCplex cplex(model);

	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::Threads, 1);
	cplex.setParam(IloCplex::PreInd,0);
	cplex.setParam(IloCplex::TiLim,t);

	double dstart = cplex.getDetTime();
	double tstart = clock();

	bool result = cplex.solve();

	double dend = cplex.getDetTime();
	double tend = clock();

	// printing statistics
	cout<<"ALGTIME;"<<(dend-dstart)/CLOCKS_PER_SEC<<"\n";
	cout<<"./TIME;"<<(tend-tstart)/CLOCKS_PER_SEC<<"\n";
	cout<<"NODES;"<<cplex.getNnodes()<<"\n";

	for (int i=0; i<n; ++i)
        c[i] = cplex.getValue(cplexc[i]);

	env.end();

}

void  Select::generate_hard_d_1(int _n, int _p, int _gamma, int _R, float _b, int _t)
{
    n = _n;
    p = _p;
    gamma = _gamma;
    R = _R;
    b = _b;
    t = _t;

   	if (R == 0)
        generate_RU_dis2();
    else if (R == 1)
        generate_R1_dis2();
    else if (R == 5)
        generate_R5_dis2();


//    if (R == 1)
//        generate_RU_dis1();
//    else if (R == 2)
//        generate_RU_con1();
//    else if (R == 3)
//        generate_RU_dis2();
//    else if (R == 4)
//        generate_RU_con2();
//    else if (R == 5)
//        generate_R1_dis1();
//    else if (R == 6)
//        generate_R1_con1();
//    else if (R == 7)
//        generate_R1_dis2();
//    else if (R == 8)
//        generate_R1_con2();
//    else if (R == 9)
//        generate_R5_dis1();
//    else if (R == 10)
//        generate_R5_con1();
//    else if (R == 11)
//        generate_R5_dis2();
//    else if (R == 12)
//        generate_R5_con2();

    //////////////////////////////////////

    vector<pair<double,int> > nomdindex(n);
    for (int i=0; i<n; ++i)
        nomdindex[i] = make_pair(nomd[i],i);

    sort(nomdindex.begin(), nomdindex.end());

    vector<double> sortc(n);
	for (int i=0; i<n; ++i)
		sortc[i] = nomc[nomdindex[i].second];

    vector<double> sortd(n);
	for (int i=0; i<n; ++i)
		sortd[i] = nomd[nomdindex[i].second];

	vector<double> mind(n);
	for (int i=0; i<n; ++i)
		mind[i] = max (0.0 , sortd[i] - b);

	vector<double> maxd(n);
	for (int i=0; i<n; ++i)
		maxd[i] = min (100.0 , sortd[i] + b);

    int K = n + 1;

    ///////////////////////////////////////

    //solve problem

	IloEnv env;
	IloModel model(env);

    //variables

    vector<IloNumVar> cplexd(n);
	for (int i=0; i<n; ++i)
		cplexd[i] = IloNumVar(env, mind[i], maxd[i], ILOFLOAT);


	vector<IloNumVar> alpha(K);
	for (int k=0; k<K; ++k)
		alpha[k] = IloNumVar(env, 0, IloInfinity, ILOFLOAT);


    vector<vector<IloNumVar> > beta(K);
	for (int k=0; k<K; ++k)
    {
        beta[k].resize(n);
        for (int i=0; i<n; ++i)
            beta[k][i] = IloNumVar(env, 0, IloInfinity, ILOFLOAT);
    }

	IloNumVar cplext(env,0,IloInfinity,ILOFLOAT);

	//constraints

    for (int k=0; k<n; ++k)
    {
        IloExpr con(env);
            con = p * alpha[k];
        for (int i=0; i<n; ++i)
            con -= beta[k][i];

        model.add(cplext <= con + (gamma * cplexd[k]));
	}

	for (int k=0; k<n; ++k)
    {
        for (int i=0; i<n; ++i)
        {
            IloExpr con(env);
                con = sortc[i] + (cplexd[i] - cplexd[k]);
            IloExpr lhs(env);
                lhs = alpha[k] - beta[k][i] ;

            if ((k) < i)
                model.add(lhs <= con);
            else
                model.add(lhs <= sortc[i]);
        }
    }

    {
        IloExpr con(env);
            con = p * alpha[0];
        for (int i=0; i<n; ++i)
            con -= beta [0][i];

        model.add(cplext <= con);
    }

    for (int i=0; i<n; ++i)
    {
        IloExpr lhs(env);
            lhs = (alpha[0] - beta[0][i]);
        IloExpr con(env);
            con = (sortc[i] + cplexd[i]);

        model.add(lhs <= con);
    }


//    for (int i=0; i<n; ++i)
//        model.add(max(0.0 , sortd[i]-b) <= cplexd[i]);

//    model.add(15 <= cplexd[0]);

//    for (int i=0; i<n; ++i)
//        model.add(max(0.0 , sortd[i]-b) <= cplexd[i]);

//	for (int i=0; i<n; ++i)
//        model.add(cplexd[i] <= min(100.0 , sortd[i]+b));

	{
        IloExpr con(env);
        for (int i=0; i<n; ++i)
            con += nomd[i];
        IloExpr lhs(env);
        for (int i=0; i<n; ++i)
            lhs += cplexd[i];

        model.add(lhs <= con);
	}

	model.add(IloMaximize(env, cplext));

	/////////////////////////////////////////////////

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

	// printing statistics

	cout<<"ALGTIME;"<<(dend-dstart)/CLOCKS_PER_SEC<<"\n";
	cout<<"./TIME;"<<(tend-tstart)/CLOCKS_PER_SEC<<"\n";
	cout<<"NODES;"<<cplex.getNnodes()<<"\n";

//	cout << sortc[0] << endl;
//	cout << sortd[0] << endl;
//	cout << mind[0] << endl;
//	cout << maxd[0] << endl;
//	cout << cplex.getValue(cplexd[0]) << endl;

	for (int i=0; i<n; ++i)
        d[nomdindex[i].second] = cplex.getValue(cplexd[i]);

	env.end();

}

void  Select::generate_hard_c_d_1(int _n, int _p, int _gamma, int _R, float _b, int _t)
{
    n = _n;
    p = _p;
    gamma = _gamma;
    R = _R;
    b = _b;
    t = _t;


	if (R == 0)
        generate_RU_dis2();
    else if (R == 1)
        generate_R1_dis2();
    else if (R == 5)
        generate_R5_dis2();


//    if (R == 1)
//        generate_RU_dis1();
//    else if (R == 2)
//        generate_RU_con1();
//    else if (R == 3)
//        generate_RU_dis2();
//    else if (R == 4)
//        generate_RU_con2();
//    else if (R == 5)
//        generate_R1_dis1();
//    else if (R == 6)
//        generate_R1_con1();
//    else if (R == 7)
//        generate_R1_dis2();
//    else if (R == 8)
//        generate_R1_con2();
//    else if (R == 9)
//        generate_R5_dis1();
//    else if (R == 10)
//        generate_R5_con1();
//    else if (R == 11)
//        generate_R5_dis2();
//    else if (R == 12)
//        generate_R5_con2();

    //////////////////////////////////////

    vector<pair<double,int> > nomdindex(n);
    for (int i=0; i<n; ++i)
        nomdindex[i] = make_pair(nomd[i],i);

    sort(nomdindex.begin(), nomdindex.end());

    vector<double> sortc(n);
	for (int i=0; i<n; ++i)
		sortc[i] = nomc[nomdindex[i].second];

    vector<double> sortd(n);
	for (int i=0; i<n; ++i)
		sortd[i] = nomd[nomdindex[i].second];

    vector<double> minc(n);
	for (int i=0; i<n; ++i)
		minc[i] = max (0.0 , sortc[i] - b);

	vector<double> maxc(n);
	for (int i=0; i<n; ++i)
		maxc[i] = min (100.0 , sortc[i] + b);

	vector<double> mind(n);
	for (int i=0; i<n; ++i)
		mind[i] = max (0.0 , sortd[i] - b);

	vector<double> maxd(n);
	for (int i=0; i<n; ++i)
		maxd[i] = min (100.0 , sortd[i] + b);

    int K = n + 1;

    ///////////////////////////////////////

    //solve problem

	IloEnv env;
	IloModel model(env);

    //variables

    vector<IloNumVar> cplexc(n);
	for (int i=0; i<n; ++i)
		cplexc[i] = IloNumVar(env, minc[i], maxc[i], ILOFLOAT);

    vector<IloNumVar> cplexd(n);
	for (int i=0; i<n; ++i)
		cplexd[i] = IloNumVar(env, mind[i], maxd[i], ILOFLOAT);


	vector<IloNumVar> alpha(K);
	for (int k=0; k<K; ++k)
		alpha[k] = IloNumVar(env, 0, IloInfinity, ILOFLOAT);


    vector<vector<IloNumVar> > beta(K);
	for (int k=0; k<K; ++k)
    {
        beta[k].resize(n);
        for (int i=0; i<n; ++i)
            beta[k][i] = IloNumVar(env, 0, IloInfinity, ILOFLOAT);
    }

	IloNumVar cplext(env,0,IloInfinity,ILOFLOAT);

    //constraints

    for (int k=0; k<n; ++k)
    {
        IloExpr con(env);
            con = p * alpha[k];
        for (int i=0; i<n; ++i)
            con += (-beta[k][i]);

        model.add(cplext <= con + (gamma * cplexd[k]));
	}

	for (int k=0; k<n; ++k)
    {
        for (int i=0; i<n; ++i)
        {
            IloExpr con(env);
                con = cplexc[i] + (cplexd[i] - cplexd[k]);
            IloExpr lhs(env);
                lhs = alpha[k] - beta[k][i] ;

            if (k < i)
                model.add(lhs <= con);
            else
                model.add(lhs <= cplexc[i]);
        }
    }

    {
        IloExpr con(env);
            con = p * alpha[0];
        for (int i=0; i<n; ++i)
            con += (-beta [0][i]);

        model.add(cplext <= con);
    }

    for (int i=0; i<n; ++i)
    {
        IloExpr lhs(env);
            lhs = (alpha[0] - beta[0][i]);
        IloExpr con(env);
            con = (cplexc[i] + cplexd[i]);

        model.add(lhs <= con);
    }

//    for (int i=0; i<n; ++i)
//        model.add((sortc[i] - b) <= cplexc[i]);
//
//	for (int i=0; i<n; ++i)
//        model.add(cplexc[i] <= (sortc[i] + b));

//    for (int i=0; i<n; ++i)
//        model.add(((1 - b) * sortc[i]) <= cplexc[i]);
//
//	for (int i=0; i<n; ++i)
//        model.add(cplexc[i] <= ((1 + b) * sortc[i]));

    {
        IloExpr con(env);
        for (int i=0; i<n; ++i)
            con += sortc[i];
        IloExpr lhs(env);
        for (int i=0; i<n; ++i)
            lhs += cplexc[i];

        model.add(lhs <= con);
	}

//    for (int i=0; i<n; ++i)
//        model.add((sortd[i] - b) <= cplexd[i]);
//
//	for (int i=0; i<n; ++i)
//        model.add(cplexd[i] <= (sortd[i] + b));
//
//    for (int i=0; i<n; ++i)
//        model.add(((1 - b) * sortd[i]) <= cplexd[i]);
//
//	for (int i=0; i<n; ++i)
//        model.add(cplexd[i] <= ((1 + b) * sortd[i]));

	{
        IloExpr con(env);
        for (int i=0; i<n; ++i)
            con += sortd[i];
        IloExpr lhs(env);
        for (int i=0; i<n; ++i)
            lhs += cplexd[i];

        model.add(lhs <= con);
	}

	model.add(IloMaximize(env, cplext));

	/////////////////////////////////////////////////

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

	// printing statistics

	cout<<"ALGTIME;"<<(dend-dstart)/CLOCKS_PER_SEC<<"\n";
	cout<<"./TIME;"<<(tend-tstart)/CLOCKS_PER_SEC<<"\n";
	cout<<"NODES;"<<cplex.getNnodes()<<"\n";

//	cout << sortc[0] << endl;
//	cout << sortd[0] << endl;
//	cout << mind[0] << endl;
//	cout << maxd[0] << endl;
//	cout << cplex.getValue(cplexd[0]) << endl;

    for (int i=0; i<n; ++i)
        c[nomdindex[i].second] = cplex.getValue(cplexc[i]);

	for (int i=0; i<n; ++i)
        d[nomdindex[i].second] = cplex.getValue(cplexd[i]);

	env.end();

}

void  Select::generate_hard_c_2(int _n, int _p, int _gamma, int _R, float _b, int _t)
{
    n = _n;
    p = _p;
    gamma = _gamma;
    R = _R;
    b = _b;
    t = _t;


	if (R == 0)
        generate_RU_dis2();
    else if (R == 1)
        generate_R1_dis2();
    else if (R == 5)
        generate_R5_dis2();



//	if (R == 1)
//        generate_RU_dis1();
//    else if (R == 2)
//        generate_RU_con1();
//    else if (R == 3)
//        generate_RU_dis2();
//    else if (R == 4)
//        generate_RU_con2();
//    else if (R == 5)
//        generate_R1_dis1();
//    else if (R == 6)
//        generate_R1_con1();
//    else if (R == 7)
//        generate_R1_dis2();
//    else if (R == 8)
//        generate_R1_con2();
//    else if (R == 9)
//        generate_R5_dis1();
//    else if (R == 10)
//        generate_R5_con1();
//    else if (R == 11)
//        generate_R5_dis2();
//    else if (R == 12)
//        generate_R5_con2();

    //////////////////////////////////////

    int K = n + 1;

    vector<double> pi(K,0);
	for (int k=0; k<(K-1); ++k)
		pi[k] = nomd[k];

//    pi[K-1] = 0;

    vector<double> minc(n);
	for (int i=0; i<n; ++i)
		minc[i] = max (0.0 , ((1 - b) * nomc[i]));

	vector<double> maxc(n);
	for (int i=0; i<n; ++i)
		maxc[i] = min (100.0 , ((1 + b) * nomc[i]));

    //////////////////////////////////////

//    set<double> tempi;
//    tempi.insert(0);
//    for (int j=0; j<n; ++j)
//        tempi.insert(nomd[j]);
//
//
//    vector<double> pi;
//    for (std::set<double>::iterator it=tempi.begin(); it!=tempi.end(); ++it)
//    {
//        pi.push_back(*it);
//    }
//
//    int K = pi.size();

    //////////////////////////////////////

	IloEnv env;
	IloModel model(env);


    vector<IloNumVar> cplexc(n);
	for (int i=0; i<n; ++i)
		cplexc[i] = IloNumVar(env, minc[i], maxc[i], ILOFLOAT);


	vector<IloNumVar> alpha(K);
	for (int k=0; k<K; ++k)
		alpha[k] = IloNumVar(env, 0, IloInfinity, ILOFLOAT);


    vector<vector<IloNumVar> > beta(K);
	for (int k=0; k<K; ++k)
    {
        beta[k].resize(n);
        for (int i=0; i<n; ++i)
            beta[k][i] = IloNumVar(env, 0, IloInfinity, ILOFLOAT);
    }

	IloNumVar cplext(env,0,IloInfinity,ILOFLOAT);

	//constraints

    for (int k=0; k<K; ++k)
    {
        IloExpr con(env);
            con = p * alpha[k];
        for (int i=0; i<n; ++i)
            con += (-beta[k][i]);

        model.add(cplext <= (con + (gamma * pi[k])));
	}


	for (int k=0; k<K; ++k)
    {
        for (int i=0; i<n; ++i)
        {
            IloExpr con(env);
                con = cplexc[i] + max (0.0 , nomd[i] - pi[k]);
            IloExpr lhs(env);
                lhs = alpha[k] - beta[k][i] ;

            model.add(lhs <= con);

        }
    }

//    for (int i=0; i<n; ++i)
//        model.add(((1 - b) * nomc[i]) <= cplexc[i]);
//
//	for (int i=0; i<n; ++i)
//        model.add(cplexc[i] <= ((1 + b) * nomc[i]));
//
//    for (int i=0; i<n; ++i)
//        model.add(((1 - b) * nomc[i]) <= cplexc[i]);
//
//	for (int i=0; i<n; ++i)
//        model.add(cplexc[i] <= ((1 + b) * nomc[i]));

	{
        IloExpr con(env);
        for (int i=0; i<n; ++i)
            con += nomc[i];
        IloExpr lhs(env);
        for (int i=0; i<n; ++i)
            lhs += cplexc[i];

        model.add(lhs <= con);
	}

	model.add(IloMaximize(env, cplext));

	IloCplex cplex(model);

	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::Threads, 1);
	cplex.setParam(IloCplex::PreInd,0);
	cplex.setParam(IloCplex::TiLim,t);

	double dstart = cplex.getDetTime();
	double tstart = clock();

	bool result = cplex.solve();

	double dend = cplex.getDetTime();
	double tend = clock();

	// printing statistics
	cout<<"ALGTIME;"<<(dend-dstart)/CLOCKS_PER_SEC<<"\n";
	cout<<"./TIME;"<<(tend-tstart)/CLOCKS_PER_SEC<<"\n";
	cout<<"NODES;"<<cplex.getNnodes()<<"\n";

	for (int i=0; i<n; ++i)
        c[i] = cplex.getValue(cplexc[i]);

	env.end();

}

void  Select::generate_hard_d_2(int _n, int _p, int _gamma, int _R, float _b, int _t)
{
    n = _n;
    p = _p;
    gamma = _gamma;
    R = _R;
    b = _b;
    t = _t;


	if (R == 0)
        generate_RU_dis2();
    else if (R == 1)
        generate_R1_dis2();
    else if (R == 5)
        generate_R5_dis2();


//    if (R == 1)
//        generate_RU_dis1();
//    else if (R == 2)
//        generate_RU_con1();
//    else if (R == 3)
//        generate_RU_dis2();
//    else if (R == 4)
//        generate_RU_con2();
//    else if (R == 5)
//        generate_R1_dis1();
//    else if (R == 6)
//        generate_R1_con1();
//    else if (R == 7)
//        generate_R1_dis2();
//    else if (R == 8)
//        generate_R1_con2();
//    else if (R == 9)
//        generate_R5_dis1();
//    else if (R == 10)
//        generate_R5_con1();
//    else if (R == 11)
//        generate_R5_dis2();
//    else if (R == 12)
//        generate_R5_con2();

    //////////////////////////////////////

    vector<pair<double,int> > nomdindex(n);
    for (int i=0; i<n; ++i)
        nomdindex[i] = make_pair(nomd[i],i);

    sort(nomdindex.begin(), nomdindex.end());

    vector<double> sortc(n);
	for (int i=0; i<n; ++i)
		sortc[i] = nomc[nomdindex[i].second];

    vector<double> sortd(n);
	for (int i=0; i<n; ++i)
		sortd[i] = nomd[nomdindex[i].second];

	vector<double> mind(n);
	for (int i=0; i<n; ++i)
		mind[i] = max (0.0 , ((1 - b) * sortd[i]));

	vector<double> maxd(n);
	for (int i=0; i<n; ++i)
		maxd[i] = min (100.0 , ((1 + b) * sortd[i]));

    int K = n + 1;

    ///////////////////////////////////////

    //solve problem

	IloEnv env;
	IloModel model(env);

    //variables

    vector<IloNumVar> cplexd(n);
	for (int i=0; i<n; ++i)
		cplexd[i] = IloNumVar(env, mind[i], maxd[i], ILOFLOAT);


	vector<IloNumVar> alpha(K);
	for (int k=0; k<K; ++k)
		alpha[k] = IloNumVar(env, 0, IloInfinity, ILOFLOAT);


    vector<vector<IloNumVar> > beta(K);
	for (int k=0; k<K; ++k)
    {
        beta[k].resize(n);
        for (int i=0; i<n; ++i)
            beta[k][i] = IloNumVar(env, 0, IloInfinity, ILOFLOAT);
    }

	IloNumVar cplext(env,0,IloInfinity,ILOFLOAT);

	//constraints

    for (int k=0; k<n; ++k)
    {
        IloExpr con(env);
            con = p * alpha[k];
        for (int i=0; i<n; ++i)
            con -= beta[k][i];

        model.add(cplext <= con + (gamma * cplexd[k]));
	}

	for (int k=0; k<n; ++k)
    {
        for (int i=0; i<n; ++i)
        {
            IloExpr con(env);
                con = sortc[i] + (cplexd[i] - cplexd[k]);
            IloExpr lhs(env);
                lhs = alpha[k] - beta[k][i] ;

            if ((k) < i)
                model.add(lhs <= con);
            else
                model.add(lhs <= sortc[i]);
        }
    }

    {
        IloExpr con(env);
            con = p * alpha[0];
        for (int i=0; i<n; ++i)
            con -= beta [0][i];

        model.add(cplext <= con);
    }

    for (int i=0; i<n; ++i)
    {
        IloExpr lhs(env);
            lhs = (alpha[0] - beta[0][i]);
        IloExpr con(env);
            con = (sortc[i] + cplexd[i]);

        model.add(lhs <= con);
    }


//    for (int i=0; i<n; ++i)
//        model.add(max(0.0 , sortd[i]-b) <= cplexd[i]);

//    model.add(15 <= cplexd[0]);

//    for (int i=0; i<n; ++i)
//        model.add(max(0.0 , sortd[i]-b) <= cplexd[i]);

//	for (int i=0; i<n; ++i)
//        model.add(cplexd[i] <= min(100.0 , sortd[i]+b));

	{
        IloExpr con(env);
        for (int i=0; i<n; ++i)
            con += nomd[i];
        IloExpr lhs(env);
        for (int i=0; i<n; ++i)
            lhs += cplexd[i];

        model.add(lhs <= con);
	}

	model.add(IloMaximize(env, cplext));

	/////////////////////////////////////////////////

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

	// printing statistics

	cout<<"ALGTIME;"<<(dend-dstart)/CLOCKS_PER_SEC<<"\n";
	cout<<"./TIME;"<<(tend-tstart)/CLOCKS_PER_SEC<<"\n";
	cout<<"NODES;"<<cplex.getNnodes()<<"\n";

//	cout << sortc[0] << endl;
//	cout << sortd[0] << endl;
//	cout << mind[0] << endl;
//	cout << maxd[0] << endl;
//	cout << cplex.getValue(cplexd[0]) << endl;

	for (int i=0; i<n; ++i)
        d[nomdindex[i].second] = cplex.getValue(cplexd[i]);

	env.end();

}

void  Select::generate_hard_c_d_2(int _n, int _p, int _gamma, int _R, float _b, int _t)
{
    n = _n;
    p = _p;
    gamma = _gamma;
    R = _R;
    b = _b;
    t = _t;

	if (R == 0)
        generate_RU_dis2();
    else if (R == 1)
        generate_R1_dis2();
    else if (R == 5)
        generate_R5_dis2();



//    if (R == 1)
//        generate_RU_dis1();
//    else if (R == 2)
//        generate_RU_con1();
//    else if (R == 3)
//        generate_RU_dis2();
//    else if (R == 4)
//        generate_RU_con2();
//    else if (R == 5)
//        generate_R1_dis1();
//    else if (R == 6)
//        generate_R1_con1();
//    else if (R == 7)
//        generate_R1_dis2();
//    else if (R == 8)
//        generate_R1_con2();
//    else if (R == 9)
//        generate_R5_dis1();
//    else if (R == 10)
//        generate_R5_con1();
//    else if (R == 11)
//        generate_R5_dis2();
//    else if (R == 12)
//        generate_R5_con2();

    //////////////////////////////////////

    vector<pair<double,int> > nomdindex(n);
    for (int i=0; i<n; ++i)
        nomdindex[i] = make_pair(nomd[i],i);

    sort(nomdindex.begin(), nomdindex.end());

    vector<double> sortc(n);
	for (int i=0; i<n; ++i)
		sortc[i] = nomc[nomdindex[i].second];

    vector<double> sortd(n);
	for (int i=0; i<n; ++i)
		sortd[i] = nomd[nomdindex[i].second];

    vector<double> minc(n);
	for (int i=0; i<n; ++i)
		minc[i] = max (0.0 , ((1 - b) * sortc[i]));

	vector<double> maxc(n);
	for (int i=0; i<n; ++i)
		maxc[i] = min (100.0 , ((1 + b) * sortc[i]));

	vector<double> mind(n);
	for (int i=0; i<n; ++i)
		mind[i] = max (0.0 , ((1 - b) * sortd[i]));

	vector<double> maxd(n);
	for (int i=0; i<n; ++i)
		maxd[i] = min (100.0 , ((1 + b) * sortd[i]));

    int K = n + 1;

    ///////////////////////////////////////

    //solve problem

	IloEnv env;
	IloModel model(env);

    //variables

    vector<IloNumVar> cplexc(n);
	for (int i=0; i<n; ++i)
		cplexc[i] = IloNumVar(env, minc[i], maxc[i], ILOFLOAT);

    vector<IloNumVar> cplexd(n);
	for (int i=0; i<n; ++i)
		cplexd[i] = IloNumVar(env, mind[i], maxd[i], ILOFLOAT);


	vector<IloNumVar> alpha(K);
	for (int k=0; k<K; ++k)
		alpha[k] = IloNumVar(env, 0, IloInfinity, ILOFLOAT);


    vector<vector<IloNumVar> > beta(K);
	for (int k=0; k<K; ++k)
    {
        beta[k].resize(n);
        for (int i=0; i<n; ++i)
            beta[k][i] = IloNumVar(env, 0, IloInfinity, ILOFLOAT);
    }

	IloNumVar cplext(env,0,IloInfinity,ILOFLOAT);

    //constraints

    for (int k=0; k<n; ++k)
    {
        IloExpr con(env);
            con = p * alpha[k];
        for (int i=0; i<n; ++i)
            con += (-beta[k][i]);

        model.add(cplext <= con + (gamma * cplexd[k]));
	}

	for (int k=0; k<n; ++k)
    {
        for (int i=0; i<n; ++i)
        {
            IloExpr con(env);
                con = cplexc[i] + (cplexd[i] - cplexd[k]);
            IloExpr lhs(env);
                lhs = alpha[k] - beta[k][i] ;

            if (k < i)
                model.add(lhs <= con);
            else
                model.add(lhs <= cplexc[i]);
        }
    }

    {
        IloExpr con(env);
            con = p * alpha[0];
        for (int i=0; i<n; ++i)
            con += (-beta [0][i]);

        model.add(cplext <= con);
    }

    for (int i=0; i<n; ++i)
    {
        IloExpr lhs(env);
            lhs = (alpha[0] - beta[0][i]);
        IloExpr con(env);
            con = (cplexc[i] + cplexd[i]);

        model.add(lhs <= con);
    }

//    for (int i=0; i<n; ++i)
//        model.add((sortc[i] - b) <= cplexc[i]);
//
//	for (int i=0; i<n; ++i)
//        model.add(cplexc[i] <= (sortc[i] + b));

//    for (int i=0; i<n; ++i)
//        model.add(((1 - b) * sortc[i]) <= cplexc[i]);
//
//	for (int i=0; i<n; ++i)
//        model.add(cplexc[i] <= ((1 + b) * sortc[i]));

    {
        IloExpr con(env);
        for (int i=0; i<n; ++i)
            con += sortc[i];
        IloExpr lhs(env);
        for (int i=0; i<n; ++i)
            lhs += cplexc[i];

        model.add(lhs <= con);
	}

//    for (int i=0; i<n; ++i)
//        model.add((sortd[i] - b) <= cplexd[i]);
//
//	for (int i=0; i<n; ++i)
//        model.add(cplexd[i] <= (sortd[i] + b));
//
//    for (int i=0; i<n; ++i)
//        model.add(((1 - b) * sortd[i]) <= cplexd[i]);
//
//	for (int i=0; i<n; ++i)
//        model.add(cplexd[i] <= ((1 + b) * sortd[i]));

	{
        IloExpr con(env);
        for (int i=0; i<n; ++i)
            con += sortd[i];
        IloExpr lhs(env);
        for (int i=0; i<n; ++i)
            lhs += cplexd[i];

        model.add(lhs <= con);
	}

	model.add(IloMaximize(env, cplext));

	/////////////////////////////////////////////////

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

	// printing statistics

	cout<<"ALGTIME;"<<(dend-dstart)/CLOCKS_PER_SEC<<"\n";
	cout<<"./TIME;"<<(tend-tstart)/CLOCKS_PER_SEC<<"\n";
	cout<<"NODES;"<<cplex.getNnodes()<<"\n";

//	cout << sortc[0] << endl;
//	cout << sortd[0] << endl;
//	cout << mind[0] << endl;
//	cout << maxd[0] << endl;
//	cout << cplex.getValue(cplexd[0]) << endl;

    for (int i=0; i<n; ++i)
        c[nomdindex[i].second] = cplex.getValue(cplexc[i]);

	for (int i=0; i<n; ++i)
        d[nomdindex[i].second] = cplex.getValue(cplexd[i]);

	env.end();

}





//
//void  Select::generate_hard_c_new(int _n, int _p, int _gamma, int _R, double _b, int _t)
//{
//    n = _n;
//    p = _p;
//    gamma = _gamma;
//    R = _R;
//    b = _b;
//    t = _t;
//
//    if (R == 0)
//        gen_RU();
//    else if (R == 1)
//        gen_R1();
//    else if (R == 5)
//        gen_R5();
//
//    //////////////////////////////////////
//
////    int K = n + 1;
////
////    vector<double> pi(K);
////	for (int k=0; k<n; ++k)
////		pi[k] = nomd[k];
////
////    pi[K-1] = 0;
//
//
//    set<double> tempi;
//    tempi.insert(0);
//    for (int j=0; j<n; ++j)
//        tempi.insert(nomd[j]);
//
//
//    vector<double> pi;
//    for (std::set<double>::iterator it=tempi.begin(); it!=tempi.end(); ++it)
//    {
//        pi.push_back(*it);
//    }
//
//    int K = pi.size();
//
//    //////////////////////////////////////
//
//	IloEnv env;
//	IloModel model(env);
//
//
//    vector<IloNumVar> cplexc(n);
//	for (int i=0; i<n; ++i)
//		cplexc[i] = IloNumVar(env, 0, IloInfinity, ILOFLOAT);
//
//
//	vector<IloNumVar> alpha(K);
//	for (int k=0; k<K; ++k)
//		alpha[k] = IloNumVar(env, 0, IloInfinity, ILOFLOAT);
//
//
//    vector<vector<IloNumVar> > beta(K);
//	for (int k=0; k<K; ++k)
//    {
//        beta[k].resize(n);
//        for (int i=0; i<n; ++i)
//            beta[k][i] = IloNumVar(env, 0, IloInfinity, ILOFLOAT);
//    }
//
//    vector<vector<IloNumVar> > cplexs(K);
//    for (int k=0; k<K; ++k)
//    {
//        cplexs[k].resize(n);
//        for (int i=0; i<n; ++i)
//            cplexs[k][i] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
//    }
//
//	IloNumVar cplext(env,0,IloInfinity,ILOFLOAT);
//
//	//constraints
//
//    for (int k=0; k<K; ++k)
//    {
//        IloExpr con(env);
//            con = p * alpha[k];
//        for (int i=0; i<n; ++i)
//            con -= (beta[k][i]);
//
//        model.add(cplext <= (con + (gamma * pi[k])));
//	}
//
//
//	for (int k=0; k<K; ++k)
//    {
//        for (int i=0; i<n; ++i)
//        {
//            IloExpr con(env);
//                con = cplexc[i] + cplexs[k][i];
//            IloExpr lhs(env);
//                lhs = alpha[k] - beta[k][i] ;
//
//            model.add(lhs <= con);
//
//        }
//    }
//
//    for (int k=0; k<K; ++k)
//        for (int i=0; i<n; ++i)
//            model.add(cplexs[k][i] <= nomd[i] - pi[k]);
//
//    for (int i=0; i<n; ++i)
//        model.add((nomc[i] - b) <= cplexc[i]);
//
//	for (int i=0; i<n; ++i)
//        model.add(cplexc[i] <= (nomc[i] + b));
//
//	{
//        IloExpr con(env);
//        for (int i=0; i<n; ++i)
//            con += nomc[i];
//        IloExpr lhs(env);
//        for (int i=0; i<n; ++i)
//            lhs += cplexc[i];
//
//        model.add(lhs <= con);
//	}
//
//	model.add(IloMaximize(env, cplext));
//
//
//    IloCplex cplex(model);
//
//	cplex.setOut(env.getNullStream());
//	cplex.setParam(IloCplex::Threads,1);
//	cplex.setParam(IloCplex::PreInd,0);
//    cplex.setParam(IloCplex::TiLim,t);
//
//	double dstart = cplex.getDetTime();
//	double tstart = clock();
//
//	bool result = cplex.solve();
//
//	double dend = cplex.getDetTime();
//	double tend = clock();
//
//
//    cout<<"ALG-TI;"<<(dend-dstart)/CLOCKS_PER_SEC<<"\n";
//	cout<<"TIME;"<<(tend-tstart)/CLOCKS_PER_SEC<<"\n";
//	cout<<"NODES;"<<cplex.getNnodes()<<"\n";
//	cout<<"OBJ;"<<cplex.getValue(cplext)<<"\n";
//
//	for (int i=0; i<n; ++i)
//        c[i] = cplex.getValue(cplexc[i]);
//
//	env.end();
//
//}
//
//void  Select::generate_hard_d_new(int _n, int _p, int _gamma, int _R, double _b, int _t)
//{
//    n = _n;
//    p = _p;
//    gamma = _gamma;
//    R = _R;
//    b = _b;
//    t = _t;
//
//    if (R == 0)
//        gen_RU();
//    else if (R == 1)
//        gen_R1();
//    else if (R == 5)
//        gen_R5();
//
//    ///////////////////////////////////////
//
////    vector<pair<double,int> > nomdindex(n);
////    for (int i=0; i<n; ++i)
////        nomdindex[i] = make_pair(nomd[i],i);
////
////    sort(nomdindex.begin(), nomdindex.end());
////
////    vector<int> sortc(n);
////	for (int i=1; i<n; ++i)
////		sortc[i] = nomc[nomdindex[i].second];
////
////    vector<int> sortd(n);
////	for (int i=1; i<n; ++i)
////		sortd[i] = nomd[nomdindex[i].second];
////
////    int K = n+1;
//
//    set<double> tempi;
//    tempi.insert(0);
//    for (int j=0; j<n; ++j)
//        tempi.insert(nomd[j]);
//
//
//    vector<double> pi;
//    for (std::set<double>::iterator it=tempi.begin(); it!=tempi.end(); ++it)
//    {
//        pi.push_back(*it);
//    }
//
//    int K = pi.size();
//
//
//    ///////////////////////////////////////
//
//	IloEnv env;
//	IloModel model(env);
//
//    vector<IloNumVar> cplexd(n);
//	for (int i=0; i<n; ++i)
//		cplexd[i] = IloNumVar(env, 0, IloInfinity, ILOFLOAT);
//
//	vector<IloNumVar> alpha(K);
//	for (int k=0; k<K; ++k)
//		alpha[k] = IloNumVar(env, 0, IloInfinity, ILOFLOAT);
//
//    vector<vector<IloNumVar> > beta(K);
//	for (int k=0; k<K; ++k)
//    {
//        beta[k].resize(n);
//        for (int i=0; i<n; ++i)
//            beta[k][i] = IloNumVar(env, 0, IloInfinity, ILOFLOAT);
//    }
//
//    vector<vector<IloNumVar> > cplexs(K);
//    for (int k=0; k<K; ++k)
//    {
//        cplexs[k].resize(n);
//        for (int i=0; i<n; ++i)
//            cplexs[k][i] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
//    }
//
//
//	IloNumVar cplext(env,0,IloInfinity,ILOFLOAT);
//
//	//constraints
//
//    for (int k=0; k<K; ++k)
//    {
//        IloExpr con(env);
//            con = p * alpha[k];
//        for (int i=0; i<n; ++i)
//            con -= (beta[k][i]);
//
//        model.add(cplext <= con + (gamma * cplexd[k]));
//	}
//
//	for (int k=0; k<K; ++k)
//    {
//        for (int i=0; i<n; ++i)
//        {
//            IloExpr con(env);
//                con = nomc[i] + cplexs[k][i];
//            IloExpr lhs(env);
//                lhs = alpha[k] - beta[k][i] ;
//
//            model.add(lhs <= con);
//
//        }
//    }
//
//    for (int k=0; k<K; ++k)
//        for (int i=0; i<n; ++i)
//            model.add(cplexs[k][i] <= cplexd[i] - cplexd[k]);
//
////    {
////        IloExpr con(env);
////            con = p * alpha[0];
////        for (int i=0; i<n; ++i)
////            con -= (beta [0][i]);
////
////        model.add(cplext <= con);
////    }
////
////    for (int i=0; i<n; ++i)
////    {
////        IloExpr lhs(env);
////            lhs = (alpha[0] - beta[0][i]);
////        IloExpr con(env);
////            con = (sortc[i] + cplexd[i]);
////
////        model.add(lhs <= con);
////    }
//
//    for (int i=0; i<n; ++i)
//        model.add((nomd[i] - b) <= cplexd[i]);
//
//	for (int i=0; i<n; ++i)
//        model.add(cplexd[i] <= (nomd[i] + b));
//
//	{
//        IloExpr con(env);
//        for (int i=0; i<n; ++i)
//            con += nomd[i];
//        IloExpr lhs(env);
//        for (int i=0; i<n; ++i)
//            lhs += cplexd[i];
//
//        model.add(lhs <= con);
//	}
//
//	model.add(IloMaximize(env, cplext));
//
//    IloCplex cplex(model);
//
//	cplex.setOut(env.getNullStream());
//	cplex.setParam(IloCplex::Threads,1);
//	cplex.setParam(IloCplex::PreInd,0);
//    cplex.setParam(IloCplex::TiLim,t);
//
//	double dstart = cplex.getDetTime();
//	double tstart = clock();
//
//	bool result = cplex.solve();
//
//	double dend = cplex.getDetTime();
//	double tend = clock();
//
//
//    cout<<"ALG-TI;"<<(dend-dstart)/CLOCKS_PER_SEC<<"\n";
//	cout<<"TIME;"<<(tend-tstart)/CLOCKS_PER_SEC<<"\n";
//	cout<<"NODES;"<<cplex.getNnodes()<<"\n";
//	cout<<"OBJ;"<<cplex.getValue(cplext)<<"\n";
//
//	for (int i=0; i<n; ++i)
//        d[i] = cplex.getValue(cplexd[i]);
//
//	env.end();
//
//}
//
//
//
//void  Select::generate_hard_c_d_new(int _n, int _p, int _gamma, int _R, double _b, int _t)
//{
//    n = _n;
//    p = _p;
//    gamma = _gamma;
//    R = _R;
//    b = _b;
//    t = _t;
//
//    if (R == 0)
//        gen_RU();
//    else if (R == 1)
//        gen_R1();
//    else if (R == 5)
//        gen_R5();
//
//    //////////////////////////////////////
//
////    vector<pair<double,int> > nomdindex(n);
////    for (int i=0; i<n; ++i)
////        nomdindex[i] = make_pair(nomd[i],i);
////
////    sort(nomdindex.begin(), nomdindex.end());
////
////    vector<int> sortc(n);
////	for (int i=1; i<n; ++i)
////		sortc[i] = nomc[nomdindex[i].second];
////
////    vector<int> sortd(n);
////	for (int i=1; i<n; ++i)
////		sortd[i] = nomd[nomdindex[i].second];
////
////    int K = n+1;
//
//    set<double> tempi;
//    tempi.insert(0);
//    for (int j=0; j<n; ++j)
//        tempi.insert(nomd[j]);
//
//
//    vector<double> pi;
//    for (std::set<double>::iterator it=tempi.begin(); it!=tempi.end(); ++it)
//    {
//        pi.push_back(*it);
//    }
//
//    int K = pi.size();
//
//
//    ///////////////////////////////////////
//
//	IloEnv env;
//	IloModel model(env);
//
//    vector<IloNumVar> cplexc(n);
//	for (int i=0; i<n; ++i)
//		cplexc[i] = IloNumVar(env, 0, IloInfinity, ILOFLOAT);
//
//    vector<IloNumVar> cplexd(n);
//	for (int i=0; i<n; ++i)
//		cplexd[i] = IloNumVar(env, 0, IloInfinity, ILOFLOAT);
//
//	vector<IloNumVar> alpha(K);
//	for (int k=0; k<K; ++k)
//		alpha[k] = IloNumVar(env, 0, IloInfinity, ILOFLOAT);
//
//    vector<vector<IloNumVar> > beta(K);
//	for (int k=0; k<K; ++k)
//    {
//        beta[k].resize(n);
//        for (int i=0; i<n; ++i)
//            beta[k][i] = IloNumVar(env, 0, IloInfinity, ILOFLOAT);
//    }
//
//    vector<vector<IloNumVar> > cplexs(K);
//    for (int k=0; k<K; ++k)
//    {
//        cplexs[k].resize(n);
//        for (int i=0; i<n; ++i)
//            cplexs[k][i] = IloNumVar(env,0,IloInfinity,ILOFLOAT);
//    }
//
//	IloNumVar cplext(env,0,IloInfinity,ILOFLOAT);
//
//	//constraints
//
//    for (int k=0; k<K; ++k)
//    {
//        IloExpr con(env);
//            con = p * alpha[k];
//        for (int i=0; i<n; ++i)
//            con -= (beta[k][i]);
//
//        model.add(cplext <= con + (gamma * cplexd[k]));
//	}
//
//	for (int k=0; k<K; ++k)
//    {
//        for (int i=0; i<K; ++i)
//        {
//            IloExpr con(env);
//                con = cplexc[i] + cplexs[k][i];
//            IloExpr lhs(env);
//                lhs = alpha[k] - beta[k][i] ;
//
//            model.add(lhs <= con);
//        }
//    }
//
//    for (int k=0; k<K; ++k)
//        for (int i=0; i<n; ++i)
//            model.add(cplexs[k][i] <= cplexd[i] - cplexd[k]);
//
////    {
////        IloExpr con(env);
////            con = p * alpha[0];
////        for (int i=0; i<n; ++i)
////            con -= (beta [0][i]);
////
////        model.add(cplext <= con);
////    }
////
////    for (int i=0; i<n; ++i)
////    {
////        IloExpr lhs(env);
////            lhs = (alpha[0] - beta[0][i]);
////        IloExpr con(env);
////            con = (cplexc[i] + cplexd[i]);
////
////        model.add(lhs <= con);
////    }
//
//    for (int i=0; i<n; ++i)
//        model.add((nomc[i] - b) <= cplexc[i]);
//
//	for (int i=0; i<n; ++i)
//        model.add(cplexc[i] <= (nomc[i] + b));
//
//    {
//        IloExpr con(env);
//        for (int i=0; i<n; ++i)
//            con += nomc[i];
//        IloExpr lhs(env);
//        for (int i=0; i<n; ++i)
//            lhs += cplexc[i];
//
//        model.add(lhs <= con);
//	}
//
//    for (int i=0; i<n; ++i)
//        model.add((nomd[i] - b) <= cplexd[i]);
//
//	for (int i=0; i<n; ++i)
//        model.add(cplexd[i] <= (nomd[i] + b));
//
//	{
//        IloExpr con(env);
//        for (int i=0; i<n; ++i)
//            con += nomd[i];
//        IloExpr lhs(env);
//        for (int i=0; i<n; ++i)
//            lhs += cplexd[i];
//
//        model.add(lhs <= con);
//	}
//
//	model.add(IloMaximize(env, cplext));
//
//    IloCplex cplex(model);
//
//	cplex.setOut(env.getNullStream());
//	cplex.setParam(IloCplex::Threads,1);
//	cplex.setParam(IloCplex::PreInd,0);
//    cplex.setParam(IloCplex::TiLim,t);
//
//	double dstart = cplex.getDetTime();
//	double tstart = clock();
//
//	bool result = cplex.solve();
//
//	double dend = cplex.getDetTime();
//	double tend = clock();
//
//
//    cout<<"ALG-TI;"<<(dend-dstart)/CLOCKS_PER_SEC<<"\n";
//	cout<<"TIME;"<<(tend-tstart)/CLOCKS_PER_SEC<<"\n";
//	cout<<"NODES;"<<cplex.getNnodes()<<"\n";
//	cout<<"OBJ;"<<cplex.getValue(cplext)<<"\n";
//
//	for (int i=0; i<n; ++i)
//	{
//        c[i] = cplex.getValue(cplexc[i]);
//        d[i] = cplex.getValue(cplexd[i]);
//    }
//
//	env.end();
//
//}




vector<double> Select::get_c()
{
	return c;
}

vector<double> Select::get_d()
{
	return d;
}



//double Select::fRand(double fMin, double fMax)
//{
//    double f = (double)rand() / RAND_MAX;
//    return fMin + f * (fMax -fMin);
//}
