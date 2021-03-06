#ifndef _H_SELECT_H
#define _H_SELECT_H

#include <vector>
#include "ilcplex/ilocplex.h"


struct Solution
{
	std::vector<double> x;
	double ub;
	int nodes;
};

class Select
{
	public:
		void generate_rand(int _n, int _p, int _N, int _R);
		void generate_hard_c(int _n, int _p, int _N, int _R);

		Solution solve_ip(double tlim = -1);

		void set_budget(int _scenbudget);
		std::vector<std::vector<double> > get_c();
		//std::vector<double> get_C();
		void set_timelimit(double _timelimit);
		void print();

	private:
		int n,p,N,R;
		std::vector<std::vector<double> > c;
		std::vector<std::vector<double> > nomc;

        //std::vector<double> C;
        //std::vector<double> nomC;

		std::vector<std::vector<double> > x;

		//bounds on c
		std::vector<std::vector<double> > cl;
		std::vector<std::vector<double> > cu;

        //std::vector<double> Cl;
        //std::vector<double> Cu;

		double solve_master();
        //double solve_master_C_and_c();

		void generate_MMR_D_U();
		void generate_MMR_D_1();
		void generate_MMR_D_2();

		int scenbudget;

        bool status;
};

#endif
