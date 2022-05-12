#ifndef _H_SELECT_H
#define _H_SELECT_H

#include <vector>
#include <set>
#include "ilcplex/ilocplex.h"
//#include "ilcplex/ilocplexi.h"

class Select
{
	public:
		void generate_rand(int _n, int _p, int _R);
		void generate_hard_c_d(int _n, int _p, int _R, double _b, int _t);

		std::vector<double> get_c();
		std::vector<double> get_d();

	private:
		int n,p,t,R;
		float b;

		std::set<double> tempi();
		std::vector<double> pi;

		std::vector<double> c;
		std::vector<double> nomc;

        std::vector<double> d;
        std::vector<double> nomd;

		void generate_MMR_I_U();
		void generate_MMR_I_1();
		void generate_MMR_I_2();

        bool status;
};

#endif

