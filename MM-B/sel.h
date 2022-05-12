#ifndef _H_SELECT_H
#define _H_SELECT_H

#include <vector>
#include <set>
#include "ilcplex/ilocplex.h"
//#include "ilcplex/ilocplexi.h"

class Select
{
	public:
		void generate_rand(int _n, int _p, int _gamma, int _R);
		void generate_hard_c_1(int _n, int _p, int _gamma, int _R, float _b, int _t);
		void generate_hard_d_1(int _n, int _p, int _gamma, int _R, float _b, int _t);
		void generate_hard_c_d_1(int _n, int _p, int _gamma, int _R, float _b, int _t);
		void generate_hard_c_2(int _n, int _p, int _gamma, int _R, float _b, int _t);
		void generate_hard_d_2(int _n, int _p, int _gamma, int _R, float _b, int _t);
		void generate_hard_c_d_2(int _n, int _p, int _gamma, int _R, float _b, int _t);

		double fRand(double fMin, double fMax);

		std::vector<double> get_c();
		std::vector<double> get_d();

	private:
		int n,p;
		int gamma;
		int t,R,K;
		float b;

		std::vector<double> c;
		std::vector<double> nomc;

        std::vector<double> d;
        std::vector<double> nomd;

        std::vector<double> minc;
		std::vector<double> maxc;

        std::vector<double> mind;
		std::vector<double> maxd;

        std::vector<double> pi;

        std::set<double> tempi;

        void generate_MM_B_U();
        void generate_MM_B_1();
        void generate_MM_B_2();

        bool status;
};

#endif

