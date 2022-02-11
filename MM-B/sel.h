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

		//std::set<double> tempi();
		//std::vector<double> pi;

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

//        void generate_RU_dis1();
//        void generate_RU_con1();
        void generate_RU_dis2();
//        void generate_RU_con2();
//
//        void generate_R1_dis1();
//        void generate_R1_con1();
        void generate_R1_dis2();
//        void generate_R1_con2();
//
//        void generate_R5_dis1();
//        void generate_R5_con1();
        void generate_R5_dis2();
//        void generate_R5_con2();

//		void gen_RU();
//		void gen_R1();
//		void gen_R5();

        bool status;
};

#endif

