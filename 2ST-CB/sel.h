#ifndef _SEL_H_
#define _SEL_H_

#include <vector>

class Sel
{
    public:
        void generate_2ST_CB_U(int _n, int _p, int _gamma);
        void generate_2ST_CB_1(int _n, int _p, int _gamma);
        void generate_2ST_CB_2(int _n, int _p, int _gamma);

        std::vector<int> get_C();
        std::vector<int> get_c();
        std::vector<int> get_d();

    private:

        int n, p;

        int gamma;

        std::vector<int> C;
        std::vector<int> c;
        std::vector<int> d;
};

#endif
