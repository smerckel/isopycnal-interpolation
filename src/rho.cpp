#include "rho.hpp"
#include <iostream>
#include <gswteos-10.h>

double Rho::density(double S, double T, double P)
{
    double R;


    //std::cout << "S:" << S  << "T:" << T << "P:" <<P << std::endl;
    R = gsw_pot_rho_t_exact(S, T, P, 0.);
    return R;
}
