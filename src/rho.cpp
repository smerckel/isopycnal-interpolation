#include "rho.hpp"
#include <iostream>
#include <gswteos-10.h>

double Rho::density(double sp, double pt, double p)
{


    double sa = gsw_sa_from_sp(sp, p, -31.0, 16.0);
    double ct = gsw_ct_from_pt(sa, pt);
    double R = gsw_sigma0(sa, ct) + 1000.0;
    return R;
}
