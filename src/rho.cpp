#include "rho.hpp"
#include <iostream>
#include <gswteos-10.h>

double Rho::densityTEOS10(double sp, double pt, double p)
{


    double sa = gsw_sa_from_sp(sp, p, -31.0, 16.0);
    double ct = gsw_ct_from_pt(sa, pt);
    double R = gsw_sigma0(sa, ct) + 1000.0;
    return R;
}

double Rho::density(double sp, double pt)
{
    double sqrt_sp=sqrt(sp);

    double rho = QR + pt*(Q01+pt*(Q02+pt*(Q03+pt*(Q04+pt*Q05)))) + sp*(Q10+pt*(Q11+pt*(Q12+pt*(Q13+pt*Q14)))+sqrt_sp*(QS0+pt*(QS1+pt*QS2))+sp*Q20);

    return rho;
}
