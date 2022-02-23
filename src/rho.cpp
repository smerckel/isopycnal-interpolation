#include "rho.hpp"
#include <iostream>
#include <cmath>


double Rho::density(double sp, double pt)
{
    double sqrt_sp=sqrt(sp);

    double rho = QR + pt*(Q01+pt*(Q02+pt*(Q03+pt*(Q04+pt*Q05)))) + sp*(Q10+pt*(Q11+pt*(Q12+pt*(Q13+pt*Q14)))+sqrt_sp*(QS0+pt*(QS1+pt*QS2))+sp*Q20);

    return rho;
}
