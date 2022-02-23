#ifndef RHO_HPP_INCLUDED
#define RHO_HPP_INCLUDED


class Rho

{
private:

double QR=+999.842594;
double Q01=+6.793952e-2;
double Q02=-9.095290e-3;
double Q03=+1.001685e-4;
double Q04=-1.120083e-6;
double Q05=+6.536332e-9;
double Q10=+0.824493;
double Q11=-4.08990e-3;
double Q12=+7.64380e-5;
double Q13=-8.24670e-7;
double Q14=+5.38750e-9;
double QS0=-5.72466e-3;
double QS1=+1.02270e-4;
double QS2=-1.65460e-6;
double Q20=+4.8314e-4;


public:


    double density(double sp, double tp);

};

#endif // RHO_HPP_INCLUDED
