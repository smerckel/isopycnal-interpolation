#ifndef RHO_HPP_INCLUDED
#define RHO_HPP_INCLUDED


class Rho

{
private:

public:

    double density(double S, double T, double P);

    double density_from_depth(double S, double T, double z)
    {
        return density(S, T, pressure_from_depth(z));
    }


    double pressure_from_depth(double z)
    {
        return z*1025*9.81*1e-4; //Pressure in dbar.
    }


};

#endif // RHO_HPP_INCLUDED
