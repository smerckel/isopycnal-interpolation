#ifndef INTERPOLATION_HPP_INCLUDED
#define INTERPOLATION_HPP_INCLUDED

#include <vector>
#include <cmath>

#include "nc_data.hpp"
#include "rho.hpp"

class Interpolation
{
private:

    Rho    density_calculations;
    size_t nx, ny, nz;

    double interpolate_linear(const double rho, const double sigma0, const double sigma1, const double v0, const double v1);
    int interpolate_at_ji(size_t & k, double & sigma0, double & sigma1, const double rho, const size_t j, const size_t i,
                          const std::vector<double> & salt, const std::vector<double> & temp,
                          const std::vector<double> & z);

    // get the index for 2 and 3 D ravelled arrays.
    size_t index(const size_t j, const size_t i);
    size_t index(const size_t k, const size_t j, const size_t i);
public:

    Interpolation() : density_calculations{}, nx{0}, ny{0}, nz{0}
    {
    };


    //desctructor
    ~Interpolation()
    {
    };

    int interpolate_onto_surface(std::vector<double> & s_field, std::vector<double> & z_field, DataNC & data, const double rho);



};

#endif // INTERPOLATION_HPP_INCLUDED
