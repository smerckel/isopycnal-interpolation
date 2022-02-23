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
    const size_t masked=99999;

    double interpolate_linear(const double rho, const double sigma0, const double sigma1, const double v0, const double v1);
    int interpolate_at_ji(size_t & k, double & sigma0, double & sigma1, const double rho, const size_t j, const size_t i,
                          const std::vector<double> & salt, const std::vector<double> & temp);

    void interpolate_on_rho_points(std::vector<double> &f, const size_t offset,
                                   const size_t nx, const size_t ny, const double rho,
                                   const std::vector<double> &iv,
                                   const std::vector<size_t> &kvec,
                                   const std::vector<double> &sigma0,
                                   const std::vector<double> &sigma1);

    void interpolate_on_u_points(std::vector<double> &f, const size_t offset,
                                   const size_t nx, const size_t ny, const double rho,
                                   const std::vector<double> &iv,
                                   const std::vector<size_t> &kvec,
                                   const std::vector<double> &sigma0,
                                   const std::vector<double> &sigma1);

    void interpolate_on_v_points(std::vector<double> &f, const size_t offset,
                                   const size_t nx, const size_t ny, const double rho,
                                   const std::vector<double> &iv,
                                   const std::vector<size_t> &kvec,
                                   const std::vector<double> &sigma0,
                                   const std::vector<double> &sigma1);

    void resolve_pycnocline_indices(const double rho, const size_t ny, const size_t nx,
                                    DataNC &data,
                                    std::vector<double> &sigma0,
                                    std::vector<double> &sigma1,
                                    std::vector<size_t> &kvec);

    // get the index for 2 and 3 D ravelled arrays.
    size_t index2(const size_t j, const size_t i);
    size_t index2(const size_t j, const size_t i, const size_t dx);
    size_t index3(const size_t k, const size_t j, const size_t i);
    size_t index3(const size_t k, const size_t j, const size_t i, const size_t dy, const size_t dx);

public:

    Interpolation() : density_calculations{}, nx{0}, ny{0}, nz{0}
    {
    };


    //desctructor
    ~Interpolation()
    {
    };
    void interpolate_onto_surface(std::map<std::string, std::vector<double>> & interpolation_variables, DataNC & data, const size_t i, const double rho);



};

#endif // INTERPOLATION_HPP_INCLUDED
