#ifndef INTERPOLATION_HPP_INCLUDED
#define INTERPOLATION_HPP_INCLUDED

#include <vector>
#include <cmath>

#include "ndarray.hpp"
#include "nc_data.hpp"
#include "rho.hpp"

class Interpolation
{
private:
    Rho    density_calculations;
    size_t nx, ny, nz;
    const size_t masked=99999;

    double interpolate_linear(const double rho, const double sigma0, const double sigma1, const double v0, const double v1);

    double integrate(const NDarray<double> &v, const NDarray<double> &z, const size_t k0, const size_t k1, const size_t j, const size_t i,
                     const size_t layer, const NDarray<double> &pycnocline_fields_v, const NDarray<double> &pycnocline_fields_z);

    int interpolate_at_ji(size_t & k, double & sigma0, double & sigma1, const double rho, const size_t j, const size_t i,
                          const NDarray<double> & salt, const NDarray<double> & temp);

    void interpolate_on_rho_points(NDarray<double> &f, const size_t offset,
                                   const size_t nx, const size_t ny, const double rho,
                                   const NDarray<double> &v,
                                   const NDarray<size_t> &kvec,
                                   const NDarray<double> &sigma0,
                                   const NDarray<double> &sigma1);

    void interpolate_on_u_points(NDarray<double> &f, const size_t offset,
                                   const size_t nx, const size_t ny, const double rho,
                                   const NDarray<double> &v,
                                   const NDarray<size_t> &kvec,
                                   const NDarray<double> &sigma0,
                                   const NDarray<double> &sigma1);

    void interpolate_on_v_points(NDarray<double> &f, const size_t offset,
                                   const size_t nx, const size_t ny, const double rho,
                                   const NDarray<double> &iv,
                                   const NDarray<size_t> &kvec,
                                   const NDarray<double> &sigma0,
                                   const NDarray<double> &sigma1);

    void compute_avg_on_rho_points(NDarray<double> &f,
                                   const size_t nx, const size_t ny, const std::vector<double> &pycnoclines,
                                   const NDarray<double> &v,
                                   const NDarray<double> &z,
                                   const std::vector<NDarray<size_t>> &kvec,
                                   const NDarray<double> &pycnocline_fields_z,
                                   const std::vector<NDarray<double>> &sigma0,
                                   const std::vector<NDarray<double>> &sigma1);

    void compute_avg_on_u_points(NDarray<double> &f,
                                   const size_t nx, const size_t ny, const std::vector<double> &pycnoclines,
                                   const NDarray<double> &v,
                                   const NDarray<double> &z,
                                   const std::vector<NDarray<size_t>> &kvec,
                                   const NDarray<double> &pycnocline_fields_z,
                                   const std::vector<NDarray<double>> &sigma0,
                                   const std::vector<NDarray<double>> &sigma1);

    void compute_avg_on_v_points(NDarray<double> &f,
                                   const size_t nx, const size_t ny, const std::vector<double> &pycnoclines,
                                   const NDarray<double> &v,
                                   const NDarray<double> &z,
                                   const std::vector<NDarray<size_t>> &kvec,
                                   const NDarray<double> &pycnocline_fields_z,
                                   const std::vector<NDarray<double>> &sigma0,
                                   const std::vector<NDarray<double>> &sigma1);

    void resolve_pycnocline_indices(const double rho, const size_t ny, const size_t nx,
                                    DataNC &data,
                                    NDarray<double> &sigma0,
                                    NDarray<double> &sigma1,
                                    NDarray<size_t> &kvec);

public:

    Interpolation() : density_calculations{}, nx{0}, ny{0}, nz{0}
    {
    };


    //desctructor
    ~Interpolation()
    {
    };

    void interpolate_onto_surface(std::map<std::string, NDarray<double>> & interpolation_variables, DataNC & data, const size_t i, const double rho);

    /* computes average values between two isopycnals with density rho0 and rho1 */
    void compute_avg_between_isopycnals(std::map<std::string, NDarray<double>> & interpolation_variables,
                                        DataNC & data,
                                        const std::vector<double> pycnoclines);


};

#endif // INTERPOLATION_HPP_INCLUDED
