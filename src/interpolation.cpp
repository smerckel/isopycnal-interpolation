#include <iostream>

#include "interpolation.hpp"
#include<functional>

/* Interpolation class
 *
 *
 *
 *
 */

 /*
 *
 * P U B L I C    I N T E R F A C E
 *
 */

 /*
    interpolate_onto_surface

    Parameters
    ----------
    a dictionary (map) of vectors representing the surfaces of the variables to interpolate (output)
    DataNC: data (input)
    const double rho (input) the density of the pycnocline

    This methods loops over all variables in the dictionary that are to be interpolated.

 */

void Interpolation::resolve_pycnocline_indices(const double rho, const size_t ny, const size_t nx,
                                               DataNC &data,
                                               std::vector<double> &sigma0,
                                               std::vector<double> &sigma1,
                                               std::vector<size_t> &kvec)
{
    int interpolation_result=0;
    size_t n=0, k=0;
    // read the prerequisites...
    std::vector<double> & salt{data.get("salt")};
    std::vector<double> & temp{data.get("temp")};
    std::vector<double> & mask_rho{data.get("mask_rho")};
    double _sigma0, _sigma1;

    for(size_t j=0; j<ny; ++j)
        for(size_t i=0; i<nx; ++i)
        {
            if (mask_rho[n] < 0.5) //mask is a double 0. or 1. For some reason 1 means NOT masked...
                kvec[n]=masked;
            else
            {
                interpolation_result = interpolate_at_ji(k, _sigma0, _sigma1, rho, j, i, salt, temp);
                if (interpolation_result==0)
                {
                    sigma0[n] = _sigma0;
                    sigma1[n] = _sigma1;
                    kvec[n] = k;
                }
                else
                    kvec[n]=masked;
            }
            n++;
        }
}

void Interpolation::interpolate_on_rho_points(std::vector<double> &f,
                                              const size_t nx, const size_t ny, const double rho,
                                              const std::vector<double> &iv,
                                              const std::vector<size_t> &kvec,
                                              const std::vector<double> &sigma0,
                                              const std::vector<double> &sigma1)
{
    size_t k, n{0};
    double s;
    for(size_t j=0; j<ny; ++j)
        for(size_t i=0; i<nx; ++i)
        {
            k = kvec[n];
            if (k!=masked)
            {
                s = interpolate_linear(rho, sigma0[n], sigma1[n], iv[index(k, j, i)], iv[index(k+1, j, i)]);
                f[n] = s;
            }
            n++;
        }
}

void Interpolation::interpolate_on_u_points(std::vector<double> &f,
                                              const size_t nx, const size_t ny, const double rho,
                                              const std::vector<double> &iv,
                                              const std::vector<size_t> &kvec,
                                              const std::vector<double> &sigma0,
                                              const std::vector<double> &sigma1)
{
    size_t k, k_right, n, n_right, m=0;
    double s;
    for(size_t j=0; j<ny; ++j)
        for(size_t i=0; i<nx-1; ++i)
        {
            n = index(j,i);
            n_right = n + 1;
            k = kvec[n];
            k_right = kvec[n_right];
            if ((k!=masked) && k_right!=masked)
            {
                s = interpolate_linear(rho,
                                       0.5*(sigma0[n]+sigma0[n_right]),
                                       0.5*(sigma1[n]+sigma1[n_right]),
                                       iv[index(k, j, i, 0, 1)], iv[index(k+1, j, i, 0, 1)]);
                f[m] = s;
            }
            else
                f[m] = NAN;
            m++;
        }
}

void Interpolation::interpolate_on_v_points(std::vector<double> &f,
                                              const size_t nx, const size_t ny, const double rho,
                                              const std::vector<double> &iv,
                                              const std::vector<size_t> &kvec,
                                              const std::vector<double> &sigma0,
                                              const std::vector<double> &sigma1)
{
    size_t k, k_up, n, n_up, m=0;
    double s;
    for(size_t j=0; j<ny-1; ++j)
        for(size_t i=0; i<nx; ++i)
        {
            n = index(j,i);
            n_up = index(j+1, i);
            k = kvec[n];
            k_up = kvec[n_up];
            if ((k!=masked) && k_up!=masked)
            {
                s = interpolate_linear(rho,
                                       0.5*(sigma0[n]+sigma0[n_up]),
                                       0.5*(sigma1[n]+sigma1[n_up]),
                                       iv[index(k, j, i, 1, 0)], iv[index(k+1, j, i, 1, 0)]);
                f[m] = s;
            }
            else
                f[m] = NAN;
            m++;
        }
}

void Interpolation::interpolate_onto_surface(std::map<std::string, std::vector<double>> & interpolation_variables, DataNC & data, const double rho)
{
    nx = data.get_nx(); // sizes of the rho-variables
    ny = data.get_ny();
    nz = data.get_nz();
    // storage for sigma0 and sigma1, which we calculate only once.
    std::vector<double> sigma0(ny*nx), sigma1(ny*nx);
    std::vector<size_t> kvec(ny*nx);

    resolve_pycnocline_indices(rho, ny,  nx, data, sigma0, sigma1, kvec);


   // Ensure the surfaces have the right size
    for(auto it=interpolation_variables.begin(); it != interpolation_variables.end(); ++it)
    {
        size_t dy = (size_t) (data.get_variable_coordinates(it->first) == data.v_coordinates);
        size_t dx = (size_t) (data.get_variable_coordinates(it->first) == data.u_coordinates);
        it->second.resize((ny-dy)  * (nx-dx));
    }

    for(auto it=interpolation_variables.begin(); it != interpolation_variables.end(); ++it)
    {
        std::vector<double> & iv{data.get(it->first)};
        int cv=data.get_variable_coordinates(it->first);
        switch (cv)
        {
            case data.rho_coordinates:
                interpolate_on_rho_points(it->second, nx, ny, rho, iv, kvec, sigma0, sigma1);
                break;
            case data.u_coordinates:
                interpolate_on_u_points(it->second, nx, ny, rho, iv, kvec, sigma0, sigma1);
                break;
            case data.v_coordinates:
                interpolate_on_v_points(it->second, nx, ny, rho, iv, kvec, sigma0, sigma1);
                break;
        }
    }
}


/*
 *
 * P R I V A T E    I N T E R F A C E
 *
 */


size_t Interpolation::index(const size_t j, const size_t i)
{
    return j*nx + i;
}

size_t Interpolation::index(const size_t j, const size_t i, const size_t dy, const size_t dx)
{
    return j*(nx-dx) + i;
}

size_t Interpolation::index(const size_t k, const size_t j, const size_t i)
{
    return k*ny*nx + j*nx + i;
}

size_t Interpolation::index(const size_t k, const size_t j, const size_t i, const size_t dy, const size_t dx)
{
    return k*(ny-dy)*(nx-dx) + j*(nx-dx) + i;
}

/*
  Note: Depths are negative, k=0 is at the bottom, k=nz is the surface.

*/


int Interpolation::interpolate_at_ji(size_t& k, double& sigma0, double& sigma1, const double rho,
                                     const size_t j, const size_t i,
                                     const std::vector<double>& salt, const std::vector<double>& temp)
{
    int result=0; // all good, until proven otherwise
    if((j==499) && (i==51))
        result=0;
    size_t index0 = index(k, j, i), index1 = index(k+1, j, i);
    sigma0 = density_calculations.density(salt[index0], temp[index0]);
    sigma1 = density_calculations.density(salt[index1], temp[index1]);

    while(1)
    {
        if ((rho<=sigma0) && (rho > sigma1))
            break;
        else if (sigma0 < rho)
        {
            if(k==0)
            {
                result = 1; // We are asked to decrease k, but we can't. So we won't find a value :-(
                break;
            }
            k--;
            index0 = index(k, j, i);
            index1 = index(k+1, j, i);
            sigma1=sigma0; //reuse previously computed densities.
            sigma0 = density_calculations.density(salt[index0], temp[index0]);
        }
        else if (sigma1 > rho)
        {
            if (k==nz-2) // -2 because we need k+1 as well.
            {
                result = 2; // We are asked to increase k beyond the bound. So we won't find a value either.
                break;
            }
            k++;
            index0 = index(k, j, i);
            index1 = index(k+1, j, i);
            sigma0=sigma1; //reuse previously computed densities.
            sigma1 = density_calculations.density(salt[index1], temp[index1]);
        }
        else
        {
            std::cout << "k     : " << k << std::endl;
            std::cout << "sigma0: " << sigma0 << std::endl;
            std::cout << "sigma1: " << sigma1 << std::endl;
            std::cout << "rho   : " << rho << std::endl;
            throw std::runtime_error("Unhandled condition.");
        }
    }
    return result;
}



double Interpolation::interpolate_linear(const double rho, const double sigma0, const double sigma1, const double v0, const double v1)
{
    double ds, dv, vi;
    ds = sigma1 - sigma0;
    dv = v1 - v0;
    vi = v0 + dv/ds * (rho-sigma0);
    return vi;
}

