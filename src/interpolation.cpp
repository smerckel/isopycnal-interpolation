#include <iostream>

#include "interpolation.hpp"


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


int Interpolation::interpolate_onto_surface(std::vector<double> & s_field, std::vector<double> & z_field, DataNC & data, const double rho)
{
    int return_value = 1; // no successful interpolation
    int interpolation_result=0;
    size_t n=0,k=0, s_offset=0;
    double s=0;

    // Ensure the surface has the right size

    nx = data.get_nx();
    ny = data.get_ny();
    nz = data.get_nz();
    s_offset = s_field.size();
    s_field.resize(s_offset + ny * nx);
    z_field.resize(s_offset + ny * nx);

    std::vector<double> salt, temp, z, O2, mask_rho;
    salt = data.get_salt();
    temp = data.get_temp();
    z = data.get_z();
    O2 = data.get_O2();
    mask_rho = data.get_mask_rho();
    double sigma0=0, sigma1=0; // density values to be computed for k and k+1

    size_t k0=-1, k1 = -1;
    for(size_t j=0; j<ny; ++j)
    {
        for(size_t i=0; i<nx; ++i)
        {
            n = index(j, i); // Gives the index in a 2d field
            if (mask_rho[n] < 0.5) //mask is a double 0. or 1. For some reason 1 means NOT masked...
            {
                s_field[s_offset + n] = NAN; // Or set it to 0?
                continue;
            }
            interpolation_result = interpolate_at_ji(k, sigma0, sigma1, rho, j, i, salt, temp, z);
            if (interpolation_result==0)
            {
                // interpolate for oxygen
                s = interpolate_linear(rho, sigma0, sigma1, O2[index(k, j, i)], O2[index(k+1, j, i)]);
                s_field[s_offset + n] = s;
                // interpolate for z
                s = interpolate_linear(rho, sigma0, sigma1, z[index(k, j, i)], z[index(k+1, j, i)]);
                z_field[s_offset + n] = s;
                return_value = 0; // at least on valid point found
            }
            else
                s_field[s_offset + n] = NAN;
        }
        if ( (k1=j*100/ny/5) != k0)
        {
            std::cout << k1*5 <<  std::endl;
            k0 = k1;
        }
    }
    return return_value;
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

size_t Interpolation::index(const size_t k, const size_t j, const size_t i)
{
    return k*ny*nx + j*nx + i;
}


/*
  Note: Depths are negative, k=0 is at the bottom, k=nz is the surface.

*/


int Interpolation::interpolate_at_ji(size_t &k, double & sigma0, double & sigma1, const double rho, const size_t j, const size_t i,
                                     const std::vector<double> & salt, const std::vector<double> & temp,
                                     const std::vector<double> & z)
{
    int result=0; // all good, until proven otherwise

    size_t index0 = index(k, j, i), index1 = index(k+1, j, i);
    sigma0 = density_calculations.density_from_depth(salt[index0], temp[index0], -z[index0]);
    sigma1 = density_calculations.density_from_depth(salt[index1], temp[index1], -z[index1]);

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
            sigma0 = density_calculations.density_from_depth(salt[index0], temp[index0], -z[index0]);
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
            sigma1 = density_calculations.density_from_depth(salt[index1], temp[index1], -z[index1]);
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

