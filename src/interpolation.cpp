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

int Interpolation::interpolate_onto_surface(std::map<std::string, std::vector<double>> & interpolation_variables, DataNC & data, const double rho)
{
    int return_value = 1; // no successful interpolation
    int interpolation_result=0;
    size_t n=0,k=0, s_offset=0;
    double s=0;

    // Ensure the surfaces have the right size

    nx = data.get_nx();
    ny = data.get_ny();
    nz = data.get_nz();
    for(auto it=interpolation_variables.begin(); it != interpolation_variables.end(); ++it)
    {
        s_offset = it->second.size();
        it->second.resize(s_offset + ny * nx);
    }
    // read the prerequisites...
    std::vector<double> & salt{data.get("salt")};
    std::vector<double> & temp{data.get("temp")};
    std::vector<double> & mask_rho{data.get("mask_rho")};

    double sigma0=0, sigma1=0; // density values to be computed for k and k+1

    for(size_t j=0; j<ny; ++j)
    {
        for(size_t i=0; i<nx; ++i)
        {
            n = index(j, i); // Gives the index in a 2d field
            if (mask_rho[n] < 0.5) //mask is a double 0. or 1. For some reason 1 means NOT masked...
            {
                for(auto it=interpolation_variables.begin(); it != interpolation_variables.end(); ++it)
                    it->second[s_offset + n] = NAN;
                continue;
            }
            interpolation_result = interpolate_at_ji(k, sigma0, sigma1, rho, j, i, salt, temp, data.get("z"));
            if (interpolation_result==0)
            {
                for(auto it=interpolation_variables.begin(); it != interpolation_variables.end(); ++it)
                {
                    std::vector<double> & iv{data.get(it->first)};
                    s = interpolate_linear(rho, sigma0, sigma1, iv[index(k, j, i)], iv[index(k+1, j, i)]);
                    it->second[s_offset + n] = s;
                }

                if ((i==705) && (j==350) && 0==1)
                {
                    std::cout << "rho: " << rho << std::endl;
                    std::cout << "k: " << k << std::endl;
                    std::cout << "z: " << data.get("z")[s_offset + n] << std::endl;
                    std::cout << "nz: " << nz << std::endl;
                    for(size_t kk=0; kk<nz;kk++)
                    {

                        Rho rho_fun;
                        double rho = rho_fun.density(salt[index(kk, j, i)], temp[index(kk, j, i)], -data.get("z")[index(kk, j, i)]);
                        std::cout << "rho(z): " << data.get("z")[index(kk, j, i)]<< " : " << " : "<< rho << " : " << salt[index(kk, j, i)] << " : " << temp[index(kk, j, i)] << " ";
                        std::cout << kk << ":" << salt[index(kk, j, i)] << " " << temp[index(kk, j, i)] << "  " << data.get("z")[index(kk, j, i)] << std::endl;
                    }
                }
                return_value = 0; // at least on valid point found
            }
            else
                for(auto it=interpolation_variables.begin(); it != interpolation_variables.end(); ++it)
                        it->second[s_offset + n] = NAN;
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

