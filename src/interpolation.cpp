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


void Interpolation::compute_avg_between_isopycnals(std::map<std::string, NDarray<double>> & interpolation_variables,
                                                   DataNC & data, const std::vector<double> pycnoclines)
{

    nx = data.get_nx(); // sizes of the rho-variables
    ny = data.get_ny();
    nz = data.get_nz();
    size_t np = pycnoclines.size(); // number of pycnoclines we have.
    // storage for sigma0 and sigma1, which we calculate only once. We'll set it up for all pycnoclines
    std::vector<NDarray<double>> sigma0(np, NDarray<double>(ny, nx));
    std::vector<NDarray<double>> sigma1(np, NDarray<double>(ny, nx));
    std::vector<NDarray<size_t>> kvec(np, NDarray<size_t>(ny, nx));
    NDarray<double> &z = data.get("z");
    for(size_t i=0; i<pycnoclines.size(); ++i)
        resolve_pycnocline_indices(pycnoclines[i], ny,  nx, data, sigma0[i], sigma1[i], kvec[i]);

    /* compute the fields for the pycnocline depths. */
    NDarray<double> pycnocline_fields_z(np, ny, nx);
    for (size_t i=0; i<np; ++i)
    {
      interpolate_on_rho_points(pycnocline_fields_z, i, nx, ny, pycnoclines[i], z, kvec[i], sigma0[i], sigma1[i]);
    }

    for(auto it=interpolation_variables.begin(); it != interpolation_variables.end(); ++it)
    {
        NDarray<double> & v{data.get(it->first)};
        int cv=data.get_variable_coordinates(it->first);
        switch (cv)
        {
            case data.rho_coordinates:
                compute_avg_on_rho_points(it->second, nx, ny, pycnoclines, v, z, kvec,
                                          pycnocline_fields_z, sigma0, sigma1);
                break;
            case data.u_coordinates:
                compute_avg_on_u_points(it->second, nx, ny, pycnoclines, v, z, kvec,
                                        pycnocline_fields_z, sigma0, sigma1);
                break;
            case data.v_coordinates:
                compute_avg_on_v_points(it->second, nx, ny, pycnoclines, v, z, kvec,
                                        pycnocline_fields_z, sigma0, sigma1);

                break;
        }
    }
}


void Interpolation::interpolate_onto_surface(std::map<std::string, NDarray<double>> & interpolation_variables, DataNC & data,
                                             const size_t i_pycnocline, const double rho)
{
    nx = data.get_nx(); // sizes of the rho-variables
    ny = data.get_ny();
    nz = data.get_nz();
    // storage for sigma0 and sigma1, which we calculate only once.
    NDarray<double> sigma0(ny, nx), sigma1(ny, nx);
    NDarray<size_t> kvec(ny, nx);

    resolve_pycnocline_indices(rho, ny,  nx, data, sigma0, sigma1, kvec);

    for(auto it=interpolation_variables.begin(); it != interpolation_variables.end(); ++it)
    {
        NDarray<double> & v{data.get(it->first)};
        int cv=data.get_variable_coordinates(it->first);
        switch (cv)
        {
            case data.rho_coordinates:
                interpolate_on_rho_points(it->second, i_pycnocline, nx, ny, rho, v, kvec, sigma0, sigma1);
                break;
            case data.u_coordinates:
                interpolate_on_u_points(it->second, i_pycnocline, nx, ny, rho, v, kvec, sigma0, sigma1);
                break;
            case data.v_coordinates:
                interpolate_on_v_points(it->second, i_pycnocline, nx, ny, rho, v, kvec, sigma0, sigma1);
                break;
        }
    }
}


///*
// *
// * P R I V A T E    I N T E R F A C E
// *
// */
//
// /* resolve_pycnocline_indices
//
//  Finds vectors for kvec, sigma0 and sima1 representing 2D arrays for the whole domain, where
//
//  kvec is the k-index below the pycnocline,
//  sigma0 the density level for k
//  sigma1 the density level for k+1
//
//  Parameters
//  ----------
//  rho: density of the pycnocline
//  ny,ny: domain size
//  data: DataNC object
//  sigma0, sigma, kvec, see above.
//
//*/
void Interpolation::resolve_pycnocline_indices(const double rho, const size_t ny, const size_t nx,
                                               DataNC &data,
                                               NDarray<double> &sigma0,
                                               NDarray<double> &sigma1,
                                               NDarray<size_t> &kvec)
{
    int interpolation_result=0;
    size_t n=0, k=0;
    // read the prerequisites...
    NDarray<double> & salt{data.get("salt")};
    NDarray<double> & temp{data.get("temp")};
    NDarray<double> & mask_rho{data.get("mask_rho")};
    double _sigma0, _sigma1;

    for(size_t j=0; j<ny; ++j)
        for(size_t i=0; i<nx; ++i)
        {
            if (mask_rho(n) < 0.5) //mask is a double 0. or 1. For some reason 1 means NOT masked...
                kvec(n)=masked;
            else
            {
                interpolation_result = interpolate_at_ji(k, _sigma0, _sigma1, rho, j, i, salt, temp);
                if (interpolation_result==0)
                {
                    sigma0(n) = _sigma0;
                    sigma1(n) = _sigma1;
                    kvec(n) = k;
                }
                else
                    kvec(n)=masked;
            }
            n++;
        }
}

void Interpolation::compute_avg_on_rho_points(NDarray<double> &f,
                                              const size_t nx, const size_t ny,
                                              const std::vector<double> &pycnoclines,
                                              const NDarray<double> &v,
                                              const NDarray<double> &z,
                                              const std::vector<NDarray<size_t>> &kvec,
                                              const NDarray<double> &pycnocline_fields_z,
                                              const std::vector<NDarray<double>> &sigma0,
                                              const std::vector<NDarray<double>> &sigma1)
{
    size_t k0, k1;
    size_t np = pycnoclines.size();
    // compute the values of iv matrix at the pycnoclines
    NDarray<double> pycnocline_fields_v(np, ny, nx);
    for (size_t layer=0; layer<np; layer++)
    {
        interpolate_on_rho_points(pycnocline_fields_v, layer, nx, ny, pycnoclines[layer],
                                  v, kvec[layer], sigma0[layer], sigma1[layer]);
    }
    for (size_t layer=0; layer<np-1; layer++)
    {
        // compute the averages for each layer.
        for(size_t j=0; j<ny; ++j)
            for(size_t i=0; i<nx; ++i)
            {
                k0 = kvec[layer+1](j,i); // k index just below the bottom pycnocline
                k1 = kvec[layer](j,i);       // k index just below the top pycnocline

                if ((k0!=masked) && (k1!=masked))
                {
                    f(layer, j, i) = integrate(v, z, k0, k1, j, i, layer,
                                               pycnocline_fields_v, pycnocline_fields_z);
                }
                else
                    f(layer, j, i) = NAN;
            }
    }
}


void Interpolation::compute_avg_on_u_points(NDarray<double> &f,
                                              const size_t nx, const size_t ny,
                                              const std::vector<double> &pycnoclines,
                                              const NDarray<double> &v,
                                              const NDarray<double> &z,
                                              const std::vector<NDarray<size_t>> &kvec,
                                              const NDarray<double> &pycnocline_fields_z,
                                              const std::vector<NDarray<double>> &sigma0,
                                              const std::vector<NDarray<double>> &sigma1)
{
    size_t k0, k1, kp0, kp1;
    size_t np = pycnoclines.size();
    // compute the values of iv matrix at the pycnoclines
    NDarray<double> pycnocline_fields_v(np, ny, nx-1);
    NDarray<double> pycnocline_fields_z_shifted(np, ny, nx-1);
    for (size_t layer=0; layer<np; layer++)
    {
        interpolate_on_u_points(pycnocline_fields_v, layer, nx, ny, pycnoclines[layer],
                                v, kvec[layer], sigma0[layer], sigma1[layer]);
        // interpolated pycnocline_fields z to velocity points.
        for(size_t j=0; j<ny; ++j)
            for(size_t i=0; i<nx-1; ++i)
            {
                pycnocline_fields_z_shifted(layer, j, i) = 0.5 * (pycnocline_fields_z(layer, j, i) + pycnocline_fields_z(layer, j, i+1));
            }
    }
    for (size_t layer=0; layer<np-1; layer++)
    {
        // compute the averages for each layer.
        for(size_t j=0; j<ny; ++j)
            for(size_t i=0; i<nx-1; ++i)
            {
                k0 = kvec[layer+1](j,i); // k index just below the bottom pycnocline
                k1 = kvec[layer](j,i);       // k index just below the top pycnocline
                kp0 = kvec[layer+1](j,i+1);
                kp1 = kvec[layer](j,i+1);

                if ((k0!=masked) && (k1!=masked) && (kp0!=masked) && (kp1 != masked))
                {
                    f(layer, j, i) = integrate(v, z, k0, k1, j, i, layer,
                                               pycnocline_fields_v, pycnocline_fields_z_shifted);
                }
                else
                    f(layer, j, i) = NAN;
            }
    }
}

void Interpolation::compute_avg_on_v_points(NDarray<double> &f,
                                              const size_t nx, const size_t ny,
                                              const std::vector<double> &pycnoclines,
                                              const NDarray<double> &v,
                                              const NDarray<double> &z,
                                              const std::vector<NDarray<size_t>> &kvec,
                                              const NDarray<double> &pycnocline_fields_z,
                                              const std::vector<NDarray<double>> &sigma0,
                                              const std::vector<NDarray<double>> &sigma1)
{
    size_t k0, k1, kp0, kp1;
    size_t np = pycnoclines.size();
    // compute the values of iv matrix at the pycnoclines
    NDarray<double> pycnocline_fields_v(np, ny-1, nx);
    NDarray<double> pycnocline_fields_z_shifted(np, ny-1, nx);
    for (size_t layer=0; layer<np; layer++)
    {
        interpolate_on_v_points(pycnocline_fields_v, layer, nx, ny, pycnoclines[layer],
                                v, kvec[layer], sigma0[layer], sigma1[layer]);
        // interpolated pycnocline_fields z to velocity points.
        for(size_t j=0; j<ny-1; ++j)
            for(size_t i=0; i<nx; ++i)
            {
                pycnocline_fields_z_shifted(layer, j, i) = 0.5 * (pycnocline_fields_z(layer, j, i) + pycnocline_fields_z(layer, j+1, i));
            }
    }
    for (size_t layer=0; layer<np-1; layer++)
    {
        // compute the averages for each layer.
        for(size_t j=0; j<ny-1; ++j)
            for(size_t i=0; i<nx; ++i)
            {
                k0 = kvec[layer+1](j,i); // k index just below the bottom pycnocline
                k1 = kvec[layer](j,i);       // k index just below the top pycnocline
                kp0 = kvec[layer+1](j+1,i);
                kp1 = kvec[layer](j+1,i);

                if ((k0!=masked) && (k1!=masked) && (kp0!=masked) && (kp1 != masked))
                {
                    f(layer, j, i) = integrate(v, z, k0, k1, j, i, layer,
                                               pycnocline_fields_v, pycnocline_fields_z_shifted);
                }
                else
                    f(layer, j, i) = NAN;
            }
    }
}


void Interpolation::interpolate_on_rho_points(NDarray<double> &f, const size_t i_pycnocline,
                                              const size_t nx, const size_t ny, const double rho,
                                              const NDarray<double> &v,
                                              const NDarray<size_t> &kvec,
                                              const NDarray<double> &sigma0,
                                              const NDarray<double> &sigma1)
{
    size_t k;
    double s;
    for(size_t j=0; j<ny; ++j)
        for(size_t i=0; i<nx; ++i)
        {
            k = kvec(j,i);
            if (k!=masked)
            {
                s = interpolate_linear(rho, sigma0(j,i), sigma1(j,i), v(k, j, i), v(k+1, j, i));
                f(i_pycnocline, j, i) = s;
            }
            else
                f(i_pycnocline, j, i) = NAN;
        }
}

void Interpolation::interpolate_on_u_points(NDarray<double> &f, const size_t i_pycnocline,
                                            const size_t nx, const size_t ny, const double rho,
                                            const NDarray<double> &v,
                                            const NDarray<size_t> &kvec,
                                            const NDarray<double> &sigma0,
                                            const NDarray<double> &sigma1)
{
    size_t k, kp;
    double s;
    for(size_t j=0; j<ny; ++j)
        for(size_t i=0; i<nx-1; ++i)
        {
            k = kvec(j,i);
            kp = kvec(j,i+1);
            if (k!=masked && kp!=masked)
            {
                s = interpolate_linear(rho,
                                       (sigma0(j,i) + sigma0(j, i+1))*0.5,
                                       (sigma1(j,i) + sigma1(j, i+1))*0.5, v(k, j, i), v(k+1, j, i));
                f(i_pycnocline, j, i) = s;
            }
            else
                f(i_pycnocline, j, i) = NAN;
        }
}

void Interpolation::interpolate_on_v_points(NDarray<double> &f, const size_t i_pycnocline,
                                            const size_t nx, const size_t ny, const double rho,
                                            const NDarray<double> &v,
                                            const NDarray<size_t> &kvec,
                                            const NDarray<double> &sigma0,
                                            const NDarray<double> &sigma1)
{
    size_t k, kp;
    double s;
    for(size_t j=0; j<ny-1; ++j)
        for(size_t i=0; i<nx; ++i)
        {
            k = kvec(j,i);
            kp = kvec(j+1,i);
            if (k!=masked && kp!=masked)
            {
                s = interpolate_linear(rho,
                                       (sigma0(j,i) + sigma0(j+1, i))*0.5,
                                       (sigma1(j,i) + sigma1(j+1, i))*0.5, v(k, j, i), v(k+1, j, i));
                f(i_pycnocline, j, i) = s;
            }
            else
                f(i_pycnocline, j, i) = NAN;
        }
}

/*
  Note: Depths are negative, k=0 is at the bottom, k=nz is the surface.

*/


int Interpolation::interpolate_at_ji(size_t& k, double& sigma0, double& sigma1, const double rho,
                                     const size_t j, const size_t i,
                                     const NDarray<double>& salt, const NDarray<double>& temp)
{
    int result=0; // all good, until proven otherwise
    sigma0 = density_calculations.density(salt(k,j,i), temp(k,j,i));
    sigma1 = density_calculations.density(salt(k+1,j,i), temp(k+1,j,i));

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
            sigma1=sigma0; //reuse previously computed densities.
            sigma0 = density_calculations.density(salt(k,j,i), temp(k,j,i));
        }
        else if (sigma1 > rho)
        {
            if (k==nz-2) // -2 because we need k+1 as well.
            {
                result = 2; // We are asked to increase k beyond the bound. So we won't find a value either.
                break;
            }
            k++;
            sigma0=sigma1; //reuse previously computed densities.
            sigma1 = density_calculations.density(salt(k+1, j, i), temp(k+1, j, i));
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

double Interpolation::integrate(const NDarray<double> &v, const NDarray<double> &z,
                                const size_t k0, const size_t k1, const size_t j, const size_t i,
                                const size_t layer,
                                const NDarray<double> &pycnocline_fields_v,
                                const NDarray<double> &pycnocline_fields_z)

/*

k==0 : bottom
k==nz-1 : surface
k0 : index just below the bottom pycnocline
k1 : index just below the top pycnocline

*/

{
    double sum{0};
    double sum_z{0};
    double dz;
    double f;
    for (size_t k = k0+1; k<k1; ++k)
    {
        dz = z(k+1, j, i) - z(k, j, i);
        f = 0.5 * (v(k+1, j, i) + v(k, j, i));
        sum += f*dz;
        sum_z += dz;
    }
    // corrections at the boundaries
    // Note that because of the isopycnal values are increasing, we have
    // thtat k0 corresponds to layer+1
    // and k1 to layer.
    dz = z(k0+1, j, i) - pycnocline_fields_z(layer+1, j,i);
    f = 0.5 * (v(k0+1, j, i) + pycnocline_fields_v(layer+1, j,i));
    sum += f*dz;
    sum_z += dz;

    dz = -z(k1, j, i) + pycnocline_fields_z(layer, j, i);
    f = 0.5 * (v(k1, j, i) + pycnocline_fields_v(layer, j, i));
    sum += f*dz;
    sum_z += dz;

    // sanity check
    if (sum_z <= 0)
    {
        std::cout << "Sum of dz's <=0 : "<< sum_z << std::endl;
        throw std::runtime_error("Negative dz");
    }
    sum /= sum_z;
    return sum;
}
