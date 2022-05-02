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

void Interpolation::compute_avg_between_isopycnals(std::map<std::string, std::vector<double>> & interpolation_variables,
                                                   DataNC & data, const std::vector<double> pycnoclines)
{

    nx = data.get_nx(); // sizes of the rho-variables
    ny = data.get_ny();
    nz = data.get_nz();
    size_t np = pycnoclines.size(); // number of pycnoclines we have.
    // storage for sigma0 and sigma1, which we calculate only once. We'll set it up for all pycnoclines
    std::vector<std::vector<double>> sigma0(np, std::vector<double>(ny*nx)), sigma1(np, std::vector<double>(ny*nx));
    std::vector<std::vector<size_t>> kvec(np, std::vector<size_t>(ny*nx));

    std::vector<size_t> s_offset(interpolation_variables.size());
    std::vector<double> &z = data.get("z");

    for(size_t i=0; i<pycnoclines.size(); ++i)
        resolve_pycnocline_indices(pycnoclines[i], ny,  nx, data, sigma0[i], sigma1[i], kvec[i]);

    // Ensure the average fields have the right size
    size_t m=0;
    for(auto it=interpolation_variables.begin(); it != interpolation_variables.end(); ++it)
    {
        size_t dy = (size_t) (data.get_variable_coordinates(it->first) == data.v_coordinates);
        size_t dx = (size_t) (data.get_variable_coordinates(it->first) == data.u_coordinates);
        s_offset[m++] = (ny-dy)  * (nx-dx); // Size required for one layer
        it->second.resize( (ny-dy)  * (nx-dx) * (np-1) ); // Make space for all layers
    }
    /* compute the fields for the pycnocline depths. */
    std::vector<std::vector<double>> pycnocline_fields_z(np, std::vector<double>(ny*nx));
    for (size_t i=0; i<np; ++i)
    {
      interpolate_on_rho_points(pycnocline_fields_z[i], 0, nx, ny, pycnoclines[i], data.get("z"), kvec[i], sigma0[i], sigma1[i]);
    }

    m=0;
    for(auto it=interpolation_variables.begin(); it != interpolation_variables.end(); ++it)
    {
        std::vector<double> & iv{data.get(it->first)};
        int cv=data.get_variable_coordinates(it->first);
        switch (cv)
        {
            case data.rho_coordinates:
                compute_avg_on_rho_points(it->second, s_offset[m++], nx, ny, pycnoclines, iv, z, kvec,
                                          pycnocline_fields_z, sigma0, sigma1);
                break;
            case data.u_coordinates:
                compute_avg_on_u_points(it->second, s_offset[m++], nx, ny, pycnoclines, iv, z, kvec,
                                        pycnocline_fields_z, sigma0, sigma1);
                break;
            case data.v_coordinates:
                compute_avg_on_v_points(it->second, s_offset[m++], nx, ny, pycnoclines, iv, z, kvec,
                                        pycnocline_fields_z, sigma0, sigma1);
                break;
        }
    }
    std::cout << "Compute average values completed." << std::endl;
}


void Interpolation::interpolate_onto_surface(std::map<std::string, std::vector<double>> & interpolation_variables, DataNC & data,
                                             const size_t i, const double rho)
{
    nx = data.get_nx(); // sizes of the rho-variables
    ny = data.get_ny();
    nz = data.get_nz();
    // storage for sigma0 and sigma1, which we calculate only once.
    std::vector<double> sigma0(ny*nx), sigma1(ny*nx);
    std::vector<size_t> kvec(ny*nx);
    std::vector<size_t> s_offset(interpolation_variables.size());

    resolve_pycnocline_indices(rho, ny,  nx, data, sigma0, sigma1, kvec);


    // Ensure the surfaces have the right size
    size_t m=0;
    for(auto it=interpolation_variables.begin(); it != interpolation_variables.end(); ++it)
    {
        size_t dy = (size_t) (data.get_variable_coordinates(it->first) == data.v_coordinates);
        size_t dx = (size_t) (data.get_variable_coordinates(it->first) == data.u_coordinates);
        s_offset[m++] = i * (ny-dy)  * (nx-dx);
        it->second.resize((i+1)*(ny-dy)  * (nx-dx));
    }

    m=0;
    for(auto it=interpolation_variables.begin(); it != interpolation_variables.end(); ++it)
    {
        std::vector<double> & iv{data.get(it->first)};
        int cv=data.get_variable_coordinates(it->first);
        switch (cv)
        {
            case data.rho_coordinates:
                interpolate_on_rho_points(it->second, s_offset[m++], nx, ny, rho, iv, kvec, sigma0, sigma1);
                break;
            case data.u_coordinates:
                interpolate_on_u_points(it->second, s_offset[m++], nx, ny, rho, iv, kvec, sigma0, sigma1);
                break;
            case data.v_coordinates:
                interpolate_on_v_points(it->second, s_offset[m++], nx, ny, rho, iv, kvec, sigma0, sigma1);
                break;
        }
    }
}


/*
 *
 * P R I V A T E    I N T E R F A C E
 *
 */

 /* resolve_pycnocline_indices

  Finds vectors for kvec, sigma0 and sima1 representing 2D arrays for the whole domain, where

  kvec is the k-index below the pycnocline,
  sigma0 the density level for k
  sigma1 the density level for k+1

  Parameters
  ----------
  rho: density of the pycnocline
  ny,ny: domain size
  data: DataNC object
  sigma0, sigma, kvec, see above.

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

void Interpolation::compute_avg_on_rho_points(std::vector<double> &f,
                                              const size_t offset,
                                              const size_t nx, const size_t ny,
                                              const std::vector<double> &pycnoclines,
                                              const std::vector<double> &iv,
                                              const std::vector<double> &z,
                                              const std::vector<std::vector<size_t>> &kvec,
                                              const std::vector<std::vector<double>> &pycnocline_fields_z,
                                              const std::vector<std::vector<double>> &sigma0,
                                              const std::vector<std::vector<double>> &sigma1)
{
    size_t k0, k1, n;
    size_t np = pycnoclines.size();
    // compute the values of iv matrix at the pycnoclines
    std::vector<std::vector<double>> pycnocline_fields_iv(np, std::vector<double>(ny*nx));
    for (size_t layer=0; layer<np; layer++)
    {
        interpolate_on_rho_points(pycnocline_fields_iv[layer], 0, nx, ny, pycnoclines[layer], iv, kvec[layer], sigma0[layer], sigma1[layer]);
    }

    for (size_t layer=0; layer<np-1; layer++)
    {
        // compute the averages for each layer.
        n=0;
        for(size_t j=0; j<ny; ++j)
            for(size_t i=0; i<nx; ++i)
            {
                k0 = kvec[layer+1][n]; // k index just below the bottom pycnocline
                k1 = kvec[layer][n];       // k index just below the top pycnocline
                if ((k0!=masked) && (k1!=masked))
                {
                    //f[offset*layer + n] = integrate(iv, z, k0, k1, j, i); // without boundary correction
                    f[offset*layer + n] = integrate(iv, z, k0, k1, j, i,
                                                    pycnocline_fields_iv[layer+1], pycnocline_fields_iv[layer],
                                                    pycnocline_fields_z[layer+1], pycnocline_fields_z[layer]);
                }
                else
                    f[offset*layer + n] = NAN;
                n++;
            }
    }
}

void Interpolation::compute_avg_on_u_points(std::vector<double> &f,
                                              const size_t offset,
                                              const size_t nx, const size_t ny,
                                              const std::vector<double> &pycnoclines,
                                              const std::vector<double> &iv,
                                              const std::vector<double> &z,
                                              const std::vector<std::vector<size_t>> &kvec,
                                              const std::vector<std::vector<double>> &pycnocline_fields_z,
                                              const std::vector<std::vector<double>> &sigma0,
                                              const std::vector<std::vector<double>> &sigma1)
{
    size_t k0, k1, n{0}, n_rho_grid, n_rho_grid_right;
    size_t np = pycnoclines.size();

    // compute the values of iv matrix at the pycnoclines
    std::vector<std::vector<double>> pycnocline_fields_iv(np, std::vector<double>(ny*nx));
    for (size_t layer=0; layer<np; layer++)
    {
        interpolate_on_u_points(pycnocline_fields_iv[layer], 0, nx, ny, pycnoclines[layer], iv, kvec[layer], sigma0[layer], sigma1[layer]);
    }
    for (size_t layer=0; layer<np-1; layer++)
    {
        // compute the averages for each layer.
        for(size_t j=0; j<ny; ++j)
            for(size_t i=0; i<nx-1; ++i)
            {
                n_rho_grid = index2(j, i);
                n_rho_grid_right = index2(j, i+1);
                k0 = kvec[layer+1][n_rho_grid];     // k index just below the bottom pycnocline
                k1 = kvec[layer][n_rho_grid];       // k index just below the top pycnocline
                if ((k0!=masked) && (k1!=masked) && (kvec[layer+1][n_rho_grid_right]!=masked)
                    && (kvec[layer][n_rho_grid_right]!=masked))
                {
                    f[offset*layer + n] = integrate(iv, z, k0, k1, j, i,
                                                    pycnocline_fields_iv[layer+1], pycnocline_fields_iv[layer],
                                                    pycnocline_fields_z[layer+1], pycnocline_fields_z[layer], 0, 1);
                    //f[offset*layer + n] = integrate(iv, z, k0, k1, j, i, 0, 1); // No so precise
                }
                else
                    f[offset*layer + n] = NAN;
                n++;
            }
    }
}

void Interpolation::compute_avg_on_v_points(std::vector<double> &f,
                                              const size_t offset,
                                              const size_t nx, const size_t ny,
                                              const std::vector<double> &pycnoclines,
                                              const std::vector<double> &iv,
                                              const std::vector<double> &z,
                                              const std::vector<std::vector<size_t>> &kvec,
                                              const std::vector<std::vector<double>> &pycnocline_fields_z,
                                              const std::vector<std::vector<double>> &sigma0,
                                              const std::vector<std::vector<double>> &sigma1)
{
    size_t k0, k1, n{0}, n_rho_grid_up;
    size_t np = pycnoclines.size();
   // compute the values of iv matrix at the pycnoclines
    std::vector<std::vector<double>> pycnocline_fields_iv(np, std::vector<double>(ny*nx));
    for (size_t layer=0; layer<np; layer++)
    {
        interpolate_on_v_points(pycnocline_fields_iv[layer], 0, nx, ny, pycnoclines[layer], iv, kvec[layer], sigma0[layer], sigma1[layer]);
    }
    for (size_t layer=0; layer<np-1; layer++)
    {
        // compute the averages for each layer.
        for(size_t j=0; j<ny-1; ++j)
            for(size_t i=0; i<nx; ++i)
            {
                // Note that we can use n here, as it is the same as index2(j,i), because the x-dimension is the same.
                // This is in contrast with the method for u velocities.
                k0 = kvec[layer+1][n];     // k index just below the bottom pycnocline
                k1 = kvec[layer][n];       // k index just below the top pycnocline
                n_rho_grid_up = index2(j+1, i);
                if ((k0!=masked) && (k1!=masked) && (kvec[layer+1][n_rho_grid_up]!=masked)
                    && (kvec[layer][n_rho_grid_up]!=masked))
                {
                    f[offset*layer + n] = integrate(iv, z, k0, k1, j, i,
                                                    pycnocline_fields_iv[layer+1], pycnocline_fields_iv[layer],
                                                    pycnocline_fields_z[layer+1], pycnocline_fields_z[layer], 1, 0);

                }
                else
                    f[offset*layer + n] = NAN + 0*pycnocline_fields_z[0][0] +sigma0[layer][0]+sigma1[layer][0]*0; //keeps compiler happy for now.
                n++;
            }
    }
}

void Interpolation::interpolate_on_rho_points(std::vector<double> &f, const size_t offset,
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
                s = interpolate_linear(rho, sigma0[n], sigma1[n], iv[index3(k, j, i)], iv[index3(k+1, j, i)]);
                f[offset + n] = s;
            }
            else
                f[offset + n] = NAN;
            n++;
        }
}

void Interpolation::interpolate_on_u_points(std::vector<double> &f, const size_t offset,
                                              const size_t nx, const size_t ny, const double rho,
                                              const std::vector<double> &iv,
                                              const std::vector<size_t> &kvec,
                                              const std::vector<double> &sigma0,
                                              const std::vector<double> &sigma1)
{
    size_t k, k_right, n, n_right, m=offset;
    double s;
    for(size_t j=0; j<ny; ++j)
        for(size_t i=0; i<nx-1; ++i)
        {
            n = index2(j,i);
            n_right = n + 1;
            k = kvec[n];
            k_right = kvec[n_right];
            if ((k!=masked) && k_right!=masked)
            {
                s = interpolate_linear(rho,
                                       0.5*(sigma0[n]+sigma0[n_right]),
                                       0.5*(sigma1[n]+sigma1[n_right]),
                                       iv[index3(k, j, i, 0, 1)], iv[index3(k+1, j, i, 0, 1)]);
                f[m] = s;
            }
            else
                f[m] = NAN;
            m++;
        }
}

void Interpolation::interpolate_on_v_points(std::vector<double> &f, const size_t offset,
                                              const size_t nx, const size_t ny, const double rho,
                                              const std::vector<double> &iv,
                                              const std::vector<size_t> &kvec,
                                              const std::vector<double> &sigma0,
                                              const std::vector<double> &sigma1)
{
    size_t k, k_up, n, n_up, m=offset;
    double s;
    for(size_t j=0; j<ny-1; ++j)
        for(size_t i=0; i<nx; ++i)
        {
            n = index2(j,i);
            n_up = index2(j+1, i);
            k = kvec[n];
            k_up = kvec[n_up];
            if ((k!=masked) && k_up!=masked)
            {
                s = interpolate_linear(rho,
                                       0.5*(sigma0[n]+sigma0[n_up]),
                                       0.5*(sigma1[n]+sigma1[n_up]),
                                       iv[index3(k, j, i, 1, 0)], iv[index3(k+1, j, i, 1, 0)]);
                f[m] = s;
            }
            else
                f[m] = NAN;
            m++;
        }
}

size_t Interpolation::index2(const size_t j, const size_t i)
{
    return j*nx + i;
}

size_t Interpolation::index2(const size_t j, const size_t i, const size_t dx)
{
    return j*(nx-dx) + i;
}

size_t Interpolation::index3(const size_t k, const size_t j, const size_t i)
{
    return k*ny*nx + j*nx + i;
}

size_t Interpolation::index3(const size_t k, const size_t j, const size_t i, const size_t dy, const size_t dx)
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
    size_t index0 = index3(k, j, i), index1 = index3(k+1, j, i);
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
            index0 = index3(k, j, i);
            index1 = index3(k+1, j, i);
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
            index0 = index3(k, j, i);
            index1 = index3(k+1, j, i);
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

double Interpolation::integrate(const std::vector<double> &iv, const std::vector<double> &z,
                                const size_t k0, const size_t k1, const size_t j, const size_t i,
                                const size_t dy, const size_t dx)

/*

k==0 : bottom
k==nz-1 : surface
k0 : index just below the bottom pycnocline
k1 : index just below the top pycnocline

*/

{
    double sum{0}, sum_z{0};
    double dz{0};
    size_t idx0, idx1; // index numbers for the variable to integrate.
    idx0 = index3(k0+1, j, i, dy, dx); // because we want to start just above the bottom pycnocline, we also need to
                                       // account for the grid used.
    for(size_t k=k0+1; k<k1; ++k)
    {
        idx1 = index3(k+1, j, i, dy , dx);
        if ((dx == 0) && (dy==0))
        {
            // rho-grid
            dz = z[idx1] - z[idx0];
        }
        else if (dx==1)
        {
            //u - grid: evaluate dz at u point.
            dz = 0.5*(z[index3(k+1, j, i)] + z[index3(k+1, j, i+1)] - z[index3(k, j, i)] - z[index3(k, j, i+1)]);
        }
        else
        {
            //v - grid: evaluate dz at v point.
            dz = 0.5*(z[index3(k+1, j, i)] + z[index3(k+1, j+1, i)] - z[index3(k, j, i)] - z[index3(k, j+1, i)]);
        }
        sum_z += dz;
        sum += (iv[idx0] + iv[idx1])*0.5*dz;
        idx0 = idx1;
    }
    if (sum_z>0)
        sum /= sum_z;
    else
    {
        std::cout<< "Current position in grid: i: " << i << "  j: " << j << std::endl;
        throw std::runtime_error("Integration distance < 0!!!\nPossible cause: pycnoclines entered are not monotonously increasing.");
    }
    return sum;
}

double Interpolation::integrate(const std::vector<double> &iv, const std::vector<double> &z,
                                const size_t k0, const size_t k1, const size_t j, const size_t i,
                                const std::vector<double> &pycnocline_fields_iv0,
                                const std::vector<double> &pycnocline_fields_iv1,
                                const std::vector<double> &pycnocline_fields_z0,
                                const std::vector<double> &pycnocline_fields_z1,
                                const size_t dy, const size_t dx)

/*

k==0 : bottom
k==nz-1 : surface
k0 : index just below the bottom pycnocline
k1 : index just below the top pycnocline

*/

{
    double sum{0}, sum_z{0};
    double dz{0};
    size_t idx0, idx1, idx2; // index numbers for the variable to integrate.
    idx0 = index3(k0+1, j, i, dy, dx); // because we want to start just above the bottom pycnocline, we also need to
                                       // account for the grid used.
    idx2 = index2(j, i, dx);
    // The contribution from pycnocline to index k0+1 (just above bottom pycnocline)
    // Here we need to evaluate z and pycnocline_fields_z0 on a u and v grids if necessary.
    if ((dx==0) && (dy==0))
    {
        // rho-coordinates
        dz = z[idx0] - pycnocline_fields_z0[idx2];
    }
    else if (dx==1)
    {
        // u-coordinates
        dz = 0.5* (z[index3(k0+1, j, i)] + z[index3(k0+1, j, i+1)]
                   - pycnocline_fields_z0[index2(j,i)] - pycnocline_fields_z0[index2(j,i+1)]) ;
    }
    else
    {
        // v-coordinates
        dz = 0.5* (z[index3(k0+1, j, i)] + z[index3(k0+1, j+1, i)]
                   - pycnocline_fields_z0[index2(j,i)] - pycnocline_fields_z0[index2(j+1,i)]) ;
    }
    sum += (iv[idx0] + pycnocline_fields_iv0[idx2]) * 0.5 * dz;
    sum_z += dz;

    // Contribution between the pycnoclines:
    for(size_t k=k0+1; k<k1; ++k)
    {
        idx1 = index3(k+1, j, i, dy , dx);
        if ((dx == 0) && (dy==0))
        {
            // rho-grid
            dz = z[idx1] - z[idx0];
        }
        else if (dx==1)
        {
            //u - grid: evaluate dz at u point.
            dz = 0.5*(z[index3(k+1, j, i)] + z[index3(k+1, j, i+1)] - z[index3(k, j, i)] - z[index3(k, j, i+1)]);
        }
        else
        {
            //v - grid: evaluate dz at v point.
            dz = 0.5*(z[index3(k+1, j, i)] + z[index3(k+1, j+1, i)] - z[index3(k, j, i)] - z[index3(k, j+1, i)]);
        }
        sum_z += dz;
        sum += (iv[idx0] + iv[idx1])*0.5*dz;
        idx0 = idx1;
    }
    // The contribution k1 -- upper pycnocline
    if ((dx==0) && (dy==0))
    {
        // rho-coordinates
        dz = pycnocline_fields_z1[idx2] - z[idx1];
    }
    else if (dx==1)
    {
        // u-coordinates
        dz = -0.5* (z[index3(k1, j, i)] + z[index3(k1, j, i+1)]
                    - pycnocline_fields_z1[index2(j,i)] - pycnocline_fields_z1[index2(j,i+1)]);
    }
    else
    {
        // v-coordinates
        dz = -0.5* (z[index3(k1, j, i)] + z[index3(k1, j+1, i)]
                    - pycnocline_fields_z1[index2(j,i)] - pycnocline_fields_z1[index2(j+1,i)]) ;
    }
    sum_z += dz;
    sum += (iv[idx1] + pycnocline_fields_iv1[idx2]) * 0.5 * dz;

    if (sum_z>0)
        sum /= sum_z;
    else
    {
        std::cout<< "Current position in grid: i: " << i << "  j: " << j << std::endl;
        throw std::runtime_error("Integration distance < 0!!!\nPossible cause: pycnoclines entered are not monotonously increasing.");
    }
    return sum;
}
