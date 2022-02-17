#ifndef NC_DATA_HPP_INCLUDED
#define NC_DATA_HPP_INCLUDED

#include <iostream>
#include <netcdf>
#include <string>
#include <vector>
#include <map>

class DataNC
{
private:

    const int DIM1=1;
    const int DIM2=2;
    const int DIM3=3;
    const int DIM4=4;
    const int MASK_RHO_MATRIX=0;
    const int O2_MATRIX=1;
    const int SALT_MATRIX=2;
    const int TEMP_MATRIX=3;
    const int H_MATRIX=4;
    const int ZETA_MATRIX=5;
    const int S_RHO_MATRIX=6;

    // private variables
    netCDF::NcFile datafile;
    std::vector<double> mask_rho, O2, salt, temp, h, zeta, s_rho, z, eta_rho, xi_rho;
    size_t nt, nx, ny, nz;

    // private methods
    int compute_z_levels();
    int set_dimensions();

public:
    DataNC() : datafile(), mask_rho{}, O2{}, salt{}, temp{}, h{}, zeta{}, s_rho{}, z{}, eta_rho{}, xi_rho{}, nt{0}, nx{0}, ny{0}, nz{0}
    {
    }

    int open(const std::string filename);

    int get_data_const(const std::string var_name);
    int get_data(const std::string var_name, const size_t time_level);
    int get_data(const size_t time_level);
    int get_data();

    size_t & get_nt();
    size_t & get_nx();
    size_t & get_ny();
    size_t & get_nz();
    std::vector<double> & get_mask_rho();
    std::vector<double> & get_z();
    std::vector<double> & get_O2();
    std::vector<double> & get_salt();
    std::vector<double> & get_temp();
    std::vector<double> & get_xi_rho();
    std::vector<double> & get_eta_rho();
};

class PycnoNC
{
    private:

    //private attributes
    netCDF::NcFile datafile;

    std::vector<double> time, iso_pycnal, eta_rho, xi_rho, O2;

    size_t nt, nz, ny, nx;

    std::map<std::string, netCDF::NcVar> surface_variables;
    std::vector<netCDF::NcDim> dimVector;

    //private methods
    void write_vector_variable(netCDF::NcVar ncv, const std::vector<double> & v);

    public:

    PycnoNC() : datafile(), time{}, iso_pycnal{}, eta_rho{}, xi_rho{}, O2{}, nt{0}, nz{0}, ny{0}, nx{0}, surface_variables{}, dimVector{}
    {
    }


    void open(const std::string filename);
    int create_dimensions(const std::vector<double> pycnoclines,
                          const std::vector<double> eta_rho,
                          const std::vector<double> xi_rho);
    void create_surface_variable(std::string variable_name, std::string units);

    int write_parameter(const std::vector<double> s, const std::string variable_name, size_t rec);
};

#endif // NC_DATA_HPP_INCLUDED
