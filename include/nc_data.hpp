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

    // private variables
    netCDF::NcFile datafile;
    std::map<std::string, std::vector<double>> variables_dict;
    std::map<std::string, std::string> units_dict;

    size_t nt, nx, ny, nz;

    // private methods
    int compute_z_levels();
    int set_dimensions();

public:
    DataNC() : datafile(), variables_dict{}, units_dict{}, nt{0}, nx{0}, ny{0}, nz{0}
    {
    }

    int open(const std::string filename);

    int get_data_const(const std::string var_name);
    int get_data(const std::string var_name, const size_t time_level);
    int get_data(const size_t time_level);
    int get_data();
    int set_unit(const std::string var_name);
    std::string & get_unit(const std::string var_name);

    size_t & get_nt();
    size_t & get_nx();
    size_t & get_ny();
    size_t & get_nz();
    std::vector<double> & get(std::string var_name);
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
