#ifndef NC_DATA_HPP_INCLUDED
#define NC_DATA_HPP_INCLUDED

#include <iostream>
#include <netcdf>
#include <string>
#include <vector>
#include <map>
#include "ndarray.hpp"


class DataNC
{
private:

    enum {DIM0, DIM1, DIM2, DIM3, DIM4};
/*    const int DIM1=1;
    const int DIM2=2;
    const int DIM3=3;
    const int DIM4=4;
*/

    // private variables
    netCDF::NcFile datafile;
    std::map<std::string, NDarray<double>> variables_dict;
    std::map<std::string, std::string> units_dict;
    std::map<std::string, int> variable_coordinates_dict;
    std::vector<std::string> interpolation_variables;

    size_t nt, nx, ny, nz;
    // variables that are to be read from NC file to compute z from s_rho
    double hc;
    std::vector<double> sc_r, Cs_r;

    // private methods
    int compute_z_levels();
    int set_dimensions();

    int get_attributes(const std::string att_name, double &f);
    int get_attributes(const std::string att_name, std::vector<double> &f);

public:
    enum {rho_coordinates, u_coordinates, v_coordinates};

    DataNC() : datafile(), variables_dict{}, units_dict{}, variable_coordinates_dict{}, \
               interpolation_variables{}, nt{0}, nx{0}, ny{0}, nz{0}, \
               hc{0}, sc_r{}, Cs_r{}
    {
    }

    int open(const std::string filename, const std::vector<std::string>);

    int get_data_const(const std::string var_name);
    int get_data(const std::string var_name, const size_t time_level);
    int get_data(const size_t time_level);
    int get_data();
    int set_unit(const std::string var_name);
    std::string & get_unit(const std::string var_name);

    int get_variable_coordinates(const std::string var_name);
    void set_variable_coordinates(const std::string var_name);

    size_t & get_nt();
    size_t & get_nx();
    size_t & get_ny();
    size_t & get_nz();
    NDarray<double> & get(std::string var_name);
};

class PycnoNC
{
    private:
    enum {rho_coordinates, u_coordinates, v_coordinates};

    //private attributes
    netCDF::NcFile datafile;

    std::vector<double> time, iso_pycnal, eta_rho, xi_rho;

    size_t nt, nz, nlayers, ny, nx;

    std::map<std::string, netCDF::NcVar> surface_variables;
    std::vector<netCDF::NcDim> dimVector_rho, dimVector_u, dimVector_v;

    //private methods
    void write_vector_variable(netCDF::NcVar ncv, const NDarray<double> & v);
    void write_vector_variable(netCDF::NcVar ncv, const NDarray<size_t> & v);
    void write_vector_variable(netCDF::NcVar ncv, const std::vector<double> & v);
    void write_vector_variable(netCDF::NcVar ncv, const std::vector<size_t> & v);


    public:

    PycnoNC() : datafile(), time{}, iso_pycnal{}, eta_rho{}, xi_rho{}, nt{0}, nz{0}, nlayers{0}, ny{0}, nx{0},
                surface_variables{}, dimVector_rho{}, dimVector_u{}, dimVector_v{}
    {
    }


    void open(const std::string filename);
    int create_dimensions(const std::vector<double> & pycnoclines,
                          const NDarray<double> & eta_rho,
                          const NDarray<double> & xi_rho,
                          const NDarray<double> & eta_v,
                          const NDarray<double> & xi_u,
                          const bool compute_averages);

    void create_surface_variable(const std::string variable_name, const std::string units,
                                 const int variable_coordinates);

    int write_parameter(const NDarray<double> s, const std::string variable_name,
                        const int variable_coordinates, const size_t rec);
};

#endif // NC_DATA_HPP_INCLUDED
