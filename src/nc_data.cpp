#include "nc_data.hpp"


int DataNC::open(const std::string filename)
{
    datafile.open(filename, netCDF::NcFile::read);
    set_dimensions();
return 0;
}

int DataNC::set_unit(const std::string var_name)
{
    netCDF::NcVar nc_var = datafile.getVar(var_name);
    int return_value=0;
    try
    {
        netCDF::NcVarAtt nc_att = nc_var.getAtt("units");
        std::string unit;
        nc_att.getValues(unit);
        units_dict[var_name] = unit;
    }
    catch (netCDF::exceptions::NcException & e)
    {
        return_value=-1;
    }
    return return_value;
}

std::string & DataNC::get_unit(const std::string var_name)
{
    return units_dict[var_name]; // no error checking...
}

int DataNC::get_data_const(const std::string var_name)
{
    netCDF::NcVar nc_var = datafile.getVar(var_name);
    std::vector<netCDF::NcDim> dims = nc_var.getDims();
    size_t n = 0;
    if (dims.size()==(size_t) DIM2)
    {
        size_t ny, nx;
        ny = dims[0].getSize();
        nx = dims[1].getSize();
        n = ny * nx;
    }
    else if (dims.size() == (size_t) DIM1)
    {
        n = dims[0].getSize();
    }
    else
        throw std::runtime_error("Unhandled dimension.");

    double *pvalues = new double [n];
    nc_var.getVar(pvalues);
    std::vector<double> v(n);
    for(size_t j=0; j<n; j++)
        v[j] = pvalues[j];
    variables_dict[var_name] = v;
    delete [] pvalues;
    return 0;
}

int DataNC::get_data(const std::string var_name, const size_t time_level)
{
    netCDF::NcVar nc_var = datafile.getVar(var_name);
    std::vector<netCDF::NcDim> dims = nc_var.getDims();
    std::vector<size_t> startp,countp;
    size_t nt = 0, n = 0;
    size_t ny=1, nx=1, nz=1;

    startp.push_back(time_level);
    countp.push_back(1);

    if (dims.size()==(size_t) DIM3)
    {
        nt = dims[0].getSize();
        ny = dims[1].getSize();
        nx = dims[2].getSize();
        startp.push_back(0);
        startp.push_back(0);
        countp.push_back(ny);
        countp.push_back(nx);
    }
    else if (dims.size() == (size_t) DIM4)
    {
        nt = dims[0].getSize();
        nz = dims[1].getSize();
        ny = dims[2].getSize();
        nx = dims[3].getSize();
        startp.push_back(0);
        startp.push_back(0);
        startp.push_back(0);
        countp.push_back(nz);
        countp.push_back(ny);
        countp.push_back(nx);
    }
    else
        throw std::runtime_error("Unhandled variable dimension.");

    n = nz * ny * nx; //total elements in this variable

    if (time_level>nt)
    {
        throw std::runtime_error("Requesting a time level that does not exist.");
    }
    double *pvalues = new double [n];

    nc_var.getVar(startp, countp, pvalues);
    std::vector<double> v(n);
    for(size_t j=0; j<n; j++)
        v[j] = pvalues[j];
    variables_dict[var_name] = v;
    delete [] pvalues;
    return 0;
}

int DataNC::get_data()
{
    get_data_const("mask_rho");
    get_data_const("h");
    get_data_const("s_rho");
    get_data_const("eta_rho");
    get_data_const("xi_rho");
    return 0;
}

int DataNC::get_data(const size_t time_level)
{
    get_data("zeta", time_level);
    get_data("O2", time_level);
    get_data("temp", time_level);
    get_data("salt", time_level);
    compute_z_levels();
    return 0;
}

size_t & DataNC::get_nt()
{
    return nt;
}

size_t & DataNC::get_nx()
{
    return nx;
}

size_t & DataNC::get_ny()
{
    return ny;
}

size_t & DataNC::get_nz()
{
    return nz;
}

std::vector<double> & DataNC::get(std::string var_name)
{
    return variables_dict[var_name];
}

// Private methods:

int DataNC::compute_z_levels()
{
    // check if we have all the ingredients ready
    std::vector<double> h = variables_dict["h"];
    std::vector<double> zeta =variables_dict["zeta"];
    std::vector<double> s_rho = variables_dict["s_rho"];
    //if ( (h.size()==0) || (zeta.size()==0) )
    //    throw std::runtime_error("No data available for h and/or zeta.");

    size_t n = nz * ny * nx;
    std::vector<double> z(n);
    n=0;
    for (size_t kdx=0; kdx<nz; kdx++)
        for (size_t idx=0; idx<ny * nx; idx++)
            z[n++] = s_rho[kdx] * (h[idx] + zeta[idx]);
    variables_dict["z"] = z;
    /* z is not a netcdf variable as such, so we cannot find
       the unit from the nc file, so we set it here.
    */
    units_dict["z"] = "meter";
    return 0;
}

int DataNC::set_dimensions()
{
    // this method relies on the variable "temp" to be present...
    netCDF::NcVar nc_var = datafile.getVar("temp");
    std::vector<netCDF::NcDim> dims = nc_var.getDims();
    // sanity check
    if (dims.size() != (size_t) DIM4)
        throw std::runtime_error("Not all dimension information provided. Failed to determine all dimensions.");
    nt = dims[0].getSize();
    nz = dims[1].getSize();
    ny = dims[2].getSize();
    nx = dims[3].getSize();
    return 0;
}


// PycnoNC data class

void PycnoNC::open(std::string filename)
{
    datafile.open(filename, netCDF::NcFile::replace); // overwrites existing one.
}

int PycnoNC::create_dimensions(const std::vector<double> pycnoclines,
                               const std::vector<double> eta_rho,
                               const std::vector<double> xi_rho)
{
    nz = pycnoclines.size();
    ny = eta_rho.size();
    nx = xi_rho.size();
    // create the dimensions
    netCDF::NcDim tDim = datafile.addDim("time"); // unlimited dimension
    netCDF::NcDim zDim = datafile.addDim("pycnocline", nz);
    netCDF::NcDim yDim = datafile.addDim("eta_rho", ny);
    netCDF::NcDim xDim = datafile.addDim("xi_rho", nx);
    // set the dimensions in the private attribute dimVector for later use.
    dimVector.push_back(tDim);
    dimVector.push_back(zDim);
    dimVector.push_back(yDim);
    dimVector.push_back(xDim);

    //write the coordinates
    netCDF::NcVar xi_rho_var = datafile.addVar("xi_rho", netCDF::ncDouble, xDim);
    netCDF::NcVar eta_rho_var = datafile.addVar("eta_rho", netCDF::ncDouble, yDim);
    netCDF::NcVar pycnoclines_var = datafile.addVar("pycnoclines", netCDF::ncDouble, zDim);

    //write the coordinate values
    write_vector_variable(pycnoclines_var, pycnoclines);
    write_vector_variable(eta_rho_var, eta_rho);
    write_vector_variable(xi_rho_var, xi_rho);

    // we could put some more attributes here, if needed.
    eta_rho_var.putAtt("standard_name", "y_grid_index");
    xi_rho_var.putAtt("standard_name", "x_grid_index");
    pycnoclines_var.putAtt("standard_name", "pycnocline levels");
    pycnoclines_var.putAtt("units", "kilogram meter^-3");
    return 0;
}

void PycnoNC::write_vector_variable(netCDF::NcVar ncv, const std::vector<double> & v)
{
    size_t n = v.size();
    double *pvalues = new double [n];
    for (size_t i=0; i<n; i++)
        pvalues[i] = v[i];
    ncv.putVar(pvalues);
    delete [] pvalues;
}

int PycnoNC::write_parameter(const std::vector<double> s, const std::string variable_name, size_t rec)
{
    // get the nc variable first
    netCDF::NcVar & v = surface_variables[variable_name];

    std::vector<size_t> startp, countp;
    startp.push_back(rec); // first time
    startp.push_back(0); // pycnocline start index
    startp.push_back(0); // eta start index
    startp.push_back(0); // xi start index
    countp.push_back(1); // write one time level at a time.
    countp.push_back(nz);
    countp.push_back(ny);
    countp.push_back(nx);

    // convert vector to pointer and write the variable.
    size_t n = s.size();
    double *pvalues = new double [n];
    for (size_t i=0; i<n; i++)
    {
        pvalues[i] = s[i];
    }
    v.putVar(startp, countp, pvalues);
    delete [] pvalues;
    return 0;
}

void PycnoNC::create_surface_variable(std::string variable_name, std::string units)
{

    netCDF::NcVar v = datafile.addVar(variable_name, netCDF::ncDouble, dimVector);
    v.putAtt("units", units);
    surface_variables[variable_name] = v;
}
