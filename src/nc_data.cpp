#include "nc_data.hpp"


int DataNC::open(const std::string filename, const std::vector<std::string> v)
{
    datafile.open(filename, netCDF::NcFile::read);
    set_dimensions();
    interpolation_variables=v;
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

int DataNC::get_variable_coordinates(const std::string var_name)
{
    return variable_coordinates_dict[var_name];
}

void DataNC::set_variable_coordinates(const std::string var_name)
{
    // z will be a default variable that gets interpolated, but it is computed and not looked up
    // in the netcdf file.
    if (var_name=="z")
        variable_coordinates_dict[var_name] = rho_coordinates;
    else
    {
        netCDF::NcVar nc_var = datafile.getVar(var_name);
        std::vector<netCDF::NcDim> dims = nc_var.getDims();
        size_t _ny=1, _nx=1;
        _ny = dims[2].getSize();
        _nx = dims[3].getSize();
        if ((ny==_ny) && (nx==_nx))
            variable_coordinates_dict[var_name] = rho_coordinates;
        else if ((ny==_ny) && (nx-1==_nx))
            variable_coordinates_dict[var_name] = u_coordinates;
        else if ((ny-1==_ny) && (nx==_nx))
            variable_coordinates_dict[var_name] = v_coordinates;
        else
            throw(std::runtime_error("Unknown coordinate type."));
    }
}

int DataNC::get_attributes(const std::string att_name, double & f)
{
    std::vector<double> _f{0};
    get_attributes(att_name, _f);
    f=_f[0];
    return 0;
}

int DataNC::get_attributes(const std::string att_name, std::vector<double> &f)
{
    netCDF::NcGroupAtt nc_att = datafile.getAtt(att_name);
    size_t n = nc_att.getAttLength();
    f.resize(n);
    double* pvalues = new double [n];
    nc_att.getValues(pvalues);
    for(size_t i=0; i< n; i++)
        f[i] = pvalues[i];
    delete [] pvalues;
    return 0;
}


int DataNC::get_data_const(const std::string var_name)
{
    // get attributes for z calculation
    get_attributes("hc", hc);
    get_attributes("Cs_r", Cs_r);
    get_attributes("sc_r", sc_r);

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
    size_t _ny=1, _nx=1, _nz=1;

    startp.push_back(time_level);
    countp.push_back(1);

    if (dims.size()==(size_t) DIM3)
    {
        nt = dims[0].getSize();
        _ny = dims[1].getSize();
        _nx = dims[2].getSize();
        startp.push_back(0);
        startp.push_back(0);
        countp.push_back(_ny);
        countp.push_back(_nx);
    }
    else if (dims.size() == (size_t) DIM4)
    {
        nt = dims[0].getSize();
        _nz = dims[1].getSize();
        _ny = dims[2].getSize();
        _nx = dims[3].getSize();
        startp.push_back(0);
        startp.push_back(0);
        startp.push_back(0);
        countp.push_back(_nz);
        countp.push_back(_ny);
        countp.push_back(_nx);
    }
    else
        throw std::runtime_error("Unhandled variable dimension.");

    n = _nz * _ny * _nx; //total elements in this variable

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
    // set the coordinate system used for this variable
    if ((ny==_ny) && (nx==_nx))
        variable_coordinates_dict[var_name] = rho_coordinates;
    else if ((ny==_ny) && (nx-1==_nx))
        variable_coordinates_dict[var_name] = u_coordinates;
    else if ((ny-1==_ny) && (nx==_nx))
        variable_coordinates_dict[var_name] = v_coordinates;
    else
        throw(std::runtime_error("Unknown coordinate type."));
    return 0;
}

int DataNC::get_data()
{
    get_data_const("mask_rho");
    get_data_const("h");
    get_data_const("eta_rho");
    get_data_const("xi_rho");
    get_data_const("eta_v");
    get_data_const("xi_u");
    return 0;
}

int DataNC::get_data(const size_t time_level)
{
    get_data("zeta", time_level);
    get_data("temp", time_level);
    get_data("salt", time_level);
    // Read the user requested variables (except z which we compute)
    for(auto it = interpolation_variables.begin(); it != interpolation_variables.end(); ++it)
    {
        if (*it == "z")
            continue;
        get_data(*it, time_level);
    }
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
    //std::vector<double> s_rho = variables_dict["s_rho"];
    //if ( (h.size()==0) || (zeta.size()==0) )
    //    throw std::runtime_error("No data available for h and/or zeta.");

    // compute depth from s_rho using the matlab code in ../doc as template
    size_t n = nz * ny * nx;
    std::vector<double> z(n);
    double cff, cff1;
    n=0;
    for (size_t kdx=0; kdx<nz; kdx++)
    {
        cff = hc *sc_r[kdx];
        cff1 = Cs_r[kdx];
        for (size_t idx=0; idx<ny * nx; idx++)
            z[n++] = zeta[idx] + (zeta[idx] + h[idx]) * (cff + cff1*h[idx]) / (h[idx] + hc);
    }
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
                               const std::vector<double> xi_rho,
                               const std::vector<double> eta_v,
                               const std::vector<double> xi_u)
{
    nz = pycnoclines.size();
    ny = eta_rho.size();
    nx = xi_rho.size();
    // create the dimensions
    netCDF::NcDim tDim = datafile.addDim("time"); // unlimited dimension
    netCDF::NcDim zDim = datafile.addDim("pycnocline", nz);
    netCDF::NcDim yDim = datafile.addDim("eta_rho", ny);
    netCDF::NcDim xDim = datafile.addDim("xi_rho", nx);
    netCDF::NcDim yvDim = datafile.addDim("eta_v", ny-1);
    netCDF::NcDim xuDim = datafile.addDim("xi_u", nx-1);
    // set the dimensions in the private attribute dimVector for later use.
    dimVector_rho.push_back(tDim);
    dimVector_rho.push_back(zDim);
    dimVector_rho.push_back(yDim);
    dimVector_rho.push_back(xDim);

    dimVector_u.push_back(tDim);
    dimVector_u.push_back(zDim);
    dimVector_u.push_back(yDim);
    dimVector_u.push_back(xuDim);

    dimVector_v.push_back(tDim);
    dimVector_v.push_back(zDim);
    dimVector_v.push_back(yvDim);
    dimVector_v.push_back(xDim);

    //write the coordinates
    netCDF::NcVar xi_rho_var = datafile.addVar("xi_rho", netCDF::ncDouble, xDim);
    netCDF::NcVar eta_rho_var = datafile.addVar("eta_rho", netCDF::ncDouble, yDim);
    netCDF::NcVar pycnoclines_var = datafile.addVar("pycnoclines", netCDF::ncDouble, zDim);
    netCDF::NcVar xi_u_var = datafile.addVar("xi_u", netCDF::ncDouble, xuDim);
    netCDF::NcVar eta_v_var = datafile.addVar("eta_v", netCDF::ncDouble, yvDim);

    //write the coordinate values
    write_vector_variable(pycnoclines_var, pycnoclines);
    write_vector_variable(eta_rho_var, eta_rho);
    write_vector_variable(xi_rho_var, xi_rho);
    write_vector_variable(eta_v_var, eta_v);
    write_vector_variable(xi_u_var, xi_u);

    // we could put some more attributes here, if needed.
    eta_rho_var.putAtt("standard_name", "y_grid_index");
    xi_rho_var.putAtt("standard_name", "x_grid_index");
    eta_v_var.putAtt("standard_name", "y_grid_index v-velocity");
    xi_u_var.putAtt("standard_name", "x_grid_index u-velocity");
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

int PycnoNC::write_parameter(const std::vector<double> s, const std::string variable_name,
                             const int variable_coordinates, const size_t rec)
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
    countp.push_back(ny - (size_t) (variable_coordinates==v_coordinates));
    countp.push_back(nx - (size_t) (variable_coordinates==u_coordinates));

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

void PycnoNC::create_surface_variable(const std::string variable_name, const std::string units,
                                      const int variable_coordinates)
{
    std::vector<netCDF::NcDim> _dimVector;
    switch (variable_coordinates)
    {
        case rho_coordinates:
            _dimVector=dimVector_rho; break;
        case u_coordinates:
            _dimVector=dimVector_u; break;
        case v_coordinates:
            _dimVector=dimVector_v; break;
    }
    netCDF::NcVar v = datafile.addVar(variable_name, netCDF::ncDouble, _dimVector);
    v.putAtt("units", units);
    surface_variables[variable_name] = v;
}
