#include <iostream>
#include <string>

#include "interpolation.hpp"
#include "rho.hpp"
#include "nc_data.hpp"
#include "cxxopts.hpp"
#include "ndarray.hpp"

const std::string VERSION{"0.3"};

int main(int argc, char** argv)
{
    std::cout << "+----------------------------------------------------------------+" << std::endl;
    std::cout << "| ROMS model netcdf isopyncal interpolation (ncii) version " << VERSION << "   |" <<std::endl;
    std::cout << "| Copyright Lucas Merckelbach (lucas.merckelbach@hereon.de) 2022 |" << std::endl;
    std::cout << "+----------------------------------------------------------------+" << std::endl << std::endl;

    std::string input_filename;
    std::string output_filename;

    /* parse command line options */
    cxxopts::Options options("ncii",
    "\nROMS model netCDF isopycnal interpolation (ncii)\n\nNCII is a C++ program that interpolates 3D variables onto one or more isopycnal\nlevels. The program reads from a netCDF file and writes the results into a new netCDF\nfile.\n");

    options.add_options()
    ("input", "Input netCDF file", cxxopts::value<std::string>())
    ("output", "Output netCDF file", cxxopts::value<std::string>())
    ("a,avg", "Average variables between isopycnal levels")
    ("v,variables", "(list of) variable(s) to interpolate", cxxopts::value<std::vector<std::string>>())
    ("p,isopycnal_levels", "(list of) isopycnal density or densities", cxxopts::value<std::vector<double>>())
    ("h,help", "Print usage")
    ("version", "Print version number")
    ;
    options.parse_positional({"input", "output"});
    int exit_value=0;
    auto result = options.parse(argc, argv);
    if (result.count("help"))
    {
      std::cout << options.help() << std::endl;
      exit(0);
    }
    if (result.count("version"))
    {
        // just exit, header already shown says it all.
        exit(0);
    }
    if (result.count("variables")==0)
    {
        std::cout << "It is mandatory to supply a list of variables. (Use -v variable0,<varibale1>,...)" << std::endl;
        exit_value |= 1;
    }
    if (result.count("isopycnal_levels")==0)
    {
        std::cout << "It is mandatory to supply a list of isopycnal levels. (Use -p 1026.3,...)" << std::endl;
        exit_value |= 2;
    }
    if ( (result["input"].count()!=1) || (result["output"].count()!=1))
    {
        std::cout << "Please supply at both <input netCDF file> and <output netCDF file>." << std::endl;
        exit_value |= 4;
    }
    if (result.unmatched().size()!=0)
    {
        std::cout << "Failed to parse the command line arguments properly." << std::endl;
        std::cout << "Note that lists are entered without spaces (-p 1026.1,1027.1)" << std::endl;
        exit_value |= 8;
    }

    if ( (result["avg"].count()==1) && ((exit_value &2) || (result["isopycnal_levels"].as<std::vector<double>>().size()<2)) )
    {
        std::cout << "Cannot compute averages if less than two isopycnal values are given." << std::endl;
        exit_value |= 16;
    }

    if (exit_value) // we had one or more command line issues. Exit so they can fix it.
    {
        if(exit_value>=7)
            std::cout << "Hint: try 'ncii --help' for more information..." << std::endl;
        exit(exit_value);
    }

    // All input should be ok.

    //Set the appropriate variables:
    input_filename = result["input"].as<std::string>();
    output_filename = result["output"].as<std::string>();

    //isopycnals holds a vector with the density values of the isopycnals:
    std::vector<double> isopycnals = result["isopycnal_levels"].as<std::vector<double>>();

    // interpolation_variables is a dictionary <variable name : 2D array>
    std::map<std::string, NDarray<double>> interpolation_variables{};

    //set the variable names to a vector of strings
    std::vector<std::string> v = result["variables"].as<std::vector<std::string>>();
    interpolation_variables["z"] = {}; // added by default.
    for (auto it = v.begin(); it!=v.end(); ++it)
        interpolation_variables[*it] = {};
    // Tell 'em what we are going to do...
    std::cout << "Going to read from: " << input_filename << ";" << std::endl;
    std::cout << "Going to write to : " << output_filename << ";" << std::endl;
    std::cout << "Going to process the following variables:" << std::endl;
    for (auto it = interpolation_variables.begin(); it!=interpolation_variables.end(); ++it)
        std::cout << "\t" << it->first << std::endl;
    std::cout << "and interpolate them on the the following isopycnals:" << std::endl;
    for (auto it = isopycnals.begin(); it!=isopycnals.end(); ++it)
        std::cout << "\t" << *it << std::endl;

    // The work starts here.
    DataNC data;
    data.open(input_filename, v);

    PycnoNC data_out;
    data_out.open(output_filename);

    Interpolation interp;

    std::cerr << "Read constant data..." << std::endl;
    data.get_data(); // gets constant fields.
    std::cerr << "Store units and coordinate system used..." << std::endl;
    // make a note of the units of all interpolation variables:
    for(auto it = interpolation_variables.begin(); it != interpolation_variables.end(); ++ it)
    {
        data.set_unit(it->first);
        data.set_variable_coordinates(it->first);
    }
    std::cerr << "Create dimensions for output file..." << std::endl;
    data_out.create_dimensions(isopycnals, data.get("eta_rho"), data.get("xi_rho"),
                               data.get("eta_v"), data.get("xi_u"), result["avg"].count()==1);

    std::cerr << "Create isopycnal-interpolated variables for output file..." << std::endl;
    for(auto it = interpolation_variables.begin(); it != interpolation_variables.end(); ++ it)
        data_out.create_surface_variable(it->first, data.get_unit(it->first),
                                         data.get_variable_coordinates(it->first));
    // Set the size of all interpolation_variables
    for(auto it = interpolation_variables.begin(); it != interpolation_variables.end(); ++ it)
    {
        size_t ny = data.get_ny() - (size_t) (data.get_variable_coordinates(it->first) == data.v_coordinates);
        size_t nx = data.get_nx() - (size_t) (data.get_variable_coordinates(it->first) == data.u_coordinates);
        size_t nlayers = isopycnals.size() - (size_t) (result["avg"].count() == 1);
        it->second.resize(nlayers, ny, nx);
    }

    std::cerr << "Start interpolation..." << std::endl;
    for (size_t j=0; j<data.get_nt(); j++)
    {
        std::cerr << "Get data for time level " << j << "..." << std::endl;
        data.get_data(j); // gets the time fields for level j

        if (result["avg"].count() == 0)
        {
            std::cerr << "Interpolate fields..." << std::endl;
            for (size_t i=0; i<isopycnals.size(); ++i)
                interp.interpolate_onto_surface(interpolation_variables, data, i, isopycnals[i]);
        }
        else
        {
            std::cerr << "Computing averages..." << std::endl;
            interp.compute_avg_between_isopycnals(interpolation_variables, data, isopycnals);
        }
        std::cerr << "Write fields..." << std::endl;
        //write the fields
        for(auto it = interpolation_variables.begin(); it != interpolation_variables.end(); ++ it)
            data_out.write_parameter(it->second, it->first, data.get_variable_coordinates(it->first), j);
    }
    return 0;
}
