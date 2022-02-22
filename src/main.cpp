#include <iostream>
#include <string>

#include "interpolation.hpp"
#include "rho.hpp"
#include "nc_data.hpp"
#include "cxxopts.hpp"

using namespace std;

int main(int argc, char** argv)
{
    std::string input_filename;
    std::string output_filename;

    /* parse command line options */
    cxxopts::Options options("pycnocline",
    "\nPycnocline is a C++ program that interpolates 3D variables onto one or more pycnocline\nlevels. The program reads from a netCDF file and writes the results into a new netCDF\nfile.\n");

    options.add_options()
    ("input", "Input netCDF file", cxxopts::value<std::string>())
    ("output", "Output netCDF file", cxxopts::value<std::string>())
    ("v,variables", "(list of) variable(s) to interpolate", cxxopts::value<std::vector<std::string>>())
    ("p,pycnocline_levels", "(list of) pycnocline density or densities", cxxopts::value<std::vector<double>>())
    ("h,help", "Print usage")
    ;
    options.parse_positional({"input", "output"});

    auto result = options.parse(argc, argv);
    if (result.count("help"))
    {
      std::cout << options.help() << std::endl;
      exit(0);
    }
    if (result.count("variables")==0)
    {
        std::cout << "It is mandatory to supply a list of variables. (Use -v variable0,<varibale1>,...)" << std::endl;
        exit(1);
    }
    if (result.count("pycnocline_levels")==0)
    {
        std::cout << "It is mandatory to supply a list of pycnocline levels. (Use -p 1026.3,...)" << std::endl;
        exit(2);
    }
    if ( (result["input"].count()!=1) || (result["output"].count()!=1))
    {
        std::cout << "Please supply at both <input netCDF file> and <output netCDF file>." << std::endl;
        exit(3);
    }
    if (result.unmatched().size()!=0)
    {
        std::cout << "Failed to parse the command line arguments properly." << std::endl;
        std::cout << "Note that lists are entered without spaces (-p 1026.1,1027.1)" << std::endl;
        exit(4);
    }
    // All input should be ok. Set the appropriate variables:
    input_filename = result["input"].as<std::string>();
    output_filename = result["output"].as<std::string>();
    std::vector<double> pycnoclines = result["pycnocline_levels"].as<std::vector<double>>();
    std::map<std::string, std::vector<double>> interpolation_variables{};
    std::vector<string> v = result["variables"].as<std::vector<std::string>>();
    interpolation_variables["z"] = {};
    for (auto it = v.begin(); it!=v.end(); ++it)
        interpolation_variables[*it] = {};
    // Tell 'em what we are going to do...
    std::cout << "Going to read from: " << input_filename << ";" << std::endl;
    std::cout << "Going to write to : " << output_filename << ";" << std::endl;
    std::cout << "Going to process the following variables:" << std::endl;
    for (auto it = interpolation_variables.begin(); it!=interpolation_variables.end(); ++it)
        std::cout << "\t" << it->first << std::endl;
    std::cout << "and interpolate them on the the following pycnoclines:" << std::endl;
    for (auto it = pycnoclines.begin(); it!=pycnoclines.end(); ++it)
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
    data_out.create_dimensions(pycnoclines, data.get("eta_rho"), data.get("xi_rho"),
                               data.get("eta_v"), data.get("xi_u"));;

    std::cerr << "Create pycnocline-interpolated variables for output file..." << std::endl;
    for(auto it = interpolation_variables.begin(); it != interpolation_variables.end(); ++ it)
        data_out.create_surface_variable(it->first, data.get_unit(it->first),
                                         data.get_variable_coordinates(it->first));

    std::cerr << "Start interpolation..." << std::endl;

    for (size_t j=0; j<data.get_nt(); j++)
    {
        std::cerr << "Get data for time level " << j << "..." << std::endl;
        data.get_data(j); // gets the time fields for level 0

        std::cerr << "Interpolate fields..." << std::endl;
        for (size_t i=0; i<pycnoclines.size(); ++i)
            interp.interpolate_onto_surface(interpolation_variables, data, pycnoclines[i]);
        std::cerr << "Write fields..." << std::endl;
        //write the fields
        for(auto it = interpolation_variables.begin(); it != interpolation_variables.end(); ++ it)
            data_out.write_parameter(it->second, it->first, data.get_variable_coordinates(it->first), j);

        if (j==0)
            break;
    }
    return 0;
}
