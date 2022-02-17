#include <iostream>
#include <string>

#include "interpolation.hpp"
#include "rho.hpp"
#include "nc_data.hpp"

using namespace std;

int main(int narg, char *argv[])
{
    std::string input_filename;
    std::string output_filename;

    if (narg != 3)
    {
        std::cout << "Supply input NC and output NC file.\n" << narg;
        std::exit(0);
    }
    else
    {
        input_filename = argv[1];
        output_filename = argv[2];
    }
    std::vector<double> pycnoclines{1026.5, 1027.1};

    DataNC data;
    data.open(input_filename);

    PycnoNC data_out;
    data_out.open(output_filename);

    Interpolation interp;

    data.get_data(); // gets constant fields.
    data_out.create_dimensions(pycnoclines, data.get_eta_rho(), data.get_xi_rho());
    data_out.create_surface_variable("O2","umol L-1");
    data_out.create_surface_variable("z_pycnocline","m");
    for (size_t j=0; j<data.get_nt(); j++)
    {
        data.get_data(j); // gets the time fields for level 0

        std::vector<double> s, z;

        for (size_t i=0; i<pycnoclines.size(); ++i)
            interp.interpolate_onto_surface(s, z, data, pycnoclines[i]);

        data_out.write_parameter(s, "O2", j);
        data_out.write_parameter(z, "z_pycnocline", j);
    }

    return 0;
}
