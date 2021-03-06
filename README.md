# NetCDF isopycnal interpolation for ROMS model data (ncii)


## Changelog

v.0.1: Initial release
v.0.2: Adds option to calculate averages between isopycnals

## Synopsis


C++ utility program to interpolate variables on a pycnocline, reading netcdf files

This utility program takes a netCDF file generated by the ROMS model as input and interpolates selected variables 
on user-defined pycnocline levels. The resulting 2D fields are ouput to a separate netCDF file for further data
analysis.

## Installing

Ncii depends on external libraries for the netcdf interface, and assumes these to be installed in the usual places.
These libraries are:
1) libnetcdf, and
2) libnetcdf_c++4.

On a linux these can usually be installed using your system's package manager. Ensure that the development packages 
are also installed (header files). On Fedora, for example, this can be achieved by

`$ sudo dnf install netcdf netcdf-devel netcdf-cxx4 netcdf-cxx4-devel`

To compile the code, you can either use the `codeblocks` IDE, or using `make` on makefile generated with `cmake`.

Using `make` you would typically make a `build` directory in the top directory of ncii, change directory to build 
and run cmake from there:

`$ mkdir build`
`$ cd build`
`$ cmake ..`
`$ make`

Optionally you can run `make install`, which installs `ncii` in `~/.local/bin` by default. The setting of the installation path
can be overriden by supplying the option `-DCMAKE_INSTALL_PATH=<your/custom/install/path` to the `cmake` command.

## Using ncii

Ncii can be called with the option --help to show a brief description of the program.

To do something useful, however, ncii needs to be called with two positional arguments: a filename pointing to the 
input (netCDF) file and a filename of a filename where the output is written to. If the output file exists, it will
be overwritten. It is manadatory to specify the option parameters -v and -p to define a list of variables that need 
to be interpolated, and a list of pycnocline densities, respectively. It is important that if a multiple values are
specified, they are entered as comma-separated list, without spaces. 

For example:

`$ ncii -v O2,temp,salt -p 1026.5,1027.1 input.nc output.nc`

would interpolate the 3D fields of O2, temp and salt of the input.nc on to the surfaces defined by the pycnoclines
with densities of 1026.5 and 1027.1 kg m<sup>-3</sup>, which are written to output.nc.

## Disclaimer

This software is provided "as is" and any expressed or implied
warranties, including, but not limited to, the implied warranties of
merchantability and fitness for a particular purpose are
disclaimed. In no event shall the regents or contributors be liable
for any direct, indirect, incidental, special, exemplary, or
consequential damages (including, but not limited to, procurement of
substitute goods or services; loss of use, data, or profits; or
business interruption) however caused and on any theory of liability,
whether in contract, strict liability, or tort (including negligence
or otherwise) arising in any way out of the use of this software, even
if advised of the possibility of such damage.

