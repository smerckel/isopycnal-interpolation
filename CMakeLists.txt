cmake_minimum_required(VERSION 3.16)

project("isopycnal_interpolation")

set(warnings "-Wall -Wextra -Werror")
set(compiler_flags "-O3")
set(CMAKE_CXX_FLAGS "${warnings} ${compiler_flags}")

set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local") # Can also be set on command line cmake -DCMAKE_INSTALL_PATH=path etc.

add_executable(ncii src/interpolation.cpp src/main.cpp src/nc_data.cpp src/rho.cpp)
include_directories(BEFORE include)
target_link_libraries(ncii netcdf netcdf_c++4)

install(TARGETS ncii
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib)
