cmake_minimum_required(VERSION 3.16)

project("pycnocline")
set(warnings "-Wall -Wextra -Werror")
set(CMAKE_CXX_FLAGS "${warnings}")
#set(CMAKE_INSTALL_PREFIX "path") # Can also be set on command line cmake -DCMAKE_INSTALL_PATH=path etc.

add_executable(pycnocline src/interpolation.cpp src/main.cpp src/nc_data.cpp src/rho.cpp)
include_directories(BEFORE include gsw_c_v3.05)
target_link_libraries(pycnocline netcdf netcdf_c++4 gswteos-10)

add_subdirectory("gsw_c_v3.05")

install(TARGETS pycnocline
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib)
