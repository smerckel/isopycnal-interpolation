add_library(gswteos-10 STATIC gsw_oceanographic_toolbox.c gsw_saar.c)
include_directories(BEFORE ".")
set(CMAKE_C_FLAGS "-O3")
