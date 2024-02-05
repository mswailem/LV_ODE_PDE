cmake_minimum_required(VERSION 3.10)

# Set the project name and version
project(LV_edt2)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find the FFTW library
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW3 REQUIRED fftw3)

# Add executable
add_executable(LV_edt2 main.cpp)
include_directories(${FFTW3_INCLUDE_DIRS})

# Link FFTW to the executable
target_link_libraries(LV_edt2 ${FFTW3_LIBRARIES})

