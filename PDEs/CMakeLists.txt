cmake_minimum_required(VERSION 3.10)

# Set the project name and version
project(LV_edt2)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Find the FFTW library
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW3 REQUIRED fftw3)

# Add executable
add_executable(solve src/solve.cpp src/ETD2_solver.cpp)
add_executable(fft_test src/fft_test.cpp)
include_directories(${FFTW3_INCLUDE_DIRS} include/)

add_custom_target(plot
    COMMAND gnuplot ${CMAKE_SOURCE_DIR}/scripts/plot.gp
    DEPENDS solve
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/output
    COMMENT "Generating plots with Gnuplot"
)
# Link FFTW to the executable
target_link_libraries(solve ${FFTW3_LIBRARIES})
target_link_libraries(fft_test ${FFTW3_LIBRARIES})

