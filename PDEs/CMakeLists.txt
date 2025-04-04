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

# Custom make command to plot the density fields
add_custom_target(plot
    COMMAND python3 ${CMAKE_SOURCE_DIR}/scripts/plot.py
    DEPENDS solve
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/output
    COMMENT "Generating plots with Gnuplot"
)

# Custom make command to plot the Fourier transform of the predator density
add_custom_target(plot_k_u
    COMMAND python3 ${CMAKE_SOURCE_DIR}/scripts/plot_k_u.py
    DEPENDS solve
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/output
    COMMENT "Generating plots with Gnuplot"
)

# Custom make command to plot the Fourier transform of the prey density
add_custom_target(plot_k_v
    COMMAND python3 ${CMAKE_SOURCE_DIR}/scripts/plot_k_v.py
    DEPENDS solve
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/output
    COMMENT "Generating plots with Gnuplot"
)

# Custom make command to create the data for the spatial fluctuations
add_custom_target(compute_std
    COMMAND python3 ${CMAKE_SOURCE_DIR}/scripts/compute_std.py
    DEPENDS solve
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/output
	COMMENT "Computing std as a function of time"
)


# Custom make command to animate the images procuced from "make plot"
add_custom_target(animate
    COMMAND python3 ${CMAKE_SOURCE_DIR}/scripts/animate.py
    DEPENDS solve
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/images
    COMMENT "Generating gif with Gnuplot"
)

# Link FFTW to the executable
target_link_libraries(solve ${FFTW3_LIBRARIES})
target_link_libraries(fft_test ${FFTW3_LIBRARIES})
