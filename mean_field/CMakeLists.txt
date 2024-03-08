cmake_minimum_required(VERSION 3.0)
project(LV_mean_field)


set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
include_directories(include)
# Find the GSL package
find_package(GSL REQUIRED)

set(SOURCE_FILES_MAIN src/main.cpp src/DESolver.cpp include/DESolver.h)
# Add executable target with source files
add_executable(main ${SOURCE_FILES_MAIN} ${HEADER_FILES_MAIN})

# Link GSL libraries to your application
target_link_libraries(main GSL::gsl GSL::gslcblas)

