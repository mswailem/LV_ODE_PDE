cmake_minimum_required(VERSION 3.0)
project(LV_mean_field)


set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
include_directories(include)
# Find the GSL package
find_package(GSL REQUIRED)
find_package(OpenMP)

if(NOT OpenMP_CXX_FOUND)
    message(SEND_ERROR "OpenMP could not be found. Please ensure it is installed and try again.")
else()
    message(STATUS "Found OpenMP.")
endif()

set(SOURCE_FILES_MAIN src/main.cpp src/DESolver.cpp)
set(SOURCE_FILES_TEST src/testing.cpp src/DESolver.cpp)
# Add executable target with source files
add_executable(main ${SOURCE_FILES_MAIN})
add_executable(test ${SOURCE_FILES_TEST})

# Link GSL libraries to your application
target_link_libraries(main GSL::gsl GSL::gslcblas)

target_link_libraries(test GSL::gsl GSL::gslcblas)

if(OpenMP_CXX_FOUND)
    target_link_libraries(main OpenMP::OpenMP_CXX)
    target_link_libraries(test OpenMP::OpenMP_CXX)
endif()
