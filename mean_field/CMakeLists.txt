cmake_minimum_required(VERSION 3.0)
project(LV_mean_field)

# Find the GSL package
find_package(GSL REQUIRED)

# Add executable target with source files
add_executable(solve src/ODE_solver.cpp)
add_executable(fm src/fm.cpp)
add_executable(homotopy src/homotopy.cpp)
add_executable(chaos_diagram src/chaos_diagram.cpp)
add_executable(debugging_test src/debugging_test.cpp)

# Link GSL libraries to your application
target_link_libraries(solve GSL::gsl GSL::gslcblas)
target_link_libraries(fm GSL::gsl GSL::gslcblas)
target_link_libraries(homotopy GSL::gsl GSL::gslcblas)
target_link_libraries(chaos_diagram GSL::gsl GSL::gslcblas)
target_link_libraries(debugging_test GSL::gsl GSL::gslcblas)

