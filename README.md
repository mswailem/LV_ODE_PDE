In this project, I am using GSL to solve the Lotka-Volterra Predator-prey system of coupled differential equations with a temporally varying parameter (the carrying capacity).

I am also solving the same system as a Partial differential equation to include the effects of space.

the mean_field directory contains the C++ code for the non-spatial model. Running it requires GSL to be installed.

To run this code, first install GSL, and then:

``` mkdir build
cd build
cmake ..
make
```

There are several exectuables that this generates in the build directory:
ODE_solver - solves the ODE system and generates a files with the densities of the prey and predator populations
fm: Solves for the Floquet Multipliers
homotopy: Solves the ODE system with a homotopy parameter, that continually connect the system with a time-independent fixed point and a time-dependent fixed point.
chaos_diagram: Generates bifrucation diagrams for the system as the homotopy parameter is varied
debugging_test: Used for debugging purposes only

I might eventually turn this into a software application that users can more easily interact with

For the PDE folder, this is still a work in progress. I have a python implemenmtation but it leads to divergences that should not be there, so instead I am using C++ in the cpp_test folder to solve the system using the FFTW library.

This code using exponential time differencing techniques to solve the system by utilizing Fourier Transforms. The C++ code is still a work in progress.

To build it, you follow the same steps as above for the mean_field program


-- Need to add more details here about how to fetch the output of the code, and add a TODO list
