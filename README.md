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

The PDE foler cotains the C++ code for the spatial model.

To run this code, first install FFTW, and then:

``` mkdir build
cd build
cmake ..
make
```

Then you can use the terminal command ./solve to run the code, with the parameters of the system specified as arguments.

There is a gnuplot script that can be used to visualize the output of the code.

This can be run using the command:

``` make plot ```


The raw output of the code is in the output directory 

This code using exponential time differencing techniques to solve the system by utilizing Fourier Transforms. 

The method can be solved directly in "regular" space, or first transformed into log-coodinates and then solved in "log" space.

The log space method ensures that the densities stay positive, however, it seems that both methods yield similar results.

I have implemented two initial conditions:

type 0: discontinous initial conditions where the initial predator and prey are laid out in a checkboard pattern (this leads to divergences in the solutions for some parameter values due to the derivative becoming too large)

type 1: smooth initial conditions where the initial predator and prey are distributed using sines and cosines (this seems to work without divergences)

example run:

``` ./solve 0 90 0.1 20 64 0.2 1 1 1 1 1 ```


TODO:
	- Implement a more user-friendly interface
	- Implement a GUI
	- Implement more initial conditions
	- Implement ETDRK4 stepping

