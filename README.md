In this project, I am using GSL to solve the Lotka-Volterra Predator-prey system of coupled differential equations with a temporally varying parameter (the carrying capacity).

I am also solving the same system as a Partial differential equation to include the effects of space.

the mean_field directory contains the C++ code for the non-spatial model. Running it requires GSL to be installed.

To run this code, first install GSL, and then:

``` mkdir build
cd build
cmake ..
make
```

Builiding this project produces a main executable that can run perform different computations beased on the first command line argument provided to it.

Here is a list of the different computations that could be performed:
| Name              | Description | Output |
|-------------------|:-----------:|-------:|
| time_series       |    Solves the system for a specific time interval     |  t u(t) v(t) |
| stability_(var_type) |    Produces the stability diagram of the system as a function of var_type (for example: k1_vs_n)     |  var1 var2 (for example: k1 n) |
| bifurcation_(var_type) |    Produces the bifurcation diagram of the system as a function of var_type (for example: alpha)     | var ustar vstar (for example: alpha ustar vstar) |

I might eventually turn this into a software application that users can more easily interact with

The PDE folder cotains the C++ code for the spatial model.

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
    - Implement more variables for the bifurcation and stability diagrams
	- Implement a GUI
	- Implement more initial conditions
	- Implement ETDRK4 stepping

