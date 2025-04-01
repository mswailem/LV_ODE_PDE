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

Fourier transform plots can be produced using:


``` make plot_k_u ```

for predators and 

``` make plot_k_v ```

for prey.

Spatial fluctuation data can be produced using:

``` make compute_std ```

Note, the above make commands require data to be present in the output directory, which should automatically be generated after running "./solve"

The raw output of the code is in the output directory 

This code using exponential time differencing techniques to solve the system by utilizing Fourier Transforms. 

The method can be solved directly in "regular" space (outdated and haven't been updated or checked in a while, not recommended to use), or first transformed into log-coodinates and then solved in "log" space.

The log space method ensures that the densities stay positive, however, it both methods yield similar results.

I have implemented the following initial conditions:

checkboard: Discontinous initial conditions where the initial predator and prey are laid out in a checkboard pattern

continous: Smooth initial conditions where the initial predator and prey are distributed using sines and cosines

constant: Constant initial conditions

random: Each grid point is assigned the fixed point density value of the average environment + a random number between 0 and 1e-4 (this parameter can change)

Initial conditions are specified by setting the corresponding initial condition in "ETD2_solver.cpp" in the solve_in_log() function. (Not the best way to implement this but it works for now)

example run:

``` ./solve 0 90 0.1 20 64 0.2 1 1 1 1 1 ```

You can just run ``` ./solve ``` to see a terminal prompt that will let you know what the required arguments are.


TODO:
	- Implement a more user-friendly interface
    - Implement more variables for the bifurcation and stability diagrams
	- Implement a GUI
	- Implement more initial conditions
	- Implement ETDRK4 stepping

