# Project 4 - Evolution of Electron Spin and Analysis of Thermodynamic Properties

**Formal write up is located in this repository, named "Rafizadeh_WRITEUP.pdf"** 

**IMPORTANT NOTE: Simulating spin evolution in crystal lattices is complex and can be done many different ways. This script uses Monte Carlo simulation, which is computationally taxing, combined with the fact that I am not great at writing efficient code makes these scripts incredibly inefficient. I have timed the scripts to take the following times to run:**

`rng_electron_spin.py`: 90 seconds

`plot_electron_spin.py`: 40 seconds

**I recommend getting a drink and putting on some smooth jazz before running the scripts.**

## Prior Generation

The simulation for eletron spin interactions is going to be from a Monte Carlo simulation to evaluate thermodynamic properties at different temperature values. It's easiest to use a uniform distribution for these values, so the data generation file includes generation of the uniform distribution for some random values that can be used in the computation. Physically, the data generated from this uniform distribution is the initial spin states of the eletrons in the lattice. A demonstration of this is shown in the plot at the very bottom of this file. The plot of the distribution is saved as `UniformDistribution.png`, and looks as follows:

![UniformDistribution](https://github.com/rafizadehn/PHSX815_Project4/assets/76142511/5d02cbe8-45a6-4b06-9ebc-81306ea7d96a)

How to generate this plot is discussed below, in Data Generation. 

## Data Generation

The thermodynamic parameters are calculated by the `rng_electron_spin.py` python file. This file requires python3 to run, and includes the following packages listed at the top of the script:

```
  from __future__ import division
  import sys
  import numpy as np
  from numpy.random import rand
  import matplotlib.pyplot as plt
```

To run this script from the terminal in linux, run:

> $ python3 rng_electron_spin.py

This runs the file with the default parameters, which are: 50 temperature points, a 10x10 atom lattice, and 555 as the default seed for the random number generator.

These values can be altered from the command line in the terminal by simply adding an argument after the file name. The arguments to change these include `-Npoints`, `-Natoms`, and `-seed` for those values, respectively. Keep in mind, the larger the lattice becomes, the longer it will take for the script to compile. 

For example, it may looks something like this in linux:

> $ python3 rng_electron_spin.py -Npoints 5 -Natoms 5 -seed607

The parameters chosen for this file are written to a text file called `parameters.txt`, and thermodynamic property information is stored in respective text files. This allows for the analysis script to read the parameters that were used in the generation file for more efficient analysis.

## Data Analysis

The generated data is plotted and analyzed by the `plot_electron_spin.py` python file. This file requires python3 to run, and includes the following packages listed at the top of the script:

```
  from __future__ import division
  import sys
  import numpy as np
  from numpy.random import rand
  import matplotlib.pyplot as plt
```

To run this script from the terminal in linux, run:

> $ python3 plot_eletron_spin.py

This creates plots from the input data sets and analyzes it. 

The first set of data analyzed is the energy data set, of which the plot with default values looks like:

![EnergyPlot](https://github.com/rafizadehn/PHSX815_Project4/assets/76142511/d8edbdb2-b4ec-44cc-b488-6915bbc6abef)

then the magnetization,

![MagnetizationPlot](https://github.com/rafizadehn/PHSX815_Project4/assets/76142511/31b35305-c3da-4ab6-8e1a-3a903e9025b7)

then the specific heat,

![SpecificHeatPlot](https://github.com/rafizadehn/PHSX815_Project4/assets/76142511/ffeb9de0-a47f-41b1-b6d1-284328d6e1a4)

then the susceptibility,

![SusceptibilityPlot](https://github.com/rafizadehn/PHSX815_Project4/assets/76142511/b409597a-11ee-4fb9-96c4-8c6d42726f46)

which can all be used to quantitatively determine the critical, or Curie, temperature for the simulated lattice. 

In addition, I found some helpful scripts online that are able to simulate the evolution of the eletron spins in the lattice graphically. Combining some of them gives:

![ElectronSpinSimulation](https://github.com/rafizadehn/PHSX815_Project4/assets/76142511/014e871a-f502-4b3a-8ca2-eb02dcada906)

which is a wonderful demonstration of how these thermodynamic properties could change with time, or temperature. The very first plot, at T=0, shows the randomly generated uniform distribution used at the beginning for the Monte Carlo simulation.  

## Sources

Code was adapted and frankenseined from:

[Using Monte Carlo methods to study lattice structures.](https://towardsdatascience.com/monte-carlo-method-applied-on-a-2d-binary-alloy-using-an-ising-model-on-python-70afa03b172b)

[Extracting thermodynamic properties from Monte Carlo configuration simulations.](https://github.com/prtkm/ising-monte-carlo/blob/master/ising-monte-carlo.org)

Sources of plots in writeup:

[Experimentally measured magnetic properties.](https://doi.org/10.1007/s10751-019-1571-1)

More information regarding magnetization:

[Magnetic Susceptibility.](https://en.wikipedia.org/wiki/Magnetic_susceptibility)

[Magnetization.](https://en.wikipedia.org/wiki/Magnetization)

