# Lorenz-Cycle

## What is the Lorenz Energy Cycle?

The Lorenz Energy Cycle (LEC) is a methodology for estimating the energy in the Atmosphere designed by the mathematician and meteorologist Edward Lorenz, in 1965. On it, the energy is partitioned into zonal and eddy components of Kinectic Energy and  Available Potential Energy. It is also computed the terms related to the conversion between each form of energy and the terms related to generation/dissipation of the energy contents. Lorenz initial work computed the energetics for the whole globe, but later developments adapted the methodology for computing the energetics of a limited area in the atmosphere, then requiring the computation of terms related to the transport of energy thorugh the boundaries of the limited area.

## What does the program do?

The Lorenz-Cycle program is designed for computing the LEC for a specific region on the atmosphere and can be run using two distinct frameworks:

1. Eulerian framewrok, where the domain is fixed in time by the file: inputs/box_limits; 

2. Lagrangian framework, where the domain can follow a pertubation on the atmosphere, which track is defined in the file: inputs/track. 

More details on running the program are provided bellow.

# Usage

## Eulerian framework

The first step is to delimite the computational domain in the csv text file [boxlimits](inputs/box_limits) that should look like this:

![image](https://user-images.githubusercontent.com/56005607/206709581-34ebe0a7-ff45-4bd4-86e0-8cce8dde91ea.png)

Then, it is required to specify in the [fvars](inputs/fvars) file the how the variables are named in the NetCDF file and which units are being used:  

![image](https://user-images.githubusercontent.com/56005607/206710982-37e2c5a3-4379-41c2-840a-ba2af0cd28a9.png)

Note that you should only change the colums corresponding to "Variable" and "Units". Modying the first column will make the program unable to look for the variables in the NetCDF file.

