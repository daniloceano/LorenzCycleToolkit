---
title: 'LorenzCycleToolkit: A Tool for Calculating the Lorenz Energy Cycle'
tags:
  - Python
  - meteorology
  - atmospheric dynamics
  - diagnostic
  - cyclones
authors:
  - name: Danilo Couto de Souza
    orcid: 0000-0003-4121-7583
    affiliation: 1
  - name: Pedro Leite da Silva Dias
    orcid: 0000-0002-4051-2962
    affiliation: 1
  - name: Matheus Bonjour Laviola da Silva
    orcid: 0000-0002-8629-8762
    affiliation: 1
    affiliation: 2
  - name: Ricardo de Camargo
    orcid: 0000-0002-9425-5391
    affiliation: 1
affiliations:
  - name: Institute of Astronomy, Geophysics and Atmospheric Sciences of the São Paulo University, Rua do Matão, 226, Cidade Universitária, 05508-090, São Paulo, Brazil
    index: 1
  - name: OceanPact, Rua da Glória, 122, Rio de Janeiro, 20241-180, Rio de Janeiro, Brazil
    index: 2
date: 2024-07-16
bibliography: paper.bib
---

# Summary

The LorenzCycleToolkit is a Python package for calculating the Lorenz Energy Cycle (LEC) in limited atmospheric regions. Introduced by Edward Lorenz in 1965, the LEC estimates atmospheric energy, including zonal and eddy components of Kinetic and Available Potential Energy, and conversions between these forms. This toolkit allows researchers to analyze energy transformations in various atmospheric phenomena, enhancing our understanding of atmospheric dynamics and aiding in weather and climate prediction.

# Background

The Lorenz Energy Cycle (LEC), introduced by [@lorenz1955], is a framework to describe the processes of energy conversion and transfer within the Earth's atmosphere. It serves as a diagnostic tool for the general circulation of the atmosphere, providing insights into its dynamic and thermodynamic structures, and information about atmospheric flows and the transport of heat, moisture, and angular momentum. The LEC divides atmospheric energy into four main reservoirs: available potential energy (APE), kinetic energy (K), and their respective zonal mean (Z) and eddy (E) components. The LEC can indicate how energy is generated, converted, and dissipated within these reservoirs.

![Figure 1: Schematic representation of the Lorenz Energy Cycle. The figure shows the four main energy reservoirs in the atmosphere: zonal available potential energy ($A_Z$), eddy available potential energy ($A_E$), zonal kinetic energy ($K_Z$), and eddy kinetic energy ($K_E$). The arrows indicate the energy transformations between these reservoirs: $C_Z$ (conversion from $A_Z$ to $K_Z$), $C_E$ (conversion from $A_E$ to $K_E$), $C_A$ (conversion between $A_Z$ and $A_E$), and $C_K$ (conversion between $K_Z$ and $K_E$). Additional terms include $B$ (boundary fluxes), $R$ (residuals), $G$ (generation terms), and $D$ (dissipation terms).](figures/LEC_example.png)

Applications of the LEC are diverse and significant in meteorology and climate science. It aids in diagnosing the performance of numerical weather prediction models by comparing simulated and observed energy cycles [@ulbrich1991global]. Moreover, it is used to study climate change impacts, such as shifts in energy transfer processes due to global warming [@li2007lorenz; @veiga2013global]. Understanding the LEC is also crucial for improving the representation of atmospheric processes in general circulation models (GCMs), which are essential for accurate climate projections [@lorenz1967nature; @boer2008energy].

## Statement of Need

`LorenzCycleToolkit` is a Python package designed for calculating the Lorenz Energy Cycle (LEC) in limited atmospheric regions, supporting both Eulerian and Semi-Lagrangian frameworks. The formulations used are the same as in @michaelides1999quasi. The Eulerian framework computes the LEC for a fixed region in the atmosphere [@dias2011energy], requiring users to provide a "box limits" file containing the minimum and maximum latitudes and longitudes for the desired computational domain. The Semi-Lagrangian framework, on the other hand, uses a domain that changes position with time [@michaelides1999quasi], which is especially useful for studying cyclonic systems. For this option, users must specify the central position of the domain for each timestep in a "track" CSV file, which can also specify the domain width and length. If these specifications are not provided, a default 15°x15° domain is used. The `LorenzCycleToolkit` also offers an option for users to interactively choose the domain position and size for each data timestep.

![Figure 2: Schematic representation of the Eulerian and Semi-Lagrangian frameworks. On the left (A), the Eulerian framework is shown with a fixed domain capturing the evolution of a cyclone as it moves through the grid. On the right (B), the Semi-Lagrangian framework is illustrated, which tracks the cyclone with a moving domain that follows the cyclone's path. For each time step, the domain shifts to remain centered on the cyclone. The latitude and longitude axes are shown for spatial reference.](figures/frameworks.pdf)


The program can handle any reanalysis or model data provided in isobaric levels. Users must adjust the "namelist" file, specifying the names for each variable and coordinate in the input file. The requirements for the input file are that it must be in NetCDF format [@rew1990netcdf], the grid must be structured, and it must have isobaric levels as vertical coordinates. These specifications are common in reanalysis, climate, and weather datasets.

The program outputs the results for each term in comma-separated values (CSV) files, generating one CSV file for each term. For energy, conversion, and generation terms, it also provides results for each vertical level, enabling two-dimensional analysis. Additionally, the program can be configured by the user to automatically generate plots. These plots include spatial representations of the computational domain defined for the analysis, time series of the system's central minimum pressure and relative vorticity at 850 hPa, time series of each LEC term, Hovmöller diagrams, and complete energy flux boxes for the LEC.


# External Libraries Used

The NetCDF data is handled by Xarray [@hoyer2017xarray], with Dask optimizing the handling of large datasets [@daniel2019data]. This combination allows for efficient slicing of large data chunks, better memory management, and improved performance when dealing with extensive data in NetCDF files. Pandas DataFrames [@reback2020pandas] are used for storing output results and saving them as CSV files. Numpy is utilized for angle conversion and computation procedures [@harris2020array]. For handling meteorological constants and variable units, the program uses the MetPy package [@may2022metpy].


# References