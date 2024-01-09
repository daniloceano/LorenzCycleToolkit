# Lorenz Energy Cycle (LEC)

## Overview
The Lorenz Energy Cycle (LEC), introduced by Edward Lorenz in 1965, is an analytical framework used to estimate atmospheric energy. It categorizes energy into zonal and eddy components of Kinetic Energy (Kz and Ke, respectively) and Available Potential Energy (Az and Ae, respectively). The LEC also quantifies conversions between these forms (Ca, Ce, Cz, and Ck), along with generation and dissipation terms (Gz, Ge, Dz, and De). Originally developed for global energetics, the framework has been adapted for regional studies, incorporating calculations for energy transport across boundaries (BAz, BAe, BKz, BKe).

The LEC budget is described by the following equations:

$\frac{\partial A_Z}{\partial t} = -C_Z - C_A + G_Z + B A_Z$

$\frac{\partial K_Z}{\partial t} = -C_Z + C_K - D_Z + B K_Z + B \Phi_Z$

$\frac{\partial A_E}{\partial t} = C_A - C_E + G_E + B A_E$

$\frac{\partial K_E}{\partial t} = C_E - C_K - D_E + B K_E + B \Phi_E$

Due to the difficulty in measuring friction terms for dissipation, both dissipation and generation are often computed as residuals from the budget equations:

$RG_Z = G_Z + \varepsilon_{AZ}$

$RG_E = G_E + \varepsilon_{AE}$

$RK_Z = B \Phi_Z - D_Z + \varepsilon_{KZ}$

$RK_E = B \Phi_E - D_E + \varepsilon_{KE}$

Where ε represents numerical errors. The complete cycle, assuming all terms are positive, is depicted below:

<img src="https://github.com/daniloceano/lorenz-cycle/assets/56005607/d59eeb31-5cef-46ac-a841-1ba4170fafbd" width="350">

## Program Description

The Lorenz-Cycle program calculates the LEC for specific atmospheric regions. It requires a netCDF file containing wind components, air temperature, and geopotential or geopotential height data. The program offers two frameworks:

1. Fixed Framework: The domain is static, defined by inputs/box_limits.
2. Moving Framework: The domain moves with an atmospheric perturbation.

# Usage

## Fixed framework


**Prerequisites:** Check the [Flags](#flags) section.

- Ensure you have a netCDF file with necessary variables (u, v, omega, air temperature, geopotential or geopotential height) on pressure-level vertical levels.
- Define the domain in 'inputs/box_limits' as a CSV:

![image](https://user-images.githubusercontent.com/56005607/206709581-34ebe0a7-ff45-4bd4-86e0-8cce8dde91ea.png)

- Specify variable names and units in 'inputs/fvars':

![image](https://user-images.githubusercontent.com/56005607/210861069-1c899cc8-860a-4212-bd44-118e308db9bd.png)

- Execute from the source code directory:

```
python lorenz-cycle.py path/to/infile.nc -r -e
```

## Moving freamework

**Prerequisites:** Check the [Flags](#flags) section.

#### Pre-defined domain

- Follow the initial steps for the Fixed Framework.
- Use a [track file](inputs/track_file) to define the system's center over time. Optionally, add length and width columns to adjust the domain size.

![image](https://user-images.githubusercontent.com/56005607/206721056-61fa32ce-aa5d-4f16-af28-c46ac2a9bf88.png)

- Execute:

```
python lorenz-cycle.py path/to/infile.nc -r -t
```

#### Interactive Domain Selection

- For interactive map-based domain selection, run:

```
python lorenz-cycle.py path/to/infile.nc -r -c
```

- First, define a slice of global data to enhance visualization and processing speed.
- Select the computational area for each timestep:

![image](https://user-images.githubusercontent.com/56005607/214921907-e19d0024-08dc-4475-ab65-c953e04e7859.png)

- Each timestep prompts the user to select a computational area:

![image](https://user-images.githubusercontent.com/56005607/214922008-5b7c094f-c160-4415-a528-07cc58730827.png)

## File Naming Convention

It's recommended to use a naming convention like "subject_source.nc" for input files. For compound names, use a hyphen, e.g., "Cyclone-20100101_NCEP-R2.nc".

## Flags

- `-r`, `--residuals`: Compute the Dissipation and Generation terms as residuals.
- `-f`, `--fixed`: Compute the energetics for a fixed domain specified by the 'box_limits' file.
- `-t`, `--track`: Define the domain using a track file.
- `-c`, `--choose`: Interactively select the domain for each time step.
- `-o`, `--outname`: Specify an output name for the results.
- `-z`, `--zeta`: Use the vorticity from the track file instead of computing it at 850 hPa.
- `-g`, `--geopotential`: Use geopotential data instead of geopotential height.
- `-m`, `--mpas`: Specify this flag if working with MPAS-A data processed with MPAS-BR routines.
- `-p`, `--plots`: Generate plots.

Ensure that the provided NetCDF file and the `fvars` configuration align with the selected flags.

**Example of 'fvars' file with geopotential:**

![image](https://user-images.githubusercontent.com/56005607/210860966-713243c8-7447-4661-a33d-a988ab1055cf.png)




