# Lorenz-Cycle

## What is the Lorenz Energy Cycle?

The Lorenz Energy Cycle (LEC) is a methodology for estimating the energy in the Atmosphere designed by the mathematician and meteorologist Edward Lorenz, in 1965. On it, the energy is partitioned into zonal and eddy components of Kinectic Energy (Kz and Ke, respectively) and Available Potential Energy (Az and Ae, respectively). It is also computed the terms related to the conversion between each form of energy (Ca, Ce, Cz and Ck) and the terms related to generation/dissipation of the energy contents (Gz, Ge, Dz and De). Lorenz initial work computed the energetics for the whole globe, but later developments adapted the methodology for computing the energetics of a limited area in the atmosphere, thus requiring the computation of terms related to the transport of energy thorugh the boundaries of the limited area (BAe, BAz, BKe and BKz). Due to difficulties in accessing the frictions terms required for computing the dissipation terms, both disspiation and generation terms are often computed as residuals from the budgets equations:

![image](https://user-images.githubusercontent.com/56005607/210858922-1d29f3b9-2446-422e-87c4-2dc5c4ced361.png)

Where:

![image](https://user-images.githubusercontent.com/56005607/210859000-07af27c1-0295-4c48-bf7c-1432100bdf60.png)

Where ε represents errors numeric errors. The complete cycle for a given period of time can be ilustred as follows:


![LEC_example](https://user-images.githubusercontent.com/56005607/210855570-d3272989-8871-4a20-996f-e373f73934c5.png)


## What does the program do?

The Lorenz-Cycle program is designed for computing the LEC for a specific region on the atmosphere. It requires a netCDF file containing wind (zonal, meridional and vertical components), air temperature and geopotential data and can be run using two distinct frameworks:

1. Fixed framewrok, where the domain is fixed in time by the file: inputs/box_limits; 

2. Moving framework, where the domain can follow a pertubation on the atmosphere. 

More details on running the program are provided bellow.

# Usage

## Fixed framework

**Important!** Before running the program be sure to check the [Flags](#flags) section.

First of all, you'll need a NetCDF file containing the following variables: zonal (u) and meridional (v) wind components, vertical wind speed (omega), air temperature and geopotential or geopotential height. The data is required to follow a pressure-level vertical levels.

Then, it is to necessary to delimite the computational domain using the csv text file [box_limits](inputs/box_limits), that should look like this:

![image](https://user-images.githubusercontent.com/56005607/206709581-34ebe0a7-ff45-4bd4-86e0-8cce8dde91ea.png)

Afterwrds, it is required to specify in the [fvars](inputs/fvars) file the how the variables are named in the NetCDF file and which units are being used. It will look like this:  

![image](https://user-images.githubusercontent.com/56005607/210861069-1c899cc8-860a-4212-bd44-118e308db9bd.png)

Note that you should only change the colums corresponding to "Variable" and "Units". Modying the first column will make the program unable to look for the variables in the NetCDF file.

Then, from the [source code folder](src) you can run, for example:

```
python lorenz-cycle.py path/to/infile.nc -r -e
```

## Moving freamework

**Important!** Before running the program be sure to check the [Flags](#flags) section.

#### Using a pre-defined domain

As in the fixed framework, the first step is to have a NetCDF file containing the following variables: zonal (u) and meridional (v) wind components, vertical wind speed (omega), air temperature and geopotential or geopotential height. The data is required to follow a pressure-level vertical levels.

Now, instead of delimiting the computational domain, it is required a [track file](inputs/track) containing the central position of the system of interest in different time steps. The program will then create, for each time step, a box with 15°x15° around this central point to compute the energetics terms. Optionally, the user can add the length and width columns to the track file, for using a distinct domain size than the default. The track file should look like this:

![image](https://user-images.githubusercontent.com/56005607/206721056-61fa32ce-aa5d-4f16-af28-c46ac2a9bf88.png)

As in the fixed framework, it is required to specify in the [fvars](inputs/fvars) file the how the variables are named in the NetCDF file and which units are being used. See above. 

Then, from the [source code folder](src) you can run, for example:

```
python lorenz-cycle.py path/to/infile.nc -r -t
```

#### Interactively choosing the domain

Instead of using a pre-defined domain, in this method a pop-up window will appear for the user displaying the world map and the vorticity field. For using such method, the user can run, for example:

```
python lorenz-cycle.py path/to/infile.nc -r -c
```

Firstly, the user must select an area for slicing the global data, which improves visualization and processing time. To select the area, the user needs to click twice on the screen. The first click selects the top-left corner of the area, and the second click selects the bottom-right corner. The order of these clicks does not matter.

![image](https://user-images.githubusercontent.com/56005607/214921907-e19d0024-08dc-4475-ab65-c953e04e7859.png)


After this, the pop-up window will display the vorticity data for the selected area only. For each time step, the user will be prompted to select a computational area by clicking on the screen.

![image](https://user-images.githubusercontent.com/56005607/214922008-5b7c094f-c160-4415-a528-07cc58730827.png)

After finishing the computations, it will be created a [track](#auxiliary files)  file containing the domain central latitude and longitude, its length and width and the minimum vorticity and maximum windspeed inside the domain.

## File naming system

Alhtough it is not necessary to, it is advisible to use a file naming convenction such as "object-of-study_file-source.nc". For example, if you wish to compute the energetics of the Katrina Hurricane, using data from the ERA5 model, the file name should be "Katrina_ERA5.nc". May you wish to use compund names for the system, separete then using a "-", for example: "Cyclone-20100101_NCEP-R2.nc".


## Flags

- -r, --residuals

The default behaviour for computing the energetics, and the one intended by the works of Lorenz, was to compute the dissipation terms directly by using wind stress fields. However, for many of the reanalysis data, those variables are not available for all levels of the atmosphere. Therefore, the dissipation terms might be estimated as residuals from the budget equations (see the [first section](## What is the Lorenz Energy Cycle?)).

- -m, --moving

Compute the LEC using a fixed domain (see [fixed framework](#fixed framework))

- -t, --track

Compute the LEC using a domain that follows the system (see [moving framework](#moving framework))

- -c, --choose

For each time step, interactively select the domain (see [Moving domain](#moving domain))

- -o, --outname

Choose a name for saving the results (optional)

- -g, --geopotential

Use this flag when instead of Geopotential Height, Geopotential data is provided. It is required that the FVars file is adjusted accordingly, for example:

![image](https://user-images.githubusercontent.com/56005607/210860966-713243c8-7447-4661-a33d-a988ab1055cf.png)




