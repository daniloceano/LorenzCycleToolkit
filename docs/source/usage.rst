Usage
=====

Fixed Framework
---------------

**Prerequisites**
- NetCDF file with necessary variables (u, v, omega, air temperature, geopotential or geopotential height).
- `inputs/box_limits` file defining the domain.
- `inputs/namelist` file specifying variable names and units.

**Command**
::

   python lorenzcycletoolkit.py path/to/infile.nc -r -f

Moving Framework
----------------

**Pre-defined Domain**
- Follow initial steps for the Fixed Framework.
- Use a `inputs/track_file` to define the system's center over time.

**Command**
::

   python lorenzcycletoolkit.py path/to/infile.nc -r -t

Interactive Domain Selection
----------------------------
- For interactive domain selection, run::

   python lorenzcycletoolkit.py path/to/infile.nc -r -c

Flags
=====

- `-r`, `--residuals`: Compute the Dissipation and Generation terms as residuals.
- `-f`, `--fixed`: Compute the energetics for a fixed domain specified by the 'box_limits' file.
- `-t`, `--track`: Define the domain using a track file.
- `-c`, `--choose`: Interactively select the domain for each time step.
- `-o`, `--outname`: Specify an output name for the results.
- `-z`, `--zeta`: Use the vorticity from the track file instead of computing it at 850 hPa.
- `-m`, `--mpas`: Specify this flag if working with MPAS-A data processed with MPAS-BR routines.
- `-p`, `--plots`: Generate plots.

Ensure that the provided NetCDF file and the `namelist` configuration align with the selected flags.
