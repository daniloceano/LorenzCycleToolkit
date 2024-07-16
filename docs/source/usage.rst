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
