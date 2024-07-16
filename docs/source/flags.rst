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
