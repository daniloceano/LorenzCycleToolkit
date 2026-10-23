Results
=======

The `LorenzCycleToolkit` generates a directory structure to store the results of the Lorenz Energy Cycle (LEC) calculations. This structure includes CSV files containing the calculated terms, a log file, and optional figures if the plotting flag is enabled. Here's an example of the results directory structure:

.. code-block:: bash

    LEC_Results/
    ├── testdata_ERA5_fixed
    │   ├── Figures
    │   │   ├── LEC
    │   │   ├── boxplots
    │   │   ├── hovmollers
    │   │   ├── timeseries
    │   ├── log.testdata_ERA5
    │   ├── results_vertical_levels
    │   │   ├── Ae_level.csv
    │   │   ├── Az_level.csv
    │   │   ├── Ca_1_level.csv
    │   │   ├── Ck_1_level.csv
    │   │   ├── C_overturning_level.csv
    │   │   ├── Cz_1_level.csv
    │   │   ├── Ge_level.csv
    │   │   ├── Gz_level.csv
    │   │   ├── Ke_level.csv
    │   │   ├── Kz_level.csv
    │   │   ├── M_level.csv
    │   ├── domain_selection.csv
    │   ├── lec_metadata.json
    │   ├── run_status.json
    │   ├── testdata_ERA5_fixed_results.csv

Each subdirectory within the `LEC_Results` directory corresponds to a different dataset or analysis run. The structure of each subdirectory includes:

- **Figures/**: Contains figures generated when the ``-p`` flag is used.
  
  - **LEC/**: Contains LEC diagrams for each time step.
  - **boxplots/**: Contains boxplots of the LEC terms for each time step.
  - **hovmollers/**: Contains Hovmöller diagrams of the LEC terms.
  - **timeseries/**: Contains time series plots of the LEC terms.
  
- **log.*** (log file): Contains detailed information about the processing and calculations performed during the analysis.

- **run_status.json**: Status of the run: whether the computation completed and whether each requested product was generated (see `Run Status File`_).

- **lec_metadata.json**: What the computation used: input file and its SHA-256, namelist, times, the retained pressure layer, the units of every column of the main results file and, for every per-level file, its units, its relation to the column value and its layout (see `Scientific Metadata File`_).

- **domain_selection.csv**: The box actually selected at each time step, compared with the box requested and with the extent of the input data (see `Domain Selection File`_).


- **results_vertical_levels/**: One CSV file per term with its vertical profile (rows: times; columns: pressure levels in Pa). The header lists the levels actually used (the levels common to all time steps, in every framework) and each value is written under its own level. The files do **not** all hold the same kind of quantity and do **not** share the units of the main results file: see `Per-Level Files: Quantities and Units`_ before plotting, comparing or integrating them.

  - ``Az_level.csv``, ``Ae_level.csv``: zonal and eddy available potential energy per unit pressure.
  - ``Kz_level.csv``, ``Ke_level.csv``: area means of the squared zonal-mean and eddy wind components, ⟨[u]²+[v]²⟩ and ⟨u*²+v*²⟩. They equal **twice** the kinetic energy per unit mass and are written **without** the factor 1/(2g) of the column integral.
  - ``Cz_level.csv``, ``Ca_level.csv``, ``Ce_level.csv``: conversions per unit pressure.
  - ``Ck_level.csv``: the Ck conversion per unit pressure multiplied by g (written **without** the factor 1/g of the column integral).
  - ``Gz_level.csv``, ``Ge_level.csv``: generation of Az and Ae per unit pressure.
  - ``C_overturning_level.csv``: domain-mean overturning pressure work per unit pressure. A diagnostic; it does not enter the budget.
  - ``M_level.csv``: mass-continuity residual M(p), in s⁻¹, which sets the numerical noise floor of the budget.
  - Split terms: ``Ca_1_level.csv`` and ``Ca_2_level.csv`` (horizontal and vertical parts of the Ca integrand), ``Ck_1_level.csv`` … ``Ck_5_level.csv`` (the five parts of the Ck integrand, multiplied by g), ``Cz_1_level.csv``/``Ce_1_level.csv`` (the pressure factor R/(p g)) and ``Cz_2_level.csv``/``Ce_2_level.csv`` (the omega–temperature covariances). In the fixed framework ``Cz_1`` and ``Ce_1`` do not depend on time and are written as two rows (the levels, then the factor).

  **Note**: Hovmöller diagrams and boxplots are produced for the energy, conversion and generation groups only, so ``C_overturning_level.csv`` and ``M_level.csv`` are archived for analysis but are not plotted at vertical levels. Their column integrals do appear in the time series figures.

Main Results File
~~~~~~~~~~~~~~~~~

The main CSV file (e.g., ``testdata_ERA5_fixed_results.csv``) contains, for each time step, the LEC terms integrated over the **retained pressure layer**: the pressure levels of the input that have complete data inside the box, from the top level to the bottom level of the input (for the ERA5 sample, 700 to 1000 hPa). It is not an integral over the entire atmospheric column. The retained levels, the top and bottom pressures and any excluded level are recorded in ``lec_metadata.json`` (``vertical``). This file provides a time series of the energy cycle and is the primary output for analysis.

Levels with a pressure lower than 10 hPa (1000 Pa) are discarded when the data are read; 10 hPa itself is kept. When the input has such levels, the terminal and the log give a warning with the discarded levels, the input grid and the grid after this cut, and ``vertical.top_cut`` in ``lec_metadata.json`` records the same.

In every framework (``-f``, ``-t``, ``--choose``) the same levels are used at every time step. A level of the input grid whose data are incomplete inside the box of **any** time step (in the fixed framework the box is the same at every time step) is excluded from **all** of them, so every column value, tendency and residual covers the same layer. Nothing is interpolated and the input file is not modified. Only levels of the input grid can be excluded: a deliberately sparse grid (for example 1000, 850, 500 and 250 hPa) with complete data is used as it is, whereas a level that is in the grid but whose data are incomplete, even entirely missing, is excluded. The terminal and the log name the framework and show the input grid, the levels complete at each time step, the excluded levels with the time steps that caused each exclusion, and the layer used. Excluding levels at the top or at the bottom only makes the layer thinner; the run stops with an error, before any result is written, when fewer than two levels remain or when an excluded level lies between the top and the bottom of the levels kept (an internal gap): the toolkit does not integrate across missing data. Two levels are only the numerical minimum for the column integral and the vertical derivatives to exist: they do not make the result representative of the atmospheric column, and with few levels the static stability and the vertical derivatives at the top and bottom of the layer are coarse. With ``--choose`` the boxes of all time steps are drawn first and the computation starts after the last one; with an internal gap the run stops after the drawing. Check ``vertical.harmonization`` in ``lec_metadata.json`` to see which layer was used.

The first column holds the time; its header cell is empty (the header line starts with a comma).

**File Columns**:

- **Time**: Timestamp of each analysis step
- **Az, Ae, Kz, Ke**: Energy reservoirs (J/m²), integrated over all pressure levels:
  
  - ``Az``: Zonal Available Potential Energy
  - ``Ae``: Eddy Available Potential Energy
  - ``Kz``: Zonal Kinetic Energy
  - ``Ke``: Eddy Kinetic Energy

- **Ca, Ce, Cz, Ck**: Conversion terms (W/m²), representing energy transformations:
  
  - ``Ca``: Conversion between Az and Ae
  - ``Ce``: Conversion from Ae to Ke
  - ``Cz``: Conversion from Az to Kz
  - ``Ck``: Conversion from Ke to Kz

- **BAz, BAe, BKz, BKe**: Boundary terms (W/m²), energy transport across domain boundaries (only for regional analyses)

- **BΦZ, BΦE**: Boundary pressure work terms (W/m²), always exported

- **C_overturning**: Domain-mean overturning pressure work (W/m²). A diagnostic; it does **not** enter ``RGz`` or ``RKz``.

- **M**: Mass-continuity residual (kg·m⁻²·s⁻¹). Sets the numerical noise floor of the budget; it does **not** enter the residuals.

- **∂Az/∂t, ∂Ae/∂t, ∂Kz/∂t, ∂Ke/∂t (finite diff.)**: Tendencies of the energy reservoirs (W/m², i.e. J m⁻² s⁻¹), from ``numpy.gradient`` over the elapsed time in seconds: second-order centred differences at interior times (non-uniform spacing allowed) and first-order one-sided differences at the first and last times.

- **Generation and Dissipation terms** (W/m²):
  
  When using ``-r/--residuals`` flag:
  
  - ``RGz, RGe``: Residual generation terms (include generation + numerical errors)
  - ``RKz, RKe``: Residual dissipation terms (include boundary flux + dissipation + numerical errors)
  
  When **not** using ``-r`` flag:
  
  - ``Dz, De``: Dissipation terms (computed from friction)

  ``Gz, Ge`` (generation, computed from diabatic heating) are exported in both
  modes.

**Usage Example**:

.. code-block:: python

   import pandas as pd
   import matplotlib.pyplot as plt
   
   # Load results
   df = pd.read_csv('LEC_Results/testdata_ERA5_fixed/testdata_ERA5_fixed_results.csv', 
                    parse_dates=['Time'], index_col='Time')
   
   # Plot energy reservoirs
   fig, ax = plt.subplots(figsize=(12, 6))
   df[['Az', 'Ae', 'Kz', 'Ke']].plot(ax=ax)
   ax.set_ylabel('Energy (J/m²)')
   ax.set_title('Lorenz Energy Cycle - Energy Reservoirs')
   plt.show()

.. _per-level-units:

Per-Level Files: Quantities and Units
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The column value of a term in the main results file is ``factor × ∫ profile dp``, with the trapezoidal rule over the pressure levels of the per-level file (pressure in Pa). The factor is **not the same for every term**: most profiles already contain every physical factor, but ``Kz``/``Ke`` are written without 1/(2g) and ``Ck`` without 1/g. Notation: [x] zonal mean, x* = x − [x], ⟨x⟩ area mean, [x]″ = [x] − ⟨[x]⟩, σ static stability, g gravitational acceleration, Q diabatic heating rate.

.. Machine-checked against src/utils/output_spec.py (tests/test_docs_units.py).

.. csv-table::
   :header: "File", "Quantity written", "Units", "Column term", "Column units", "Factor"
   :widths: 10, 40, 12, 10, 10, 8

   "Az", "⟨[T]″²⟩/(2σ)", "J/m**2/Pa", "Az", "J/m**2", "1"
   "Ae", "⟨T*²⟩/(2σ)", "J/m**2/Pa", "Ae", "J/m**2", "1"
   "Kz", "⟨[u]²+[v]²⟩ (twice the zonal kinetic energy per unit mass)", "m**2/s**2", "Kz", "J/m**2", "1/(2g)"
   "Ke", "⟨u*²+v*²⟩ (twice the eddy kinetic energy per unit mass)", "m**2/s**2", "Ke", "J/m**2", "1/(2g)"
   "Cz", "−Cz_1·Cz_2 (positive: Az to Kz)", "W/m**2/Pa", "Cz", "W/m**2", "1"
   "Ca", "−(Ca_1+Ca_2) (positive: Az to Ae)", "W/m**2/Pa", "Ca", "W/m**2", "1"
   "Ce", "−Ce_1·Ce_2 (positive: Ae to Ke)", "W/m**2/Pa", "Ce", "W/m**2", "1"
   "Ck", "Ck_1+…+Ck_5, g times Ck per unit pressure (positive: Ke to Kz)", "W/kg", "Ck", "W/m**2", "1/g"
   "C_overturning", "−⟨ω⟩α/g, α = R⟨T⟩/p", "W/m**2/Pa", "C_overturning", "W/m**2", "1"
   "Gz", "⟨[Q]″[T]″⟩/(c_p σ)", "W/m**2/Pa", "Gz", "W/m**2", "1"
   "Ge", "⟨Q*T*⟩/(c_p σ)", "W/m**2/Pa", "Ge", "W/m**2", "1"
   "M", "⟨D_EW⟩+⟨D_NS⟩+∂⟨ω⟩/∂p", "1/s", "M", "kg/m**2/s", "special"
   "Ca_1", "⟨v*T* ∂[T]″/∂φ⟩/(aσ)", "W/m**2/Pa", "", "", ""
   "Ca_2", "⟨ω*T* ∂[T]″/∂p⟩/σ", "W/m**2/Pa", "", "", ""
   "Ck_1", "⟨u*v* cosφ ∂([u]/cosφ)/∂φ⟩/a", "W/kg", "", "", ""
   "Ck_2", "⟨v*² ∂[v]/∂φ⟩/a", "W/kg", "", "", ""
   "Ck_3", "⟨tanφ u*² [v]⟩/a", "W/kg", "", "", ""
   "Ck_4", "⟨ω*u* ∂[u]/∂p⟩", "W/kg", "", "", ""
   "Ck_5", "⟨ω*v* ∂[v]/∂p⟩", "W/kg", "", "", ""
   "Cz_1", "R/(p g)", "m/Pa/K", "", "", ""
   "Ce_1", "R/(p g)", "m/Pa/K", "", "", ""
   "Cz_2", "⟨[ω]″[T]″⟩", "K*Pa/s", "", "", ""
   "Ce_2", "⟨ω*T*⟩", "K*Pa/s", "", "", ""

Units are written as Pint strings (``J/m**2/Pa`` is J m⁻² Pa⁻¹). The split terms combine as Ca = −(Ca_1 + Ca_2), Ck = Ck_1 + … + Ck_5, Cz = −Cz_1·Cz_2 and Ce = −Ce_1·Ce_2. For ``M`` the column is (∫(⟨D_EW⟩ + ⟨D_NS⟩) dp + ω_bottom − ω_top)/g: the vertical part uses the exact boundary difference of ω, so ∫ M dp / g does not reproduce the column value. The sign conventions are those of the budget residuals (e.g. RKz = ∂Kz/∂t − Cz − Ck − BKz).

Consequences for users of the per-level files:

- ``Kz_level.csv`` and ``Ke_level.csv`` are not energies per unit area or per unit mass. Integrating them in pressure and dividing by g gives **twice** the column value; divide by 2g.
- ``Ck_level.csv`` is in W/kg, not W m⁻² Pa⁻¹: divide by g before comparing it level by level with Cz, Ca or Ce.
- Profiles of different terms are comparable level by level only when their units in the table are the same.

The same information is available in Python as ``src.utils.output_spec.VERTICAL_PROFILE_SPEC``. ``Dz_level.csv`` and ``De_level.csv``, written only by the direct dissipation pathway (without ``-r``, with friction fields in the input), are **not verified**: their units, level and column relation are under review and are not documented here.

Vertical-Level Figures
~~~~~~~~~~~~~~~~~~~~~~

With ``-p``, three kinds of figures are drawn from the per-level files, each for the energy (Az, Ae, Kz, Ke), conversion (Cz, Ca, Ck, Ce) and generation (Gz, Ge) groups:

- ``boxplots/boxplot_vertical_levels_<group>.png``: one box per pressure level, summarising the values at that level over all times. The statistics are **unweighted** by the interval between times.
- ``boxplots/vertical_timeseries_<group>.png``: one box per time, summarising the values at all pressure levels at that time. The statistics are **unweighted** by layer thickness: a box is not a vertical integral or a mass-weighted vertical mean, and its median is not related to the column value.
- ``hovmollers/hovmoller_<group>.png``: time–pressure diagrams of the profiles.

Each panel (axis or colour bar) is labelled with the quantity and units of its own term, as in the table above; the panels of one figure can have different units, which a note under the figure recalls. The plotted values are those of the CSV files, unconverted. Pressure is shown in hPa on the axis labels and ticks only; the files keep it in Pa.

Figures Produced by Earlier Versions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Versions up to and including 2.0.0 drew the vertical-level figures with the units of the column values:

- **Affected files**: ``boxplots/boxplot_vertical_levels_{energy,conversion,generation_dissipation}.png``, ``boxplots/vertical_timeseries_{…}.png`` and ``hovmollers/hovmoller_{…}.png``. The time series, the boxplots of the column terms (``boxplots/boxplot_terms_*.png``), the LEC diagrams and the maps were not affected.
- **What was wrong**: the labels. Each figure carried one group label for all its panels: "Energy (J·m⁻²)" for the energy group, "Conversion (W·m⁻²)" for the conversion group and "Generation/Dissipation (W·m⁻²)" for the generation group. In the boxplots by level and by time this label was shown on the panels of the left column only; in the Hovmöller diagrams every colour bar carried it, and the pressure axis had no label. The plotted values were the per-level values described in the table above (for Kz and Ke, ⟨[u]²+[v]²⟩ and ⟨u*²+v*²⟩ in m² s⁻²; for Ck, W/kg). These values are correct as per-level quantities of the version that produced them, but they are not in the labelled units.
- **Missing level**: ``boxplot_vertical_levels_*.png`` omitted the first pressure level of the file (the level of lowest pressure, e.g. 700 hPa in the sample data). The Hovmöller diagrams and the boxplots by time included every level.
- **How to recognise an affected figure**: group labels with column units (J·m⁻² or W·m⁻²) on a vertical-level figure; pressure ticks in Pa (e.g. "85000.0"); one box fewer than the number of levels in ``results_vertical_levels``.
- **How to regenerate**: run the toolkit again with ``-p`` on the same input with this version; the CSV files of the computation are unchanged by this correction, so only the figures differ. Note that numerical results of versions before 2.0.0 differ for other reasons (see the 2.0.0 entry of the changelog).

Log File
~~~~~~~~

The log file (``log.<filename>``) provides detailed information about the processing and calculations performed during the analysis. It includes:

- Input data information and validation checks
- Configuration used for the analysis
- Processing steps and timing information
- Any errors or warnings encountered during execution
- Download progress (when using ``--cdsapi``)

Run Status File
~~~~~~~~~~~~~~~

``run_status.json`` is written by every command-line run, first with the state ``running`` before the data are read and then with the final state. It is an additional file: the CSV files are not changed. The file is replaced in one step (temporary file and rename), so a run that is killed leaves ``running``, never a complete-looking status. Schema version 1:

- ``schema_version`` (``1``), ``lectk_version``, ``command`` (arguments, method, input file), ``started_utc``, ``finished_utc``, ``log_file``.
- ``status``: ``success`` (computation complete, every requested product generated), ``partial`` (computation complete, at least one requested figure failed), ``failed`` (the computation or the run failed) or ``running``.
- ``exit_code``: the exit status of the command line (0, 3 or 1; ``null`` while running or after an interruption). ``allow_partial_plots`` and ``plots_requested`` record the options.
- ``scientific``: ``status`` of the computation and the ``error`` (type and message) when it failed. ``interrupted`` and ``error`` describe a run stopped by an exception.
- ``summary``: number of products in each state.
- ``products``: one entry per product with ``name``, ``kind`` (``scientific``, ``figure`` or ``diagnostic``), ``required``, ``status``, ``expected`` (files the product must write), ``written`` (files written by this run), ``missing``, ``reason``, ``error``, ``failures`` (one per failed call, e.g. per time step), ``depends_on``, ``consequences`` and ``notes`` (caveats that do not change the state, e.g. that the cyclophaser periods are a technical result only, or that with fewer than 20 steps cyclophaser 1.7.3 replaces the whole band-pass series with its low-pass series). ``offline_maps`` records the ``--offline-maps`` option. Paths are relative to the results directory.

Product states: ``success`` (requested and every expected file written by this run), ``failed`` (requested and not, or not completely, generated), ``unavailable`` (not applicable to the data for a known reason checked before the call; not a failure: currently the cyclophaser periods of series shorter than the tested minimum, and the products that need them), ``not_requested`` (optional product without ``-p``). ``pending`` appears only while the run is in progress.

Required products: ``results_csv``, ``vertical_level_csvs``, in the moving framework ``trackfile``, and ``lec_metadata`` and ``domain_selection``. Optional products with ``-p`` — fixed framework: ``timeseries``, ``box_limits_map``, ``boxplots``, ``hovmoller``, ``lec_cycle``; moving framework: ``domain_maps`` (one map per time step), ``track_map``, ``zeta_hgt_timeseries``, ``timeseries``, ``hovmoller``, ``boxplots``, ``periods`` (cyclophaser life-cycle periods and ``periods.csv``), ``lec_cycle``, ``lec_periods`` and ``lps`` (these two need ``periods``).

Example (moving framework, ERA5 sample: 5 hourly steps, fewer than the 8 needed for the periods):

.. code-block:: json

   {
     "schema_version": 1,
     "status": "success",
     "exit_code": 0,
     "offline_maps": false,
     "scientific": {"status": "success", "error": null},
     "products": [
       {"name": "results_csv", "kind": "scientific", "required": true, "status": "success",
        "written": ["testdata_ERA5_track_results.csv"]},
       {"name": "periods", "kind": "diagnostic", "required": false, "status": "unavailable",
        "reason": "5 usable time steps < minimum 8 for the high_res options of cyclophaser 1.7.3 (Savitzky-Golay window not larger than the polynomial order); cyclophaser not called",
        "consequences": ["lec_periods", "lps"]},
       {"name": "lps", "kind": "figure", "required": false, "status": "unavailable",
        "depends_on": ["periods"],
        "reason": "not generated: requires 'periods', which is unavailable (5 usable time steps < minimum 8 ...)"}
     ]
   }

(Abridged: the file lists every product and every field.)

``previous_results`` is ``null`` unless the run reused, with ``--overwrite``, a results directory that already held files (see :doc:`usage`). Then it gives ``n_files_before`` and ``stale_files``: files of the earlier run that this run did not rewrite. They are left in place, untouched, and are **not** results of this run (for example the figures of an earlier run with ``-p`` when the new run has no ``-p``).

Scientific Metadata File
~~~~~~~~~~~~~~~~~~~~~~~~

``lec_metadata.json`` is written by every run once the computation is complete, before the figures. It is an additional file: no CSV file changes. Schema version 1:

- ``lectk_version``, ``generated_utc``, ``command`` (method, ``-r``, ``--outname``, box-limits or track file, ``-z``, ``-m``).
- ``input``: file name, ``sha256`` and size of the input file, and the namelist (variable name and declared units of each row).
- ``time``: number of times, first and last time, and the distinct time steps in seconds.
- ``horizontal``: the extent of the input data, the name of ``domain_selection.csv`` and ``boxes_cut_at_input_edge``, the time steps (``all`` in the fixed framework) whose requested box was cut at an edge of the data, with the edges.

- ``vertical``: the input levels, the levels used (``levels_used_Pa``), the excluded levels (``dropped_levels_Pa``), ``p_top_Pa`` and ``p_bottom_Pa`` of the retained layer, and the levels in the header of the per-level files. In the moving framework ``by_time`` lists, for each time step, the levels used and the excluded levels. In every framework the levels are common to all time steps (``levels_used_vary_with_time`` is ``false``). ``harmonization`` records the framework (``mode``), the rule, the internal-gap policy (``internal_gap_policy``), the minimum number of levels (``minimum_levels``), the number of time steps (``time_steps``), the levels available, used and excluded, the time steps that caused each exclusion (``excluded_because_incomplete_at``), the levels complete at each time step (``complete_levels_by_time``) and the input and effective layers (``input_layer_Pa``, ``effective_layer_Pa``). ``top_cut`` (DEC-09) records the rule of the cut at 10 hPa, the grid of the input before the cut (``input_grid_before_cut_Pa``) and the levels it discarded (``discarded_levels_Pa``); ``input_levels_Pa`` is the grid after the cut.
- ``main_results``: the file, the empty header of the time column and the units of every column (``J/m**2`` for Az, Ae, Kz, Ke; ``kg/m**2/s`` for M; ``W/m**2`` for the other terms, the tendencies and the residuals). ``Dz`` and ``De`` are marked ``NOT VERIFIED``.
- ``vertical_level_files``: for each of the 23 per-level files, its ``units``, ``role``, ``column_term``, ``column_units`` and ``column_factor`` (as in `Per-Level Files: Quantities and Units`_) and its ``layout``: ``time_by_level`` (one row per time, one column per level) or ``time_invariant_profile`` (``Cz_1`` and ``Ce_1`` in the fixed framework: a row ``level`` with the levels, then a row with the values, which do not depend on time).

Reading a per-level file in both layouts as a time × level table (example; not part of the toolkit's API):

.. code-block:: python

   import json
   import pandas as pd

   meta = json.load(open("LEC_Results/testdata_ERA5_fixed/lec_metadata.json"))
   times = pd.read_csv("LEC_Results/testdata_ERA5_fixed/testdata_ERA5_fixed_results.csv",
                       index_col=0, parse_dates=True).index
   entry = {e["term"]: e for e in meta["vertical_level_files"]}["Cz_1"]
   path = "LEC_Results/testdata_ERA5_fixed/" + entry["file"]
   if entry["layout"] == "time_invariant_profile":
       rows = pd.read_csv(path, index_col=0)          # rows "level" and "Cz_1"
       profile = rows.iloc[-1].astype(float)          # values per level (Pa)
       table = pd.DataFrame([profile.values] * len(times), index=times,
                            columns=[float(c) for c in rows.columns])
   else:
       table = pd.read_csv(path, index_col=0, parse_dates=True)

Domain Selection File
~~~~~~~~~~~~~~~~~~~~~

``domain_selection.csv`` has one row per time step (moving framework) or one row with time ``all`` (fixed framework). Limits are in degrees, in the longitude convention of the input data after preprocessing (for example −60 rather than 300):

- ``requested_west/east/south/north``: the box limits given to the computation (the box-limits file, or the track centre ± width/2 and ± length/2).
- ``selected_west/east/south/north``: the first and last longitude and latitude of the fields actually used. Each limit is the grid point nearest to the request.
- ``offset_*``: selected minus requested.
- ``n_lon``, ``n_lat``: number of grid points of the box.
- ``input_west/east/south/north``: extent of the input data; ``*_at_input_edge``: whether the selected edge is the edge of the data.
- ``requested_outside_input``: the edges requested beyond the data (``none`` otherwise). Such a box is cut at the edge of the data without an error, in every framework, and the computation goes on with the box actually used; the terminal and the log give a warning, for each time step concerned, with the limits requested, the limits of the fields actually used and the extent of the data.
- ``fields_share_selection``: whether temperature, wind, omega and geopotential use the same grid points.

Use the ``-v/--verbosity`` flag for more detailed debug-level logging.