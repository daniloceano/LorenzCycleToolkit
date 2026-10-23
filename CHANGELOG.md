# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Added

- **Warning for the levels above 10 hPa** (DEC-09): the levels with a
  pressure lower than 1000 Pa, which are discarded when the data are read, as
  before, are now reported in the terminal and in the log, with the input grid
  and the grid after the cut; `vertical.top_cut` in `lec_metadata.json`
  records the same. 10 hPa itself is kept. No warning when nothing is
  discarded; no change in the computation.
- **Warning for a box cut at the data edge** (DEC-11): a box requested beyond
  the input data is still cut at the data edge and the run goes on; the
  terminal and the log now warn, for each time step concerned and in every
  framework, with the limits requested, the limits of the fields actually
  used and the extent of the data; `horizontal.boxes_cut_at_input_edge` in
  `lec_metadata.json` lists the time steps and edges. No change in the
  computation.

- **`vertical.harmonization`** in `lec_metadata.json` (DEC-10 B, every
  framework): the framework, the internal-gap policy, the minimum number of
  levels, the levels available, used and excluded, the time steps that
  caused each exclusion, the levels complete at each time step and the input
  and effective layers.

- **`lec_metadata.json`** (E0-17, E0-18, E0-19, DEC-23): every run records the
  input file and its SHA-256, the namelist, the time axis, the retained
  pressure layer (levels used, excluded levels, top and bottom pressure), the
  units of every column of the main results file and, for every per-level
  file, its units, column relation and layout (`Cz_1`/`Ce_1` of the fixed
  framework are time-invariant profiles). Additional file; no CSV changes.
- **`domain_selection.csv`** (E0-22): the box actually selected at each time
  step (grid points nearest to the request), the offset from the requested
  box, the number of points and whether each edge is the edge of the input
  data. A box requested beyond the data, which is cut at the data edge
  without an error, now shows here. Additional file.
- **`--overwrite`** (DEC-08): reuse a results directory that already holds
  files; see the change below.

- **`run_status.json`** (DEC-02): every command-line run writes, next to its
  results, the state of the computation and of each requested product
  (`success`, `failed`, `unavailable`, `not_requested`), the cause of each
  failure and the exit status. It is written atomically, first as `running`.
  The CSV files are unchanged. Schema version 1, described in the results
  documentation.
- **`--offline-maps`**: draw maps from local Natural Earth data only and never
  download them; a map whose data are missing is reported as failed (partial
  output). Cartopy's default behaviour is unchanged without the option.
- **`--allow-partial-plots`** (DEC-02): exit with status 0 when the
  computation is complete but some requested figures failed; the run is still
  reported as partial.
- **`src/utils/output_spec.py`** (E0-14): `VERTICAL_PROFILE_SPEC` records,
  for each of the 23 per-level CSV files, the quantity written, its units,
  the column term and the factor in
  `column = factor * integral(profile dp)`; `SPLIT_TERM_RELATIONS` records how
  the split terms combine. The computation does not read it; no value
  changes. `Dz`/`De` are listed apart as not verified.

### Changed

- **Box plots of the column terms without notches**: `boxplot_terms_*.png`
  draw conventional boxes. With few time steps the notch interval (median
  ± 1.57 IQR/√n) reached beyond the quartiles and the outline of the box
  folded outward. Same data, quartiles, medians, whiskers, outliers, colours,
  order of the terms, units, axis limits and file names; only the outline of
  the boxes changes, and the median line spans the whole box. The box plots
  by level and by time were not notched and are unchanged.

- **One vertical control volume for every time step, in every framework
  (`-f`, `-t`, `--choose`)** (DEC-10, option B; DEC-10-L, option B; behaviour
  change): a pressure level of the input grid whose data are incomplete
  inside the box of any time step is excluded from all time steps. In the
  track and `--choose` modes it used to be dropped at that time step only, so
  consecutive column integrals covered different layers and the
  finite-difference tendencies and the residuals carried the change of the
  integration domain (with one level missing at one time step in the ERA5
  sample, ∂Kz/∂t changed by about −1000% and +550% at the neighbouring
  times). Runs without incomplete levels are unchanged. Fewer than two common
  levels, or an excluded level between the top and the bottom of the levels
  kept, stop the run with an error in every framework (two levels are a
  numerical minimum, not a guarantee of a representative column); the fixed
  framework used to integrate across such a gap, and failed with an
  `IndexError` with fewer than two levels. A level that is not in the input
  grid is not a gap: sparse grids with complete data are unchanged. With
  `--choose` the boxes of all time steps are drawn before the computation.
  The terminal and the log name the framework and show the input grid, the
  levels complete at each time step, the exclusions and the layer used.
- **Per-level files of the moving framework: every value sits under its own
  level** (E0-21a, defect fixed): when a level was dropped at one time step,
  the values of that row were written one column to the left, under the
  header of the level above (850 hPa under 700 hPa). The header lists the
  levels used, in every framework.
- **Runtime dependencies declared** (P-A1): `cmocean`, `cyclophaser==1.7.3`,
  `dask`, `lorenz-phase-space` and `shapely` are now in `pyproject.toml`;
  `pip install` of the package alone used to give a command line that failed
  with `No module named 'cmocean'`.

- **`-o/--outname` names the directory and every file** (E0-12, DEC-07,
  behaviour change): `LEC_Results/NAME_<method>/` with
  `NAME_<method>_results.csv`, `NAME_<method>_trackfile` and `log.NAME`, as
  the documentation already said. It used to rename only the main CSV of the
  fixed framework (`NAME.csv`, in the directory named after the input file)
  and was ignored by the moving framework. Invalid names (empty, `.`, `..`,
  with `/` or `\`) exit with status 2. Without `--outname` nothing changes.
- **Earlier results are protected** (DEC-08, behaviour change): a run whose
  results directory already contains files stops with exit status 2 before
  writing anything. With `--overwrite` it rewrites its files in place and
  deletes nothing; files of the earlier run that it does not produce are
  listed in `run_status.json` (`previous_results`).
- **Main results are integrals over the retained pressure layer** (E0-19,
  documentation): the documentation said "entire atmospheric column"; the
  values have always been integrals from the top to the bottom retained
  level of the input. No value changes.

- **Exit status 3 for partial output** (DEC-02, behaviour change): when the
  computation is complete but a requested figure could not be generated, the
  command now exits with status 3. It used to exit with status 0 when the
  error was only logged (fixed framework; cyclophaser in the moving
  framework) and with status 1 and a traceback when a moving-framework
  figure raised. Errors of the computation still exit with status 1 and are
  never reported as partial.
- **Figure failures are reported, not hidden** (E0-07): a figure counts as
  generated only if its files were written by the run. "🎉 Analysis
  complete!" and "🎉 All plots generated successfully!" are written only when
  every requested product was generated; otherwise the log ends with an
  "Analysis complete with PARTIAL output" warning naming the figures.
- **Moving framework: a map failure no longer stops the computation**
  (E0-05): the map drawn at each time step with `-p` failed the whole run
  before the CSV was written when Natural Earth data could not be obtained
  (e.g. offline). The failure is now recorded and the energetics are
  computed; the figures after the results are each attempted and recorded.
- **Cyclophaser errors are reported** (E0-06): `plot_periods` logged and
  hid them; the periods, and the Lorenz cycle diagrams of the periods and the
  Lorenz Phase Space that need them, are now reported as failed with the
  cyclophaser error.
- **Periods unavailable for short series** (E0-06): series with fewer usable
  time steps than cyclophaser 1.7.3 accepts with the toolkit options (8 for
  high-resolution data or vorticity from the track file, 3 otherwise) are not
  passed to cyclophaser; the periods and the products that need them are
  reported as `unavailable`, which is not a failure (the ERA5 sample with `-p`
  now exits with status 0). These are technical limits, recorded with a note
  that the phases are not physically validated. Another cyclophaser version,
  other options or a time step other than 1, 3 or 6 hours give a
  `CyclophaserCompatibilityError` instead. The cyclophaser options are
  unchanged.

- **Log handlers** (E0-09): `initialize_logging` now closes and removes the
  handlers of a previous run in the same Python process, so each run writes
  only to its own `log.<input>` file and console lines are no longer repeated.
  A single command-line run is unaffected.
- **Running without arguments** now prints the usage and exits with status 2
  (E0-11). It previously replaced the empty command line with built-in example
  arguments (`test.nc -t -r -p -v --cdsapi`), which created a results directory,
  read `inputs/track` and started a CDS API download. Explicit `--cdsapi` runs
  are unchanged.
- **Vertical-level figures: units of each term** (E0-15, N-07; figures only,
  no CSV changes). The boxplots by level, the boxplots by time
  (`vertical_timeseries_*.png`) and the Hovmöller diagrams labelled every
  panel with the units of the column values ("Energy (J·m⁻²)",
  "Conversion (W·m⁻²)", "Generation/Dissipation (W·m⁻²)"). The per-level
  files hold other quantities, and each panel and colour bar is now labelled
  with the quantity and units of its own term: J m⁻² Pa⁻¹ (Az, Ae),
  m² s⁻² (Kz, Ke: ⟨[u]²+[v]²⟩ and ⟨u*²+v*²⟩, twice the kinetic energy per
  unit mass), W m⁻² Pa⁻¹ (Cz, Ca, Ce, Gz, Ge) and W kg⁻¹ (Ck, g times its
  value per unit pressure). The plotted values are unchanged. Pressure is
  shown in hPa on the axis labels and ticks; the files keep Pa. A note under
  each figure says when panels have different units and that the boxes are
  unweighted statistics. The Hovmöller figures are wider.
- **Boxplots by level include the first pressure level** (E0-16): the level
  of lowest pressure (700 hPa in the samples) was dropped from every
  `boxplot_vertical_levels_*.png`.
- **Figures produced by earlier versions**: up to and including 2.0.0, the
  vertical-level figures carried the units of the column values and the
  boxplots by level omitted the first level. The plotted values were the
  per-level values of the CSV files, which are unchanged. Regenerating the
  figures with `-p` gives the corrected labels; see "Figures Produced by
  Earlier Versions" in the results documentation.

### Documentation

- **Per-level files** (E0-13, E0-14, E0-18): table of the quantity, units,
  column term and factor of every per-level file; `Kz_level.csv` and
  `Ke_level.csv` hold ⟨[u]²+[v]²⟩ and ⟨u*²+v*²⟩ without 1/(2g), `Ck_level.csv`
  is written without 1/g, and `M` has its own column relation. Units and
  finite-difference scheme of the tendencies ∂X/∂t (W/m²). Meaning of the
  boxplots (unweighted statistics) and how to recognise and regenerate
  figures of earlier versions.

### Tests

- `tests/test_cut_and_clipping_reports.py` (DEC-09, DEC-11): levels above
  10 hPa are reported and give the run of the input without them; no warning
  without them; 10 hPa is kept and 9.99 hPa is discarded; the same cut with
  levels in Pa; a box beyond the data is cut and reported, with the limits
  used, in the fixed, track and `--choose` modes; a box inside the data is not
  reported; a cut box computes with the box actually used.

- `tests/test_boxplot_terms.py`: every box of `boxplot_terms_*` is a
  rectangle between the quartiles, and its median, whiskers, caps and
  outliers equal those of `matplotlib.cbook.boxplot_stats`, for 5, 2 and 40
  time steps, outliers, a constant series and the two reference runs; the
  colours and the order of the terms.

- `tests/test_vertical_harmonization.py` (DEC-10 B, E0-21a, DEC-10-L B): in
  the fixed, track and `--choose` modes, a level missing at the top or at the
  bottom at one time step gives, byte for byte, the run of the input without
  that level (tendencies and residuals included); values under their own
  level; internal gap and fewer than two levels stop the run; a sparse grid
  is not a gap; the messages and the metadata name the framework; boxes of
  different sizes; `--choose` draws every box before computing.
- `tests/test_packaging.py` (P-A1): every imported third-party module is
  declared; the cyclophaser pin matches `cyclophaser_compat`.
- `tests/test_invariances.py`: homogeneity of each term in the winds
  (a = −1, 2), fixed = moving framework for one box, column = factor ×
  integral of the per-level profile.
- `tests/test_maps_offline.py` (P-L1): the runs whose maps fail are compared
  byte for byte with the same case run without figures on the same machine,
  and with the reference within its registered tolerances; the reference is
  no longer required to be byte-identical on another platform (it was
  produced on x86_64). A failed map must not write any file.

- **Offline maps** (DEC-03): the tests that draw maps require
  `CARTOPY_DATA_DIR` to point to a copy of the Natural Earth files listed,
  with their sha256, in `tests/data/naturalearth_manifest.tsv`, and fail
  explicitly otherwise; the shapefiles are not in the repository. A session
  guard blocks and records every network attempt and fails the run on an
  unexpected one. New tests cover maps with the complete copy and without
  any data (results unchanged, maps reported as failed).
- **Run status and exit status** (DEC-02, E0-05, E0-06, E0-07): unit tests of
  `run_status.json`; integration tests check the status of every product;
  injected figure failures, silent plotting functions and errors of the
  computation; the real command line in a subprocess for exit statuses 0, 3,
  0 with `--allow-partial-plots`, 1 and 2.
- **Run logs are checked after the whole session** (E0-09): the suite fails if
  any run log grows after its test has ended, holds another run, or has lines
  after "Analysis complete"; a pair of tests also fails if log handlers are
  inherited from an earlier test.
- **Integration tests check their outputs** (E0-08): required files, terms,
  time steps and finite values, and, for the ERA5 fixed and track cases, a
  cell-by-cell comparison with reference outputs (`tests/data/reference`) under
  the tolerances registered for the audit baseline. Negative controls in
  `tests/test_output_checks.py`.
- **Integration tests no longer touch `inputs/`** (E0-10): each test runs in a
  private temporary directory with copies of the templates it needs; the
  session fails if any file under `inputs/` or `samples/` changes.
- **Per-level convention and vertical-level figures** (E0-13 to E0-16,
  E0-18, N-07): units of every profile checked with Pint when it is written;
  profile-to-column factors and split-term combinations on the reference
  outputs; column relation of `M` against a closed form; labels of every
  panel, levels drawn and values handed to matplotlib against the CSV files;
  documentation table against `output_spec`.

## [2.0.0] - 2026-08-18

Numerical results change in this release. Diagnostics produced with earlier
versions are not directly comparable with diagnostics produced from 2.0.0
onwards; see the notes below for which terms are affected.

### Changed

- **`C_A`, meridional term**: the latitude derivative now acts on `T*` alone.
  It previously acted on `T* cos(phi)`, which introduced a spurious
  `cos(phi)` factor and an additional `-T* sin(phi)` contribution. The
  leading factor is now `1/(a sigma)` rather than `1/(2 a sigma)`, consistent
  with the reservoir definition `A_E = int T'^2 / (2 sigma) dp` used by the
  toolkit. `C_A` typically grows by a factor of 2.5 to 3.0.
- **`C_K`, fifth term**: the eddy momentum flux now acts on the vertical shear
  of the zonal-mean meridional wind, `d[v]/dp`, in place of `d[u]/dp`. The
  effect is case-dependent and can reverse the sign of `C_K`.
- **`BPhi_Z`**: the east and west faces are now diagnosed. The previous
  expression evaluated a product of two zonal means across those faces, which
  is identically zero. The north/south and top/bottom faces use the area
  departures `Phi*` and `omega*`.
- **`BPhi_E`**: all three faces now use eddy covariances.
- **Static stability floor**: the `sigma = 0.03` floor is retained as the
  default but is now explicit, reversible via `apply_floor=False`, and
  reported when it engages.
- **Residual plot labels**: `RG` and `RK` are labelled as composite residuals
  rather than as named physical processes.
- **Main results CSV schema**: the fixed framework now emits the columns
  `C_overturning`, `BPhi_Z`, `BPhi_E` and `M` alongside the existing terms, so
  the fixed and moving frameworks share one column set. Downstream scripts
  that read the CSV by column *position* must be updated; scripts that read by
  column *name* are unaffected. None of the four enters the residuals.
- **Input file resolution**: `inputs/namelist` and `inputs/box_limits` are no
  longer tracked by git. Only `inputs/namelist.default` and
  `inputs/box_limits.default` are distributed, and the toolkit falls back to
  the `.default` companion when the working copy is absent, so a clean clone
  runs out of the box. Copy the `.default` file to the name without the suffix
  to customise a run.
- **Zonal and area averages** are normalised by the same trapezoidal rule used
  to integrate the field, rather than by the analytic domain width. The two
  differ by `O(dphi^2)`, about `1.6e-4` on a 2.5 degree grid and `1.6e-6` on a
  0.25 degree grid, so most terms move by that amount. The average of a
  constant is now exact, which makes the departures from the zonal and area
  means satisfy their defining identity and makes every diagnosed flux
  independent of the arbitrary reference level of the geopotential. The static
  stability `sigma` and the two horizontal terms of `M` use that same
  normalisation, so a single area-averaging convention holds throughout.

### Added

- **`C_overturning`**: domain-mean overturning diagnostic. Exported for
  interpretation; it does not enter `RGz` or `RKz`.
- **`M`**: mass-continuity residual, exported as a vertical profile and as a
  column integral. It sets the numerical noise floor of the budget.
- Analytic test suite for the energy, conversion and boundary equations.
- Disabling the `sigma = 0.03` static-stability floor is now reported at INFO
  level in the run log. The floor is disabled by passing `apply_floor=False` to
  `StaticStability`, which is **experimental** and documented as such in
  `docs/source/usage.rst`: `sigma` sits in the denominator of `Az`, `Ae`, `Ca`,
  `Gz` and `Ge`, so those results are not comparable with the defaults. There
  is deliberately no command-line or environment-variable switch for it.
- Time series figures for the new outputs: `BPhi_Z`/`BPhi_E`,
  `C_overturning` and `M` now have their own groups in
  `src/plots/utils.py::TERM_DETAILS`. Hovmoller diagrams and boxplots still
  cover the energy, conversion and generation/dissipation groups only.

### Fixed

- **`Gz`, `Ge`, `Dz` and `De` were returned unconverted**: the unit conversion
  was computed and its result discarded, so the four terms carried whatever
  units the formula produced instead of `W/m^2`. The magnitudes are unchanged
  (the conversion is 1:1), but the returned values now carry the declared
  units, and a dimensional error in any of the four raises instead of passing
  silently.
- **`--box_limits` was ignored**: `slice_domain` read the hard-coded
  `inputs/box_limits`, so a run passing a different file sliced the dataset
  with one domain and computed the energetics on the intersection with
  another. The requested and realised domains are now logged.
- **Longitude conventions**: box limits are translated into the dataset's own
  convention, and wrapped domains fail explicitly instead of silently
  selecting the wrong region. Limits given in decreasing order (`min_lon`
  greater than `max_lon`) are reported as swapped values in the box-limits
  file, rather than as a dataset in the wrong longitude convention.
- **850 hPa constant**: the interactive `--choose` diagnostic used `8500`
  where the vertical coordinate is in Pa, selecting the level nearest
  85 hPa. It is now `85000`.
- **Vertical-level output**: column headers are written after the vertical
  control volume is fixed, so the per-level CSV files line up with the levels
  actually used.
- **Time tendencies**: finite differences now use elapsed seconds rather than
  assuming a uniform time step.
- **NaN policy**: missing values are interpolated along the vertical without
  dropping pressure levels, and the control volume is decided once for all
  terms instead of independently per term.
- **Unit-error handling**: `pint` raises `DimensionalityError`, a subclass of
  `TypeError`, which the previous `ValueError` guard could not catch.
- **Direct dissipation**: the non-residual dissipation pathway now raises
  `NotImplementedError` instead of returning a value built from a single
  scalar friction field. Use the residual formulation (`-r`). The `BoxData`
  attributes carrying the friction components were renamed from `ust`/`vst`
  (a conventional name for friction *velocity*, m s^-1) to `fric_u`/`fric_v`,
  since they now hold friction *forces* in m s^-2.

### Documentation

- `docs/source/math.rst` updated for the equations above, with the source of
  each adopted form recorded alongside it.
## [1.1.11] - 2025-02-26

### Changed

- **Code Organization**: Refactored `src/utils/tools.py` (1096 lines) into three focused modules:
  - `tools.py`: Utility functions (logging, CDS API, coordinate conversions) - 453 lines
  - `validation.py`: Input validation functions (track files, namelists, coordinates) - 449 lines
  - `preprocessing.py`: Data preprocessing functions (loading, processing, preparation) - 592 lines
  
- **Improved Maintainability**:
  - Separation of concerns: utilities, validation, and preprocessing logic are now separate
  - Better testability with focused modules
  - Clearer import statements showing dependencies
  - Module-level docstrings explaining each module's purpose

- **Import Updates**:
  - `lorenzcycletoolkit.py`: Now imports `prepare_data` from `src.utils.preprocessing`
  - `initialize_logging` remains in `src.utils.tools`
  - `get_cdsapi_data` remains in `src.utils.tools` (no changes to tests)
  
- **Module Structure**:
  - `tools.py` contains:
    - `initialize_logging()` - Setup application logging
    - `convert_longitude_range()` - Coordinate conversion utility
    - `find_extremum_coordinates()` - Find min/max points
    - `get_cdsapi_data()` - ERA5 data download from CDS API
    
  - `validation.py` contains:
    - `validate_track_file()` - Track file format validation
    - `validate_namelist_file()` - Namelist loading and validation
    - `validate_variable_match()` - Check namelist-dataset variable matching
    - `validate_required_coordinates()` - Validate presence of required coordinates
    
  - `preprocessing.py` contains:
    - `get_data()` - NetCDF file opening with error handling
    - `process_data()` - Data preprocessing (timestep validation, coordinate conversion, unit conversion)
    - `prepare_data()` - Main orchestration function

### Developer Notes

- No breaking changes for end users - all functionality remains the same
- Test suite remains unchanged (all tests use functions that stayed in `tools.py`)
- Error messages and logging behavior unchanged
- This refactoring improves code navigation and makes future maintenance easier

## [1.1.10] - 2026-02-26

### Improved

- **Comprehensive Error Messages**: All error messages now include:
  - Clear description of the problem
  - Specific solutions for users
  - Developer information (file paths, exception details, context)
  - Visual formatting with emojis and separators for better readability

- **File Opening Errors**: Enhanced error handling when opening NetCDF files:
  - FileNotFoundError: Shows absolute path, current directory, existence check
  - OSError: Displays file permissions, size, readability status
  - Generic errors: Includes full traceback and diagnostic information

- **Namelist Validation**: Improved namelist file error messages:
  - FileNotFoundError: Shows available preset namelists and copy commands
  - EmptyDataError: Explains the issue with clear solution
  - Parsing errors: Indicates format requirements and encoding issues

- **Coordinate Validation**: Added validation for required coordinates:
  - Checks for longitude, latitude, vertical level, and time coordinates
  - Shows which coordinates are missing vs available in dataset
  - Provides examples of common coordinate names
  - Lists dataset dimensions and coordinates for debugging

- **Unit Conversion Errors**: Enhanced pressure level unit conversion handling:
  - Assumes hPa if units attribute is missing (with warning)
  - Detailed error messages showing current and target units
  - Suggestions for fixing units with ncatted command
  - Developer info includes MetPy error details

- **Box Limits Validation** (Fixed Framework):
  - FileNotFoundError: Shows absolute path and example format
  - Missing fields: Lists required vs found fields
  - Invalid ranges: Explains min/max relationship requirements
  - All errors include file path and clear solutions

- **CDS API Error Handling**: Comprehensive error messages for:
  - Authentication errors: Links to registration and API key setup
  - Network errors: Suggests checking CDS API status page
  - File creation errors: Checks disk space and permissions
  - Download errors: Shows request parameters and common causes

### Added

- **Developer Information**: All error messages now include technical details:
  - Exception types and messages
  - File paths (both relative and absolute)
  - Variable values at point of failure
  - Stack traces where relevant
  - Configuration parameters

## [1.1.9] - 2026-02-26

### Added

- **Track File Format Validation**: Added comprehensive validation for track file format:
  - Automatic delimiter detection (`;` or `,`)
  - Validation of required columns (`time`, `Lat`, `Lon`)
  - Validation of date format (`YYYY-MM-DD-HHMM`)
  - Clear error messages showing expected format with examples
  - Support for both semicolon and comma delimiters

### Fixed

- **Track/Data Timestep Comparison**: Fixed inverted logic in timestep validation:
  - Now correctly allows track with larger timesteps than data (e.g., 6h track on 1h data)
  - Properly errors when data has larger timesteps than track (missing timesteps)
  - Added detailed error messages showing both data and track temporal information

### Improved

- **Enhanced Error Messages for Timestep Issues**: When timestep or timestamp mismatches occur, error messages now display:
  - Data time resolution (e.g., 1 hour, 6 hours)
  - Track time resolution
  - First and last timestamps for both data and track
  - Clear explanation of the problem
  - Actionable solutions for resolving the issue
  - Visual indicators (arrows, emojis) highlighting problematic timestamps

## [1.1.8] - 2026-02-26

### Changed

- **Namelist File Management**: The `inputs/namelist` file is now excluded from version control to prevent conflicts, as it should be customized by each user for their specific dataset.

### Improved

- **Enhanced Error Messages**: When the namelist doesn't match the dataset, the error message now displays:
  - Available coordinates in the dataset (with units and long names when available)
  - Available variables in the dataset (with units, long names, and standard names when available)
  - Helpful suggestions for which preset namelist to use
  - This makes it much easier for users to configure the correct namelist for their data

### Documentation

- **Configuration Guide**: Clarified that users must create an `inputs/namelist` file from one of the provided presets before running the toolkit
- **Usage Guide**: Added comprehensive sections on:
  - Prerequisites for each framework (fixed, moving, interactive)
  - Input data requirements (regular grid, isobaric levels, required variables)
  - Data preparation best practices to avoid memory issues
  - Examples of data preprocessing using both Python (xarray) and CDO
  - Links to the Configuration page for detailed setup instructions

## [1.1.7] - 2026-01-25

### Fixed

- **CDS API Last Timestep Issue**: Fixed a bug where the final timestamp from the track file was excluded when downloading ERA5 data via CDS API. The download logic now correctly includes the last hour of the track period.

### Added

- **Automatic Namelist Selection**: When using the `--cdsapi` flag, the program now automatically uses the ERA5-compatible namelist (`inputs/namelist_ERA5-cdsapi`). Users no longer need to manually copy or configure the namelist file for CDS API downloads.

### Improved

- **Enhanced Logging with Emojis**: Completely revamped logging throughout the application for better user experience:
  - 🚀 Clear startup messages
  - 📥 Download progress indicators
  - ✅ Success confirmations
  - ❌ Error messages with clear explanations
  - 🎨 Plot generation status
  - ⚡ Energy calculation progress
  - 🔄 Conversion term computation
  - 🏁 Boundary term computation
  - 🔥 Generation/dissipation term computation
  - 📈 Budget estimation progress
  - 💾 File save confirmations
  - 🗑️ Cleanup operations
  - 🎉 Completion messages

- **More Informative Error Messages**: Error messages now include specific timestamps and values to help users diagnose issues (e.g., track vs data timestamp mismatches).

### Documentation

- Updated `usage.rst` to reflect automatic namelist selection when using `--cdsapi` flag.

## [1.1.6] - 2026-01-25

### Changed

 - **`src/utils/tools.py` improvements**: Improved helper logging and error messages in the tools module.

## [1.1.5] - 2026-01-21

### Improved

- **CDS API Download Logging**: Significantly enhanced logging for ERA5 data downloads:
  - Shows list of expected files before download starts
  - Displays progress with checkmarks (✓) for successful downloads
  - Shows time ranges being downloaded for each day (start time to end time)
  - Provides clear success summary with file size and time coverage
  - More informative error messages if downloads fail

- **Smart Time Range Detection**: Fixed issue where downloads always requested full days (00Z-23Z):
  - Now detects actual start and end times from track file
  - For first day: downloads only from track start time onwards
  - For last day: downloads only until track end time
  - For middle days: downloads complete day as before
  - Respects the `--time-resolution` parameter for all calculations
  - Example: Track from 18Z to 23Z on same day will only download those hours

### Technical Details

- Track files can now start at any hour (e.g., 18Z) and end at any hour (e.g., 15Z)
- Time ranges are automatically rounded to the nearest `time_resolution` interval
- Reduces unnecessary data downloads and API costs
- Better handles single-day tracks with partial day coverage

## [1.1.4] - 2026-01-21

### Changed

- **Dependencies Update**: Updated CDS API and multiurl dependencies to use flexible version constraints:
  - `cdsapi`: Changed from `==0.7.0` to `>=0.7.0` (now compatible with version 0.7.7)
  - `multiurl`: Changed from `==0.3.1` to `>=0.3.1` (now compatible with version 0.3.7)
  - This resolves dependency conflicts when users upgrade cdsapi to the latest version
  - The package now automatically adapts to newer versions of these dependencies

### Notes

- Users can now safely upgrade cdsapi without encountering dependency conflicts
- The deprecated CDS API endpoint warning will persist until Copernicus fully migrates to the new API

## [1.1.3] - 2026-01-21

### Fixed

- **CDS API Download Issues**: Modified `get_cdsapi_data` function to download ERA5 data day-by-day instead of all at once. This resolves issues with large file downloads that were causing failures. The function now:
  - Downloads data separately for each day
  - Stores temporary daily files in a system temporary directory
  - Concatenates all daily files into the final output file
  - Automatically cleans up temporary files after concatenation
  - Provides better progress logging for multi-day downloads
- **Test Suite Warnings**: Fixed FutureWarnings in test suite by updating pandas `freq` parameter from `'H'` to `'h'`. Added `pytest.ini` to filter deprecation warnings from external libraries (cdsapi, pkg_resources, metpy).

### Changed

- **CDS API Workflow**: The CDS API download process is now more robust and reliable for large datasets spanning multiple days. Users will see progress updates for each day being downloaded.

## [1.1.2] - 2026-01-21

### Fixed / Improved

- **Environment & Installation**: Added `environment.yml` and clarified recommended Conda-based installation method for developers. Ensured `pip install -e .` is documented and included in developer install steps.
- **Package installation**: Confirmed `setup.py` includes the root module (`py_modules=['lorenzcycletoolkit']`) and bumped package version to `1.1.2` to reflect this release.
- **Documentation**: Rewrote `installation.rst` and updated `README.md` to keep the README concise while directing full install and usage instructions to the documentation. Fixed RST formatting issues and rebuilt docs to resolve related warnings.
- **Docs version**: Updated Sphinx `release` to `1.1.2` in `docs/source/conf.py`.
- 
## [1.1.1] - 2026-01-20

### Added

- **--time-resolution Flag**: Added `--time-resolution` command-line flag to control the temporal resolution of ERA5 data downloads via CDS API. Default is 3 hours to avoid exceeding API cost limits.
- **environment.yml**: Added `environment.yml` file to simplify conda environment creation for users.

### Fixed

- **CDS API Cost Limits**: Fixed issue where automatic 1-hour temporal resolution detection often exceeded CDS API cost limits. Now defaults to 3-hour resolution with user-configurable options (3, 6, 12, or 24 hours).
- **Test Suite**: Updated `test_cdsapi.py` to include `time_resolution` parameter in test fixtures and adjusted time step calculation test expectations.
- **Dependencies**: Fixed `docutils` version conflict between `sphinx` and `sphinx_rtd_theme` by downgrading to `docutils==0.20.1`.
- **Package Installation**: Fixed `setup.py` to properly include `lorenzcycletoolkit.py` module by adding `py_modules=['lorenzcycletoolkit']` parameter. This resolves the `ModuleNotFoundError` when installing the package.

### Changed

- **CDS API Temporal Resolution**: Changed from automatically detecting temporal resolution from track file to using a configurable `--time-resolution` parameter with a safe default of 3 hours.
- **Documentation Theme**: Updated documentation theme from Alabaster to ReadTheDocs (`sphinx_rtd_theme`) for better readability and navigation.
- **Repository URL**: Updated installation documentation to reflect new repository URL: `https://github.com/daniloceano/LorenzCycleToolkit`.
- **Documentation Version**: Updated documentation version to 1.1.1 in `conf.py`.

### Documentation

- Completely rewrote `installation.rst` with clearer instructions for:
  - Installing from PyPI for regular users
  - Installing from source with Conda (recommended for developers)
  - Installing from source with venv
  - Added troubleshooting section for common installation issues
- Updated `usage.rst` with information about the new `--time-resolution` flag and cost limit guidelines.
- Updated `cdsapi_example.rst` with detailed examples and guidelines for choosing appropriate temporal resolutions.
- Updated `CDSAPI_USAGE.md` with troubleshooting tips for cost limit errors and temporal resolution examples.
- Updated `README.md` with comprehensive installation instructions and examples for both users and developers.
- Updated `installation.rst` with corrected repository URLs.
- Added `sphinx_rtd_theme` to `requirements.txt`.


## [1.1.0] - 2025-11-18

### Updated

- **get_cdsapi_data Function**: Modernized the `get_cdsapi_data` function to align with the latest CDS API syntax (v0.7.x). This includes improved error handling, logging, and support for new request parameters such as `retry_max`.

### Added

- **Test Suite for get_cdsapi_data**: Created a comprehensive test suite (`tests/test_cdsapi.py`) to validate the `get_cdsapi_data` function. The suite includes 13 test cases covering parameter validation, error handling, and edge cases.



## [1.0.9] - 2024-10-21

### Fixed

- **LEC Fixed Framework**: Fixed a bug where the vertical levels output file was not saving the correct pressure levels as columns. Now, the output correctly reflects the pressure levels.
  
- **LEC Moving Framework**: Resolved a bug that prevented the framework from correctly handling maximum vorticity for the Northern Hemisphere track files. The framework now properly processes maximum vorticity values for the Northern Hemisphere.

- **Plot Periods**: Enhanced the functionality to check if the track is for the Northern Hemisphere. If so, the vorticity values are multiplied by `-1` to ensure correct behavior in the CycloPhaser program.

### Added

- **Tools.py**: Added a preprocessing stage to remove the new dimensions (`expver` and `number`) introduced in recent ERA5 datasets, ensuring compatibility with previous workflows.

## [1.0.8] - 2024-09-05

### Added
- Added badges for JOSS, Zenodo, documentation status, and license (GPL v3) to the README.

## [1.0.7] - 2024-09-05

### Fixed
- Fixed the broken link to the **Contributing** section in the `README.md` file.

## [1.0.6] - 2024-09-05

### Added
- Expanded the **Contributing** section to include detailed guidelines on branch naming conventions, code style (autopep8, isort, flake8), testing, and PR submission.
- Added a suggestion for running tests (`pytest`) and code formatting/linting commands (`autopep8`, `isort`, `flake8`) before submitting pull requests.

### Changed
- Improved the **Contributing** section in the documentation to provide better clarity on contribution processes and coding standards.
- Updated the **Contributing** section to include information about the Continuous Integration (CI) pipeline and the required coding practices.
- Removed "flags" from `index.rst` as the content was moved to the **Usage** section in version 1.0.5 but not reflected in the index.

### Fixed
- Resolved a merge conflict in the `README.md` file caused during cherry-picking of changes from the `joss-submission` branch to the `main` branch.


## [1.0.5] - 2024-09-04

### Added
- Updated README to clarify that installation instructions are included in the documentation.
- Revised README and documentation to reflect the GNU license, correcting the previously stated MIT license.
- Moved the "Flags" section to the "Usage" section to improve tutorial clarity and make the process more linear for users.
- Added a note in the documentation explaining that results generated using the LorenzCycleToolkit will be referenced in a future journal publication to aid users in interpreting outputs.

### Fixed
- Corrected links to the Contributing Guide and License in the README and documentation.

## [1.0.4] - 2024-08-28

### Added
- Added a `support` section in the documentation to provide guidelines for users seeking help or reporting issues. Users are encouraged to use the GitHub repository for reporting, with the option to contact via email if needed.

### Updated
- Updated the "Moving (Semi-Lagrangian) Framework Example" section in the documentation to include information about using the default track file for testing and the requirement to overwrite the `inputs/track` file with custom track data.
- Clarified that the preferred method for running the LorenzCycleToolkit is through command line arguments executed from the top-level directory of the project.

### Fixed
- Improved clarity in the documentation by specifying the appropriate directories and command usage for various examples and tutorials.

## [1.0.3] - 2024-07-27
### Bug Fixes
- Deploying documentation

## [1.0.2] - 2024-07-27
### Added
- Improved CI/CD pipeline by ensuring up-to-date gh-pages branch before deployment.
- Updated installation instructions to include pip installation from PyPI.

## [1.0.1] - 2024-07-27
### Added
- Added `CHANGELOG.md` file to document changes.
- Published project to PyPI.
- Integrated CircleCI for continuous integration and deployment.

## [1.0.0] - 2024-07-16
### Initial Release
- Initial release of LorenzCycleToolkit.
- Implemented core functionalities for calculating the Lorenz Energy Cycle (LEC) in specific atmospheric regions.
- Added support for both Eulerian and Semi-Lagrangian frameworks.
- Provided options for fixed and moving computational domains.
- Generated detailed results in CSV format for energy, conversion, and generation terms.
- Included automatic plot generation for spatial and temporal analysis.
- Utilized Xarray, Dask, Pandas, and MetPy for data handling and computations.
- Comprehensive documentation and usage examples.
