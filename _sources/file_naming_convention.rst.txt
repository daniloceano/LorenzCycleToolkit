File Naming Convention
======================

Input Files
-----------

For input NetCDF files, we recommend using a descriptive naming convention that follows this pattern:

**Format**: ``<subject>_<datasource>.nc``

**Examples**:

- ``Catarina_NCEP-R2.nc`` - Hurricane Catarina using NCEP Reanalysis 2 data
- ``testdata_ERA5.nc`` - Test dataset using ERA5 reanalysis
- ``Cyclone-20100101_NCEP-R2.nc`` - Cyclone from January 1, 2010 (compound names use hyphens)

**Guidelines**:

- Use underscores (``_``) to separate the subject from the data source
- For compound subject names, use hyphens (``-``) within the subject part
- Include the data source (e.g., ERA5, NCEP-R2, MPAS-A) for clarity
- Use the ``.nc`` extension for NetCDF files

Track Files
-----------

For track files, use a similar convention:

**Format**: ``track_<subject>`` or ``track_<subject>_<date>``

**Examples**:

- ``track_Catarina``
- ``track_20100101``
- ``track_testdata_ERA5``

Box Limits Files
----------------

For box limits files:

**Format**: ``box_limits_<region>`` or ``box_limits_<study>``

**Examples**:

- ``box_limits_SouthAtlantic``
- ``box_limits_Reg1``
- ``box_limits_testcase``

**Note**: These naming conventions are recommendations to help organize your work. The toolkit will work with any valid filename.
