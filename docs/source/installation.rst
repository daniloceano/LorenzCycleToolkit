Installation
============

Prerequisites
-------------
- Python 3.12 or later
- Conda (recommended) or pip

Installation Methods
--------------------

Method 1: Install from Source with Conda (Recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to contribute to the project or modify the code, follow these steps:

1. Clone the repository::

      git clone https://github.com/daniloceano/LorenzCycleToolkit.git
      cd LorenzCycleToolkit

2. Create a Conda environment with Python 3.12::

      conda create --name lorenz python=3.12 -y
      conda activate lorenz

3. Install dependencies::

      pip install -r requirements.txt

4. Install the package in editable mode::

      pip install -e .

5. Verify the installation::

      lorenzcycletoolkit --help

Method 2: Install from Source with venv
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you prefer using Python's built-in virtual environment:

1. Clone the repository::

      git clone https://github.com/daniloceano/LorenzCycleToolkit.git
      cd LorenzCycleToolkit

2. Create a virtual environment and activate it::

      python -m venv venv
      source venv/bin/activate  # On Windows: venv\Scripts\activate

3. Install dependencies::

      pip install -r requirements.txt

4. Install the package in editable mode::

      pip install -e .

5. Verify the installation::

      lorenzcycletoolkit --help

Optional Packages for the Download Routes
-----------------------------------------

The CDS API (``--cdsapi``) and NCAR GDEX (``--ncar``) routes need nothing beyond ``requirements.txt``. Three routes
need optional packages, installed as extras of the toolkit (from the cloned repository):

.. list-table::
   :header-rows: 1
   :widths: 18 22 60

   * - Route
     - Command
     - Packages
   * - ``--wekeo``
     - ``pip install -e ".[wekeo]"``
     - ``hda`` (2.40 or later) and ``requests``
   * - ``--arco``
     - ``pip install -e ".[arco]"``
     - ``zarr`` 2 (zarr 3 is not supported), ``numcodecs`` < 0.16, ``fsspec`` and ``aiohttp``
   * - ``--destine``
     - ``pip install -e ".[destine]"``
     - the same packages as ``[arco]``

Without them, the route stops at once and prints the installation command. Accounts and credentials are described in
:doc:`remote_data_acquisition`.

Troubleshooting
---------------

**Issue**: Command ``lorenzcycletoolkit`` not found after installation.

**Solution**: Make sure you've activated the correct environment:

- For Conda: ``conda activate lorenz``
- For venv: ``source venv/bin/activate``

**Issue**: ``ModuleNotFoundError: No module named 'lorenzcycletoolkit'``

**Solution**: Install the package using ``pip install -e .`` (for development) or ``pip install LorenzCycleToolkit`` (for regular use).

**Issue**: Tests failing with import errors.

**Solution**: Ensure all dependencies are installed by running ``pip install -r requirements.txt`` in your activated environment

Next Steps
----------

After successful installation, proceed to:

- :doc:`configuration` - Set up your namelist and domain files before running analyses
- :doc:`usage` - Learn about command-line options and frameworks
- :doc:`examples_and_tutorials` - Follow step-by-step examples with sample data
