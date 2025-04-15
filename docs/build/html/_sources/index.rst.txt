MESH-Scripts-PyLib IO Documentation
=====================================

**MESH-Scripts-PyLib** is a Python library for preprocessing hydrometric and soil data, performing spatial analysis, and generating NetCDF files for use in hydrological modeling. The package includes utilities for:

- Streamflow file preparation
- Soil data processing
- Spatial analysis
- NetCDF file generation
- Vector Aggregation
- Generate MESH_CLASS.ini file
- Generate MESH_HYDROLOGY.ini file
- Plot Variables from Drainage Database and Parameters file

.. note::
   This project is under active development.

Installation
------------

You can install the package using pip directly from GitHub:

.. code-block:: bash

   pip install git+https://github.com/MESH-Model/MESH-Scripts-PyLib.git

Overview
--------

This library provides several utilities that streamline data preprocessing for hydrological models like MESH. The key functionalities include:

- **Streamflow Data Processing:** Extraction and processing of streamflow data from Canadian and US hydrometric sources.
- **Soil Data Handling:** Loading, merging, and cleaning soil data from CSV files.
- **Spatial Analysis:** Flagging non-contributing areas using GeoDataFrames.
- **Basin & River Aggregation:** Aggregation of basin and river networks with custom thresholds.
- **NetCDF Generation:** Writing processed soil data to NetCDF format for hydrological modeling.

Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: File preparation:

   Aggregation_vector.rst
   NetCDFWriter.rst
   convert_ddbnetcdf.rst
   gdf_edit.rst
   gen_streamflow_file.rst
   gsde_soil.rst
   generate_mesh_class_ini_from_excel.rst
   generate_mesh_hydrology_ini_from_excel.rst
   remap_climate_to_ddb.rst

Visualization
---------------

.. toctree::
   :maxdepth: 2
   :caption: Visualization:
   
   plt_var_vector_setup.rst

Jupyter Notebooks
-----------------

.. toctree::
   :maxdepth: 2
   :caption: Jupyter Notebooks:

   MESH_StreamflowFilePrep.rst

Indices and Tables
------------------

- :ref:`genindex`
- :ref:`modindex`
- :ref:`search`
