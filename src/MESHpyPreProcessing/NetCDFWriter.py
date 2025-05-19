"""

Overview
========

The ``NetCDFWriter`` class is designed to generate model-ready NetCDF files (e.g., MESH_parameters.nc) containing 
soil and other geophysical subbasin data integrated from a vector shapefile and a NetCDF drainage database. 
This class is typically used in workflows that prepare input parameters for land surface models like MESH.

It supports flexible handling of both **layer-dependent** (e.g., soil properties per depth layer) and 
**layer-independent** (e.g., slope, contributing area) variables. The output conforms to CF conventions 
and includes appropriate coordinate reference metadata for spatial consistency.

Function Descriptions
=====================

.. py:class:: NetCDFWriter(nc_filename, shapefile_path, input_ddb_path)

   Initializes the NetCDF writer with paths to the output file, input shapefile, and NetCDF drainage database.

   :param nc_filename: Path to the NetCDF output file to be created.
   :type nc_filename: str
   :param shapefile_path: Path to the input shapefile containing the attributes.
   :type shapefile_path: str
   :param input_ddb_path: Path to the NetCDF drainage database used to extract coordinates.
   :type input_ddb_path: str

.. py:method:: read_shapefile()

   Reads the input shapefile and converts it into a GeoDataFrame. The file is automatically reprojected to EPSG:4326 (WGS 84).

.. py:method:: set_coordinates()

   Extracts `lon`, `lat`, and `subbasin` values from the NetCDF drainage database file. These values serve as the spatial base for NetCDF output.

.. py:method:: set_num_soil_layers(num_layers)

   Sets the number of vertical soil layers that will be written into the NetCDF file.

   :param num_layers: The number of soil layers (e.g., 4 for a 4-layer soil profile).
   :type num_layers: int

.. py:method:: add_var_attrs(var, attrs)

   Adds metadata attributes to a NetCDF variable, such as units, standard name, and axis designation.

   :param var: The NetCDF variable to modify.
   :type var: netCDF4.Variable
   :param attrs: Dictionary of attributes to apply.
   :type attrs: dict

.. py:method:: write_netcdf(properties, variable_info)

   Writes the actual NetCDF file using the specified properties and metadata.

   :param properties: Dictionary specifying which variables are layer-dependent vs. layer-independent.
   :type properties: dict
   :param variable_info: Dictionary mapping each variable to a tuple of (NetCDF name, data type, unit).
   :type variable_info: dict

Example Usage
=============

.. code-block:: python

    from MESHpyPreProcessing.NetCDFWriter import NetCDFWriter

    # Paths for NetCDFWriter
    nc_filename = 'MESH_parameters3.nc'
    output_shapefile = 'merged_soil_data_shapefile4.shp'
    input_ddb = '/scratch/fuaday/sras-agg-model/MESH-sras-agg/MESH_drainage_database.nc'
    mesh_intervals = [(0, 0.1), (0.1, 0.35), (0.35, 1.2), (1.2, 4.1)]

    # Initialize NetCDFWriter with the necessary paths
    nc_writer = NetCDFWriter(
        nc_filename=nc_filename,
        shapefile_path=output_shapefile,
        input_ddb_path=input_ddb
    )

    # Step 1: Read the attribute shapefile and extract spatial coordinates from the drainage database
    nc_writer.read_shapefile()
    nc_writer.set_coordinates()

    # Step 2: Specify the number of vertical soil layers to include in the output
    nc_writer.set_num_soil_layers(num_layers=len(mesh_intervals))

    # Step 3: Define which variables are layer-dependent vs. layer-independent
    properties = {
        'layer_dependent': ['CLAY', 'SAND', 'OC'],  # Varies by soil layer and subbasin
        'layer_independent': ['ncontr', 'meanBDRICM', 'meanBDTICM', 'xslp', 'dd']  # Varies only by subbasin
    }

    # Step 4: Provide metadata for each variable to be written to NetCDF
    variable_info = {
        'CLAY': ('CLAY', 'f4', 'Percentage'),
        'SAND': ('SAND', 'f4', 'Percentage'),
        'OC': ('ORGM', 'f4', 'Percentage'),
        'ncontr': ('IWF', 'i4', '1'),
        'meanBDRICM': ('BDRICM', 'f4', 'Meters'),
        'meanBDTICM': ('BDTICM', 'f4', 'Meters'),
        'xslp': ('xslp', 'f4', 'degree'),
        'dd': ('dd', 'f4', 'm_per_km2')
    }

    # Step 5: Write the final NetCDF file with structured metadata and spatial consistency
    nc_writer.write_netcdf(properties=properties, variable_info=variable_info)
"""
import os
import numpy as np
import geopandas as gpd
import netCDF4 as nc
from datetime import datetime
import tempfile
import xarray as xs

class NetCDFWriter:
    """
    A class to generate NetCDF files with soil data merged from shapefiles and NetCDF drainage databases.

    Attributes:
    -----------
    nc_filename : str
        Path to the output NetCDF file.
    shapefile_path : str
        Path to the input shapefile.
    input_ddb_path : str
        Path to the NetCDF drainage database.
    merged_gdf : geopandas.GeoDataFrame
        GeoDataFrame containing merged shapefile data.
    lon : list
        List of longitude values from the NetCDF drainage database.
    lat : list
        List of latitude values from the NetCDF drainage database.
    segid : list
        List of subbasin identifiers.
    num_soil_lyrs : int
        Number of soil layers in the dataset.
    """
    def __init__(self, nc_filename, shapefile_path, input_ddb_path):
        self.nc_filename = nc_filename
        self.shapefile_path = shapefile_path
        self.input_ddb_path = input_ddb_path
        self.merged_gdf = gpd.GeoDataFrame()
        self.lon = []
        self.lat = []
        self.segid = []
        self.num_soil_lyrs = 0

    def read_shapefile(self):
        """
        Reads the shapefile and converts it into a GeoDataFrame.
        
        This function reads the shapefile, reprojects it to EPSG:4326 (WGS 84),
        and stores the result in the merged GeoDataFrame.
        """
        self.merged_gdf = gpd.read_file(self.shapefile_path).to_crs(epsg=4326)

    def set_coordinates(self):
        """
        Extracts longitude, latitude, and subbasin IDs from the NetCDF drainage database.
        """
        db = xs.open_dataset(self.input_ddb_path)
        self.lon = db.variables['lon'].values
        self.lat = db.variables['lat'].values
        self.segid = db.variables['subbasin'].values
        db.close()

    def set_num_soil_layers(self, num_layers):
        """
        Sets the number of soil layers for the NetCDF file.
        
        Parameters:
        -----------
        num_layers : int
            Number of soil layers to be included in the NetCDF file.
        """
        self.num_soil_lyrs = num_layers

    def add_var_attrs(self, var, attrs):
        """
        Adds attributes to a NetCDF variable.
        
        Parameters:
        -----------
        var : netCDF4.Variable
            The NetCDF variable to which attributes will be added.
        attrs : dict
            A dictionary of attribute names and values.
        """
        for attr, value in attrs.items():
            var.setncattr(attr, value)

    def write_netcdf(self, properties, variable_info):
        """
        Creates a NetCDF file with processed soil data.
        
        Parameters:
        -----------
        properties : dict
            A dictionary with two keys:
            - 'layer_dependent': List of property names tied to the number of soil layers.
            - 'layer_independent': List of property names dependent only on the subbasin.
        variable_info : dict
            A dictionary mapping property names to tuples containing:
            (new variable name in NetCDF, data type code, unit).
        """
        try:
            rootgrp = nc.Dataset(self.nc_filename, "w", format="NETCDF4")
        except PermissionError:
            temp_dir = tempfile.gettempdir()
            self.nc_filename = os.path.join(temp_dir, "MESH_parameters.nc")
            rootgrp = nc.Dataset(self.nc_filename, "w", format="NETCDF4")

        subbasin_dim = rootgrp.createDimension("subbasin", len(self.lon))
        nsol_dim = rootgrp.createDimension("nsol", self.num_soil_lyrs)

        lon_var = rootgrp.createVariable("lon", "f4", ("subbasin",), fill_value=-1.0)
        lat_var = rootgrp.createVariable("lat", "f4", ("subbasin",), fill_value=-1.0)
        subbasin_var = rootgrp.createVariable("subbasin", "i4", ("subbasin",))

        lon_var.units = "degrees_east"
        lat_var.units = "degrees_north"
        subbasin_var.units = "1"

        lon_var[:] = np.array(self.lon)
        lat_var[:] = np.array(self.lat)
        subbasin_var[:] = np.array(self.segid, dtype=np.int32)

        rootgrp.setncattr("Conventions", "CF-1.0")
        rootgrp.setncattr("institution", "ECCC")
        rootgrp.setncattr("history", f"Fuad Yassin, {datetime.now().strftime('%Y-%m-%d')}")

        proj = rootgrp.createVariable("crs", "i4", ())
        self.add_var_attrs(proj, {
            "grid_mapping_name": "latitude_longitude",
            "longitude_of_prime_meridian": 0,
            "semi_major_axis": 6378137.0,
            "inverse_flattening": 298.257223563
        })

        rootgrp.close()

