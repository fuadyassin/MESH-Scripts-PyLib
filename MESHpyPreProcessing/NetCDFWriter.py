"""
NetCDF File Generation
======================

NetCDFWriter.py contains a class NetCDFWriter that creates NetCDF files with processed soil data merged from shapefiles and NetCDF drainage databases.

Example Usage:
--------------
>>> from MESHpyPreProcessing.NetCDFWriter import NetCDFWriter
>>> writer = NetCDFWriter(nc_filename='output.nc', shapefile_path='path/to/shapefile.shp', input_ddb_path='path/to/input_ddb.nc')
>>> writer.read_shapefile()
>>> writer.set_coordinates()
>>> writer.set_num_soil_layers(3)
>>> properties = {'layer_dependent': ['OC', 'Sand'], 'layer_independent': ['Drainage_Area']}
>>> variable_info = {'OC': ('OrganicCarbon', 'f4', 'kg/m^2'), 'Sand': ('SandContent', 'f4', '%'), 'Drainage_Area': ('DrainageArea', 'f4', 'km^2')}
>>> writer.write_netcdf(properties, variable_info)

Functions:
----------
- read_shapefile: Reads the shapefile and prepares the GeoDataFrame.
- set_coordinates: Extracts longitude and latitude from the NetCDF drainage database.
- set_num_soil_layers: Sets the number of soil layers for the NetCDF file.
- write_netcdf: Writes processed data to a NetCDF file.
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


# import os
# import numpy as np
# import geopandas as gpd
# import netCDF4 as nc
# from datetime import datetime
# import tempfile
# import xarray as xs

# class NetCDFWriter:
#     def __init__(self, nc_filename, shapefile_path, input_ddb_path):
#         self.nc_filename = nc_filename
#         self.shapefile_path = shapefile_path
#         self.input_ddb_path = input_ddb_path
#         self.merged_gdf = gpd.GeoDataFrame()
#         self.lon = []
#         self.lat = []
#         self.segid = []
#         self.num_soil_lyrs = 0

#     def read_shapefile(self):
#         """
#         Read the shapefile and set the merged GeoDataFrame.
#         """
#         self.merged_gdf = gpd.read_file(self.shapefile_path).to_crs(epsg=4326)

#     def set_coordinates(self):
#         """
#         Set longitude and latitude values from a NetCDF drainage database.
#         """
#         db = xs.open_dataset(self.input_ddb_path)
#         self.lon = db.variables['lon'].values
#         self.lat = db.variables['lat'].values
#         self.segid = db.variables['subbasin'].values
#         db.close()

#     def set_num_soil_layers(self, num_layers):
#         """
#         Set the number of soil layers.
#         """
#         self.num_soil_lyrs = num_layers

#     def add_var_attrs(self, var, attrs):
#         """
#         Add attributes to a variable.
#         """
#         for attr, value in attrs.items():
#             var.setncattr(attr, value)

#     def write_netcdf(self, properties, variable_info):
#         """
#         Create a NetCDF file with the processed soil data.
#         properties: dict
#             A dictionary with two keys:
#             - 'layer_dependent': List of property names tied to the number of soil layers.
#             - 'layer_independent': List of property names dependent only on the subbasin.
#         variable_info: dict
#             A dictionary with keys as property names and values as tuples containing 
#             (new variable name in NetCDF, data type code).
#         """
#         try:
#             rootgrp = nc.Dataset(self.nc_filename, "w", format="NETCDF4")
#         except PermissionError:
#             temp_dir = tempfile.gettempdir()
#             self.nc_filename = os.path.join(temp_dir, f"MESH_parameters.nc")
#             rootgrp = nc.Dataset(self.nc_filename, "w", format="NETCDF4")

#         # Calculate the indices of the matching COMID values
#         ind = []
#         for i in range(len(self.segid)):
#             fid = np.where(np.int32(self.merged_gdf['COMID'].values) == self.segid[i])[0]
#             ind = np.append(ind, fid)
#         ind = np.int32(ind)

#         subbasin_dim = rootgrp.createDimension("subbasin", len(self.lon))
#         nsol_dim = rootgrp.createDimension("nsol", self.num_soil_lyrs)

#         lon_var = rootgrp.createVariable("lon", "f4", ("subbasin",), fill_value=-1.0)
#         lat_var = rootgrp.createVariable("lat", "f4", ("subbasin",), fill_value=-1.0)
#         time_var = rootgrp.createVariable("time", "f4", ("subbasin",), fill_value=-1.0)
#         subbasin_var = rootgrp.createVariable("subbasin", "i4", ("subbasin",))

#         lon_var.units = "degrees_east"
#         lat_var.units = "degrees_north"
#         time_var.units = "days since 1980-10-01 00:00:00.0 -0:00"

#         lon_var[:] = np.array(self.lon)
#         lat_var[:] = np.array(self.lat)
#         time_var[:] = np.zeros(len(self.lon))
#         subbasin_var[:] = np.array(self.segid, dtype=np.int32)

#         # Set variable attributes
#         self.add_var_attrs(lon_var, {"standard_name": "longitude", "axis": "X", "grid_mapping": "crs"})
#         self.add_var_attrs(lat_var, {"standard_name": "latitude", "axis": "Y", "grid_mapping": "crs"})
#         self.add_var_attrs(time_var, {"standard_name": "time", "axis": "T", "grid_mapping": "crs"})
#         self.add_var_attrs(subbasin_var, {"standard_name": "subbasin_id", "long_name": "Subbasin ID"})

#         # Handle properties tied to number of soil layers
#         if 'layer_dependent' in properties:
#             for prop in properties['layer_dependent']:
#                 new_name, dtype = variable_info[prop]
#                 data_var = rootgrp.createVariable(new_name, dtype, ("subbasin", "nsol"), fill_value=-1.0)
#                 data_var.long_name = f"{prop} Content of Soil Layer"
#                 self.add_var_attrs(data_var, {"grid_mapping": "crs", "standard_name": prop.lower(), "coordinates": "lat lon time"})
#                 for i in range(self.num_soil_lyrs):
#                     data_var[:, i] = np.array(self.merged_gdf[f'mesh{prop}{i+1}'].values[ind])

#         # Handle properties dependent only on subbasin
#         if 'layer_independent' in properties:
#             for prop in properties['layer_independent']:
#                 new_name, dtype = variable_info[prop]
#                 data_var = rootgrp.createVariable(new_name, dtype, ("subbasin",), fill_value=-1.0)
#                 data_var.long_name = f"{prop} Content per Subbasin"
#                 self.add_var_attrs(data_var, {"grid_mapping": "crs", "standard_name": prop.lower(), "coordinates": "lat lon time"})
#                 data_var[:] = np.array(self.merged_gdf[prop].values[ind])

#         # Global attributes
#         rootgrp.setncattr("Conventions", "CF-1.0")
#         rootgrp.setncattr("source", "MERIT geogabrics and GSDE soil")
#         rootgrp.setncattr("institution", "ECCC")
#         rootgrp.setncattr("references", "xx et al. (xxxx) journal xx:xx-xx")
#         rootgrp.setncattr("history", f"Fuad Yassin, {datetime.now().strftime('%Y-%m-%d')}")
#         rootgrp.setncattr("featureType", "point")

#         # CRS variable
#         proj = rootgrp.createVariable("crs", "i4", ())
#         self.add_var_attrs(proj, {
#             "grid_mapping_name": "latitude_longitude",
#             "longitude_of_prime_meridian": 0,
#             "semi_major_axis": 6378137.0,
#             "inverse_flattening": 298.257223563
#         })

#         rootgrp.close()
