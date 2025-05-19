"""
Spatial Variable Visualization from Vector and NetCDF Files
============================================================

This module provides a unified plotting interface to visualize spatial variables from NetCDF datasets 
(e.g., MESH drainage database or parameters file) overlaid on a vector shapefile (e.g., subbasin polygons).
It supports:
- Single variable plotting per subbasin
- Multi-class land use variable visualization (e.g., GRUs)
- Layer-wise soil variable visualization (e.g., SAND, CLAY, OC)

It automatically interprets the variable's dimensions to choose the appropriate visualization layout 
(single map or subplots).

----------------------------------------------
Function: plt_var_from_vector_ddb_netcdf
---------------------------------------------

Description:
-------------
Plots spatial variables stored in NetCDF files on top of vector shapefiles using GeoPandas and Matplotlib.
This function is useful for visualizing hydrological model inputs such as Grouped Respone Units (GRUs), 
soil layers (e.g., SAND, CLAY), or other scalar variables (e.g., Elevation, Drainage Area) in spatial context.

Parameters:
------------
- output_basin_path : str
    Full path to the basin shapefile (.shp). Must contain a unique subbasin identifier column (e.g., 'COMID').
- ddbnetcdf_path : str
    Path to the input NetCDF file containing the spatial variable (e.g., MESH_drainage_database.nc or MESH_parameters.nc).
- variable_name : str
    Name of the variable to extract and visualize from the NetCDF file.
- save_path : str, optional (default='plot.png')
    Full path to save the output image file (.png).
- text_location : tuple, optional (default=(0.55, 0.95))
    (x, y) coordinates in axes fraction (0–1) to place text annotations inside each subplot.
- font_size : int, optional (default=10)
    Font size for subplot titles and annotations.
- cmap : str, optional (default='viridis')
    Matplotlib colormap to use for shading data values.
- cbar_location : list of float, optional (default=[0.96, 0.15, 0.02, 0.7])
    Location of colorbar [left, bottom, width, height] in normalized figure coordinates.
- subplot_adjustments : dict, optional
    Dictionary of subplot spacing and padding adjustments passed to `fig.subplots_adjust()`.
    Example: {'left': 0.1, 'right': 0.9, 'bottom': 0.1, 'top': 0.9, 'wspace': 0.1, 'hspace': 0.2}
- subbasin_var : str, optional (default='subbasin')
    Name of the NetCDF variable that contains subbasin indices.
- comid_var : str, optional (default='COMID')
    Name of the shapefile column used to join NetCDF data to the shapefile geometry.
- landuse_classes : list of str, optional
    Used when `variable_name` has GRU land use classes. Overrides land use names from the NetCDF file.
- grudim : str, optional (default='NGRU')
    Name of the NetCDF dimension for land use classes.
- grunames_var : str, optional (default='LandUse')
    Name of the NetCDF variable containing land use class names.
- soldim : str, optional (default='nsol')
    Name of the NetCDF dimension representing soil layers.

Processing Logic:
------------------
- If `variable_name` has **1D dimension** (e.g., `[subbasin]`), a single map is plotted.
- If it has **2D shape with land use dim** (e.g., `[subbasin, NGRU]`), a grid of subplots is created per GRU.
- If it has **2D shape with soil layer dim** (e.g., `[subbasin, nsol]`), a grid of subplots is created per soil layer.

Input Format:
--------------
1. **Shapefile (.shp)** with a subbasin identifier column (e.g., 'COMID')
2. **NetCDF file (.nc)** with variables such as:
   - `lat`, `lon`
   - 1D or 2D variable of interest (e.g., GRU, SAND)
   - Dimension definitions for subbasin, land use, or soil layers

Output:
--------
- A static plot (.png) saved at `save_path` location
- Subplots with dissolved basin outline and percentage annotations
- Shared colorbar with variable name and units

Example Usage:
---------------

>>> import os
>>> from MESHpyPostProcessing.plt_var_vector_setup import plt_var_from_vector_ddb_netcdf

>>> # Example 1: GRU (Grouped Land Use Classes) Visualization
>>> base_dir = 'D:/Coding/GitHub/Repos/MESH-Scripts-PyLib/MESHpyPostProcessing/ExampleFiles'
>>> shapefile_path = os.path.join(base_dir, 'sras_subbasins_MAF_Agg.shp')
>>> netcdf_path = os.path.join(base_dir, 'MESH_drainage_database.nc')
>>> output_path = os.path.join(base_dir, 'Outputs', 'GRU.png')
>>> lclass = [
...     'Temperate or sub-polar needleleaf forest', 'Sub-polar taiga needleleaf forest',
...     'Temperate or sub-polar broadleaf deciduous forest', 'Mixed forest',
...     'Temperate or sub-polar shrubland', 'Temperate or sub-polar grassland',
...     'Sub-polar or polar grassland-lichen-moss', 'Wetland', 'Cropland',
...     'Barren lands', 'Urban', 'Water', 'Dump'
... ]
>>> plt_var_from_vector_ddb_netcdf(
...     output_basin_path=shapefile_path,
...     ddbnetcdf_path=netcdf_path,
...     variable_name='GRU',
...     save_path=output_path,
...     text_location=(0.55, 0.95),
...     font_size=10,
...     cmap='gnuplot2_r',
...     cbar_location=[0.91, 0.15, 0.02, 0.7],
...     subplot_adjustments={'left': 0.1, 'right': 0.9, 'bottom': 0.1, 'top': 0.9, 'wspace': 0.1, 'hspace': 0.2},
...     subbasin_var='subbasin',
...     comid_var='COMID',
...     landuse_classes=lclass,
...     grudim='NGRU',
...     grunames_var='LandUse'
... )

.. image:: Figures/GRU.png
   :width: 600
   :alt: GRU Output Plot
   :align: center

>>> # Example 2: Soil Variable Visualization (e.g., SAND across layers)
>>> netcdf_path = os.path.join(base_dir, 'MESH_parameters.nc')
>>> output_path = os.path.join(base_dir, 'Outputs', 'SAND.png')
>>> plt_var_from_vector_ddb_netcdf(
...     output_basin_path=shapefile_path,
...     ddbnetcdf_path=netcdf_path,
...     variable_name='SAND',
...     save_path=output_path,
...     text_location=(0.55, 0.95),
...     font_size=10,
...     cmap='gnuplot2_r',
...     cbar_location=[0.91, 0.15, 0.02, 0.7],
...     subplot_adjustments={'left': 0.1, 'right': 0.9, 'bottom': 0.1, 'top': 0.9, 'wspace': 0.1, 'hspace': 0.2},
...     subbasin_var='subbasin',
...     comid_var='COMID'
... )

.. image:: Figures/SAND.png
   :width: 600
   :alt: SAND Layer Output Plot
   :align: center


Notes:
-------
- All spatial plotting is done in EPSG:4326 (WGS84)
- The function automatically merges the NetCDF data with shapefile geometry using the specified ID fields
- Color scaling (`vmin`, `vmax`) and transparency for `0` values is handled internally

Dependencies:
-------------
- geopandas
- pandas
- matplotlib
- numpy
- netCDF4

"""
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import netCDF4 as nc
import matplotlib.colors as mcolors
import math
import numpy as np

def plt_var_from_vector_ddb_netcdf(
    output_basin_path, 
    ddbnetcdf_path, 
    variable_name,
    save_path='plot.png',
    text_location=(0.55, 0.95), 
    font_size=10,
    cmap='viridis',  # Set the default colormap to 'viridis'
    cbar_location=[0.96, 0.15, 0.02, 0.7],  # Set default colorbar location
    subplot_adjustments=None,  # Optional: Custom subplot adjustments
    subbasin_var='subbasin',  # Default to 'subbasin'
    comid_var='COMID',  # Default to 'COMID'
    landuse_classes=None,  # Optional land use classes as a numpy array or list
    grudim='NGRU',    # from netcdf ddb
    grunames_var='LandUse',  # from netcdf ddb
    soldim='nsol'    # from netcdf parameter file
):
    # Load and dissolve the shapefile to get a single boundary sub_agg_gdf
    sub_agg_gdf = gpd.read_file(output_basin_path).to_crs(epsg=4326)
    dissolved_basin = sub_agg_gdf.dissolve().boundary

    # Load NetCDF data
    with nc.Dataset(ddbnetcdf_path) as dataset:
        variable_data = dataset.variables[variable_name][:]
        subbasin_ids = dataset.variables[subbasin_var][:].astype(int)
        variable_units = dataset.variables[variable_name].units
        
        # Extract latitude and longitude data
        latitudes = dataset.variables['lat'][:]
        longitudes = dataset.variables['lon'][:]
        
        # Check the dimensions of the variable
        dims = dataset.variables[variable_name].dimensions

        if len(dims) == 1:  # Case where the variable only depends on subbasin
            df = pd.DataFrame({
                variable_name: variable_data, 
                'lat': latitudes, 
                'lon': longitudes
            }, index=subbasin_ids)

            # Check for duplicate columns before merging
            duplicate_columns = set(df.columns).intersection(set(sub_agg_gdf.columns))
            if duplicate_columns:
                raise ValueError(f"Duplicate columns found during merge: {duplicate_columns}")

            sub_agg_ddb_merged_gdf = sub_agg_gdf.merge(df, left_on=comid_var, right_index=True, how='left')

            vmin, vmax = sub_agg_ddb_merged_gdf[variable_name].min(), sub_agg_ddb_merged_gdf[variable_name].max()

            # Plot the single variable with lat/lon as x and y axes
            fig, ax = plt.subplots(figsize=(8, 6))  #This needs to be inputs
            sub_agg_ddb_merged_gdf.plot(column=variable_name, ax=ax, cmap=cmap, legend=False, vmin=vmin, vmax=vmax)
            ax.set_title(variable_name, fontsize=font_size)
            ax.set_xlabel('Longitude', fontsize=font_size)
            ax.set_ylabel('Latitude', fontsize=font_size)
            dissolved_basin.plot(ax=ax, edgecolor='black', linewidth=1, facecolor='none')

        elif len(dims) == 2 and dims[1] == grudim:  # Case where the variable depends on subbasin and NGRU
            print ("len(dims) == 2 and dims[1] == grudim")
            landuse_classes = np.array(landuse_classes) if landuse_classes else dataset.variables[grunames_var][:]
            df = pd.DataFrame({landuse: variable_data[:, i] for i, landuse in enumerate(landuse_classes)}, index=subbasin_ids)

            # Check for duplicate columns before merging
            duplicate_columns = set(df.columns).intersection(set(sub_agg_gdf.columns))
            if duplicate_columns:
                raise ValueError(f"Duplicate columns found during merge: {duplicate_columns}")

            sub_agg_ddb_merged_gdf = sub_agg_gdf.merge(df, left_on=comid_var, right_index=True, how='left')

            percentage_pairs = [(landuse, round(sub_agg_ddb_merged_gdf[landuse].mean() * 100, 2)) for landuse in landuse_classes]
            sorted_landuse_columns, sorted_percentages = zip(*sorted(percentage_pairs, key=lambda x: x[1], reverse=True))

            num_plots = len(sorted_landuse_columns)
            num_rows = math.floor(math.sqrt(num_plots))
            num_cols = math.ceil(num_plots / num_rows) if num_rows > 0 else 1  # Ensure num_cols is defined even if num_rows is 0
            fig, axes = plt.subplots(num_rows, num_cols, figsize=(2 * num_cols, 2 * num_rows), sharex=True, sharey=True)

            cmap = plt.cm.get_cmap(cmap).copy()
            cmap.set_under('white', alpha=0)
            
            # Function to split title into two lines
            def split_title(title):
                words = title.split()
                split_point = len(words) // 2
                return ' '.join(words[:split_point]) + '\n' + ' '.join(words[split_point:])

            for i, (col, percent) in enumerate(zip(sorted_landuse_columns, sorted_percentages)):
                ax = axes.flatten()[i]
                sub_agg_ddb_merged_gdf.plot(column=col, ax=ax, cmap=cmap, vmin=0.01, vmax=1)
                dissolved_basin.plot(ax=ax, edgecolor='black', linewidth=1, facecolor='none')
                # Convert title to two lines
                title = split_title(col)
                ax.set_title(title, fontsize=font_size, ha='center')
                #ax.set_title(col, fontsize=font_size, ha='center')
                ax.text(text_location[0], text_location[1], f"{percent}%", transform=ax.transAxes, fontsize=font_size,
                        verticalalignment='top', bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.5))
                ax.set_xticks([])
                ax.set_yticks([])

            for i in range(num_plots, axes.flatten().shape[0]):
                axes.flatten()[i].set_visible(False)

        elif len(dims) == 2 and dims[1] == soldim:  # Case where the variable depends on subbasin and NSOL
            print ("len(dims) == 2 and dims[1] ==  soldim")
            nsol = dataset.dimensions['nsol'].size
            soil_layer_classes = [f'{variable_name} Layer {i+1}' for i in range(nsol)]
            df = pd.DataFrame({soil_layer: variable_data[:, i] for i, soil_layer in enumerate(soil_layer_classes)}, index=subbasin_ids)

            # Check for duplicate columns before merging
            duplicate_columns = set(df.columns).intersection(set(sub_agg_gdf.columns))
            if duplicate_columns:
                raise ValueError(f"Duplicate columns found during merge: {duplicate_columns}")

            sub_agg_ddb_merged_gdf = sub_agg_gdf.merge(df, left_on=comid_var, right_index=True, how='left')

            num_plots = len(soil_layer_classes)
            num_rows = math.floor(math.sqrt(num_plots))
            num_cols = math.ceil(num_plots / num_rows) if num_rows > 0 else 1  # Ensure num_cols is defined even if num_rows is 0
            fig, axes = plt.subplots(num_rows, num_cols, figsize=(3 * num_cols, 3 * num_rows), sharex=True, sharey=True)

            cmap = plt.cm.get_cmap(cmap).copy()
            cmap.set_under('white', alpha=0)

            for i, col in enumerate(soil_layer_classes):
                ax = axes.flatten()[i]
                sub_agg_ddb_merged_gdf.plot(column=col, ax=ax, cmap=cmap, vmin=0.01, vmax=100)
                dissolved_basin.plot(ax=ax, edgecolor='black', linewidth=1)
                ax.set_title(col, fontsize=font_size, ha='center')
                ax.set_xticks([])
                ax.set_yticks([])
                percentage = variable_data[:, i].mean()
                ax.text(text_location[0], text_location[1], f"{percentage:.2f}%", transform=ax.transAxes, fontsize=font_size,
                        verticalalignment='top', bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.5))

            for i in range(num_plots, axes.flatten().shape[0]):
                axes.flatten()[i].set_visible(False)

        else:
            raise ValueError(f"Unsupported variable dimensions for plotting: {dims}")

    if subplot_adjustments:
        fig.subplots_adjust(**subplot_adjustments)

    # Set default values for vmin and vmax
    if len(dims) == 1:
        vmin, vmax = sub_agg_ddb_merged_gdf[variable_name].min(), sub_agg_ddb_merged_gdf[variable_name].max()
    else:
        vmin = 0.01
        vmax = 1 if len(dims) == 2 and dims[1] == grudim else 100  # Set vmax based on the condition

    # Colorbar setup
    cbar_ax = fig.add_axes(cbar_location)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=mcolors.Normalize(vmin=vmin, vmax=vmax))
    sm._A = []
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=10)
    cbar.set_label(f'{variable_name} ({variable_units})', fontsize=font_size)

    plt.savefig(save_path)
    plt.show()

