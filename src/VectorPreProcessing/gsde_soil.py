"""

Overview
========

The ``GSDESoil`` class provides a pipeline to process, clean, interpolate, and integrate soil property data
into hydrological model inputs, such as those required by MESH. It is designed to handle GSDE-derived statistics
stored in CSV files, convert them to model-ready format using weighted depth-averaging, and merge them into a 
basin shapefile based on unique identifiers (e.g., ``COMID``).

Function Descriptions
=====================

.. py:class:: GSDESoil(directory, input_basin, output_shapefile)

   Initializes the processor with input/output paths.

   :param directory: Directory containing input CSV files.
   :type directory: str
   :param input_basin: Path to the input shapefile with a COMID field.
   :type input_basin: str
   :param output_shapefile: Path where the merged shapefile will be saved.
   :type output_shapefile: str

.. py:method:: load_data(file_names, search_replace_dict=None, suffix_dict=None)

   Loads and merges soil data from multiple CSV files. Columns can be renamed using 
   search/replace rules and optionally suffixed to avoid name collisions.

   :param file_names: List of CSV file names.
   :type file_names: list
   :param search_replace_dict: Dictionary with filename as key and (search_list, replace_list) as value.
   :type search_replace_dict: dict, optional
   :param suffix_dict: Dictionary with filename as key and string suffix as value.
   :type suffix_dict: dict, optional

.. py:method:: fill_and_clean_data(exclude_cols=['COMID'], exclude_patterns=['OC', 'BD', 'BDRICM', 'BDTICM'], max_val=100)

   Cleans soil data by removing outliers, rescaling BDRICM/BDTICM, and filling missing values via forward/backward fill.

   :param exclude_cols: Columns to ignore during cleaning.
   :type exclude_cols: list
   :param exclude_patterns: Substrings used to skip certain columns during range checks.
   :type exclude_patterns: list
   :param max_val: Maximum threshold for valid data (values above this become NaN).
   :type max_val: float

.. py:method:: calculate_weights(gsde_intervals, mesh_intervals)

   Computes weights to map GSDE soil depth intervals to model mesh layers.

   :param gsde_intervals: List of tuples representing GSDE depth layers (e.g., [(0, 0.045), ...]).
   :type gsde_intervals: list of tuple
   :param mesh_intervals: List of tuples representing target model layer depths.
   :type mesh_intervals: list of tuple

.. py:method:: calculate_mesh_values(column_names)

   Applies weights to calculate layer-averaged MESH-compatible soil properties.

   :param column_names: Dictionary mapping each property (e.g., "CLAY", "OC") to its source columns.
   :type column_names: dict

.. py:method:: merge_and_save_shapefile()

   Merges the processed soil data with the input basin shapefile using ``COMID`` and saves the final output.

.. py:method:: set_coordinates(input_ddb)

   Optionally reads spatial reference (lon, lat, subbasin) from a NetCDF drainage database.

   :param input_ddb: Path to the NetCDF drainage database file.
   :type input_ddb: str

Example Usage
=============

.. code-block:: python

    from gsde_soil import GSDESoil

    # Step 1: Initialize the soil processor with paths to your directories and files
    gsde = GSDESoil(
        directory='/home/fuaday/scratch/sras-agg-model/gistool-outputs',
        input_basin='/home/fuaday/scratch/sras-agg-model/geofabric-outputs/sras_subbasins_MAF_Agg2.shp',
        output_shapefile='merged_soil_data_shapefile4.shp'
    )

    # Step 2: Define the list of input CSV files
    file_names = [
        'sras_model_stats_CLAY1.csv', 'sras_model_stats_CLAY2.csv',
        'sras_model_stats_SAND1.csv', 'sras_model_stats_SAND2.csv',
        'sras_model_stats_OC1.csv',   'sras_model_stats_OC2.csv',
        'sras_model_stats_BDRICM_M_250m_ll.csv',
        'sras_model_stats_BDTICM_M_250m_ll.csv',
        'sras_model_slope_degree.csv', 'sras_model_riv_0p1_2.csv'
    ]

    # Step 3: Prepare renaming instructions for each file (search/replace patterns)
    search_replace_dict = {
        'sras_model_stats_CLAY1.csv': (['.CLAY_depth=4.5', '.CLAY_depth=9.1000004', '.CLAY_depth=16.6', '.CLAY_depth=28.9'], ['CLAY1', 'CLAY2', 'CLAY3', 'CLAY4']),
        'sras_model_stats_CLAY2.csv': (['.CLAY_depth=49.299999', '.CLAY_depth=82.900002', '.CLAY_depth=138.3', '.CLAY_depth=229.60001'], ['CLAY5', 'CLAY6', 'CLAY7', 'CLAY8']),
        'sras_model_stats_SAND1.csv': (['.SAND_depth=4.5', '.SAND_depth=9.1000004', '.SAND_depth=16.6', '.SAND_depth=28.9'], ['SAND1', 'SAND2', 'SAND3', 'SAND4']),
        'sras_model_stats_SAND2.csv': (['.SAND_depth=49.299999', '.SAND_depth=82.900002', '.SAND_depth=138.3', '.SAND_depth=229.60001'], ['SAND5', 'SAND6', 'SAND7', 'SAND8']),
        'sras_model_stats_OC1.csv': (['.OC_depth=4.5', '.OC_depth=9.1000004', '.OC_depth=16.6', '.OC_depth=28.9'], ['OC1', 'OC2', 'OC3', 'OC4']),
        'sras_model_stats_OC2.csv': (['.OC_depth=49.299999', '.OC_depth=82.900002', '.OC_depth=138.3', '.OC_depth=229.60001'], ['OC5', 'OC6', 'OC7', 'OC8'])
    }

    # Step 4: Optionally specify suffixes to distinguish overlapping columns
    suffix_dict = {
        'sras_model_stats_BDRICM_M_250m_ll.csv': 'BDRICM',
        'sras_model_stats_BDTICM_M_250m_ll.csv': 'BDTICM'
    }

    # Step 5: Load the data, applying renaming and suffixes
    gsde.load_data(
        file_names=file_names,
        search_replace_dict=search_replace_dict,
        suffix_dict=suffix_dict
    )

    # Step 6: Clean and prepare the soil data (e.g., remove outliers, fill NaNs)
    gsde.fill_and_clean_data()

    # Step 7: Define soil profile intervals for GSDE and MESH (depths in meters)
    gsde_intervals = [(0, 0.045), (0.045, 0.091), (0.091, 0.166), (0.166, 0.289),
                      (0.289, 0.493), (0.493, 0.829), (0.829, 1.383), (1.383, 2.296)]

    mesh_intervals = [(0, 0.1), (0.1, 0.35), (0.35, 1.2), (1.2, 4.1)]

    gsde.calculate_weights(gsde_intervals, mesh_intervals)

    # Step 8: Compute mesh-compatible weighted averages of soil properties
    column_names = {
        'CLAY': ['CLAY1', 'CLAY2', 'CLAY3', 'CLAY4', 'CLAY5', 'CLAY6', 'CLAY7', 'CLAY8'],
        'SAND': ['SAND1', 'SAND2', 'SAND3', 'SAND4', 'SAND5', 'SAND6', 'SAND7', 'SAND8'],
        'OC':   ['OC1', 'OC2', 'OC3', 'OC4', 'OC5', 'OC6', 'OC7', 'OC8']
    }
    gsde.calculate_mesh_values(column_names)

    # Step 9: Merge processed soil data into the basin shapefile and save output
    gsde.merge_and_save_shapefile()
"""
import os
import numpy as np
import pandas as pd
import geopandas as gpd
from functools import reduce
import xarray as xr

class GSDESoil:
    """
    A class to process, clean, interpolate, and merge soil property data from CSV files
    with a given basin shapefile, producing model-ready soil inputs.

    Attributes
    ----------
    directory : str
        Directory containing input CSV files with soil properties.
    input_basin : str
        Path to the input basin shapefile with a 'COMID' identifier.
    output_shapefile : str
        Path to the output shapefile with processed soil attributes.
    file_paths : list
        List of full file paths for input CSVs.
    gsde_df : pandas.DataFrame
        Combined soil property table after processing.
    merged_gdf : geopandas.GeoDataFrame
        Final spatial dataset with soil properties merged to polygons.
    weights_used : list of list
        Weights used to interpolate soil layers into mesh layers.
    mesh_intervals : list of tuple
        Target depth intervals used for model input (e.g., MESH layers).
    lon : ndarray
        Longitude values loaded from a NetCDF drainage database.
    lat : ndarray
        Latitude values loaded from a NetCDF drainage database.
    segid : ndarray
        Segment IDs (e.g., subbasin or COMID) from a drainage database.
    num_soil_lyrs : int
        Number of output mesh layers.
    """
    def __init__(self, directory, input_basin, output_shapefile):
        self.directory = directory
        self.input_basin = input_basin
        self.output_shapefile = output_shapefile
        self.file_paths = []
        self.gsde_df = pd.DataFrame()
        self.merged_gdf = gpd.GeoDataFrame()
        self.weights_used = []
        self.mesh_intervals = []
        self.lon = []
        self.lat = []
        self.segid = []
        self.num_soil_lyrs = 0

    def load_data(self, file_names, search_replace_dict=None, suffix_dict=None):
        """
        Load and merge multiple CSV files into a single DataFrame. Optionally apply
        search-and-replace logic and suffixes to column names to ensure compatibility.

        Parameters
        ----------
        file_names : list of str
            List of filenames to load from the given directory.
        search_replace_dict : dict, optional
            Dictionary where keys are filenames and values are (search_list, replace_list) tuples
            used to rename columns (e.g., depth labels to CLAY1, CLAY2, etc.).
        suffix_dict : dict, optional
            Dictionary where keys are filenames and values are suffix strings
            to append to column names (useful for distinguishing overlapping variables).
        """
        self.file_paths = [os.path.join(self.directory, filename) for filename in file_names]
        self.gsde_df = self.load_and_merge_files(self.file_paths, search_replace_dict, suffix_dict)

    @staticmethod
    def load_and_merge_files(file_list, search_replace_dict=None, suffix_dict=None, key='COMID'):
        """
        Load multiple CSV files and merge them on a common key. Renames and suffixes
        column names as needed during the loading process.

        Parameters
        ----------
        file_list : list of str
            List of full CSV file paths.
        search_replace_dict : dict, optional
            Column renaming instructions for each file.
        suffix_dict : dict, optional
            Suffix strings to append to column names by file.
        key : str
            Primary key used to merge all data files (default is 'COMID').

        Returns
        -------
        pandas.DataFrame
            Merged DataFrame containing columns from all input files.
        """
        dfs = []
        for fp in file_list:
            df = pd.read_csv(fp)
            file_name = os.path.basename(fp)
            
            # Apply search and replace for column names if specified
            if search_replace_dict and file_name in search_replace_dict:
                search_list, replace_list = search_replace_dict[file_name]
                for search, replace in zip(search_list, replace_list):
                    df.columns = [col.replace(search, str(replace)) for col in df.columns]
            
            # Remove periods from column names
            df.columns = [col.replace('.', '') for col in df.columns]
            
            # Optionally append a suffix to the column names
            if suffix_dict and file_name in suffix_dict:
                suffix = suffix_dict[file_name]
                if suffix:  # Only apply if the suffix is not an empty string
                    df.columns = [f"{col}{suffix}" if col != key else col for col in df.columns]
            
            dfs.append(df)
        
        # Merge all the dataframes on the specified key
        return reduce(lambda left, right: pd.merge(left, right, on=key, how='outer'), dfs)

    def fill_and_clean_data(self, exclude_cols=['COMID'], exclude_patterns=['OC', 'BD', 'BDRICM', 'BDTICM'], max_val=100):
        """
        Clean the soil data by:
        - Replacing extreme values with NaN (based on max_val).
        - Normalizing and capping specific fields (e.g., BDRICM/BDTICM).
        - Filling missing values using forward and backward fill.

        Parameters
        ----------
        exclude_cols : list of str
            Columns to exclude from NaN replacement.
        exclude_patterns : list of str
            Column name substrings to skip when applying value caps.
        max_val : float
            Maximum valid threshold for general soil values.
        """
        for col in self.gsde_df.columns:
            if col not in exclude_cols:
                if not any(pattern in col for pattern in exclude_patterns):
                    self.gsde_df.loc[self.gsde_df[col] > max_val, col] = np.nan
            
            if 'BDRICM' in col or 'BDTICM' in col:
                self.gsde_df.loc[self.gsde_df[col] == 9999.0, col] = np.nan
                self.gsde_df[col] = self.gsde_df[col] / 100.0
                if 'BDRICM' in col:
                    self.gsde_df.loc[self.gsde_df[col] > 2.0, col] = 2.0
                    self.gsde_df.loc[self.gsde_df[col] < 0.1, col] = 0.1
                if 'BDTICM' in col:
                    self.gsde_df.loc[self.gsde_df[col] > 4.1, col] = 4.1
                    self.gsde_df.loc[self.gsde_df[col] < 0.1, col] = 0.1
                
        self.gsde_df.sort_values('COMID', inplace=True)
        self.gsde_df.fillna(method='ffill', inplace=True)
        self.gsde_df.fillna(method='bfill', inplace=True)

    def calculate_weights(self, gsde_intervals, mesh_intervals):
        """
        Calculate the contribution weights from each GSDE layer to each model-defined
        mesh layer based on depth intervals.

        Parameters:
        -----------
        gsde_intervals : list of tuple
            List of tuples representing GSDE depth layers (e.g., [(0, 0.045), ...]).
        mesh_intervals : list of tuple
            Target model layer depths (e.g., [(0, 0.1), (0.1, 0.35), ...]).
        """
        self.mesh_intervals = mesh_intervals
        self.num_soil_lyrs = len(mesh_intervals)
        weights_used = []
        for mesh_interval in mesh_intervals:
            start, end = mesh_interval
            weights = []
            for gsde_interval in gsde_intervals:
                gsde_start, gsde_end = gsde_interval
                overlap_start = max(start, gsde_start)
                overlap_end = min(end, gsde_end)
                weight = (overlap_end - overlap_start) / (end - start) if overlap_start < overlap_end else 0
                weights.append(weight)
            weights_used.append([w / sum(weights) for w in weights if sum(weights) > 0])
        self.weights_used = weights_used

    def calculate_mesh_values(self, column_names):
        """
        Apply the calculated weights to soil property columns and generate
        epth-integrated values for each mesh layer.

        Parameters:
        -----------
        column_names : dict
            Dictionary mapping each property (e.g., "CLAY", "OC") to its source columns.
            Example: {'CLAY': ['CLAY1', 'CLAY2', ...], 'OC': ['OC1', 'OC2', ...]}
        """
        for prop, cols in column_names.items():
            extracted_data = self.gsde_df[cols]
            weights_array = np.array(self.weights_used).T
            mesh_values = np.dot(extracted_data, weights_array)
            if prop == 'OC':
                mesh_values *= 0.01 * 1.72
            for i, mesh_col in enumerate([f'mesh{prop}{j+1}' for j in range(len(self.mesh_intervals))]):
                self.gsde_df[mesh_col] = mesh_values[:, i]

    def merge_and_save_shapefile(self):
        """
         Merge the processed soil data (via COMID) into the input shapefile and save the result.
         Output is a GeoDataFrame with mesh values appended as new attributes.
        """
        gdf = gpd.read_file(self.input_basin).to_crs(epsg=4326)
        gdf['COMID'] = gdf['COMID'].astype(int)
        self.merged_gdf = gdf.merge(self.gsde_df, on='COMID', how='left')
        self.merged_gdf.to_file(self.output_shapefile)

    def set_coordinates(self, input_ddb):
        """
        Set longitude and latitude values from a NetCDF drainage database.
        
        Parameters: 
        -----------
        input_ddb : str
            Path to the NetCDF drainage database file.
        """
        db = xr.open_dataset(input_ddb)
        self.lon = db.variables['lon'].values
        self.lat = db.variables['lat'].values
        self.segid = db.variables['subbasin'].values
        db.close()
