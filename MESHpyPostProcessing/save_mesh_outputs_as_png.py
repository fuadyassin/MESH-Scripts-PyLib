import os
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import netCDF4 as nc
from datetime import timedelta
import xarray as xr
import matplotlib

def save_mesh_outputs_as_png(
    shapefile_path,
    netcdf_dir,
    ddb_path,
    varnames,
    filenames,
    cbar_labels,
    outdir,
    indices_to_save,
    mode='monthly',  # can now also be 'hourly'
    domain_name='Basin',
    cmap='gnuplot2_r'
):
    font = {'family': 'DejaVu Serif', 'weight': 'bold', 'size': 24}
    matplotlib.rc('font', **font)

    os.makedirs(outdir, exist_ok=True)

    db = xr.open_dataset(ddb_path)
    segid = db['subbasin'].values
    db.close()

    df = pd.DataFrame({'ID': segid})
    shp = gpd.read_file(shapefile_path).sort_values(by='COMID').reset_index(drop=True)

    example_file = os.path.join(netcdf_dir, filenames[0])
    with nc.Dataset(example_file) as ds:
        time_var = ds.variables['time']
        time_units = time_var.units.replace('years', 'days')  # fallback fix
        times = time_var[:] * 365 if 'years' in time_var.units else time_var[:]
        calendar = getattr(time_var, 'calendar', 'standard')
        dates = nc.num2date(times, units=time_units, calendar=calendar)
        starting_date = dates[1]

    global_min_max = {}
    for i, fname in enumerate(filenames):
        with nc.Dataset(os.path.join(netcdf_dir, fname)) as ds:
            data = ds.variables[varnames[i]][1:]
            global_min_max[varnames[i]] = (np.nanmin(data), np.nanmax(data))

    for i in range(len(varnames)):
        for idx in indices_to_save:
            fig, ax = plt.subplots(figsize=(20, 20))
            with nc.Dataset(os.path.join(netcdf_dir, filenames[i])) as ds:
                values = ds.variables[varnames[i]][idx + 1, :]

            merge_df = shp.copy()
            df['value'] = values
            merged = merge_df.merge(df, left_on='COMID', right_on='ID', how='left')

            mn, mx = global_min_max[varnames[i]]

            if 'IG1' in filenames[i]: layer = 'Layer1'
            elif 'IG2' in filenames[i]: layer = 'Layer2'
            elif 'IG3' in filenames[i]: layer = 'Layer3'
            else: layer = None

            # Handle time label
            if mode == 'yearly':
                date = starting_date + timedelta(days=365 * idx)
                label = date.strftime('%Y')
            elif mode == 'monthly':
                date = starting_date + timedelta(days=30 * idx)
                label = date.strftime('%Y-%m')
            elif mode == 'daily':
                date = starting_date + timedelta(days=idx)
                label = date.strftime('%Y-%m-%d')
            elif mode == 'hourly':
                date = starting_date + timedelta(hours=idx)
                label = date.strftime('%Y-%m-%d_%H:%M')
            else:
                raise ValueError("Invalid mode. Choose from 'daily', 'monthly', 'yearly', or 'hourly'.")

            title = f'{domain_name}_MESH_{varnames[i]}{"_" + layer if layer else ""}_{label}'
            ax.set_title(title)
            ax.set_xlabel('Longitude')
            ax.set_ylabel('Latitude')

            norm = colors.Normalize(vmin=mn, vmax=mx) if mn != mx else None
            cmap_used = plt.cm.get_cmap(cmap)
            sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap_used)

            merged.plot(column='value', cmap=cmap, edgecolor='k', linewidth=0.1, ax=ax, vmin=mn, vmax=mx)

            cbar_ax = fig.add_axes([0.82, 0.2, 0.02, 0.6])
            fig.colorbar(sm, cax=cbar_ax, orientation='vertical', extend='max')
            cbar_ax.set_ylabel(cbar_labels[i])

            fig.subplots_adjust(left=0.05, right=0.80, top=0.95, bottom=0.05)
            filename = f"{varnames[i]}_{mode}_frame_{idx:03d}_{label.replace(':','-')}.png"
            fig.savefig(os.path.join(outdir, filename), dpi=300)
            plt.close(fig)



""" import os
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import netCDF4 as nc
from matplotlib.animation import FuncAnimation, PillowWriter
from datetime import timedelta
import xarray as xr
import matplotlib

def save_mesh_outputs_as_png(
    shapefile_path,
    netcdf_dir,
    ddb_path,
    varnames,
    filenames,
    cbar_labels,
    outdir,
    indices_to_save,
    mode='monthly',
    domain_name='Basin',
    cmap='gnuplot2_r'
):
    font = {'family': 'DejaVu Serif', 'weight': 'bold', 'size': 24}
    matplotlib.rc('font', **font)

    os.makedirs(outdir, exist_ok=True)

    db = xr.open_dataset(ddb_path)
    segid = db['subbasin'].values
    db.close()

    df = pd.DataFrame({'ID': segid})
    shp = gpd.read_file(shapefile_path).sort_values(by='COMID').reset_index(drop=True)

    example_file = os.path.join(netcdf_dir, filenames[0])
    with nc.Dataset(example_file) as ds:
        time_var = ds.variables['time']
        time_units = time_var.units.replace('years', 'days')
        times = time_var[:] * 365 if 'years' in time_var.units else time_var[:]
        calendar = getattr(time_var, 'calendar', 'standard')
        dates = nc.num2date(times, units=time_units, calendar=calendar)
        starting_date = dates[1]

    global_min_max = {}
    for i, fname in enumerate(filenames):
        with nc.Dataset(os.path.join(netcdf_dir, fname)) as ds:
            data = ds.variables[varnames[i]][1:]
            global_min_max[varnames[i]] = (np.nanmin(data), np.nanmax(data))

    for i in range(len(varnames)):
        for idx in indices_to_save:
            fig, ax = plt.subplots(figsize=(20, 20))
            with nc.Dataset(os.path.join(netcdf_dir, filenames[i])) as ds:
                values = ds.variables[varnames[i]][idx + 1, :]

            merge_df = shp.copy()
            df['value'] = values
            merged = merge_df.merge(df, left_on='COMID', right_on='ID', how='left')

            mn, mx = global_min_max[varnames[i]]

            if 'IG1' in filenames[i]: layer = 'Layer1'
            elif 'IG2' in filenames[i]: layer = 'Layer2'
            elif 'IG3' in filenames[i]: layer = 'Layer3'
            else: layer = None

            if mode == 'yearly':
                date = starting_date + timedelta(days=365 * idx)
                label = date.strftime('%Y')
            elif mode == 'monthly':
                date = starting_date + timedelta(days=30 * idx)
                label = date.strftime('%Y-%m')
            elif mode == 'daily':
                date = starting_date + timedelta(days=idx)
                label = date.strftime('%Y-%m-%d')
            else:
                raise ValueError("Invalid mode. Choose from 'daily', 'monthly', 'yearly'.")

            title = f'{domain_name}_MESH_{varnames[i]}{"_" + layer if layer else ""}_{label}'
            ax.set_title(title)
            ax.set_xlabel('Longitude')
            ax.set_ylabel('Latitude')

            norm = colors.Normalize(vmin=mn, vmax=mx) if mn != mx else None
            cmap_used = plt.cm.get_cmap(cmap)
            sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap_used)

            merged.plot(column='value', cmap=cmap, edgecolor='k', linewidth=0.1, ax=ax, vmin=mn, vmax=mx)

            cbar_ax = fig.add_axes([0.82, 0.2, 0.02, 0.6])
            fig.colorbar(sm, cax=cbar_ax, orientation='vertical', extend='max')
            cbar_ax.set_ylabel(cbar_labels[i])

            fig.subplots_adjust(left=0.05, right=0.80, top=0.95, bottom=0.05)
            filename = f"{varnames[i]}_{mode}_frame_{idx:03d}_{label}.png"
            fig.savefig(os.path.join(outdir, filename), dpi=300)
            plt.close(fig)
 """