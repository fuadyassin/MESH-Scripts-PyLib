"""
r2c_utils.py

Utility functions to read, process, and write EnSim‐format .r2c ASCII files,
reorder MATLAB GRU data arrays, convert to 2D grids, and clean attribute names.

Utility functions to read, process, and write EnSim-format .r2c ASCII files,
reorder MATLAB GRU data arrays, convert to 2D grids, clean attribute names,
and produce small-multiples georeferenced plots (absolute or difference).
"""

import numpy as np
import re
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
from math import ceil
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def read_r2c_file(file_path):
    """
    Read an EnSim‐format .r2c ASCII file and return header info, attribute metadata,
    and the 3D data array.

    Parameters
    ----------
    file_path : str
        Path to the .r2c file.

    Returns
    -------
    header_info : dict[str, str]
        Mapping from header keys (e.g. ':xCount', ':yOrigin') to their values.
    attributes : dict[int, dict]
        Mapping from attribute ID → {'name': str, 'type': str|None, 'units': str|None}.
    data_matrix : np.ndarray, shape (n_attributes, n_rows, n_cols)
        The numeric grid values, reshaped according to header counts.
    """
    header_info = {}
    attributes  = {}
    data_values = []

    # ─── 1) Parse header ─────────────────────────────────────────────────────────
    with open(file_path, 'r') as f:
        line = f.readline().strip()
        while ':EndHeader' not in line:
            if line.startswith(':'):
                parts = line.split()
                key   = parts[0]
                if 'Attribute' not in key:
                    header_info[key] = ' '.join(parts[1:])
                else:
                    attr_id   = int(parts[1])
                    meta_dict = attributes.setdefault(
                        attr_id, {'name': None, 'type': None, 'units': None}
                    )
                    if 'AttributeName' in key:
                        meta_dict['name'] = ' '.join(parts[2:])
                    elif 'AttributeType' in key:
                        meta_dict['type'] = ' '.join(parts[2:])
                    elif 'AttributeUnits' in key:
                        meta_dict['units'] = ' '.join(parts[2:])
            line = f.readline().strip()

    # ─── 2) Read numeric grid values ───────────────────────────────────────────────
    with open(file_path, 'r') as f:
        # skip header
        for line in f:
            if line.strip() == ':EndHeader':
                break
        # read values
        for line in f:
            data_values.extend([float(tok) for tok in line.split()])

    # ─── 3) Reshape into (attrs, rows, cols) ──────────────────────────────────────
    n_cols       = int(header_info[':xCount'])
    n_rows       = int(header_info[':yCount'])
    n_attributes = len(attributes)
    expected = n_attributes * n_rows * n_cols
    if len(data_values) != expected:
        raise ValueError(f"Expected {expected} values, got {len(data_values)}")

    data_matrix = np.array(data_values, dtype=float)
    data_matrix = data_matrix.reshape((n_attributes, n_rows, n_cols))

    return header_info, attributes, data_matrix

def reorder_attributes(attribute_data, new_order_ids, attributes):
    """
    Reorder a flat name→2D-array dict according to a list of attribute IDs.

    Parameters
    ----------
    attribute_data : dict[str, np.ndarray]
        Mapping attribute name → 2D array as read from a .r2c.
    new_order_ids : list[int]
        1-based attribute IDs in the order you want to write them.
    attributes : dict[int, dict]
        Original mapping from attribute ID → metadata dict returned by read_r2c_file.

    Returns
    -------
    ordered_attributes : dict[int, dict]
        New ID→metadata, where IDs run 1..N in the new sequence.
    ordered_attribute_data : dict[str, np.ndarray]
        Mapping attribute name → 2D array, in the same new sequence.
    """
    ordered_attributes = {
        new_idx: attributes[attr_id]
        for new_idx, attr_id in enumerate(new_order_ids, start=1)
    }
    ordered_attribute_data = {
        ordered_attributes[i]['name']: attribute_data[
            attributes[attr_id]['name']
        ]
        for i, attr_id in enumerate(new_order_ids, start=1)
    }
    return ordered_attributes, ordered_attribute_data


def reorder_matlab_data(gru_data, order_ids):
    """
    Reorder a MATLAB‐loaded GRU struct array into a specified field sequence.

    Parameters
    ----------
    gru_data : np.ndarray, dtype=object or structured
        A 1×1 struct array from scipy.io.loadmat.
    order_ids : list[int]
        1-based indices indicating the new order of fields.

    Returns
    -------
    np.ndarray of object, shape (len(order_ids),)
        Each entry is the original gru_data[field][0,0] in the new sequence.
    """
    # get original field names
    field_names = [name for name, _ in gru_data.dtype.descr]
    # map to new order
    new_fields = [field_names[i-1] for i in order_ids]
    # extract values
    reordered = np.array(
        [gru_data[field][0, 0] for field in new_fields],
        dtype=object
    )
    return reordered


def convert_to_2d(gru_data):
    """
    Convert a 1D object array of 1D arrays into a 2D float array.

    Parameters
    ----------
    gru_data : np.ndarray, dtype=object
        Object array where each element is a 1D numeric sequence.

    Returns
    -------
    np.ndarray or None
        2D array of shape (layers, length) or None on failure.
    """
    if isinstance(gru_data, np.ndarray) and gru_data.dtype == object:
        try:
            return np.vstack([np.array(layer, dtype=float) for layer in gru_data])
        except Exception as e:
            print(f"[convert_to_2d] Conversion failed: {e}")
            return None
    return gru_data


def write_new_r2c(new_file_path, header_info, ordered_attributes, attribute_data):
    """
    Write a new .r2c file from header info, attribute metadata, and grid data.

    Parameters
    ----------
    new_file_path : str
        Path to write the .r2c file.
    header_info : dict[str, str]
        Header key→value mapping (as from read_r2c_file).
    ordered_attributes : dict[int, dict]
        Attribute ID→metadata dict specifying write order.
    attribute_data : dict[str, np.ndarray]
        Mapping attribute name→2D array of shape (n_rows, n_cols).
    """
    with open(new_file_path, 'w') as out:
        # write non-attribute header lines
        for key, val in header_info.items():
            if not key.startswith(':Attribute'):
                out.write(f"{key} {val}\n")
        # write attribute metadata
        for attr_id, meta in ordered_attributes.items():
            out.write(f":AttributeName {attr_id} {meta['name']}\n")
            if meta.get('type'):
                out.write(f":AttributeType {attr_id} {meta['type']}\n")
            if meta.get('units'):
                out.write(f":AttributeUnits {attr_id} {meta['units']}\n")
        out.write(":EndHeader\n")
        # write grid rows for each attribute
        for attr_id, meta in ordered_attributes.items():
            name = meta['name']
            if name not in attribute_data:
                raise KeyError(f"Missing data for attribute '{name}'")
            grid = attribute_data[name]
            for row in grid:
                out.write(" ".join(map(str, row)) + "\n")


def clean_attribute_name(name):
    """
    Strip parentheses and trailing whitespace from an attribute name.

    Parameters
    ----------
    name : str
        Raw attribute name, possibly containing "(units)" or similar.

    Returns
    -------
    str
        Cleaned name with text in parentheses removed.
    """
    return re.sub(r"\(.*?\)", "", name).strip()

# ─── PLOTTING FUNCTIONS ──────────────────────────────────────────────────────────

def plot_diff_panel_geo(
    base_mat, scen_mat,
    header_info, attributes,
    subset_ids, skip_idxs,
    panel_title, out_file
):
    """
    Small-multiples of (scen – base) for each attribute in subset_ids.
    """
    diff_list = [scen_mat[aid-1] - base_mat[aid-1] for aid in subset_ids]
    _plot_panel_geo(
        diff_list, header_info, attributes,
        subset_ids, skip_idxs, panel_title, out_file,
        cmap='seismic_r', vmin=-1, vmax=1
    )


def plot_r2c_panel_geo(
    mat,
    header_info, attributes,
    subset_ids, skip_idxs,
    panel_title, out_file,
    cmap='viridis', vmin=None, vmax=None
):
    """
    Small-multiples of the raw .r2c layers in mat.
    """
    arr_list = [mat[aid-1] for aid in subset_ids]
    _plot_panel_geo(
        arr_list, header_info, attributes,
        subset_ids, skip_idxs, panel_title, out_file,
        cmap=cmap, vmin=vmin, vmax=vmax
    )


def _plot_panel_geo(
    arr_list, header_info, attributes,
    subset_ids, skip_idxs,
    panel_title, out_file,
    cmap, vmin, vmax
):
    """
    Core routine for both absolute and difference panels.
    """
    # grid geometry
    x0 = float(header_info[':xOrigin'])
    y0 = float(header_info[':yOrigin'])
    dx = float(header_info[':xDelta'])
    dy = float(header_info.get(':yDelta', header_info[':xDelta']))
    nx = int(header_info[':xCount'])
    ny = int(header_info[':yCount'])

    lons = x0 + np.arange(nx)*dx
    lats = y0 + np.arange(ny)*dy
    lon_edges = np.r_[lons - dx/2, lons[-1] + dx/2]
    lat_edges = np.r_[lats - dy/2, lats[-1] + dy/2]

    # shapefile (SRB sub-drainages)
    shp = gpd.read_file('/home/fuaday/matlablaptop/SRBshape/SaskRB_SubDrainage2.shp')

    # build panels
    titles = [clean_attribute_name(attributes[aid]['name']) for aid in subset_ids]
    panels = [
        (t,a)
        for idx,(t,a) in enumerate(zip(titles,arr_list))
        if idx not in skip_idxs
    ]

    n    = len(panels)
    cols = 3
    rows = ceil(n/cols)
    fig,axes = plt.subplots(rows,cols,
                           figsize=(2.5*cols, 2.3*rows),
                           constrained_layout=True,
                           facecolor='white')
    axes = axes.flatten()

    for ax,(title,arr) in zip(axes, panels):
        arrp = np.where(arr==0, np.nan, arr)
        mesh = ax.pcolormesh(lon_edges, lat_edges, arrp,
                             cmap=cmap, vmin=vmin, vmax=vmax,
                             shading='auto')
        shp.boundary.plot(ax=ax, edgecolor='black', linewidth=0.2)
        ax.set_xlim(lon_edges[0], lon_edges[-1])
        ax.set_ylim(lat_edges[0], lat_edges[-1])
        ax.grid(True, which='major',
                color='gray', linestyle='--',
                linewidth=0.5, alpha=0.7)

        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Longitude", fontsize=7)
        ax.set_ylabel("Latitude", fontsize=7)
        ax.tick_params(axis='both', labelsize=7)

        # inset boxplot
        flat = arr.flatten()
        flat = flat[~np.isnan(flat)]
        flat = flat[flat!=0]

        inset = ax.inset_axes([0.30, 0.11, 0.65, 0.11])
        if flat.size:
            inset.boxplot(
                flat, vert=False, notch=False, whis=1.5, widths=0.5,
                patch_artist=True, showcaps=True, showfliers=True,
                showmeans=False, meanline=False,
                boxprops=dict(facecolor='lightgray', edgecolor='blue'),
                whiskerprops=dict(color='black', linewidth=1.0),
                capprops=dict(color='black', linewidth=1.0),
                medianprops=dict(color='red', linewidth=1.2),
                flierprops=dict(marker='+', color='black', alpha=0.6, markersize=3),
                manage_ticks=True
            )
            lo,hi = np.percentile(flat,[1,99])
            inset.set_xlim(lo,hi)
            q1,q2,q3 = np.percentile(flat,[1,50,99])
            inset.set_xticks([q1,q2,q3])
            inset.set_xticklabels([f"{q1:.2f}",f"{q2:.2f}",f"{q3:.2f}"], fontsize=6)
            inset.tick_params(axis='x', pad=1, labelsize=6, colors='blue')
        else:
            inset.text(0.5,0.5,"no change",ha='center',va='center',fontsize=6)
            inset.set_xticks([])

        inset.patch.set_facecolor('none')
        inset.set_yticks([])
        inset.set_frame_on(False)

    # disable unused axes
    for ax in axes[n:]:
        ax.axis('off')

    # shared colorbar
    cbar = fig.colorbar(mesh, ax=axes[:n].tolist(),
                        orientation='horizontal',
                        fraction=0.02, pad=0.04)
    cbar.set_label(f"{panel_title} (Δ fractional cover)")

    fig.suptitle(panel_title, fontsize=16, y=1.02)
    fig.savefig(out_file, dpi=300)
    plt.show()