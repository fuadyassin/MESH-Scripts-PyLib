import os
import pandas as pd
def generate_mesh_class_ini_from_excel(excel_file, output_file="MESH_class_output.ini", num_cels=None, lat=None, lon=None, sheet_name='class_ini'):
    """
    -----------------------------------------------------------------------------------------
    Description:
    -----------------------------------------------------------------------------------------
    This function processes an Excel file containing land cover parameter values and 
    generates a MESH-compatible `.ini` file. The function is designed to extract, 
    format, and write the required parameters for active land covers into the 
    output file, ensuring compatibility with the MESH hydrological model.
    The function includes detailed steps for processing the Excel file, 
    extracting relevant data, and formatting it into the required structure for MESH. 
    It also validates the input data and ensures that all necessary parameters are present.
    The existing description below provides further details about the function's 
    parameters, processing steps, and output format.

    ------------
    Parameters:
    ------------
    excel_file : str
        Path to the Excel file containing parameter values.
    output_file : str
        Path to the output `.ini` file where processed values will be written.
    num_cels : int
        Number of grid cells in the model domain.
    lat : float
        Latitude of the location.
    lon : float
        Longitude of the location.

    Overview:
    ------------
    This function extracts land cover parameter values from an Excel file and 
    writes them into a MESH-compatible `.ini` file. Only active land covers are 
    included, as indicated by the 'status' row in the Excel sheet.

    Processing Steps:
    -----------------
    1. Load the Excel file and normalize column names.
         - The function reads the Excel file using `pandas` and converts all column names to lowercase for consistency.

    2. Identify active land covers (status > 0):
       - The function identifies the row labeled 'status' in the Excel file, which contains numerical values for each land cover column.
       - These numerical values indicate the activation status of each land cover. A value greater than 0 signifies that the land cover is active.
       - The function filters out inactive land covers (status <= 0) and retains only the active ones.
       - The active land covers are then sorted in ascending order based on their status values to ensure a consistent and logical ordering in the output file.

    3. Verify required rows such as 'colum'.

    4. Extract vegetation and land cover parameters.

    5. Write formatted values into an `.ini` file with the required MESH structure.

    Output Format:
    --------------
    The output file is structured to be compatible with MESH, including:
    - Header defining basin information.
    - Land cover-specific vegetation and hydrological parameters.
    - Footer containing model time initialization values.

    File Structure:
    ----------------
    The output file consists of:

    1. **Header Information**: Includes metadata such as location, author, and details.

    2. **Land Cover Blocks**: Each selected land cover is written separately, including:
       - Vegetation parameters (written in pairs)
       - One-to-One parameter assignments (written in pairs)
       - Multi-value parameter assignments (written in structured format)
       
    3. **Final Footer**: Contains three mandatory lines required for MESH processing.

    Example Usage:
    ----------------
    >>> # The excel database can be found in the MESHpyPreProcessing folder.
    >>> # Directly dowolnload the excel file from the MESHpyPreProcessing folder.
    >>> import requests
    >>> url = "https://raw.githubusercontent.com/MESH-Scripts-PyLib/MESH-Scripts-PyLib/main/MESHpyPreProcessing/meshparametersvalues2.xlsx"
    >>> # Local path to save the file
    >>> local_path = "D:/Coding/GitHub/Repos/MESH-Scripts-PyLib/MESHpyPreProcessing/meshparametersvalues2.xlsx"
    >>> # Send a GET request to download the file
    >>> response = requests.get(url)
    >>> # Check if the request was successful (status code 200)
    >>> if response.status_code == 200: 
    >>>     # Write the content to a local file
    >>>     with open(local_path, 'wb') as file:
    >>>         file.write(response.content)
    >>>     print(f"File downloaded successfully and saved to {local_path}")
    >>> else:
    >>>     print(f"Failed to download file. Status code: {response.status_code}")

    >>> # Example usage of the function
    >>> from MESHpyPreProcessing.generate_mesh_class_ini_from_excel import generate_mesh_class_ini_from_excel as gen_classini
    >>> gen_classini(
    ...     excel_file="D:/Coding/GitHub/Repos/MESH-Scripts-PyLib/MESHpyPreProcessing/meshparametersvalues2.xlsx",
    ...     output_file="MESH_output2.ini",
    ...     num_cels=7408,
    ...     lat=53.18,
    ...     lon=-99.25
    ... )

    Example Output:
    ----------------

    .. code-block::  text
    
        Basin
        Author
        Org
        53.18     -99.25     40.00     40.00     50.00   -1.0    1 7408    5
        0.000      0.000      0.000      0.000      1.000      0.000      0.000      0.000      0.000             # fcan, lamx
        0.000      0.000      0.000      0.000     -5.352      0.000      0.000      0.000      0.000             # lnz, lamn
        0.000      0.000      0.000      0.000      0.400      0.000      0.000      0.000      0.000             # alvc, cmas
        0.000      0.000      0.000      0.000      0.250      0.000      0.000      0.000      0.000             # alir, root
        0.000      0.000      0.000      0.000                 0.000      0.000      0.000      0.000             # rsmn, qa50
        0.000      0.000      0.000      0.000                 0.000      0.000      0.000      0.000             # vpda, vpdp
        0.000      0.000      0.000      0.000                 0.000      0.000      0.000      0.000             # psga, psgb
        1.000      2.000      1.000     50.000  # drn, sdep, fare, dden
        0.100      0.100      0.388      0.070      8.000     Barren  # xslp, grkf, man, wfci, mid, name
        73.000     45.000     66.000  # sand
        6.400     29.600     20.440  # clay
        7.400      0.000      0.000  # org
        4.000      2.000      1.000        4.000        0.000        4.000  # soit, cant, snot, pndt
        0.250      0.250      0.250        0.000      0.000      0.000        0.000  # soiwf, soiif, pond
        0.000        0.000        0.000        0.500      100.000        1.000  # rcan, scan, sno, albs, rho, gro

        0.000      1.000      0.000      0.000      0.000      0.000      4.148      0.000      0.000             # fcan, lamx
        0.000      0.475      0.000      0.000      0.000      0.000      0.512      0.000      0.000             # lnz, lamn
        0.000      0.049      0.000      0.000      0.000      0.000     23.583      0.000      0.000             # alvc, cmas
        0.000      0.257      0.000      0.000      0.000      0.000      1.249      0.000      0.000             # alir, root
        0.000    110.930      0.000      0.000                 0.000     37.506      0.000      0.000             # rsmn, qa50
        0.000      0.427      0.000      0.000                 0.000      0.765      0.000      0.000             # vpda, vpdp
        0.000    100.000      0.000      0.000                 0.000      5.000      0.000      0.000             # psga, psgb
        1.000      2.000      1.000     50.000  # drn, sdep, fare, dden
        0.040      0.100      0.119      0.100      2.000   MixedForest  # xslp, grkf, man, wfci, mid, name
        75.000     66.000     66.000  # sand
        5.000     30.000     17.000  # clay
        5.000      0.000      0.000  # org
        4.000      2.000      1.000        4.000        0.000        4.000  # soit, cant, snot, pndt
        0.250      0.250      0.250        0.000      0.000      0.000        0.000  # soiwf, soiif, pond
        0.000        0.000        0.000        0.500      100.000        1.000  # rcan, scan, sno, albs, rho, gro

        1.000      0.000      0.000      0.000      0.000      1.162      0.000      0.000      0.000             # fcan, lamx
        -0.187      0.000      0.000      0.000      0.000      0.509      0.000      0.000      0.000             # lnz, lamn
        0.046      0.000      0.000      0.000      0.000     13.624      0.000      0.000      0.000             # alvc, cmas
        0.154      0.000      0.000      0.000      0.000      1.768      0.000      0.000      0.000             # alir, root
        236.390      0.000      0.000      0.000                27.913      0.000      0.000      0.000             # rsmn, qa50
        0.502      0.000      0.000      0.000                 0.812      0.000      0.000      0.000             # vpda, vpdp
        100.000      0.000      0.000      0.000                 5.000      0.000      0.000      0.000             # psga, psgb
        1.000      2.000      1.000     50.000  # drn, sdep, fare, dden
        0.040      0.122      0.151      1.380      1.000     Forest  # xslp, grkf, man, wfci, mid, name
        73.200     43.300     72.903  # sand
        5.200     29.100     21.362  # clay
        5.000      0.000      0.000  # org
        4.000      2.000      1.000        4.000        0.000        4.000  # soit, cant, snot, pndt
        0.250      0.250      0.250        0.000      0.000      0.000        0.000  # soiwf, soiif, pond
        0.000        0.000        0.000        0.500      100.000        1.000  # rcan, scan, sno, albs, rho, gro

        0.000      0.000      1.000      0.000      0.000      0.000      0.000      5.807      0.000             # fcan, lamx
        0.000      0.000     -2.541      0.000      0.000      0.000      0.000      0.000      0.000             # lnz, lamn
        0.000      0.000      0.061      0.000      0.000      0.000      0.000      1.923      0.000             # alvc, cmas
        0.000      0.000      0.298      0.000      0.000      0.000      0.000      1.200      0.000             # alir, root
        0.000      0.000     63.107      0.000                 0.000      0.000     22.541      0.000             # rsmn, qa50
        0.000      0.000      0.568      0.000                 0.000      0.000      0.804      0.000             # vpda, vpdp
        0.000      0.000    100.000      0.000                 0.000      0.000      5.000      0.000             # psga, psgb
        1.000      2.000      1.000     40.000  # drn, sdep, fare, dden
        0.010      0.503      0.164      1.432      4.000       Crop  # xslp, grkf, man, wfci, mid, name
        50.000     55.000     60.000  # sand
        35.000     35.000     35.000  # clay
        15.000     10.000      5.000  # org
        4.000      2.000      1.000        4.000        0.000        4.000  # soit, cant, snot, pndt
        0.250      0.250      0.250        0.000      0.000      0.000        0.000  # soiwf, soiif, pond
        0.000        0.000        0.000        0.500      100.000        1.000  # rcan, scan, sno, albs, rho, gro

        0.000      0.000      0.000      1.000      0.000      0.000      0.000      0.000      2.000             # fcan, lamx
        0.000      0.000      0.000     -2.806      0.000      0.000      0.000      0.000      0.500             # lnz, lamn
        0.000      0.000      0.000      0.050      0.000      0.000      0.000      0.000      1.200             # alvc, cmas
        0.000      0.000      0.000      0.340      0.000      0.000      0.000      0.000      1.200             # alir, root
        0.000      0.000      0.000    100.000                 0.000      0.000      0.000     30.000             # rsmn, qa50
        0.000      0.000      0.000      0.500                 0.000      0.000      0.000      1.000             # vpda, vpdp
        0.000      0.000      0.000    100.000                 0.000      0.000      0.000      5.000             # psga, psgb
        1.000      2.000      1.000     50.000  # drn, sdep, fare, dden
        0.010      0.200      0.233      0.100      5.000      Grass  # xslp, grkf, man, wfci, mid, name
        72.100     49.700     70.020  # sand
        5.000     21.400     23.682  # clay
        7.800      0.000      0.000  # org
        4.000      2.000      1.000        4.000        0.000        4.000  # soit, cant, snot, pndt
        0.250      0.250      0.250        0.000      0.000      0.000        0.000  # soiwf, soiif, pond
        0.000        0.000        0.000        0.500      100.000        1.000  # rcan, scan, sno, albs, rho, gro

            0         0         0         0                                  20
            0         0         0         0                                  21
            0         0         0         0                                  22 IHOUR/IMINS/IJDAY/IYEAR

    """
    if not os.path.isfile(excel_file):
        raise FileNotFoundError(f"Excel file not found: {excel_file}")

    if any(v is None for v in [num_cels, lat, lon]):
        raise ValueError("You must provide values for num_cels, lat, and lon.")

    # Load Excel file
    df = pd.read_excel(excel_file, sheet_name=sheet_name)
    df.columns = df.columns.str.lower()  # Normalize column names
    df['par'] = df['par'].str.lower()    # Normalize parameter names

    # Extract and sort active land covers based on status
    status_row = pd.to_numeric(df[df['par'] == 'status'].iloc[0], errors='coerce')
    status_row = status_row.dropna()
    active_land_covers = status_row[status_row > 0].sort_values().index.tolist()

    # Ensure 'colum' row exists
    colum_row = df[df['par'] == 'colum']
    if colum_row.empty:
        raise ValueError("The 'colum' row is missing in the provided Excel file.")

    land_cover_count = len(active_land_covers)

    with open(output_file, 'w') as f:
        # Write header
        f.write("Basin\n")
        f.write("Author\n")
        f.write("Org\n")
        f.write(f"     {lat}     {lon}     40.00     40.00     50.00   -1.0    1 {num_cels}    {land_cover_count}\n")

        # Define vegetation columns
        vegetation_cols = ['v_nforest', 'v_bforest', 'v_crop', 'v_grass', 'v_bare']
        empty_space_params = {'lamx', 'lamn', 'cmas', 'root', 'qa50', 'vpdp', 'psgb', 'psga', 'vpda', 'rsmn'}
        empty_space = " " * 8  # Fixed length of empty space

        for lc in active_land_covers:
            lc_lower = lc.lower()
            if lc_lower not in df.columns:
                raise ValueError(f"Land cover '{lc}' is not found in the Excel columns.")

            # Block 1: Vegetation Parameters
            vegetation_pairs = [
                ('fcan', 'lamx'),
                ('lnz', 'lamn'),
                ('alvc', 'cmas'),
                ('alir', 'root'),
                ('rsmn', 'qa50'),
                ('vpda', 'vpdp'),
                ('psga', 'psgb')
            ]

            for pair in vegetation_pairs:
                values_pair = []
                for param in pair:
                    values = {col: f"{0.000:8.3f}" for col in vegetation_cols}  # Default is zero

                    if param in df['par'].values:
                        if param in empty_space_params:
                            values['v_bare'] = empty_space
                        assigned_col = colum_row[lc_lower].values[0].lower()
                        if assigned_col in vegetation_cols:
                            try:
                                if assigned_col == 'v_bare' and param in empty_space_params:
                                    print(values['v_bare'])
                                else:
                                    values[assigned_col] = f"{float(df[df['par'] == param][lc_lower].values[0]):8.3f}"
                            except (ValueError, IndexError):
                                values[assigned_col] = f"{0.000:8.3f}"
                    values_pair.append(values)

                f.write("  " + "   ".join(values_pair[i][col] for i in range(len(pair)) for col in vegetation_cols) + f"  # {', '.join(pair)}\n")

            # Block 2: One-to-One Parameter Assignments
            one_to_one_pairs = [
                ('drn', 'sdep', 'fare', 'dden'),
                ('xslp', 'grkf', 'man', 'wfci', 'mid', 'name')
            ]
            for pair in one_to_one_pairs:
                values_pair = []
                for param in pair:
                    if param in df['par'].values:
                        try:
                            if param == "name":
                                param_value = str(df[df['par'] == param][lc_lower].values[0]) if lc_lower in df.columns else "N/A"
                            else:
                                param_value = float(df[df['par'] == param][lc_lower].values[0]) if lc_lower in df.columns else 0.000
                        except (ValueError, IndexError):
                            param_value = "N/A" if param == "name" else 0.000
                        values_pair.append(f"{param_value:8.3f}" if param != "name" else f"{param_value:>8}")
                    else:
                        values_pair.append(f"{0.000:8.3f}")
                f.write("  " + "   ".join(values_pair) + f"  # {', '.join(pair)}\n")

            f.write("\n")  # Add a blank line to separate land cover blocks

            # Block 3: Multi-Value Assignments
            multi_value_pairs = [
                ('sand',),
                ('clay',),
                ('org',),
                ('soit', 'cant', 'snot', 'pndt'),
                ('soiwf', 'soiif', 'pond'),
                ('rcan', 'scan', 'sno', 'albs', 'rho', 'gro')
            ]
            for pair in multi_value_pairs:
                values_pair = []
                for param in pair:
                    if param in df['par'].values:
                        raw_value = df[df['par'] == param][lc_lower].values[0]
                        if isinstance(raw_value, str) and '{' in raw_value:
                            values = raw_value.strip('{}').split(',')
                            try:
                                values = [float(v) for v in values]
                            except ValueError:
                                values = [0.000] * len(values)  # Default in case of parsing error
                        else:
                            try:
                                values = [float(raw_value)]
                            except ValueError:
                                values = [0.000]
                        values_pair.append("  " + "   ".join(f"{v:8.3f}" for v in values))
                    else:
                        values_pair.append("  " + f"{0.000:8.3f}")
                f.write("  " + "   ".join(values_pair) + f"  # {', '.join(pair)}\n")

            f.write("\n")  # Add a blank line to separate land cover blocks

        # Footer
        f.write("         0         0         0         0                                  20\n")
        f.write("         0         0         0         0                                  21\n")
        f.write("         0         0         0         0                                  22 IHOUR/IMINS/IJDAY/IYEAR\n")

    print(f"MESH parameter file '{output_file}' created successfully!")

# Example usage
#from MESHpyPreProcessing.generate_mesh_class_ini_from_excel import generate_mesh_class_ini_from_excel as gen_classini
#gen_classini(excel_file="D:/Coding/GitHub/Repos/MESH-Scripts-PyLib/MESHpyPreProcessing/meshparametersvalues2.xlsx", output_file="MESH_output2.ini", num_cels=7408, lat=53.18, lon=-99.25, sheet_name='class_ini')
# generate_mesh_class_ini_from_excel(
#     excel_file="D:/Coding/GitHub/Repos/MESH-Scripts-PyLib/MESHpyPreProcessing/meshparametersvalues2.xlsx",
#     output_file="MESH_output2.ini",
#     num_cels=7408,
#     lat=53.18,
#     lon=-99.25
#     sheet_name='class_ini'
# )