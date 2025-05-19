import pandas as pd

def generate_mesh_hydrology_ini_from_excel(excel_file, output_file="MeshHydrology.ini", sheet_name="hydrology_ini"):
    """
    ---------------------------------------------------------------------------------
    Description:
    ---------------------------------------------------------------------------------
    This function reads hydrology parameter values from an Excel file and generates 
    a formatted MESH Hydrology `.ini` file. The generated file is structured into 
    various parameter categories with appropriate headers, formatting, and spacing.

    ---------------------------------------------------------------------------------
    Parameters:
    ---------------------------------------------------------------------------------
    - excel_file  : str  -> Path to the input Excel file containing parameter values.
    - sheet_name  : str  -> The specific sheet in the Excel file that holds the data.
                            Default is "hydrology_ini".
    - output_file : str  -> Path to the output `.ini` file.

    ---------------------------------------------------------------------------------
    Function Overview:
    ---------------------------------------------------------------------------------
    1. Reads the Excel file and extracts three main parameter categories:
       - **Channel Routing Parameters**
       - **GRU-Independent Parameters**
       - **GRU-Dependent Parameters**
    2. Filters parameters where `status == 1` (active parameters).
    3. Sorts and formats values for proper alignment:
       - Parameter names are **left-aligned** (12 characters).
       - Values are **right-aligned** (10 characters, formatted as `.3f` if float).
    4. Writes the output `.ini` file with section headers, formatted values, 
       and comments where necessary.

    Example Hydrology Input File (Version 2.0):
    ----------------------------------------------------------------

    .. code-block:: text

        2.0: MESH Hydrology parameters input file (Version 2.0)
        !> Any line leading with '!' is ignored
        !> Spacing, whether by tabs or spaces doesn't matter
        !> Parameter/field names are not case sensitive

        ##### Option Flags #####
        ----#
            0 # Number of option flags

        ##### Channel routing parameters #####
        -----#
        4           # Number of channel routing parameters
        R2N              0.100     0.045     0.049     0.013     0.050     0.150     0.141     0.150     0.135     0.314
        R1N              0.200     0.167     0.156     0.049     0.064     0.596     0.514     0.422     0.171     0.544
        FLZ              0.000     0.000     0.000     0.000     0.000     0.000     0.000     0.000     0.000     0.000
        PWR              1.618     2.570     2.721     2.993     2.823     1.361     1.361     1.361     1.361     1.361

        ##### GRU-independent parameters #####
        -------#
        5           # Number of GRU independent hydrologic parameters
        SOIL_POR_MAX     0.200
        SOIL_DEPTH           4
        S0                   1
        T_ICE_LENS         -10
        t0_ACC               0         0         0         0         0

        ##### GRU-dependent parameters #####
        -------#
        !> Active headers in sorted order: Forest, BroadLeafForest, MixedForest, wetland, water, Shrub
        24          # Number of GRU-dependent parameters
        IWF                  1         1         1         1         1         1
        ZSNL             0.173     0.393     0.393     0.083     0.110     0.268
        ZPLS             0.023     0.097     0.097     0.090     0.090     0.149
        ZPLG             0.190     0.282     0.282     0.260     0.090     0.045
        FREZTH               0         0         0         0         0         0
        SWELIM               0         0         0         0         0         0
        SNDENLIM             0         0         0         0         0         0
        fetch              300       300       300       300       300       600
        Ht                   3         4         4     1.500     0.030     1.200
        N_S                  1         1         1         1         1         1
        A_S                  1         1         1         1         1         1
        Distrib              6         6         6         0         0         1
        CMAX             1.500     1.500     1.500     1.500     1.500     1.500
        CMIN                 0         0         0         0         0         0
        B                    3         3         3         3         3         3
        K1                   0         0         0         0         0         0
        K2                   0         0         0         0         0         0
        irflg                0         0         0         0         0         0
        irthlmin         0.500     0.500     0.500     0.500     0.500     0.500
        irignd               0         0         0         0         0         0
        irt1                 6         6         6         6         6         6
        irt2                10        10        10        10        10        10
        irijday1             0         0         0         0         0         0
        irijday2             0         0         0         0         0         0

    ---------------------------------------------------------------------------------
    File Structure:
    ---------------------------------------------------------------------------------
    - **Header section**: Explains format rules.
    - **Option Flags**: Indicates no optional flags.
    - **Channel Routing Parameters**: Defines routing-related parameters.
    - **GRU-Independent Parameters**: Contains parameters that apply to all land units.
    - **GRU-Dependent Parameters**: Contains parameters that vary by land category.

    ---------------------------------------------------------------------------------
    GRU-Dependent Parameters Selection and Ordering:
    ---------------------------------------------------------------------------------
    - The `GRU_class_dependent_active` row in the Excel sheet is used to determine 
      which GRU-dependent headers (columns) are active. Each column corresponding 
      to a GRU class (e.g., Forest, Grass, Wetland) has a value indicating its 
      active status:
    - A value > 0 indicates the column is active.
    - A value <= 0 indicates the column is inactive.
    - Active columns are sorted in ascending order based on their values in the 
      `GRU_class_dependent_active` row. This sorted order determines the sequence 
      in which GRU-dependent parameters are written to the output file.
    - The sorted active column names are included as a comment in the output file 
      for reference.

    ---------------------------------------------------------------------------------
    Example Usage:
    ---------------------------------------------------------------------------------
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

    >>> from MESHpyPreProcessing.generate_mesh_hydrology_ini_from_excel import generate_mesh_hydrology_ini_from_excel
    >>> gen_hydini = generate_mesh_hydrology_ini_from_excel
    >>> gen_hydini(
            sheet_name="hydrology_ini",
            output_file="MeshHydrology2.ini",
            excel_file="D:/Coding/GitHub/Repos/MESH-Scripts-PyLib/MESHpyPreProcessing/meshparametersvalues2.xlsx"
        )
    
    - This will generate a properly formatted `MeshHydrology.ini` file in the specified path.
    """

    # Read the Excel file and specified sheet
    data = pd.read_excel(excel_file, sheet_name=sheet_name)

    # Filter rows for 'channel_routing_header' to get channel routing parameters
    channel_routing_data = data[data['channel_routing_header'] == 'channel_routing']

    # Further filter rows where 'status' is 1 for channel routing
    active_channel_routing = channel_routing_data[channel_routing_data['status'] == 1]

    # Count the number of active channel routing parameters
    num_active_channel_routing = len(active_channel_routing)

    # Define the header content
    header = """2.0: MESH Hydrology parameters input file (Version 2.0)
!> Any line leading with '!' is ignored
!> Spacing, whether by tabs or spaces doesn't matter
!> Parameter/field names are not case sensitive
##### Option Flags #####
----#
    0 # Number of option flags
##### Channel routing parameters #####
-----#
"""

    # Add the number of active channel routing parameters
    content = header + f"{num_active_channel_routing:<12}# Number of channel routing parameters\n"

    # Write the active channel routing parameters
    for idx, row in active_channel_routing.iterrows():
        param_name = row['par']
        param_values = row['value']
        if isinstance(param_values, str) and param_values.startswith("{"):
            values_list = param_values.strip("{}").split(",")
            formatted_values = "".join([f"{float(val):>10.3f}" if '.' in val else f"{int(val):>10}" for val in values_list])
        else:
            formatted_values = f"{float(param_values):>10.3f}" if isinstance(param_values, float) else f"{int(param_values):>10}"
        content += f"{param_name:<12}{formatted_values}"
        if idx != active_channel_routing.index[-1]:
            content += "\n"

    # Add GRU-independent parameters section
    content += "\n##### GRU-independent parameters #####\n-------#\n"

    # Filter GRU-independent parameters
    gru_independent_data = data[data['channel_routing_header'] == 'GRU_class_independent']
    active_gru_independent = gru_independent_data[gru_independent_data['status'] == 1]
    num_active_gru_independent = len(active_gru_independent)

    # Add count and parameters
    content += f"{num_active_gru_independent:<12}# Number of GRU independent hydrologic parameters\n"
    for idx, row in active_gru_independent.iterrows():
        param_name = row['par']
        param_values = row['value']
        if isinstance(param_values, str) and param_values.startswith("{"):
            values_list = param_values.strip("{}").split(",")
            formatted_values = "".join([f"{float(val):>10.3f}" if '.' in val else f"{int(val):>10}" for val in values_list])
        else:
            formatted_values = f"{float(param_values):>10.3f}" if isinstance(param_values, float) else f"{int(param_values):>10}"
        content += f"{param_name:<12}{formatted_values}"
        if idx != active_gru_independent.index[-1]:
            content += "\n"

    # Add GRU-dependent parameters section
    content += "\n##### GRU-dependent parameters #####\n-------#\n"

    # Process GRU-dependent parameters
    header_row = data[data.iloc[:, 0] == 'GRU_class_dependent_header']
    if not header_row.empty:
        new_headers = header_row.iloc[0].values
        data.columns = new_headers
        data = data[data.iloc[:, 0] != 'GRU_class_dependent_header']
        active_row = data[data.iloc[:, 0] == 'GRU_class_dependent_active']
        if not active_row.empty:
            active_row = active_row.iloc[0].apply(pd.to_numeric, errors='coerce')
            active_columns = active_row[active_row > 0].index
            sorted_active_columns = active_row[active_columns].sort_values(ascending=True)
            active_header_names = ", ".join(sorted_active_columns.index)
            content += f"!> Active headers in sorted order: {active_header_names}\n"
            dependent_params = data[data.iloc[:, 0] == 'GRU_class_dependent']
            active_params = dependent_params[dependent_params['status'] == 1]
            num_active_params = len(active_params)
            content += f"{num_active_params:<12}# Number of GRU-dependent parameters\n"
            for idx, row in active_params.iterrows():
                param_name = row['par']
                values = [row[col] for col in sorted_active_columns.index]
                formatted_values = "".join([
                    f"{float(val):>10.3f}" if isinstance(val, float) or "." in str(val) else f"{int(val):>10}" 
                    for val in values
                ])
                content += f"{param_name:<12}{formatted_values}"
                if idx != active_params.index[-1]:
                    content += "\n"

    with open(output_file, 'w') as f:
        f.write(content)
# from MESHpyPreProcessing.generate_mesh_hydrology_ini_from_excel import generate_mesh_hydrology_ini_from_excel
# gen_hydini(output_file="MeshHydrology2.ini",excel_file="D:/Coding/GitHub/Repos/MESH-Scripts-PyLib/MESHpyPreProcessing/meshparametersvalues2.xlsx")
#gen_hydini(sheet_name="hydrology_ini",output_file="MeshHydrology2.ini",excel_file="D:/Coding/GitHub/Repos/MESH-Scripts-PyLib/MESHpyPreProcessing/meshparametersvalues2.xlsx")