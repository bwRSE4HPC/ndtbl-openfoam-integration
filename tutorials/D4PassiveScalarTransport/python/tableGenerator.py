import numpy as np
import sys

# --- USER INPUT ---
# Anzahl der Stützstellen
RESOLUTION = 21

# Optionen für MODE:
# "constant":   All tables are filled with the constant value 1
# "linear":     (p1 + p2 + p3 + p4) / 4.0
# "quadratic":  (p1**2 + p2**2 + p3**2 + p4**2) / 4.0
# "cubic":      (p1**3 + p2**3 + p3**3 + p4**3) / 4.0
# "mix":        Table1 -> constant, Table2 -> linear, Table3 -> quadratic, Table4 -> cubic
MODE = "mix"

def get_value(p, mode, table_idx):
    current_mode = mode
    if mode == "mix":
        modes = ["constant", "linear", "quadratic", "cubic"]
        current_mode = modes[table_idx % 4]

    p1, p2, p3, p4 = p

    if current_mode == "constant":
        return 1.0
    elif current_mode == "linear":
        return float(p1 + p2 + p3 + p4) / 4.0
    elif current_mode == "quadratic":
        return float(p1**2 + p2**2 + p3**2 + p4**2) / 4.0
    elif current_mode == "cubic":
        return float(p1**3 + p2**3 + p3**3 + p4**3) / 4.0
    return 0.0

def format_openfoam_list(data, indent=0):
    spaces = "    " * indent
    if isinstance(data, np.ndarray) and data.ndim > 1:
        content = f"{spaces}{len(data)}\n{spaces}(\n"
        content += "\n".join([format_openfoam_list(sub, indent + 1) for sub in data])
        content += f"\n{spaces})"
        return content
    else:
        inner_spaces = "    " * (indent + 1)
        vals = "\n".join([f"{inner_spaces}{v:.6g}" for v in data])
        return f"{spaces}{len(data)}\n{spaces}(\n{vals}\n{spaces})"

def generate_tables(res):

    header_template = """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           | Version:  10
     \\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       dictionary;
    location    "constant";
    object      Tables;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

{table_name}
"""
    FILENAME_TEMPLATE = "Table{}_table"
    axis = np.linspace(0, 1, res)

    for i in range(1, 5):
        table_name = FILENAME_TEMPLATE.format(i)
        print(f"Processing {table_name}...")

        data_4d = np.zeros((res, res, res, res))
        for it in np.ndindex(data_4d.shape):
            params = [axis[index] for index in it]
            data_4d[it] = get_value(params, MODE, i-1)

        file_name = "./constant/" + table_name

        with open(file_name, "w") as f:
            f.write(header_template.format(table_name=table_name))
            f.write(format_openfoam_list(data_4d))
            f.write(";")
            f.write("\n\n// ************************************************************************* //")

if __name__ == "__main__":
    res = int(sys.argv[1]) if len(sys.argv) > 1 else RESOLUTION
    print("\nStart generating Tables")
    print("\nResolution: ", res)
    print("\nMode: ", MODE)
    generate_tables(res)
    print("\nSuccess!")
