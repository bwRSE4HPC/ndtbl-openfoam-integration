from pathlib import Path

import numpy as np

import ndtbl
import sys


# --- USER INPUT ---
# Anzahl der Stützstellen
RESOLUTION = 101
TABLE_SUFFIX = "f"


def dtype_from_suffix(suffix):
    if suffix in ("f", "float", "float32"):
        return np.float32
    if suffix in ("d", "double", "float64"):
        return np.float64

    raise ValueError("table suffix must be one of: f, d, float32, float64")


def generate_table(res, dtype, output_path):
    axis0 = ndtbl.UniformAxis(min=0.0, max=1.0, size=res)
    axis1 = ndtbl.UniformAxis(min=0.0, max=1.0, size=res)
    axis2 = ndtbl.UniformAxis(min=0.0, max=1.0, size=res)

    axis0values = axis0.coordinates()
    axis1values = axis1.coordinates()
    axis2values = axis2.coordinates()

    grid0, grid1, grid2 = np.meshgrid(
        axis0values,
        axis1values,
        axis2values,
        indexing="ij",
    )

    value0 = 0.25 * (grid0 + grid1 + grid2)
    value1 = 0.25 * (np.pow(grid0, 2) + np.pow(grid1, 2) + np.pow(grid2, 2))
    value2 = 0.25 * (np.pow(grid0, 3) + np.pow(grid1, 3) + np.pow(grid2, 3))

    values = np.stack((value0, value1, value2), axis=-1).astype(dtype)

    group = ndtbl.FieldGroup(
        axes=(axis0, axis1, axis2),
        field_names=("Table1", "Table2", "Table3"),
        values=values,
    )

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    ndtbl.write_group(output_path, group, max_size_mib=9000)

if __name__ == "__main__":
    res = int(sys.argv[1]) if len(sys.argv) > 1 else RESOLUTION
    suffix = sys.argv[2] if len(sys.argv) > 2 else TABLE_SUFFIX
    output_path = sys.argv[3] if len(sys.argv) > 3 else "constant/Tables.ndtbl"
    dtype = dtype_from_suffix(suffix)

    print("\nStart generating Tables")
    print("\nResolution: ", res)
    print("\nTable suffix: ", suffix)
    print("\nValue dtype: ", np.dtype(dtype))
    print("\nOutput path: ", output_path)
    generate_table(res, dtype, output_path)
    print("\nSuccess!")
