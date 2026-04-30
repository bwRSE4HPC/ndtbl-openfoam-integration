from pathlib import Path

import numpy as np

import ndtbl
import sys


# --- USER INPUT ---
# Anzahl der Stützstellen
RESOLUTION = 21

def generate_table(res):
    axis0 = ndtbl.UniformAxis(min=0.0, max=1.0, size=res)
    axis1 = ndtbl.UniformAxis(min=0.0, max=1.0, size=res)
    axis2 = ndtbl.UniformAxis(min=0.0, max=1.0, size=res)
    axis3 = ndtbl.UniformAxis(min=0.0, max=1.0, size=res)

    axis0values = axis0.coordinates()
    axis1values = axis1.coordinates()
    axis2values = axis2.coordinates()
    axis3values = axis3.coordinates()

    grid0, grid1, grid2, grid3 = np.meshgrid(
        axis0values,
        axis1values,
        axis2values,
        axis3values,
        indexing="ij",
    )

    value0 = 0.25 * np.ones(shape=(res, res, res, res))
    value1 = 0.25 * (grid0 + grid1 + grid2 + grid3)
    value2 = 0.25 * (np.pow(grid0, 2) + np.pow(grid1, 2) + np.pow(grid2, 2) + np.pow(grid3, 2))
    value3 = 0.25 * (np.pow(grid0, 3) + np.pow(grid1, 3) + np.pow(grid2, 3) + np.pow(grid3, 3))

    values = np.stack((value0, value1, value2, value3), axis=-1).astype(np.float64)

    group = ndtbl.FieldGroup(
        axes=(axis0, axis1, axis2, axis3),
        field_names=("Table1", "Table2", "Table3", "Table4"),
        values=values,
    )

    output_path = Path("constant/tables.ndtbl")

    ndtbl.write_group(output_path, group)

if __name__ == "__main__":
    res = int(sys.argv[1]) if len(sys.argv) > 1 else RESOLUTION
    print("\nStart generating Tables")
    print("\nResolution: ", res)
    generate_table(res)
    print("\nSuccess!")
