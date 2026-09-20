import sys
from pathlib import Path

import numpy as np

import ndtbl


# --- USER INPUT ---
# Anzahl der Stützstellen
RESOLUTION = 101
TABLE_SUFFIX = "f"
TABLE_LAYOUT = "combined"


def dtype_from_suffix(suffix):
    if suffix in ("f", "float", "float32"):
        return np.float32
    if suffix in ("d", "double", "float64"):
        return np.float64

    raise ValueError("table suffix must be one of: f, d, float32, float64")


def split_output_path(output_path, field_name):
    output_path = Path(output_path)
    suffix = output_path.suffix or ".ndtbl"
    stem = output_path.stem if output_path.suffix else output_path.name
    return output_path.with_name(f"{stem}.{field_name}{suffix}")


def field_values(coordinates, power):
    axis0values, axis1values, axis2values = coordinates
    grid0 = axis0values[:, np.newaxis, np.newaxis]
    grid1 = axis1values[np.newaxis, :, np.newaxis]
    grid2 = axis2values[np.newaxis, np.newaxis, :]
    return 0.25 * (
        np.pow(grid0, power)
        + np.pow(grid1, power)
        + np.pow(grid2, power)
    )


def write_table(axes, field_names, values, dtype, output_path):
    group = ndtbl.FieldGroup(
        axes=axes,
        field_names=field_names,
        values=np.asarray(values, dtype=dtype),
    )

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    ndtbl.write_group(output_path, group, max_size_mib=9000)


def generate_table(res, dtype, output_path, layout="combined"):
    axis0 = ndtbl.UniformAxis(min=0.0, max=1.0, size=res)
    axis1 = ndtbl.UniformAxis(min=0.0, max=1.0, size=res)
    axis2 = ndtbl.UniformAxis(min=0.0, max=1.0, size=res)
    axes = (axis0, axis1, axis2)
    coordinates = tuple(axis.coordinates() for axis in axes)
    field_names = ("Table1", "Table2", "Table3")
    output_path = Path(output_path)

    if layout == "split":
        output_paths = []
        for power, field_name in enumerate(field_names, start=1):
            field_output_path = split_output_path(output_path, field_name)
            values = field_values(coordinates, power)[..., np.newaxis]
            write_table(axes, (field_name,), values, dtype, field_output_path)
            output_paths.append(field_output_path)
        return tuple(output_paths)

    if layout != "combined":
        raise ValueError("table layout must be one of: combined, split")

    values = np.stack(
        [field_values(coordinates, power) for power in range(1, 4)],
        axis=-1,
    )
    write_table(axes, field_names, values, dtype, output_path)
    return (output_path,)

if __name__ == "__main__":
    res = int(sys.argv[1]) if len(sys.argv) > 1 else RESOLUTION
    suffix = sys.argv[2] if len(sys.argv) > 2 else TABLE_SUFFIX
    output_path = sys.argv[3] if len(sys.argv) > 3 else "constant/Tables.ndtbl"
    layout = sys.argv[4] if len(sys.argv) > 4 else TABLE_LAYOUT
    dtype = dtype_from_suffix(suffix)

    print("\nStart generating Tables")
    print("\nResolution: ", res)
    print("\nTable suffix: ", suffix)
    print("\nValue dtype: ", np.dtype(dtype))
    print("\nTable layout: ", layout)
    output_paths = generate_table(res, dtype, output_path, layout)
    print("\nOutput paths: ", ", ".join(map(str, output_paths)))
    print("\nSuccess!")
