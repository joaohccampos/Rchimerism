from pathlib import Path

import polars as pl

__all__: list[str] = ["write_rchimerism_file"]

_OUTPUT_COLUMNS: tuple[str, ...] = (
    "Dye/Sample Peak",
    "Sample File Name",
    "Marker",
    "Allele",
    "Size",
    "Height",
    "Area",
    "Data Point",
)


def write_rchimerism_file(df: pl.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    present = [col for col in _OUTPUT_COLUMNS if col in df.columns]
    ordered = df.select(present)

    with output_path.open("w", encoding="utf-8") as file_handle:
        file_handle.write("\t".join(_OUTPUT_COLUMNS) + "\t\n")
        for row in ordered.iter_rows():
            line = "\t".join("" if value is None else str(value) for value in row)
            file_handle.write(line + "\t\n")
