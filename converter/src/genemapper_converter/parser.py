from pathlib import Path
from typing import Final

import polars as pl

__all__: list[str] = ["parse_genemapper", "list_sample_names"]

_EXPECTED_COLUMNS: Final[list[str]] = [
    "Dye/Sample Peak",
    "Sample File Name",
    "Marker",
    "Allele",
    "Size",
    "Height",
    "Area",
    "Data Point",
]


def parse_genemapper(file_path: Path) -> pl.DataFrame:
    raw = pl.read_csv(
        file_path,
        separator="\t",
        truncate_ragged_lines=True,
        null_values=["", " "],
        quote_char=None,
        schema_overrides={"Allele": pl.String},
    )
    available = [col for col in _EXPECTED_COLUMNS if col in raw.columns]
    return raw.select(available)


def list_sample_names(df: pl.DataFrame) -> list[str]:
    return (
        df.select("Sample File Name").unique(maintain_order=True).to_series().to_list()
    )
