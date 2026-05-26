import polars as pl
import pytest


@pytest.fixture()
def minimal_genemapper_df() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "Dye/Sample Peak": [
                '"B,1"',
                '"B,2"',
                '"Y,1"',
                '"B,1"',
                '"B,2"',
                '"Y,1"',
            ],
            "Sample File Name": [
                "CASE01_DOADOR.fsa",
                "CASE01_DOADOR.fsa",
                "CASE01_DOADOR.fsa",
                "CASE01_PRE.fsa",
                "CASE01_PRE.fsa",
                "CASE01_PRE.fsa",
            ],
            "Marker": ["D3S1358", "D3S1358", "TH01", "D3S1358", "D3S1358", "TH01"],
            "Allele": ["14", "15", "7", "15", None, "9"],
            "Size": [120.5, 124.3, 166.0, 124.3, None, 173.9],
            "Height": [3000, 3100, 2500, 2800, None, 2200],
            "Area": [15000, 15500, 12000, 14000, None, 11000],
            "Data Point": [2086, 2123, 2456, 2123, None, 2554],
        }
    )
