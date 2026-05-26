from pathlib import Path

import polars as pl
import pytest

from genemapper_converter.parser import list_sample_names, parse_genemapper


@pytest.fixture()
def genemapper_tsv(tmp_path: Path) -> Path:
    sep = "\t"
    header = sep.join(
        [
            "Dye/Sample Peak",
            "Sample File Name",
            "Marker",
            "Allele",
            "Size",
            "Height",
            "Area",
            "Data Point",
            "",
        ]
    )
    rows = [
        sep.join(['"B,1"', "CASO_DOADOR.fsa", "", "", "", "295", "4928", "748", ""]),
        sep.join(
            [
                '"B,40"',
                "CASO_DOADOR.fsa",
                "D3S1358",
                "14",
                "120.5",
                "3000",
                "15000",
                "2086",
                "",
            ]
        ),
        sep.join(
            [
                '"B,41"',
                "CASO_DOADOR.fsa",
                "D3S1358",
                "15",
                "124.3",
                "3100",
                "15500",
                "2123",
                "",
            ]
        ),
        sep.join(
            [
                '"B,40"',
                "CASO_PRE.fsa",
                "D3S1358",
                "15",
                "124.3",
                "2800",
                "14000",
                "2123",
                "",
            ]
        ),
    ]
    content = "\n".join([header, *rows]) + "\n"
    tsv_file = tmp_path / "test.txt"
    tsv_file.write_text(content)
    return tsv_file


def test_parse_genemapper_returns_expected_columns(genemapper_tsv: Path) -> None:
    df = parse_genemapper(genemapper_tsv)
    assert "Dye/Sample Peak" in df.columns
    assert "Sample File Name" in df.columns
    assert "Marker" in df.columns
    assert "Allele" in df.columns


def test_parse_genemapper_allele_column_is_string(genemapper_tsv: Path) -> None:
    df = parse_genemapper(genemapper_tsv)
    assert df["Allele"].dtype == pl.String


def test_parse_genemapper_preserves_dye_peak_quotes(genemapper_tsv: Path) -> None:
    df = parse_genemapper(genemapper_tsv)
    first = df["Dye/Sample Peak"][0]
    assert first == '"B,1"'


def test_parse_genemapper_empty_fields_become_null(genemapper_tsv: Path) -> None:
    df = parse_genemapper(genemapper_tsv)
    row = df.filter(pl.col("Dye/Sample Peak") == '"B,1"')
    assert row["Marker"][0] is None
    assert row["Allele"][0] is None


def test_list_sample_names_returns_unique_ordered(
    minimal_genemapper_df: pl.DataFrame,
) -> None:
    names = list_sample_names(minimal_genemapper_df)
    assert names == ["CASE01_DOADOR.fsa", "CASE01_PRE.fsa"]


def test_list_sample_names_no_duplicates(minimal_genemapper_df: pl.DataFrame) -> None:
    names = list_sample_names(minimal_genemapper_df)
    assert len(names) == len(set(names))
