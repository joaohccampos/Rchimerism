from pathlib import Path

import polars as pl
import pytest

from genemapper_converter.writer import write_rchimerism_file


@pytest.fixture()
def sample_df() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "Dye/Sample Peak": ['"B,1"', '"B,2"'],
            "Sample File Name": ["CASE01_DOADOR.fsa", "CASE01_DOADOR.fsa"],
            "Marker": ["D3S1358", "D3S1358"],
            "Allele": ["14", "15"],
            "Size": [120.5, 124.3],
            "Height": [3000, 3100],
            "Area": [15000, 15500],
            "Data Point": [2086, 2123],
        }
    )


def test_write_rchimerism_file_creates_file(
    tmp_path: Path, sample_df: pl.DataFrame
) -> None:
    output = tmp_path / "ddata.txt"
    write_rchimerism_file(sample_df, output)
    assert output.exists()


def test_write_rchimerism_file_creates_parent_directories(
    tmp_path: Path, sample_df: pl.DataFrame
) -> None:
    output = tmp_path / "nested" / "dir" / "ddata.txt"
    write_rchimerism_file(sample_df, output)
    assert output.exists()


def test_write_rchimerism_file_header_matches_expected(
    tmp_path: Path, sample_df: pl.DataFrame
) -> None:
    output = tmp_path / "ddata.txt"
    write_rchimerism_file(sample_df, output)
    lines = output.read_text().splitlines()
    header_fields = lines[0].split("\t")
    assert header_fields[0] == "Dye/Sample Peak"
    assert header_fields[2] == "Marker"
    assert header_fields[3] == "Allele"


def test_write_rchimerism_file_data_row_count_matches(
    tmp_path: Path, sample_df: pl.DataFrame
) -> None:
    output = tmp_path / "ddata.txt"
    write_rchimerism_file(sample_df, output)
    lines = [line for line in output.read_text().splitlines() if line.strip()]
    assert len(lines) == 3  # 1 header + 2 data rows


def test_write_rchimerism_file_tab_separated(
    tmp_path: Path, sample_df: pl.DataFrame
) -> None:
    output = tmp_path / "ddata.txt"
    write_rchimerism_file(sample_df, output)
    first_data_line = output.read_text().splitlines()[1]
    assert "\t" in first_data_line
    assert '"B,1"' in first_data_line


def test_write_rchimerism_file_null_fields_written_as_empty(
    tmp_path: Path,
) -> None:
    df = pl.DataFrame(
        {
            "Dye/Sample Peak": ['"B,1"'],
            "Sample File Name": ["CASE01.fsa"],
            "Marker": [None],
            "Allele": [None],
            "Size": [None],
            "Height": [295],
            "Area": [4928],
            "Data Point": [748],
        }
    )
    output = tmp_path / "out.txt"
    write_rchimerism_file(df, output)
    data_line = output.read_text().splitlines()[1]
    fields = data_line.split("\t")
    assert fields[2] == ""
    assert fields[3] == ""
