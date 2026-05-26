from pathlib import Path

import pytest
from typer.testing import CliRunner

from genemapper_converter.__main__ import app

runner = CliRunner()


@pytest.fixture()
def genemapper_file(tmp_path: Path) -> Path:
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
        sep.join(
            [
                '"B,40"',
                "CASO_QUI.fsa",
                "D3S1358",
                "14",
                "120.5",
                "2900",
                "14500",
                "2086",
                "",
            ]
        ),
        sep.join(
            [
                '"B,41"',
                "CASO_QUI.fsa",
                "D3S1358",
                "15",
                "124.3",
                "2950",
                "14800",
                "2123",
                "",
            ]
        ),
    ]
    content = "\n".join([header, *rows]) + "\n"
    tsv = tmp_path / "input.txt"
    tsv.write_text(content)
    return tsv


def test_list_command_shows_samples(genemapper_file: Path) -> None:
    result = runner.invoke(app, ["list", str(genemapper_file)])
    assert result.exit_code == 0
    assert "CASO_DOADOR.fsa" in result.output
    assert "CASO_PRE.fsa" in result.output


def test_list_command_shows_auto_detected_roles(genemapper_file: Path) -> None:
    result = runner.invoke(app, ["list", str(genemapper_file)])
    assert "donor" in result.output
    assert "recipient" in result.output
    assert "sample" in result.output


def test_convert_with_explicit_args_creates_files(
    genemapper_file: Path, tmp_path: Path
) -> None:
    output_dir = tmp_path / "output"
    result = runner.invoke(
        app,
        [
            "convert",
            str(genemapper_file),
            "--output-dir",
            str(output_dir),
            "--donor",
            "CASO_DOADOR.fsa",
            "--recipient",
            "CASO_PRE.fsa",
            "--sample",
            "CASO_QUI.fsa",
        ],
    )
    assert result.exit_code == 0
    assert (output_dir / "ddata.txt").exists()
    assert (output_dir / "rdata.txt").exists()
    assert (output_dir / "sdata.txt").exists()


def test_convert_with_auto_detect_creates_files(
    genemapper_file: Path, tmp_path: Path
) -> None:
    output_dir = tmp_path / "output"
    result = runner.invoke(
        app,
        [
            "convert",
            str(genemapper_file),
            "--output-dir",
            str(output_dir),
            "--auto-detect",
        ],
    )
    assert result.exit_code == 0
    assert (output_dir / "ddata.txt").exists()


def test_convert_fails_with_unknown_sample(
    genemapper_file: Path, tmp_path: Path
) -> None:
    result = runner.invoke(
        app,
        [
            "convert",
            str(genemapper_file),
            "--output-dir",
            str(tmp_path),
            "--donor",
            "NONEXISTENT.fsa",
        ],
    )
    assert result.exit_code != 0


def test_convert_fails_without_any_role_assignment(
    genemapper_file: Path, tmp_path: Path
) -> None:
    result = runner.invoke(
        app,
        ["convert", str(genemapper_file), "--output-dir", str(tmp_path)],
    )
    assert result.exit_code != 0


def test_convert_output_file_has_correct_header(
    genemapper_file: Path, tmp_path: Path
) -> None:
    output_dir = tmp_path / "output"
    runner.invoke(
        app,
        [
            "convert",
            str(genemapper_file),
            "--output-dir",
            str(output_dir),
            "--donor",
            "CASO_DOADOR.fsa",
        ],
    )
    content = (output_dir / "ddata.txt").read_text()
    lines = content.splitlines()
    assert lines[0].startswith("Dye/Sample Peak\t")


def test_convert_output_preserves_dye_peak_quotes(
    genemapper_file: Path, tmp_path: Path
) -> None:
    output_dir = tmp_path / "output"
    runner.invoke(
        app,
        [
            "convert",
            str(genemapper_file),
            "--output-dir",
            str(output_dir),
            "--donor",
            "CASO_DOADOR.fsa",
        ],
    )
    content = (output_dir / "ddata.txt").read_text()
    assert '"B,' in content
