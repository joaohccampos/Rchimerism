from pathlib import Path
from typing import Annotated

import typer

from .autodetect import detect_role
from .converter import extract_sample
from .models import ROLE_TO_FILENAME, SampleRole
from .parser import list_sample_names, parse_genemapper
from .writer import write_rchimerism_file

__all__: list[str] = []

app = typer.Typer(help="Split GeneMapper exports into rchimerism input files.")


@app.command("list")
def list_samples(
    input_file: Annotated[Path, typer.Argument(help="GeneMapper TSV export file")],
) -> None:
    """List all sample file names found in the GeneMapper export."""
    df = parse_genemapper(input_file)
    names = list_sample_names(df)
    typer.echo(f"Found {len(names)} sample(s):")
    for name in names:
        role = detect_role(name)
        role_hint = (
            f"  → auto-detected as: {role.value}" if role else "  → no role detected"
        )
        typer.echo(f"  {name!r}{role_hint}")


@app.command("convert")
def convert(
    input_file: Annotated[Path, typer.Argument(help="GeneMapper TSV export file")],
    output_dir: Annotated[
        Path, typer.Option("--output-dir", "-o", help="Output directory")
    ] = Path("."),
    donor: Annotated[
        str | None,
        typer.Option("--donor", help="Sample File Name for single donor"),
    ] = None,
    donor1: Annotated[
        str | None,
        typer.Option("--donor1", help="Sample File Name for donor 1 (double donor)"),
    ] = None,
    donor2: Annotated[
        str | None,
        typer.Option("--donor2", help="Sample File Name for donor 2 (double donor)"),
    ] = None,
    recipient: Annotated[
        str | None,
        typer.Option("--recipient", "-r", help="Sample File Name for recipient"),
    ] = None,
    sample: Annotated[
        str | None,
        typer.Option("--sample", "-s", help="Sample File Name for chimera sample"),
    ] = None,
    auto_detect: Annotated[
        bool,
        typer.Option("--auto-detect", help="Auto-detect roles from sample file names"),
    ] = False,
) -> None:
    """Split GeneMapper export into rchimerism input files by role."""
    df = parse_genemapper(input_file)
    found_names = list_sample_names(df)

    role_map: dict[SampleRole, str] = {}

    if auto_detect:
        for name in found_names:
            detected = detect_role(name)
            if detected is not None:
                if detected in role_map:
                    typer.echo(
                        f"Warning: multiple samples for role '{detected.value}'. "
                        f"Keeping '{role_map[detected]}', ignoring '{name}'.",
                        err=True,
                    )
                else:
                    role_map[detected] = name
                    typer.echo(f"Auto-detected: {name!r} → {detected.value}")

    explicit_assignments: list[tuple[SampleRole, str]] = [
        (role, value)
        for role, value in [
            (SampleRole.DONOR, donor),
            (SampleRole.DONOR_1, donor1),
            (SampleRole.DONOR_2, donor2),
            (SampleRole.RECIPIENT, recipient),
            (SampleRole.SAMPLE, sample),
        ]
        if value is not None
    ]
    for role, name in explicit_assignments:
        role_map[role] = name

    if not role_map:
        typer.echo(
            "No role assignments found. "
            "Use --auto-detect or provide explicit sample file names.",
            err=True,
        )
        raise typer.Exit(code=1)

    for role, name in role_map.items():
        if name not in found_names:
            typer.echo(f"Error: sample '{name}' not found in {input_file}.", err=True)
            typer.echo(f"Available samples: {found_names}", err=True)
            raise typer.Exit(code=1)

        extracted = extract_sample(df, name)
        output_path = output_dir / ROLE_TO_FILENAME[role]
        write_rchimerism_file(extracted, output_path)
        typer.echo(
            f"Written: {output_path} ({len(extracted)} peaks, role: {role.value})"
        )


if __name__ == "__main__":
    app()
