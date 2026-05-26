import polars as pl

from genemapper_converter.converter import extract_sample


def test_extract_sample_returns_only_matching_rows(
    minimal_genemapper_df: pl.DataFrame,
) -> None:
    result = extract_sample(minimal_genemapper_df, "CASE01_DOADOR.fsa")
    assert all(
        name == "CASE01_DOADOR.fsa" for name in result["Sample File Name"].to_list()
    )


def test_extract_sample_row_count_matches_source(
    minimal_genemapper_df: pl.DataFrame,
) -> None:
    result = extract_sample(minimal_genemapper_df, "CASE01_DOADOR.fsa")
    assert len(result) == 3


def test_extract_sample_preserves_all_columns(
    minimal_genemapper_df: pl.DataFrame,
) -> None:
    result = extract_sample(minimal_genemapper_df, "CASE01_DOADOR.fsa")
    assert set(result.columns) == set(minimal_genemapper_df.columns)


def test_extract_sample_preserves_dye_peak_values(
    minimal_genemapper_df: pl.DataFrame,
) -> None:
    result = extract_sample(minimal_genemapper_df, "CASE01_DOADOR.fsa")
    assert '"B,1"' in result["Dye/Sample Peak"].to_list()


def test_extract_sample_unknown_name_returns_empty(
    minimal_genemapper_df: pl.DataFrame,
) -> None:
    result = extract_sample(minimal_genemapper_df, "NONEXISTENT.fsa")
    assert len(result) == 0


def test_extract_sample_does_not_mutate_input(
    minimal_genemapper_df: pl.DataFrame,
) -> None:
    original_len = len(minimal_genemapper_df)
    extract_sample(minimal_genemapper_df, "CASE01_DOADOR.fsa")
    assert len(minimal_genemapper_df) == original_len


def test_extract_sample_null_rows_are_preserved(
    minimal_genemapper_df: pl.DataFrame,
) -> None:
    result = extract_sample(minimal_genemapper_df, "CASE01_PRE.fsa")
    null_allele_rows = result.filter(pl.col("Allele").is_null())
    assert len(null_allele_rows) == 1
