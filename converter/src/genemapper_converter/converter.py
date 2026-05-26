import polars as pl

__all__: list[str] = ["extract_sample"]


def extract_sample(df: pl.DataFrame, sample_file_name: str) -> pl.DataFrame:
    return df.filter(pl.col("Sample File Name") == sample_file_name)
