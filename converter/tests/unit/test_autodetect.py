import pytest

from genemapper_converter.autodetect import detect_role
from genemapper_converter.models import SampleRole


@pytest.mark.parametrize(
    ("sample_name", "expected_role"),
    [
        ("0028MD01-071020 01D SG DOADOR", SampleRole.DONOR),
        ("SAMPLE_DONOR_001", SampleRole.DONOR),
        ("0028MR-071020 01D SG PRE", SampleRole.RECIPIENT),
        ("CASE01_RECEPTOR", SampleRole.RECIPIENT),
        ("CP_QUI353", SampleRole.SAMPLE),
        ("QUI001_POST", SampleRole.SAMPLE),
        ("UNKNOWN_SAMPLE_XYZ", None),
    ],
)
def test_detect_role_matches_expected(
    sample_name: str, expected_role: SampleRole | None
) -> None:
    assert detect_role(sample_name) == expected_role
