from .models import AUTODETECT_PATTERNS, SampleRole

__all__: list[str] = ["detect_role"]


def detect_role(sample_name: str) -> SampleRole | None:
    upper = sample_name.upper()
    for role, patterns in AUTODETECT_PATTERNS.items():
        if any(pattern in upper for pattern in patterns):
            return role
    return None
