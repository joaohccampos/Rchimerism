from enum import StrEnum
from typing import Final

__all__: list[str] = ["SampleRole", "ROLE_TO_FILENAME", "AUTODETECT_PATTERNS"]


class SampleRole(StrEnum):
    DONOR = "donor"
    DONOR_1 = "donor1"
    DONOR_2 = "donor2"
    RECIPIENT = "recipient"
    SAMPLE = "sample"


ROLE_TO_FILENAME: Final[dict[SampleRole, str]] = {
    SampleRole.DONOR: "ddata.txt",
    SampleRole.DONOR_1: "d1data.txt",
    SampleRole.DONOR_2: "d2data.txt",
    SampleRole.RECIPIENT: "rdata.txt",
    SampleRole.SAMPLE: "sdata.txt",
}

AUTODETECT_PATTERNS: Final[dict[SampleRole, list[str]]] = {
    SampleRole.DONOR: ["DOADOR", "DONOR", "DON"],
    SampleRole.RECIPIENT: ["PRE", "RECEP", "RECIPIENT", "RECEPTOR"],
    SampleRole.SAMPLE: ["QUI", "CHIM", "POS", "POST"],
}
