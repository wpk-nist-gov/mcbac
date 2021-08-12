from .core import CifHelper
from .internaldata import get_database_0, get_database_1, get_database_2

# Version info
try:
    import pkg_resources

    __version__ = pkg_resources.get_distribution("icpmsflow").version
except Exception:
    # Local copy or not installed with setuptools.
    # Disable minimum version checks on downstream libraries.
    __version__ = "999"

__author__ = """William P. Krekelberg"""
__email__ = "wpk@nist.gov"


__all__ = [
    "CifHelper",
    "get_database_0",
    "get_database_1",
    "get_database_2",
]
