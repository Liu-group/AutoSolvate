"""
AutoSolvate
Automated workflow for adding explict solvent to molecules
"""

# Add imports here
from .autosolvate import *
from .generatetrajs import *
from .clustergen import *
from .multicomponent import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
