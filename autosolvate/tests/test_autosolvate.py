"""
Unit and regression test for the autosolvate package.
"""

# Import package, test suite, and other packages as needed
import autosolvate
import pytest
import sys

def test_autosolvate_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "autosolvate" in sys.modules
