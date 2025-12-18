"""Environment detection helpers for AutoSolvate interactive wizard."""
from __future__ import annotations

import os
import shutil
from typing import Dict, List, Optional

from ..Common import ANTECHAMBER, PACKMOL, TLEAP, PARMCHK, OBABEL


def detect_amberhome() -> Optional[str]:
    amberhome = os.getenv("AMBERHOME")
    if amberhome:
        return os.path.expanduser(os.path.expandvars(amberhome))
    return None


def detect_amber_tools(amberhome: Optional[str] = None) -> Dict[str, Optional[str]]:
    amberhome = amberhome or detect_amberhome()
    tools = {"antechamber": ANTECHAMBER, "parmchk2": PARMCHK, "tleap": TLEAP, "packmol": PACKMOL, "obabel": OBABEL}
    paths: Dict[str, Optional[str]] = {}
    for label, exe in tools.items():
        path = shutil.which(exe)
        if path:
            paths[label] = path
            continue
        if amberhome:
            candidate = os.path.join(amberhome, "bin", exe)
            if os.path.exists(candidate):
                paths[label] = candidate
                continue
        paths[label] = None
    return paths


def available_qm_software() -> Dict[str, Optional[str]]:
    """Return detected QM executables by shorthand name."""
    candidates = {
        "orca": "orca",
        "gau": "gau",
        "g09": "g09",
        "g03": "g03",
        "gms": "rungms",
    }
    found: Dict[str, Optional[str]] = {}
    for short, exe in candidates.items():
        path = shutil.which(exe)
        found[short] = path
    return found


def cpu_count() -> int:
    try:
        return os.cpu_count() or 1
    except Exception:
        return 1
