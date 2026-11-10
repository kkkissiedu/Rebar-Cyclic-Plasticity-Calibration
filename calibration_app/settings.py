"""
settings.py
===========

Per-user settings and tool discovery, so nothing machine-specific is hardcoded.

Settings live in one JSON file in the user's home folder (``~/.rebar_calibrator.json``). Each value can also be given
by an environment variable, which wins over the file:

    ABAQUS_CMD       full path (or command name) of the ABAQUS launcher
    VCVARS64         full path of Visual Studio's vcvars64.bat (C++ user materials only)
    REBAR_SESSIONS   folder where run results are written

Discovery order for ABAQUS: environment variable, saved setting, ``abaqus`` on PATH, the standard SIMULIA command
folders. For MSVC: environment variable, saved setting, ``vswhere``, then the standard Visual Studio folders.
Only the FE backend and "Transfer to CAE" need ABAQUS; the surrogate pipeline needs neither ABAQUS nor MSVC.
"""

from __future__ import annotations

import glob
import json
import os
import shutil
import subprocess
import sys
from typing import Dict, Optional

_FILE = os.path.join(os.path.expanduser("~"), ".rebar_calibrator.json")
_KEYS = ("abaqus_cmd", "vcvars", "sessions_dir")


def _load() -> Dict[str, str]:
    try:
        with open(_FILE, encoding="utf-8") as f:
            data = json.load(f)
        return {k: data[k] for k in _KEYS if isinstance(data.get(k), str) and data[k]}
    except (OSError, ValueError):
        return {}


def get(key: str) -> Optional[str]:
    """Saved value of a setting (None if unset)."""
    return _load().get(key)


def set_value(key: str, value: Optional[str]) -> None:
    """Save (or clear, with None or "") one setting."""
    if key not in _KEYS:
        raise KeyError(key)
    data = _load()
    if value:
        data[key] = value
    else:
        data.pop(key, None)
    with open(_FILE, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=1)


# ---------------------------------------------------------------- ABAQUS
def _abaqus_candidates() -> list:
    found = []
    on_path = shutil.which("abaqus")
    if on_path:
        found.append(on_path)
    if sys.platform == "win32":
        for root in ("C:\\SIMULIA\\Commands", "D:\\SIMULIA\\Commands"):
            for name in ("abaqus.bat", "abaqus.exe"):
                p = os.path.join(root, name)
                if os.path.exists(p):
                    found.append(p)
    else:
        found += glob.glob("/opt/simulia/*/*/abaqus") + glob.glob("/usr/local/bin/abaqus")
    return found


def abaqus_command() -> Optional[str]:
    """The ABAQUS launcher to use, or None if none can be found."""
    env = os.environ.get("ABAQUS_CMD")
    if env:
        return env
    saved = get("abaqus_cmd")
    if saved:
        return saved
    cands = _abaqus_candidates()
    return cands[0] if cands else None


def resolve_abaqus(cmd: Optional[str] = None) -> str:
    """Turn the default 'abaqus' (or None) into the discovered launcher; keep explicit values as given."""
    if cmd and cmd != "abaqus":
        return cmd
    return abaqus_command() or "abaqus"


# ---------------------------------------------------------------- MSVC
_VS_EDITIONS = ("Community", "Professional", "Enterprise", "BuildTools")


def _vcvars_candidates() -> list:
    out = []
    pf86 = os.environ.get("ProgramFiles(x86)", r"C:\Program Files (x86)")
    pf = os.environ.get("ProgramFiles", r"C:\Program Files")
    vswhere = os.path.join(pf86, "Microsoft Visual Studio", "Installer", "vswhere.exe")
    if os.path.exists(vswhere):
        try:
            r = subprocess.run([vswhere, "-latest", "-products", "*", "-property", "installationPath"],
                               capture_output=True, text=True, timeout=20)
            for line in r.stdout.splitlines():
                if line.strip():
                    out.append(os.path.join(line.strip(), "VC", "Auxiliary", "Build", "vcvars64.bat"))
        except (OSError, subprocess.SubprocessError):
            pass
    for base, year in ((pf, "2022"), (pf86, "2019"), (pf, "2019")):
        for ed in _VS_EDITIONS:
            out.append(os.path.join(base, "Microsoft Visual Studio", year, ed, "VC", "Auxiliary", "Build",
                                    "vcvars64.bat"))
    return out


def vcvars_path() -> Optional[str]:
    """vcvars64.bat for compiling C++ user materials, or None if no MSVC is found."""
    env = os.environ.get("VCVARS64")
    if env and os.path.exists(env):
        return env
    saved = get("vcvars")
    if saved and os.path.exists(saved):
        return saved
    for p in _vcvars_candidates():
        if os.path.exists(p):
            return p
    return None


# ---------------------------------------------------------------- results folder
def sessions_root() -> str:
    """Absolute folder for run results (never depends on the launch directory)."""
    root = os.environ.get("REBAR_SESSIONS") or get("sessions_dir")
    if not root:
        root = os.path.join(os.path.expanduser("~"), "RebarCalibrator", "sessions")
    return os.path.abspath(root)


# ---------------------------------------------------------------- preflight
def preflight() -> Dict[str, str]:
    """One status line per external tool, for the settings dialog and for clear error messages."""
    cmd = abaqus_command()
    vc = vcvars_path()
    return {
        "abaqus": cmd or "NOT FOUND: set the ABAQUS command in Settings (only needed for FE verification).",
        "msvc": vc or "NOT FOUND: install Visual Studio Build Tools with the C++ workload "
                      "(only needed for the UVC and Ohno-Wang finite-element runs).",
        "sessions": sessions_root(),
    }
