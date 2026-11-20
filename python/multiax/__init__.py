"""multiax - data-oriented C++ core of the Multiax-Invariant method (Bolchoun 2018), Python interface.

The compiled extension is ``multiax._core``. Its only runtime dependency besides the CRT is
``libomp.dll`` (OpenMP), which is installed next to the extension.
"""
import os as _os
import sys as _sys

if _sys.platform == 'win32':
    # libomp.dll lies next to _core*.pyd; make sure the folder is searched (Python >= 3.8 ignores PATH)
    _here = _os.path.dirname(_os.path.abspath(__file__))
    try:
        _os.add_dll_directory(_here)
    except (AttributeError, FileNotFoundError, OSError):
        pass

from ._core import build_info, WoehlerCurve  # noqa: E402

__version__ = build_info()['version']
__all__ = ['build_info', 'WoehlerCurve', '__version__']
