from __future__ import annotations

import ctypes as _ct
from pathlib import Path
from typing import Sequence


def _load_lib() -> _ct.CDLL:
    here = Path(__file__).resolve().parent
    cand = here / "libbssm_circle_api.so"
    if not cand.exists():
        # Manylinux naming fallback if CMake used default SONAME
        for name in ("bssm_circle_api.so", "libbssm_circle_api.dylib", "bssm_circle_api.dylib"):
            alt = here / name
            if alt.exists():
                cand = alt
                break
    return _ct.CDLL(str(cand))


_lib = _load_lib()

_lib.bssm_circle_create.restype = _ct.c_int64
_lib.bssm_circle_free.argtypes = (_ct.c_int64,)
_lib.bssm_circle_prior_logpdf.argtypes = (_ct.c_int64, _ct.POINTER(_ct.c_double), _ct.c_int, _ct.POINTER(_ct.c_double))
_lib.bssm_circle_lik_logpdf.argtypes = (_ct.c_int64, _ct.POINTER(_ct.c_double), _ct.c_int, _ct.c_int, _ct.POINTER(_ct.c_double))


class CircleModel:
    def __init__(self) -> None:
        self._h = _lib.bssm_circle_create()
        if self._h <= 0:
            raise RuntimeError("failed to create BSSM circle model")

    def close(self) -> None:
        if getattr(self, "_h", 0):
            _lib.bssm_circle_free(self._h)
            self._h = 0

    def __del__(self) -> None:  # pragma: no cover
        try:
            self.close()
        except Exception:
            pass

    @staticmethod
    def _as_cvec(x: Sequence[float]):
        arr = (_ct.c_double * len(x))(*map(float, x))
        return arr

    def prior_logpdf(self, params: Sequence[float]) -> float:
        buf = self._as_cvec(params)
        out = _ct.c_double()
        _lib.bssm_circle_prior_logpdf(self._h, buf, _ct.c_int(len(params)), _ct.byref(out))
        return float(out.value)

    def lik_logpdf(self, params: Sequence[float], T: int) -> float:
        buf = self._as_cvec(params)
        out = _ct.c_double()
        _lib.bssm_circle_lik_logpdf(self._h, buf, _ct.c_int(len(params)), _ct.c_int(T), _ct.byref(out))
        return float(out.value)

