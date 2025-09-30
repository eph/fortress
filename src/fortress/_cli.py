from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path


def _find_smc_exe() -> Path | None:
    # Prefer the packaged binary installed by CMake
    pkg_dir = Path(__file__).resolve().parent
    cand = pkg_dir / "bin" / "smc_driver"
    if cand.exists():
        return cand
    # Fallback: search on PATH (useful for editable builds)
    path = shutil.which("smc_driver")
    return Path(path) if path else None


def run_smc(argv: list[str] | None = None) -> int:
    if argv is None:
        argv = sys.argv[1:]
    exe = _find_smc_exe()
    if not exe:
        print("smc_driver executable not found. Try rebuilding with BUILD_DRIVERS=ON.", file=sys.stderr)
        return 2
    # Ensure local rpath-friendly lookup for sibling libs
    env = os.environ.copy()
    # No special env needed if we linked statically to dependencies; keep hook for future
    try:
        proc = subprocess.run([str(exe), *argv], check=False)
        return proc.returncode
    except OSError as e:
        print(f"Failed to launch {exe}: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(run_smc())

