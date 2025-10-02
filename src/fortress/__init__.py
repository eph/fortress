try:
    from ._fortress import *  # type: ignore # noqa: F401,F403
except Exception:  # Extension not built; CLI-only usage still works
    pass

from .fortress import SMCDriver, load_estimates, make_model_file, make_smc

__all__ = [name for name in dir() if not name.startswith("_")]
