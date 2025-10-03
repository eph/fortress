from __future__ import annotations

import json
import os
import shlex
import shutil
import subprocess
import textwrap
from pathlib import Path
from typing import Dict, Mapping, MutableMapping, Optional, Sequence, Union
import site

import numpy as np
import pandas as p
import tqdm


_BASE_ENV: MutableMapping[str, str] = os.environ.copy()
_BASE_ENV.setdefault("OPENBLAS_NUM_THREADS", "1")


def load_estimates(
    file_pattern: Union[str, Path],
    *,
    resample: bool = True,
    paranames: Optional[Sequence[str]] = None,
    posterior: str = "final",
) -> Optional[p.DataFrame]:
    """Load posterior draws from fortress output JSON files."""

    pattern = str(file_pattern)
    output_files = sorted(Path().glob(pattern))
    files = [path.resolve() for path in output_files if path.exists()]
    if not files:
        print("No files found")
        return None

    results = []
    for output in files:
        payload = json.loads(output.read_text())
        posterior_keys = sorted(k for k in payload if k.startswith("posterior"))
        if not posterior_keys:
            continue

        key = posterior_keys[-1] if posterior == "final" else posterior
        frame = p.DataFrame(payload[key])

        if resample and "weights" in frame:
            weights = frame["weights"].to_numpy()
            indices = np.random.choice(
                frame.index.size, size=frame.index.size, p=weights
            )
            frame = frame.iloc[indices].reset_index(drop=True)

        if paranames is not None:
            var_columns = [c for c in frame.columns if c.startswith("var")]
            rename_map = {column: name for column, name in zip(var_columns, paranames)}
            frame = frame.rename(columns=rename_map)

        frame["logmdd"] = np.array(payload.get("Z_estimates", [])).sum()
        results.append(frame)

    if not results:
        return None
    if len(results) == 1:
        return results[0]
    return p.concat(results, axis=0, keys=[str(f) for f in files])


class SMCDriver:
    """Wrapper around a compiled SMC driver executable."""

    def __init__(self, executable: Union[str, Path]):
        self.executable = Path(executable).resolve()
        if not self.executable.exists():
            raise FileNotFoundError(f"SMC executable not found: {self.executable}")

        try:
            res = subprocess.run(
                [str(self.executable), "--help"],
                env=_BASE_ENV,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
                text=True,
            )
            help_text = res.stdout or res.stderr
        except OSError:
            help_text = ""
        self.help = help_text

    def run(
        self,
        *,
        nproc: int = 1,
        mpi_command: Optional[Union[str, Sequence[str]]] = None,
        progress: bool = True,
        env: Optional[Mapping[str, str]] = None,
        output_file: Union[str, Path] = "output.json",
        **cli_args,
    ) -> Dict:
        """Execute the driver and return the parsed JSON output."""

        runtime_env = _BASE_ENV.copy()
        if env:
            runtime_env.update(env)

        if mpi_command is None:
            mpi_parts = ["mpirun", "-n", str(nproc)] if nproc > 1 else []
        elif isinstance(mpi_command, str):
            mpi_parts = shlex.split(mpi_command)
        else:
            mpi_parts = list(mpi_command)

        arg_list: list[str] = []
        for key, value in cli_args.items():
            flag = f"--{key.replace('_', '-')}"
            if isinstance(value, bool):
                if value:
                    arg_list.append(flag)
                continue
            arg_list += [flag, str(value)]

        cmd = mpi_parts + [str(self.executable)] + arg_list

        proc = subprocess.Popen(
            cmd,
            env=runtime_env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )

        pbar = tqdm.tqdm(total=1.0, disable=not progress)
        try:
            assert proc.stdout is not None
            for line in proc.stdout:
                line = line.strip()
                if line.lower().startswith("iteration"):
                    parts = line.split()
                    if len(parts) >= 6:
                        try:
                            current = float(parts[1])
                            total = float(parts[5])
                            if total:
                                pbar.n = current / total
                                pbar.refresh()
                        except ValueError:
                            continue
        finally:
            pbar.close()

        proc.wait()
        if proc.returncode != 0:
            raise subprocess.CalledProcessError(proc.returncode, cmd)

        payload_path = Path(output_file)
        if not payload_path.exists():
            raise FileNotFoundError(f"Expected output file not found: {payload_path}")
        return json.loads(payload_path.read_text())


def make_model_file(
    lik_body: str,
    *,
    npara: int,
    T: int,
    other_functions: str = "",
    other_includes: str = "",
) -> str:
    """Create a minimal Fortran model module around the provided likelihood."""

    template = textwrap.dedent(
        """
        module model_t
        use, intrinsic :: iso_fortran_env, only: wp => real64

        use fortress_bayesian_model_t, only: fortress_abstract_bayesian_model
        use fortress_prior_t, only: model_prior => prior
        {other_includes}

        implicit none

        type, public, extends(fortress_abstract_bayesian_model) :: model
        contains
        procedure :: lik
        end type model

        interface model
        module procedure new_model
        end interface model

        contains

        type(model) function new_model() result(self)
        character(len=144) :: name, datafile, priorfile
        integer :: nobs, T_local, npara_local

        name = 'generated-model'
        datafile = 'data.txt'
        priorfile = 'prior.txt'
        nobs = 1
        T_local = {T}
        npara_local = {npara}

        call self%construct_model(name, datafile, npara_local, nobs, T_local)
        allocate(self%prior, source=model_prior(priorfile))
        end function new_model

        function lik(self, para, T) result(l)
        class(model), intent(inout) :: self
        real(wp), intent(in) :: para(self%npara)
        integer, intent(in), optional :: T
        real(wp) :: l
        {lik_body}
        end function lik

        {other_functions}

        end module model_t
        """
    ).strip("\n")

    indented_body = textwrap.indent(lik_body.strip(), " " * 10)
    indented_functions = (
        textwrap.indent(other_functions.strip(), " " * 8) if other_functions else ""
    )
    includes = textwrap.indent(other_includes.strip(), " " * 2) if other_includes else ""

    return template.format(
        lik_body=indented_body,
        npara=npara,
        T=T,
        other_functions=indented_functions,
        other_includes=includes,
    )


DRIVER_SOURCE = textwrap.dedent(
    """
    program smc_driver
    use iso_fortran_env, only: wp => real64
    use fortress, only: fortress_smc
    use model_t, only: model
    implicit none
    include 'mpif.h'

    type(fortress_smc) :: smc
    type(model) :: smc_model
    integer :: mpierror, rank, nproc

    call mpi_init(mpierror)
    call mpi_comm_size(MPI_COMM_WORLD, nproc, mpierror)
    call mpi_comm_rank(MPI_COMM_WORLD, rank, mpierror)

    smc_model = model()
    smc = fortress_smc(smc_model, nproc)
    call smc%estimate(rank)

    call mpi_finalize(mpierror)
    end program smc_driver
    """
).strip()


CMAKE_TEMPLATE = textwrap.dedent(
    """
    cmake_minimum_required(VERSION 3.21)
    project(fortress_smc_driver LANGUAGES Fortran)

    find_package(fortress CONFIG REQUIRED)
    find_package(MPI COMPONENTS Fortran)
    find_package(OpenMP COMPONENTS Fortran)

    add_executable(smc_driver
    smc_driver.f90
    model_t.f90
    @FORTRESS_EXTRA_SOURCES@
    )

    target_link_libraries(smc_driver PRIVATE fortress::fortress_static)

    if(MPI_Fortran_FOUND)
    target_link_libraries(smc_driver PRIVATE MPI::MPI_Fortran)
    endif()

    if(OpenMP_Fortran_FOUND)
    target_link_libraries(smc_driver PRIVATE OpenMP::OpenMP_Fortran)
    endif()

    target_compile_options(smc_driver PRIVATE
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-ffree-line-length-none>
    )

    set_target_properties(smc_driver PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod"
    )
    """
).strip()


def _resolve_cmake_executable() -> str:
    cmake = os.environ.get("CMAKE")
    if cmake and shutil.which(cmake):
        return cmake

    path = shutil.which("cmake")
    if path:
        return path

    try:
        import cmake  # type: ignore

        return str(Path(cmake.CMAKE_BIN_DIR) / "cmake")  # type: ignore[attr-defined]
    except Exception as exc:  # pragma: no cover
        raise RuntimeError(
            "Unable to locate the `cmake` executable. Please install CMake >= 3.21 "
            "or set the CMAKE environment variable."
        ) from exc





def _discover_cmake_package_dir(
    explicit_dir: Optional[Union[str, Path]] = None
) -> Path:
    """
    Discover the fortress CMake package directory.
    """
    if explicit_dir:
        return Path(explicit_dir)

    # Search in standard site-packages locations
    for site_packages_dir in site.getsitepackages():
        candidate = Path(site_packages_dir) / "fortress" / "lib" / "cmake" / "fortress"
        if (candidate / "fortressConfig.cmake").exists():
            return candidate

    # Fallback for development environments
    package_root = Path(__file__).resolve().parent
    dev_candidate = package_root / "lib" / "cmake" / "fortress"
    if (dev_candidate / "fortressConfig.cmake").exists():
        return dev_candidate

    raise RuntimeError(
        "Could not locate fortressConfig.cmake. Set FORTRESS_CMAKE_DIR or install the package."
    )


def _write_auxiliary_file(
    dest: Path, contents: Union[str, Path, np.ndarray, p.DataFrame]
) -> None:
    if isinstance(contents, (np.ndarray, p.DataFrame)):
        np.savetxt(dest, np.asarray(contents))
    else:
        candidate = Path(str(contents))
        if candidate.exists():
            shutil.copyfile(candidate, dest)
        else:
            dest.write_text(str(contents))


def make_smc(
    model_source: Union[str, Path],
    *,
    output_directory: Union[str, Path] = "_fortress_tmp",
    other_files: Optional[
        Mapping[str, Union[str, Path, np.ndarray, p.DataFrame]]
    ] = None,
    extra_sources: Optional[Mapping[str, Union[str, Path]]] = None,
    cmake_args: Optional[Sequence[str]] = None,
    cmake_generator: Optional[str] = None,
    build_type: str = "Release",
    fortress_cmake_dir: Optional[Union[str, Path]] = None,
    env: Optional[Mapping[str, str]] = None,
    check: bool = True,
) -> SMCDriver:
    """Compile an SMC driver for a generated model and return an executable wrapper."""

    output_dir = Path(output_directory).resolve()
    src_dir = output_dir / "src"
    build_dir = output_dir / "build"
    src_dir.mkdir(parents=True, exist_ok=True)
    build_dir.mkdir(parents=True, exist_ok=True)


    model_text = str(model_source)
    model_text = model_text.format(output_directory=str(output_dir))

    other_files = dict(other_files or {})
    for name, contents in list(other_files.items()):
        target = output_dir / Path(name).name
        _write_auxiliary_file(target, contents)
        model_text = model_text.replace(Path(name).name, str(target))

    (src_dir / "model_t.f90").write_text(model_text)
    (src_dir / "smc_driver.f90").write_text(DRIVER_SOURCE)

    extra_sources = dict(extra_sources or {})
    for filename, contents in extra_sources.items():
        _write_auxiliary_file(src_dir / filename, contents)

    extra_list = "\n".join(f"      {name}" for name in extra_sources.keys())
    cmake_text = CMAKE_TEMPLATE.replace("@FORTRESS_EXTRA_SOURCES@", extra_list)
    (src_dir / "CMakeLists.txt").write_text(cmake_text)

    cmake_bin = _resolve_cmake_executable()
    cmake_dir = _discover_cmake_package_dir(fortress_cmake_dir)

    configure_cmd = [cmake_bin, "-S", str(src_dir), "-B", str(build_dir)]
    configure_cmd.append(f"-DCMAKE_BUILD_TYPE={build_type}")
    configure_cmd.append(f"-Dfortress_DIR={cmake_dir}")
    if cmake_generator:
        configure_cmd.extend(["-G", cmake_generator])
    if cmake_args:
        configure_cmd.extend(cmake_args)

    build_env = _BASE_ENV.copy()
    if env:
        build_env.update(env)

    subprocess.run(configure_cmd, check=check, env=build_env)

    build_cmd = [cmake_bin, "--build", str(build_dir), "--target", "smc_driver"]
    if cmake_generator and not cmake_generator.lower().startswith("ninja"):
        build_cmd.extend(["--config", build_type])
    subprocess.run(build_cmd, check=check, env=build_env)

    candidates = [build_dir / "smc_driver", build_dir / "smc_driver.exe"]
    candidates.extend(build_dir.glob("**/smc_driver*"))
    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return SMCDriver(candidate)

    raise FileNotFoundError("Unable to locate the compiled smc_driver executable")


__all__ = [
    "SMCDriver",
    "load_estimates",
    "make_model_file",
    "make_smc",
]
