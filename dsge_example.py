from __future__ import annotations

import argparse
import os
import shlex
from pathlib import Path
from typing import Mapping, Optional, Sequence, Union


import numpy as np

from dsge import read_yaml
from dsge.translate import smc as translate_smc
from fortress import SMCDriver, make_smc


def _normalise_path(value: Union[str, Path]) -> Path:
    return Path(value).expanduser().resolve()


def _parse_passthrough(pairs: Sequence[str]) -> Mapping[str, str]:
    """
    Convert KEY=VALUE command-line pieces into a mapping that we can forward to make_smc/run.
    """
    parsed: dict[str, str] = {}
    for item in pairs:
        if "=" not in item:
            raise ValueError(f"Expected KEY=VALUE form for extra argument, got {item!r}")
        key, value = item.split("=", 1)
        parsed[key] = value
    return parsed


def build_smc_driver(
    yaml_path: Union[str, Path],
    *,
    output_directory: Union[str, Path],
    fortress_dir: Optional[Union[str, Path]] = None,
    cmake_args: Optional[Sequence[str]] = None,
    cmake_generator: Optional[str] = None,
    build_type: str = "Release",
) -> SMCDriver:

    yaml_file = _normalise_path(yaml_path)
    out_dir = _normalise_path(output_directory)
    src_dir = out_dir / "src"
    out_dir.mkdir(parents=True, exist_ok=True)


    model =read_yaml(str(yaml_file))
    model_code = translate_smc(model)

    compiled = model.compile_model()

    driver = make_smc(
        model_code,
        output_directory=out_dir,
        other_files=None,
        extra_sources=None,
        cmake_args=cmake_args,
        cmake_generator=cmake_generator,
        build_type=build_type,
        fortress_cmake_dir=fortress_dir,
    )



    print(f"SMC driver built at {driver.executable}")

    return driver





def main(argv: Optional[Sequence[str]] = None) -> None:

    parser = argparse.ArgumentParser(
           description="Compile (and optionally run) a fortress SMC driver from a DSGE YAML specification."
    )

    parser.add_argument("yaml", help="Path to the DSGE YAML file.")

    parser.add_argument(

        "--output",

        default="_fortress_tmp",

        help="Directory for generated sources/build artifacts (default: _fortress_tmp).",

    )

    parser.add_argument(

        "--fortress-dir",

        help="Override the fortress CMake package location (defaults to the installed package).",

    )

    parser.add_argument(

        "--cmake-arg",

        dest="cmake_args",

        action="append",

        default=[],

        help="Extra -D style arguments to forward to the CMake configure step.",

    )

    parser.add_argument(

        "--cmake-generator",

        help="Explicit CMake generator (e.g. Ninja, Unix Makefiles).",

    )

    parser.add_argument(

        "--build-type",

        default="Release",

        help="Build type passed to CMake (default: Release).",

    )

    parser.add_argument(

        "--run",

        action="store_true",

        help="Run the driver immediately after building it.",

    )

    parser.add_argument(

        "--nproc",

        type=int,

        default=1,

        help="Number of MPI processes when running (default: 1).",

    )

    parser.add_argument(

        "--mpi-command",

        help="Override the MPI launcher (defaults to mpirun -n <nproc>).",

    )

    parser.add_argument(

        "--output-file",

        default="output.json",

        help="SMC output filename to read after running (default: output.json).",

    )

    parser.add_argument(

        "extra",

        nargs="*",

        help="Optional KEY=VALUE pairs forwarded to the SMC driver CLI.",

    )



    args = parser.parse_args(argv)



    driver = build_smc_driver(

        args.yaml,

        output_directory=args.output,

        fortress_dir=args.fortress_dir,

        cmake_args=args.cmake_args,

        cmake_generator=args.cmake_generator,

        build_type=args.build_type,

    )



    if args.run:

        extra_cli = _parse_passthrough(args.extra)

        result = driver.run(

            nproc=args.nproc,

            mpi_command=args.mpi_command,

            output_file=args.output_file,

            **extra_cli,

        )

        print("SMC run complete.")

        print(json.dumps(result, indent=2))





if __name__ == "__main__":

    main()

