fortress
========

Fortress is a set of modern Fortran components for Bayesian time–series estimation
and the tooling required to expose those components to Python. The project serves
two audiences:

* **Python users** who install the package and call the high-level API/CLI that
  scikit-build-core produces.
* **Fortran developers** who want a reusable library (modules plus static/shared
  libraries) that can be consumed from another CMake build.

This guide covers installation, testing, and integrating Fortress into an external
Fortran project.

---

Installation
------------

Fortress uses [scikit-build-core] and honours CMake cache options. The two most
common workflows are:

### Python / pip

```
# Serial build (default)
pip install fortress

# Enable MPI/OpenMP support
ENABLE_MPI=ON ENABLE_OPENMP=ON pip install fortress
```

If you prefer pip’s configuration interface (>=23.1), you can pass the flags as
`--config-settings`, for example:

```
pip install fortress --config-settings=cmake.define.ENABLE_MPI=ON
```

The same options work with `uv sync`/`uv pip`; the repository already documents
`uv sync --config-settings=cmake.define.ENABLE_MPI=ON` as the canonical editable
workflow.

### From source (developer)

```
uv sync --config-settings=cmake.define.ENABLE_MPI=ON \
        --config-settings=cmake.define.BUILD_FORTRAN_TESTS=ON
uvx --from cmake cmake -S . -B build -DENABLE_MPI=ON
uvx --from cmake cmake --build build
```

If you installed CMake ≥3.21 via another route, swap the `uvx --from cmake` calls
for your local executable.

---

Testing
-------

The project registers its Fortran unit suite as `fortress_unit_tests` inside the
CMake tree. After configuring and building (as above) run:

```
uvx --from cmake ctest --test-dir build
```

The suite defaults to the MPI-enabled executables when `ENABLE_MPI=ON`. When MPI
is disabled the tests that rely on `mpif.h` will be skipped automatically.

---

Using Fortress as a Fortran Dependency
--------------------------------------

Fortress installs both static and shared libraries together with Fortran module
files and a CMake package configuration. After running `cmake --install .` (or
installing the built wheel), an external CMake project can integrate Fortress in
either of two ways:

### `find_package`

```cmake
# CMakeLists.txt in your project
find_package(fortress CONFIG REQUIRED)          # add -Dfortress_ROOT=/path/to/install if needed

add_executable(example main.f90)
target_link_libraries(example PRIVATE fortress::fortress_static)

# Optional: request the shared library instead
# target_link_libraries(example PRIVATE fortress::fortress_shared)
```

Provide `-Dfortress_ROOT=/install/prefix` or extend `CMAKE_PREFIX_PATH` if the
package is not in a default search location.

### `FetchContent`

```cmake
include(FetchContent)

FetchContent_Declare(
  fortress
  GIT_REPOSITORY https://github.com/eph/fortress.git   # or point at your local clone
  GIT_TAG main
)
FetchContent_MakeAvailable(fortress)

add_executable(example main.f90)
target_link_libraries(example PRIVATE fortress::fortress_static)
```

To switch between a local checkout and GitHub just update the `GIT_REPOSITORY`
(or override `FetchContent_SOURCE_DIR_fortress`). Compile-time options propagate
normally: pass `-DENABLE_MPI=ON`/`OFF` when configuring your top-level project and
Fortress will rebuild accordingly.

---

Frequently Used CMake Options
-----------------------------

* `ENABLE_MPI` – build MPI-aware components (`OFF` by default).
* `ENABLE_OPENMP` – enable OpenMP parallelism when MPI is not available (`ON` by
  default, auto-disables if the compiler lacks OpenMP support).
* `BUILD_FORTRAN_TESTS` – compile the FRUIT-based test driver and register it with
  CTest; useful when packaging or running CI from the source tree.

[scikit-build-core]: https://scikit-build-core.readthedocs.io/
