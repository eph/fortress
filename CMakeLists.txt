cmake_minimum_required(VERSION 3.21)
project(fortress VERSION 0.1.5 LANGUAGES C Fortran)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(CTest)

# Build all libraries with PIC so they can link into the Python module.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(ENABLE_MPI "Build with MPI" OFF)
option(ENABLE_OPENMP "Enable OpenMP parallelism when MPI is disabled" ON)
option(USE_SYSTEM_JSONFORTRAN "Use system-installed JSON-Fortran" OFF)
option(USE_SYSTEM_FLAP "Use system-installed FLAP" OFF)
option(BUILD_FORTRAN_TESTS "Build and register Fortran unit tests" ${BUILD_TESTING})

find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development.Module)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -c "import numpy,sys;sys.stdout.write(numpy.get_include())"
  OUTPUT_VARIABLE NUMPY_INCLUDE
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Determine the correct Python extension module suffix (e.g.,
# ".cpython-312-x86_64-linux-gnu.so"). Fallback to ".so" if unavailable.
execute_process(
  COMMAND "${Python_EXECUTABLE}" -c "import sysconfig,sys;sys.stdout.write(sysconfig.get_config_var('EXT_SUFFIX') or '')"
  OUTPUT_VARIABLE PY_EXT_SUFFIX
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT PY_EXT_SUFFIX)
  set(PY_EXT_SUFFIX ".so")
endif()

set(FYPP_ARGS "-m time -DGFORTRAN" CACHE STRING "Arguments passed to the fypp preprocessor")
find_program(FYPP_EXECUTABLE fypp REQUIRED)
separate_arguments(FYPP_ARGS_LIST NATIVE_COMMAND "${FYPP_ARGS}")

if(USE_SYSTEM_JSONFORTRAN)
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(JSONFORTRAN REQUIRED json-fortran)
  add_library(jsonfortran INTERFACE IMPORTED)
  set_target_properties(jsonfortran PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${JSONFORTRAN_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${JSONFORTRAN_LIBRARIES}"
  )
  set(JSONFORTRAN_TARGET jsonfortran)
  set(JSONFORTRAN_PROVIDER_TARGET jsonfortran)
else()
  include(cmake/FetchJSONFortran.cmake)
  # JSON-Fortran upstream uses a few different target names depending on
  # version/options; choose the first one that exists.
  # Prefer static library to avoid runtime loader issues with rpaths.
  set(_jsonfortran_candidates jsonfortran-static jsonfortran_static jsonfortran jsonfortran::jsonfortran)
  foreach(_t IN LISTS _jsonfortran_candidates)
    if(TARGET ${_t})
      set(JSONFORTRAN_TARGET ${_t})
      break()
    endif()
  endforeach()
  if(JSONFORTRAN_TARGET)
    set(JSONFORTRAN_PROVIDER_TARGET ${JSONFORTRAN_TARGET})
  endif()
  if(NOT JSONFORTRAN_TARGET)
    # Fall back to a sensible default name
    set(JSONFORTRAN_TARGET jsonfortran)
  endif()
  if(NOT JSONFORTRAN_PROVIDER_TARGET)
    set(JSONFORTRAN_PROVIDER_TARGET jsonfortran)
  endif()
endif()

if(USE_SYSTEM_FLAP)
  # System FLAP: allow user/env-provided paths to MODs and library
  # Hints: set FLAP_MOD_DIR and FLAP_LIBRARY (or rely on default search paths)
  find_path(FLAP_MOD_DIR
    NAMES flap.mod
    HINTS $ENV{FLAP_MOD_DIR}
    PATH_SUFFIXES flap modules include
    PATHS /usr/local/include /usr/include
  )
  find_library(FLAP_LIBRARY
    NAMES flap
    HINTS $ENV{FLAP_LIBRARY}
    PATHS /usr/local/lib /usr/lib
  )
  if(NOT FLAP_LIBRARY)
    message(FATAL_ERROR "Could not find system FLAP library; set FLAP_LIBRARY or adjust paths.")
  endif()
  add_library(flap UNKNOWN IMPORTED)
  set_target_properties(flap PROPERTIES
    IMPORTED_LOCATION "${FLAP_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${FLAP_MOD_DIR}"
  )
  set(FLAP_TARGET flap)
else()
  include(cmake/FetchFLAP.cmake)
  include(FetchContent)
  # Try to detect the actual FLAP target name exported by the subproject.
  # Common possibilities across versions/configs.
  set(_flap_candidates
    FLAP::flap
    FLAP::FLAP
    flap
    flap_static
    flap-shared
    flaplib
  )
  foreach(_t IN LISTS _flap_candidates)
    if(TARGET ${_t})
      set(FLAP_TARGET ${_t})
      break()
    endif()
  endforeach()
  # Do not attempt to alias namespaced targets; FLAP upstream may export an
  # ALIAS already (e.g., FLAP::FLAP). We'll link to whatever exists.
  if(NOT FLAP_TARGET)
    # As a fallback, assume a plain 'flap' will exist during generation.
    set(FLAP_TARGET flap)
  endif()
  # We will add include dirs after 'fortress_mod' is defined, using target
  # properties and FetchContent paths.
endif()

if(ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS Fortran)
endif()

if(ENABLE_OPENMP)
  find_package(OpenMP COMPONENTS Fortran)
  if(OpenMP_Fortran_FOUND)
    list(APPEND FYPP_ARGS_LIST "-DENABLE_OPENMP")
  else()
    message(WARNING "OpenMP requested but not found; continuing without OpenMP support")
    set(ENABLE_OPENMP OFF)
  endif()
endif()

set(BLA_VENDOR OpenBLAS)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

set(SRC_FORTRAN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/fortran")
file(GLOB_RECURSE fortress_src CONFIGURE_DEPENDS "${SRC_FORTRAN_DIR}/*.f90")

set(processed_sources)
if(ENABLE_MPI)
  list(APPEND FYPP_ARGS_LIST "-DENABLE_MPI")
endif()

foreach(src ${fortress_src})
  file(RELATIVE_PATH rel_src "${SRC_FORTRAN_DIR}" "${src}")
  set(out "${CMAKE_CURRENT_BINARY_DIR}/fypp/${rel_src}")
  get_filename_component(out_dir "${out}" DIRECTORY)
  add_custom_command(
    OUTPUT "${out}"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${out_dir}"
    COMMAND ${FYPP_EXECUTABLE} ${FYPP_ARGS_LIST} "${src}" "${out}"
    DEPENDS "${src}"
    COMMENT "Preprocessing ${rel_src}"
    VERBATIM
  )
  list(APPEND processed_sources "${out}")
endforeach()

# Generate all fypp outputs up-front so CMake can scan module deps reliably.
add_custom_target(fypp_gen ALL DEPENDS ${processed_sources})

# Build the Python module directly from the processed sources.
add_library(fortress_mod MODULE ${processed_sources})
set_target_properties(fortress_mod PROPERTIES
  PREFIX ""
  SUFFIX "${PY_EXT_SUFFIX}"
  OUTPUT_NAME "_fortress"
  POSITION_INDEPENDENT_CODE ON
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod_fortress"
)

# Make sure the target’s own module output dir is also on the search path.
function(fortress_configure_fortran_target target_name)
  target_include_directories(${target_name} PRIVATE
    "$<TARGET_PROPERTY:${target_name},Fortran_MODULE_DIRECTORY>"
    "${NUMPY_INCLUDE}"
  )

  target_compile_options(${target_name} PRIVATE
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-ffree-line-length-none>
  )

  if(TARGET ${FLAP_TARGET})
    get_target_property(_flap_iface_inc ${FLAP_TARGET} INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(_flap_mod_dir   ${FLAP_TARGET} Fortran_MODULE_DIRECTORY)
    if(_flap_iface_inc)
      target_include_directories(${target_name} PRIVATE "${_flap_iface_inc}")
    endif()
    if(_flap_mod_dir)
      target_include_directories(${target_name} PRIVATE "${_flap_mod_dir}")
    endif()
  endif()

  if(JSONFORTRAN_PROVIDER_TARGET AND TARGET ${JSONFORTRAN_PROVIDER_TARGET})
    get_target_property(_json_iface_inc ${JSONFORTRAN_PROVIDER_TARGET} INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(_json_mod_dir   ${JSONFORTRAN_PROVIDER_TARGET} Fortran_MODULE_DIRECTORY)
    if(_json_iface_inc)
      target_include_directories(${target_name} PRIVATE "${_json_iface_inc}")
    endif()
    if(_json_mod_dir)
      target_include_directories(${target_name} PRIVATE "${_json_mod_dir}")
    endif()
  endif()
endfunction()

fortress_configure_fortran_target(fortress_mod)

target_link_libraries(fortress_mod PRIVATE
  ${JSONFORTRAN_TARGET}
  ${FLAP_TARGET}
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
)

# Static library variant used by native executables/tests.
add_library(fortress_static STATIC ${processed_sources})
set_target_properties(fortress_static PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod_fortress_static"
)

fortress_configure_fortran_target(fortress_static)

target_link_libraries(fortress_static PRIVATE
  ${JSONFORTRAN_TARGET}
  ${FLAP_TARGET}
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
)

# Ensure dependencies are built before compiling sources that USE their modules
if(TARGET ${FLAP_TARGET})
  add_dependencies(fortress_mod ${FLAP_TARGET})
  add_dependencies(fortress_static ${FLAP_TARGET})
else()
  message(STATUS "FLAP target '${FLAP_TARGET}' not found at configure time; relying on fallback include paths.")
endif()
if(JSONFORTRAN_PROVIDER_TARGET AND TARGET ${JSONFORTRAN_PROVIDER_TARGET})
  add_dependencies(fortress_mod ${JSONFORTRAN_PROVIDER_TARGET})
  add_dependencies(fortress_static ${JSONFORTRAN_PROVIDER_TARGET})
endif()

# Ensure generated sources exist before compiling objects (ordering).
add_dependencies(fortress_mod fypp_gen)
add_dependencies(fortress_static fypp_gen)

if(TARGET ${FLAP_TARGET})
  get_target_property(FLAP_MOD_DIR ${FLAP_TARGET} Fortran_MODULE_DIRECTORY)
  if(NOT FLAP_MOD_DIR)
    # Try FetchContent's recorded binary dir
    include(FetchContent)
    FetchContent_GetProperties(flap)
    if(DEFINED flap_BINARY_DIR)
      if(EXISTS "${flap_BINARY_DIR}/include")
        set(FLAP_MOD_DIR "${flap_BINARY_DIR}/include")
      elseif(EXISTS "${flap_BINARY_DIR}/modules")
        set(FLAP_MOD_DIR "${flap_BINARY_DIR}/modules")
      else()
        set(FLAP_MOD_DIR "${flap_BINARY_DIR}")
      endif()
    endif()
  endif()
  if(FLAP_MOD_DIR)
    target_include_directories(fortress_mod PRIVATE "${FLAP_MOD_DIR}")
    install(DIRECTORY "${FLAP_MOD_DIR}/" DESTINATION fortress/modules FILES_MATCHING PATTERN "*.mod")
  endif()
endif()

# Note: avoid adding non-existent fallback include dirs to keep compilers quiet.

if(ENABLE_MPI)
  target_link_libraries(fortress_mod PRIVATE MPI::MPI_Fortran)
  target_link_libraries(fortress_static PUBLIC MPI::MPI_Fortran)
endif()

if(ENABLE_OPENMP AND OpenMP_Fortran_FOUND)
  target_link_libraries(fortress_mod PRIVATE OpenMP::OpenMP_Fortran)
  target_link_libraries(fortress_static PUBLIC OpenMP::OpenMP_Fortran)
endif()

if(APPLE)
  set_target_properties(fortress_mod PROPERTIES INSTALL_RPATH "@loader_path")
else()
  set_target_properties(fortress_mod PROPERTIES INSTALL_RPATH "$ORIGIN")
endif()

install(TARGETS fortress_mod
  LIBRARY DESTINATION fortress
  RUNTIME DESTINATION fortress
  ARCHIVE DESTINATION fortress
)



# For SKBUILD (pip install), also install CMake config and static library to the package
if(SKBUILD)
  # Install the static library to the package
  install(TARGETS fortress_static
    ARCHIVE DESTINATION fortress/lib
  )

  # Install Fortran module files to the fortress package's include directory
  install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod_fortress_static/"
    DESTINATION fortress/include
    FILES_MATCHING PATTERN "*.mod"
  )

  # For editable installs, also copy mod files after build using a custom command
  # This ensures mod files are available even when install() commands don't run
  add_custom_command(TARGET fortress_static POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory "${SKBUILD_PLATLIB_DIR}/fortress/include"
    COMMAND ${CMAKE_COMMAND} -E copy_directory
            "${CMAKE_CURRENT_BINARY_DIR}/mod_fortress_static"
            "${SKBUILD_PLATLIB_DIR}/fortress/include"
    COMMENT "Copying Fortran module files to install directory for editable installs"
  )

  # Configure and install CMake config files for the package
  # For SKBUILD installs, paths are relative to the fortress package root
  set(FORTRESS_INSTALL_INCLUDEDIR "include")
  set(FORTRESS_INSTALL_MODDIR "include")

  configure_package_config_file(
    cmake/fortressConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/fortressConfig_package.cmake
    INSTALL_DESTINATION fortress/lib/cmake/fortress
    PATH_VARS FORTRESS_INSTALL_INCLUDEDIR FORTRESS_INSTALL_MODDIR
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
  )

  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/fortressConfigVersion_package.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
  )

  install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/fortressConfig_package.cmake
    DESTINATION fortress/lib/cmake/fortress
    RENAME fortressConfig.cmake
  )

  install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/fortressConfigVersion_package.cmake
    DESTINATION fortress/lib/cmake/fortress
    RENAME fortressConfigVersion.cmake
  )
endif()

# --- Optional native drivers (CLI executables) ---
option(BUILD_DRIVERS "Build native driver executables (e.g., smc_driver)" ON)

if(BUILD_DRIVERS AND ENABLE_MPI)
  # Build smc_driver by compiling fortress sources again (reliable) plus driver sources.
  add_executable(smc_driver
    ${processed_sources}
    ${CMAKE_CURRENT_SOURCE_DIR}/smc_driver.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/test/test_model_circle_t.f90
  )

  add_dependencies(smc_driver fypp_gen)

  target_include_directories(smc_driver PRIVATE
    "${NUMPY_INCLUDE}"
    "$<TARGET_PROPERTY:smc_driver,Fortran_MODULE_DIRECTORY>"
  )
  set_target_properties(smc_driver PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod_smc")

  if(TARGET ${FLAP_TARGET})
    get_target_property(_flap_iface_inc2 ${FLAP_TARGET} INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(_flap_mod_dir2   ${FLAP_TARGET} Fortran_MODULE_DIRECTORY)
    if(_flap_iface_inc2)
      target_include_directories(smc_driver PRIVATE "${_flap_iface_inc2}")
    endif()
    if(_flap_mod_dir2)
      target_include_directories(smc_driver PRIVATE "${_flap_mod_dir2}")
    endif()
  endif()

  if(JSONFORTRAN_PROVIDER_TARGET AND TARGET ${JSONFORTRAN_PROVIDER_TARGET})
    get_target_property(_json_iface_inc2 ${JSONFORTRAN_PROVIDER_TARGET} INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(_json_mod_dir2   ${JSONFORTRAN_PROVIDER_TARGET} Fortran_MODULE_DIRECTORY)
    if(_json_iface_inc2)
      target_include_directories(smc_driver PRIVATE "${_json_iface_inc2}")
    endif()
    if(_json_mod_dir2)
      target_include_directories(smc_driver PRIVATE "${_json_mod_dir2}")
    endif()
  endif()

  target_link_libraries(smc_driver PRIVATE
    ${JSONFORTRAN_TARGET}
    ${FLAP_TARGET}
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
  )

  # Allow long free-form lines.
  target_compile_options(smc_driver PRIVATE
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-ffree-line-length-none>
  )

  if(ENABLE_MPI)
    target_link_libraries(smc_driver PRIVATE MPI::MPI_Fortran)
  endif()
  if(ENABLE_OPENMP AND OpenMP_Fortran_FOUND)
    target_link_libraries(smc_driver PRIVATE OpenMP::OpenMP_Fortran)
  endif()

  if(APPLE)
    set_target_properties(smc_driver PROPERTIES INSTALL_RPATH "@loader_path")
  else()
    set_target_properties(smc_driver PROPERTIES INSTALL_RPATH "$ORIGIN")
  endif()

  install(TARGETS smc_driver RUNTIME DESTINATION fortress/bin)
elseif(BUILD_DRIVERS AND NOT ENABLE_MPI)
  message(STATUS "BUILD_DRIVERS requested but ENABLE_MPI=OFF; skipping smc_driver (mpif.h). Set -DENABLE_MPI=ON to build it.")
endif()

# --- Optional BSSM C API (ctypes-accessible shared lib) ---
option(BUILD_BSSM_C_API "Build ctypes-friendly BSSM C API (circle demo)" ON)
if(BUILD_BSSM_C_API)
  add_library(bssm_circle_api SHARED
    ${processed_sources}
    ${CMAKE_CURRENT_SOURCE_DIR}/test/test_model_circle_t.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/bssm/bssm_circle_api.f90
  )
  add_dependencies(bssm_circle_api fypp_gen)
  set_target_properties(bssm_circle_api PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod_bssm")
  target_include_directories(bssm_circle_api PRIVATE "$<TARGET_PROPERTY:bssm_circle_api,Fortran_MODULE_DIRECTORY>")
  target_compile_options(bssm_circle_api PRIVATE
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-ffree-line-length-none>
  )
  target_link_libraries(bssm_circle_api PRIVATE
    ${JSONFORTRAN_TARGET}
    ${FLAP_TARGET}
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
  )
  if(ENABLE_MPI)
    target_link_libraries(bssm_circle_api PRIVATE MPI::MPI_Fortran)
  endif()
  if(ENABLE_OPENMP AND OpenMP_Fortran_FOUND)
    target_link_libraries(bssm_circle_api PRIVATE OpenMP::OpenMP_Fortran)
  endif()
  if(APPLE)
    set_target_properties(bssm_circle_api PROPERTIES INSTALL_RPATH "@loader_path")
  else()
    set_target_properties(bssm_circle_api PROPERTIES INSTALL_RPATH "$ORIGIN")
  endif()
  install(TARGETS bssm_circle_api
    LIBRARY DESTINATION fortress
    RUNTIME DESTINATION fortress
    ARCHIVE DESTINATION fortress)
endif()

# --- Fortran unit tests -----------------------------------------------------
if(BUILD_TESTING AND BUILD_FORTRAN_TESTS)
  if(NOT ENABLE_MPI)
    message(STATUS "ENABLE_MPI=OFF -> skipping Fortran unit tests (require MPI for mpif.h)")
  else()
    set(fortress_unit_test_sources
      ${CMAKE_CURRENT_SOURCE_DIR}/src/testing/fruit_stub.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/test/test_json.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/test/test_prior.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/test/test_model_t.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/test/test_model.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/test/test_linalg.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/test/test_util.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/test/test_gensys.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/test/test_random.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/test/test_smc.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/test/test_particles.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/test_driver.f90
    )

    add_executable(fortress_tests ${fortress_unit_test_sources})
    set_target_properties(fortress_tests PROPERTIES
      Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod_tests"
    )

    fortress_configure_fortran_target(fortress_tests)
    target_include_directories(fortress_tests PRIVATE
      "$<TARGET_PROPERTY:fortress_static,Fortran_MODULE_DIRECTORY>"
    )

    target_link_libraries(fortress_tests PRIVATE fortress_static)
    add_dependencies(fortress_tests fortress_static)

    if(ENABLE_OPENMP AND OpenMP_Fortran_FOUND)
      target_link_libraries(fortress_tests PRIVATE OpenMP::OpenMP_Fortran)
    endif()

    if(ENABLE_MPI AND MPI_Fortran_FOUND)
      add_test(NAME fortress_unit_tests
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:fortress_tests> ${MPIEXEC_POSTFLAGS})
    else()
      add_test(NAME fortress_unit_tests COMMAND $<TARGET_FILE:fortress_tests>)
    endif()
    set_tests_properties(fortress_unit_tests PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
  endif()
endif()

# --- CMake package export ---------------------------------------------------
set(FORTRESS_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}/fortress")

if(NOT SKBUILD)
  # Install fortress_static and its module files
  install(TARGETS fortress_static
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${FORTRESS_INSTALL_INCLUDEDIR}
  )

  install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod_fortress_static/"
    DESTINATION ${FORTRESS_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.mod"
  )

  set(FORTRESS_INSTALL_MODDIR "${FORTRESS_INSTALL_INCLUDEDIR}")

  # Create a simple config file that just sets paths (no target exports for now)
  configure_package_config_file(
    cmake/fortressConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/fortressConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fortress
    PATH_VARS FORTRESS_INSTALL_INCLUDEDIR FORTRESS_INSTALL_MODDIR
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
  )

  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/fortressConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
  )

  install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/fortressConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/fortressConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fortress
  )
endif()
