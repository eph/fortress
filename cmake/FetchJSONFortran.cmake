include(FetchContent)
set(ENABLE_TESTS OFF CACHE BOOL "Disable json-fortran bundled tests" FORCE)
FetchContent_Declare(jsonfortran
  GIT_REPOSITORY https://github.com/jacobwilliams/json-fortran.git
  GIT_TAG 9.0.5
)
FetchContent_MakeAvailable(jsonfortran)
