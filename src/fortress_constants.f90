 module fortress_constants
  use, intrinsic :: iso_fortran_env, only: wp => real64

  real(wp), parameter :: REALLY_NEGATIVE = -100000000.0_wp
  real(wp), parameter :: BAD_LOG_LIKELIHOOD = REALLY_NEGATIVE 
end module fortress_constants
