 module fortress_constants
  use, intrinsic :: iso_fortran_env, only: wp => real64

  real(wp), parameter :: REALLY_NEGATIVE = -100000000.0_wp
  real(wp), parameter :: BAD_LOG_LIKELIHOOD = REALLY_NEGATIVE

  ! make a constant for 2 * pi
  real(wp), parameter :: PI = 3.14159265358979323846264338327950288419716939937510_wp
  real(wp), parameter :: TWO_PI = 2.0_wp * pi

  ! make a constant for the inverse of 2 * pi
  real(wp), parameter :: PI_INV = 1.0_wp / pi
  real(wp), parameter :: TWO_PI_INV = 1.0_wp / (2.0_wp * pi)
  
end module fortress_constants
