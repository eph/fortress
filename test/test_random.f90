module test_random
  use, intrinsic :: iso_fortran_env, only: wp => real64

  use fruit
  implicit none

contains

  subroutine test_random_init
    use fortress, only: fortress_random

    type(fortress_random) :: rn
    real(wp) :: x(10000,1)

    rn = fortress_random()

    x = rn%norm_rvs(10000,1)
    call assert_equals(0.0_wp, sum(x)/10000, 0.005_wp)
    x = rn%norm_rvs(10000,1,mu=1.0_wp)
    call assert_equals(0.0_wp, sum(x)/10000, 1.005_wp)

  end subroutine test_random_init

end module test_random
