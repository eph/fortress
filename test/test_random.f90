module test_random
  use, intrinsic :: iso_fortran_env, only: wp => real64

  use fruit
  implicit none

contains

  subroutine test_random_init
    use fortress, only: fortress_random

    type(fortress_random) :: rn
    real(wp) :: x(10000,1), y(100000,1)

    rn = fortress_random()

    x = rn%norm_rvs(10000,1)
    call assert_equals(0.0_wp, sum(x)/10000, 0.005_wp)
    x = rn%norm_rvs(10000,1,mu=1.0_wp)
    call assert_equals(1.0_wp, sum(x)/10000, 0.005_wp)

    x = rn%uniform_rvs(10000,1)
    call assert_equals(0.5_wp, sum(x)/10000, 0.005_wp)

    y = rn%gamma_rvs(100000, 1, theta=0.3_wp, k=2.0_wp)
    call assert_equals(0.6_wp, sum(y)/100000, 0.005_wp)
    call assert_equals(0.3_wp**2*2.0_wp, sum( (y - 0.6_wp)**2 )/100000, 0.03_wp)

    y = rn%beta_rvs(100000, 1, a=2.0_wp, b=2.0_wp)
    call assert_equals(0.5_wp, sum(y)/100000, 0.005_wp)
    call assert_equals( 4.0_wp/(16.0_wp*5.0_wp), sum( (y - 0.5_wp)**2 )/100000, 0.002_wp)
    
    y = rn%beta_rvs(100000, 1, a=2.625_wp, b=2.625_wp)
    call assert_equals(0.5_wp, sum(y)/100000, 0.005_wp)
    call assert_equals( 0.04_wp, sum( (y - 0.5_wp)**2 )/100000, 0.002_wp)


  end subroutine test_random_init

end module test_random
