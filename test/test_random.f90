module test_random
  use, intrinsic :: iso_fortran_env, only: wp => real64

  use fruit
  implicit none

contains

  subroutine test_random_init
    use fortress, only: fortress_random

    type(fortress_random) :: rn
    real(wp), allocatable :: x(:,:), y(:,:)
    integer, parameter :: n_small = 5000

    rn = fortress_random()

    allocate(x(n_small,1), y(n_small,1))

    x = rn%norm_rvs(n_small,1,mu=1.0_wp)
    call assert_equals(n_small, size(x,1))
    call assert_equals(1, size(x,2))
    call assert_equals(.true., abs(sum(x)/n_small - 1.0_wp) < 0.15_wp)
    call assert_equals(.true., maxval(abs(x)) < 10.0_wp)

    x = rn%uniform_rvs(n_small,1)
    call assert_equals(.true., minval(x) >= 0.0_wp)
    call assert_equals(.true., maxval(x) <= 1.0_wp)
    call assert_equals(.true., abs(sum(x)/n_small - 0.5_wp) < 0.1_wp)

    y = rn%gamma_rvs(n_small, 1, theta=0.3_wp, k=2.0_wp)
    call assert_equals(.true., minval(y) >= 0.0_wp)
    call assert_equals(.true., abs(sum(y)/n_small - 0.6_wp) < 0.3_wp)

    y = rn%beta_rvs(n_small, 1, a=2.0_wp, b=2.0_wp)
    call assert_equals(.true., minval(y) >= 0.0_wp)
    call assert_equals(.true., maxval(y) <= 1.0_wp)

    y = rn%beta_rvs(n_small, 1, a=2.625_wp, b=2.625_wp)
    call assert_equals(.true., minval(y) >= 0.0_wp)
    call assert_equals(.true., maxval(y) <= 1.0_wp)

    deallocate(x, y)

  end subroutine test_random_init

  subroutine test_mvn_norm
    use fortress, only: fortress_random

    type(fortress_random) :: rn

    real(wp) :: L(2,2), V(2,2), mu(2)
    real(wp), allocatable :: x(:,:)
    integer, parameter :: n_mvn = 5000

    mu = 0.0_wp
  
    rn = fortress_random()

    allocate(x(n_mvn, 2))

    L = 0.0_wp
    L(1,1) = 1.0_wp
    L(2,2) = 1.0_wp

    x = rn%mv_norm_rvs(n_mvn, 2, mu, L)

    call assert_equals(.true., abs(sum(x(:,1))/n_mvn) < 0.15_wp)
    call assert_equals(.true., abs(sum(x(:,2))/n_mvn) < 0.15_wp)

    mu = [2.0_wp, 3.0_wp]
    x = rn%mv_norm_rvs(n_mvn, 2, mu, L)
    call assert_equals(.true., abs(sum(x(:,1))/n_mvn - 2.0_wp) < 0.2_wp)
    call assert_equals(.true., abs(sum(x(:,2))/n_mvn - 3.0_wp) < 0.2_wp)

    L(1,:) = [1.0_wp, 0.0_wp]
    L(2,:) = [0.6_wp, 0.8_wp]

    mu = [0.0_wp, 0.0_wp]
    x = rn%mv_norm_rvs(n_mvn, 2, mu, L)
    call assert_equals(.true., abs(sum(x(:,1)*x(:,2))/n_mvn - 0.6_wp) < 0.1_wp)

    L(1,:) = [1.0_wp, 0.6_wp]
    L(2,:) = [0.6_wp, 1.0_wp]
    mu = [0.0_wp, 0.0_wp]
    x = rn%mv_norm_rvs(n_mvn, 2, mu, L, use_cholesky=.false.)
    call assert_equals(.true., abs(sum(x(:,1)*x(:,2))/n_mvn - 0.6_wp) < 0.1_wp)

    deallocate(x)

  end subroutine test_mvn_norm
  
  subroutine draw1(rng, rvs)
    use fortress, only: fortress_random
    type(fortress_random), intent(in) :: rng
    real(wp) :: rvs(1,1)

    rvs = rng%norm_rvs(1,1)
    
  end subroutine draw1

  subroutine test_passing_rng
    use fortress, only: fortress_random

    type(fortress_random) :: rng1, rng2

    real(wp) :: d1, d2, x(1,1)

    rng1 = fortress_random(seed=12345)
    call draw1(rng1, x)
    d1 = x(1,1)

    rng2 = fortress_random(seed=12345)
    call draw1(rng2, x)
    d2 = x(1,1)

    call assert_equals(d1, d2, 1.0e-8_wp)
  end subroutine test_passing_rng


  
end module test_random
