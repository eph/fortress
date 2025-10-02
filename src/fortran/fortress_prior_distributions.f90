!*******************************************************************************
!>
! author: Ed Herbst
!
!# Description
!  Distribution functions for custom priors - can be inlined in generated code.
!  Provides logpdf and rvs for common distributions.

module fortress_prior_distributions

  use, intrinsic :: iso_fortran_env, only: wp => real64
  use fortress_random_t, only : fortress_random
  use logbeta, only : betaln, gamln

  implicit none

  private
  public :: beta_logpdf, gamma_logpdf, normal_logpdf, invgamma_logpdf, uniform_logpdf
  public :: beta_rvs, gamma_rvs, normal_rvs, invgamma_rvs, uniform_rvs

  real(wp), parameter :: M_PI = 3.14159265358979d0
  real(wp), parameter :: LOG_2PI = 1.83787706640935d0

contains

  !*******************************************************************************
  ! Beta distribution: Beta(mean*nu, (1-mean)*nu) where nu = mean*(1-mean)/std^2 - 1
  !*******************************************************************************
  function beta_logpdf(x, mean, std) result(lpdf)
    real(wp), intent(in) :: x, mean, std
    real(wp) :: lpdf
    real(wp) :: alpha, beta, nu

    ! Convert mean/std to alpha/beta parameters
    nu = mean * (1.0_wp - mean) / (std * std) - 1.0_wp
    alpha = mean * nu
    beta = (1.0_wp - mean) * nu

    if (x <= 0.0_wp .or. x >= 1.0_wp) then
       lpdf = -1.0d10
    else
       lpdf = (alpha - 1.0_wp) * log(x) + (beta - 1.0_wp) * log(1.0_wp - x) - betaln(alpha, beta)
    end if
  end function beta_logpdf

  function beta_rvs(mean, std, rng) result(draw)
    real(wp), intent(in) :: mean, std
    type(fortress_random), intent(inout) :: rng
    real(wp) :: draw
    real(wp) :: alpha, beta, nu
    real(wp) :: temp_array(1,1)

    nu = mean * (1.0_wp - mean) / (std * std) - 1.0_wp
    alpha = mean * nu
    beta = (1.0_wp - mean) * nu

    temp_array = rng%beta_rvs(1, 1, alpha, beta)
    draw = temp_array(1, 1)
  end function beta_rvs

  !*******************************************************************************
  ! Gamma distribution: Gamma(shape, scale) where shape = (mean/std)^2, scale = std^2/mean
  !*******************************************************************************
  function gamma_logpdf(x, mean, std) result(lpdf)
    real(wp), intent(in) :: x, mean, std
    real(wp) :: lpdf
    real(wp) :: shape, scale

    ! Convert mean/std to shape/scale
    shape = (mean / std) ** 2
    scale = (std * std) / mean

    if (x <= 0.0_wp) then
       lpdf = -1.0d10
    else
       lpdf = -gamln(shape) - shape * log(scale) + (shape - 1.0_wp) * log(x) - x / scale
    end if
  end function gamma_logpdf

  function gamma_rvs(mean, std, rng) result(draw)
    real(wp), intent(in) :: mean, std
    type(fortress_random), intent(inout) :: rng
    real(wp) :: draw
    real(wp) :: shape, scale
    real(wp) :: temp_array(1,1)

    shape = (mean / std) ** 2
    scale = (std * std) / mean

    temp_array = rng%gamma_rvs(1, 1, scale, shape)
    draw = temp_array(1, 1)
  end function gamma_rvs

  !*******************************************************************************
  ! Normal distribution
  !*******************************************************************************
  function normal_logpdf(x, mean, std) result(lpdf)
    real(wp), intent(in) :: x, mean, std
    real(wp) :: lpdf

    lpdf = -0.5_wp * LOG_2PI - log(std) - 0.5_wp * ((x - mean) / std) ** 2
  end function normal_logpdf

  function normal_rvs(mean, std, rng) result(draw)
    real(wp), intent(in) :: mean, std
    type(fortress_random), intent(inout) :: rng
    real(wp) :: draw
    real(wp) :: temp_array(1,1)

    temp_array = rng%norm_rvs(1, 1, mean, std)
    draw = temp_array(1, 1)
  end function normal_rvs

  !*******************************************************************************
  ! Inverse Gamma distribution (parameterized by shape and scale)
  ! For IG(s, nu): mean = nu/(s-1), var = nu^2/((s-1)^2(s-2))
  !*******************************************************************************
  function invgamma_logpdf(x, shape, scale) result(lpdf)
    real(wp), intent(in) :: x, shape, scale
    real(wp) :: lpdf

    if (x <= 0.0_wp) then
       lpdf = -1.0d10
    else
       lpdf = shape * log(scale) - gamln(shape) - (shape + 1.0_wp) * log(x) - scale / x
    end if
  end function invgamma_logpdf

  function invgamma_rvs(shape, scale, rng) result(draw)
    real(wp), intent(in) :: shape, scale
    type(fortress_random), intent(inout) :: rng
    real(wp) :: draw
    real(wp) :: temp_array(1,1)

    ! Draw from Gamma and take reciprocal
    temp_array = rng%gamma_rvs(1, 1, 1.0_wp / scale, shape)
    draw = 1.0_wp / temp_array(1, 1)
  end function invgamma_rvs

  !*******************************************************************************
  ! Uniform distribution
  !*******************************************************************************
  function uniform_logpdf(x, lower, upper) result(lpdf)
    real(wp), intent(in) :: x, lower, upper
    real(wp) :: lpdf

    if (x < lower .or. x > upper) then
       lpdf = -1.0d10
    else
       lpdf = -log(upper - lower)
    end if
  end function uniform_logpdf

  function uniform_rvs(lower, upper, rng) result(draw)
    real(wp), intent(in) :: lower, upper
    type(fortress_random), intent(inout) :: rng
    real(wp) :: draw
    real(wp) :: temp_array(1,1)

    temp_array = rng%uniform_rvs(1, 1, lower, upper)
    draw = temp_array(1, 1)
  end function uniform_rvs

end module fortress_prior_distributions
