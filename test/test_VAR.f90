module test_var
  use, intrinsic :: iso_fortran_env, only: wp => real64

  use fruit

  use fortress_prior_t, only: fortress_abstract_prior


  implicit none

  type, public :: prior_set
     class(fortress_abstract_prior), pointer :: pr
  end type prior_set


contains

  subroutine test_svar_init

    use fortress_VAR_t, only: SimsZhaSVARPrior

    real(wp) :: parasim(100000,3), lpdf

    type(SimsZhaSVARPrior) :: swz
    character(len=144) :: mufile, varfile

    mufile = 'test/var/swz_mu.txt'
    varfile ='test/var/swz_var.txt'
    swz = SimsZhaSVARPrior(3, mufile, varfile)

    call assert_equals(3, swz%npara)
    call assert_equals(0.0_wp, swz%mu(1))
    call assert_equals(1.0_wp, swz%var(3,3))

    parasim = swz%rvs(100000)

    call assert_equals(0.0_wp, sum(parasim)/100000, 0.005_wp)
    
    varfile = 'test/var/swz_var_corr.txt'
    swz = SimsZhaSVARPrior(3, mufile, varfile)
    lpdf = swz%logpdf([0.2_wp, 0.1_wp,-0.6_wp])
    call assert_equals(-3.225898_wp, lpdf, 0.00005_wp)

  end subroutine test_svar_init

  subroutine test_var_init

    use tvt_t, only: MinnesotaPrior

    real(wp) :: parasim(100000,3), lpdf
    type(MinnesotaPrior) :: minnpr

    minnpr = MinnesotaPrior(1,1,1,1.0_wp,1.0_wp,1.0_wp,1.0_wp,1.0_wp,1.0_wp,1.0_wp,[0.0_wp],[0.0_wp])
    
    call assert_equals(3, minnpr%npara)
    call assert_equals(1, minnpr%p)
    call assert_equals(1, minnpr%constant)
    
  end subroutine test_var_init

  subroutine test_VAR_mult_rvs

    use tvt_t, only: MinnesotaPrior
    use fortress_random_t, only: fortress_random
    
    type(MinnesotaPrior) :: minnpr

    type(fortress_random) :: rng
    class(prior_set), allocatable, dimension(:) :: coeff_prior

    real(wp) :: rvs1(3,1000), rvs2(3,1000)

    rng = fortress_random()
    minnpr = MinnesotaPrior(1,1,1,1.0_wp,1.0_wp,1.0_wp,1.0_wp,1.0_wp,1.0_wp,1.0_wp,[0.0_wp],[0.0_wp])
    call assert_equals(0.0_wp, 1.0_wp, 0.00001_wp)
    allocate(coeff_prior(2))

    allocate(coeff_prior(1)%pr,source=minnpr)
    minnpr = MinnesotaPrior(1,1,1,1.0_wp,1.0_wp,1.0_wp,1.0_wp,1.0_wp,1.0_wp,1.0_wp,[0.0_wp],[0.0_wp])

    allocate(coeff_prior(1)%pr,source=minnpr)

    rvs1 = coeff_prior(1)%pr%rvs(1000, rng=rng)
    rvs2 = coeff_prior(2)%pr%rvs(1000, rng=rng)


    deallocate(coeff_prior)
  end subroutine test_VAR_mult_rvs

  subroutine test_svar

    use test_VAR_t, only: model

    real(wp) :: parasim(100000,54), lpdf

    real(wp) :: p0(54)
    type(model) :: swz
    character(len=144) :: mufile, varfile

    ! swz = model()

    ! call assert_equals(54, swz%npara)
    ! call assert_equals(183, swz%T)
    ! call assert_equals(3, swz%nobs)
    ! call assert_equals(5, swz%p)


    ! p0 = [125.730_wp, 84.578_wp, 54.303_wp,-54.529_wp, -3.103_wp, 20.981_wp,-54.890_wp,  8.429_wp,-95.042_wp,142.540_wp, 13.495_wp, 32.189_wp, 13.586_wp, -1.255_wp, 25.950_wp,-39.483_wp, 14.809_wp, 29.036_wp, 17.699_wp, -6.200_wp,  5.008_wp, -0.107_wp,164.214_wp, 99.716_wp, 38.117_wp, 32.166_wp, 19.109_wp,-22.681_wp, 21.814_wp, -7.998_wp,-35.420_wp, 34.484_wp,-25.540_wp, -7.258_wp, 44.071_wp,  6.020_wp, -8.435_wp, -0.025_wp,-247.662_wp, 82.674_wp, 13.246_wp, 10.592_wp,-89.664_wp, -2.935_wp, 17.063_wp, -3.500_wp, -6.658_wp, 18.117_wp,  4.108_wp, 33.954_wp, 17.469_wp,-12.105_wp, -8.910_wp, -0.125_wp]

    ! call assert_equals(-249.2239_wp, swz%prior%logpdf(p0), 0.001_wp)
    ! call assert_equals(-3988.8155_wp, swz%lik(p0), 0.02_wp)

    ! call assert_equals(0.0_wp, swz%prior%mu(1))
    ! parasim = swz%rvs(100000)
    ! call assert_equals(0.0_wp, sum(parasim)/100000, 0.005_wp)
    ! varfile = 'test/var/swz_var_corr.txt'
    ! swz = SimsZhaSVARPrior(3, mufile, varfile)
    ! lpdf = swz%logpdf([0.2_wp,0.1_wp,-0.6_wp])
    ! call assert_equals(-3.225898_wp, lpdf, 0.00005_wp)

  end subroutine test_svar



end module test_var
