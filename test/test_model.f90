module test_model
  use, intrinsic :: iso_fortran_env, only: wp => real64

  use fruit
  implicit none

contains

  subroutine test_model_initialize

    use test_model_t, only : model, new_model

    type(model) :: m

    m = model()

    call assert_equals('ss', m%name)
    call assert_equals(1, m%nobs)
    call assert_equals(80, m%T)
    call assert_equals(2, m%npara)
    call assert_equals(2, m%ns)
    call assert_equals(1, m%neps)
    call assert_equals(0.2_wp, m%p0(1))
    call assert_equals(0.5_wp, m%p0(2))
    call assert_equals(.false., m%MISSING_DATA)
    
  end subroutine test_model_initialize
  
  subroutine test_model_likelihood

    use test_model_t, only : model
    use filter, only : REALLY_NEG
    type(model) :: m
    integer :: i

    m = model()
    call assert_equals(-85.82166_wp, m%lik(m%p0), 0.0001_wp)


    i = 0
    if (m%inbounds([-1.0_wp, 1.0_wp])) i = 1
    call assert_equals(0, i)

    i = 0
    if ( m%inbounds(m%p0) ) i = 1

    call assert_equals(1, i)

    call assert_equals(REALLY_NEG, m%lik([-1.0_wp, 0.3_wp]), 0.0001_wp)


  end subroutine test_model_likelihood

  subroutine test_model_prior

    use test_model_t, only : model
  
    type(model) :: m
    
    m = model()
    call assert_equals(0.0_wp, m%prior%logpdf(m%p0), 0.0001_wp)
  
  end subroutine test_model_prior


  subroutine test_model_policy

    use test_model_t, only : model

    type(model) :: m
    logical :: converged 
    real(wp) :: states_new(2)

    m = model()

    converged = m%solve(m%p0)

    states_new = m%policy_function([0.0_wp, 0.0_wp], [0.0_wp])

    call assert_equals(0.0_wp, states_new(1))
    call assert_equals(0.0_wp, states_new(2))

    converged = m%solve([0.4_wp, 0.3_wp])
    states_new = m%policy_function([0.2_wp, 0.4_wp], [1.0_wp])
    
    call assert_equals(1.032_wp, states_new(1))
    call assert_equals(0.48_wp, states_new(2))


  end subroutine test_model_policy


  subroutine test_model_policy_1

    use test_model_t, only : model

    type(model) :: m
    logical :: converged 
    real(wp) :: states_new(2)

    m = model()

    converged = m%solve(m%p0)

    states_new = m%policy_function([0.0_wp, 0.0_wp], [0.0_wp])

    converged = m%solve([0.4_wp, 0.3_wp])
    states_new = m%policy_function([0.2_wp, 0.4_wp], [1.0_wp])

  end subroutine test_model_policy_1



end module
