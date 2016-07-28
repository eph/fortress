module test_smc
  use, intrinsic :: iso_fortran_env, only: wp => real64

  use fruit
  implicit none

contains

  subroutine test_smc_initialize

    use fortress, only: fortress_smc
    use test_model_t, only : model 

    type(fortress_smc) :: test_smc
    type(model) :: test_model


    test_model = model()
    test_smc = fortress_smc(test_model)

    call assert_equals('ss', test_smc%model%name)
    call assert_equals(1, test_smc%model%nobs)
    call assert_equals(80, test_smc%model%T)
    call assert_equals(2, test_smc%model%npara)

    !call test_smc%estimate

  end subroutine test_smc_initialize
  
end module test_smc
