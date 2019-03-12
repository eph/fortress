module test_model_circle_t
  use, intrinsic :: iso_fortran_env, only: wp => real64

  use fortress_prior_t, only: fortress_abstract_prior, model_prior => prior 
  use fortress_bayesian_model_t, only: fortress_abstract_bayesian_model, test

  implicit none

  
  type, public, extends(fortress_abstract_bayesian_model) :: model
     real(wp) :: PSI = 5000.0_wp

    contains

      procedure :: lik => model_lik 
      procedure :: dlik => model_dlik
   end type model

  interface model
     module procedure new_model
  end interface model
  
  
contains

  type(model) function new_model() result(self)

    character(len=144) :: name, datafile, priorfile
    integer :: nobs, T, ns, npara, neps
  
    ! name = 'ss'
    ! datafile = '/home/eherbst/Dropbox/code/fortress/test/test_data.txt'
    ! priorfile = '/home/eherbst/Dropbox/code/fortress/test/test_prior_model.txt'

    ! nobs = 1
    ! T = 80

    self%npara = 2
    self%T = 80
    self%nobs = 1
    allocate(self%prior, source=model_prior('/home/eherbst/Dropbox/code/fortress/test/test_prior_circle.txt'))
    !call self%construct_model(name, datafile, priorfile, npara, nobs, T, ns, neps)

    self%p0 = [1.26300312_wp, 1.95483334_wp]

  end function new_model

  subroutine print_test()
    print*,test
  end subroutine print_test

  function model_lik(self, para, T) result(l)

    class(model), intent(inout) :: self
    real(wp), intent(in) :: para(self%npara)
    integer, intent(in), optional :: T
    real(wp) :: l
    integer :: error

    l = -0.5_wp*self%PSI*( (para(1) - 1.0_wp)**2 + (para(2) - 1.0_wp)**2 - 1)**2
    
  end function model_lik


  function model_dlik(self, para, T) result(dl)

    class(model), intent(inout) :: self
    real(wp), intent(in) :: para(self%npara)
    integer, intent(in), optional :: T
    real(wp) :: dl(self%npara)
    integer :: error


    dl(1) = -2.0_wp*self%PSI*(para(1) - 1.0_wp)*((para(1) - 1.0_wp)**2 + (para(2) - 1.0_wp)**2 - 1.0_wp)
    dl(2) = -2.0_wp*self%PSI*(para(2) - 1.0_wp)*((para(1) - 1.0_wp)**2 + (para(2) - 1.0_wp)**2 - 1.0_wp)


  end function model_dlik



end module test_model_circle_t
