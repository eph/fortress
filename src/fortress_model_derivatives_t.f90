module fortress_bayesian_model_derivatives_t

  use, intrinsic :: iso_fortran_env, only: wp => real64
  use fortress_bayesian_model_t, only: fortress_abstract_bayesian_model
  use fortress_prior_t, only: fortress_abstract_prior

  implicit none

  real(wp), parameter :: test = 1.0_wp

  ! type, abstract, extends(fortress_abstract_prior) :: fortress_abstract_prior_derivatives
  !  contains
  !    procedure(dlogpdf_interface), deferred :: dlogpdf

  ! end type fortress_abstract_prior_derivatives

  ! interface fortress_abstract_prior_derivatives

  !    function dlogpdf_interface(self, para) result(dlpdf)
  !      use, intrinsic :: iso_fortran_env, only: wp => real64
  !      import fortress_abstract_prior_derivatives
       
  !      class(fortress_abstract_prior), intent(inout) :: self
  !      real(wp), intent(in) :: para(self%npara)
  !      real(wp) :: dlpdf
  !    end function dlogpdf_interface

  ! end interface fortress_abstract_prior_derivatives

  type,public,extends(fortress_abstract_bayesian_model), abstract :: fortress_abstract_bayesian_model_derivatives
   contains
     procedure(dlik_func_i), public, deferred :: dlik
  end type fortress_abstract_bayesian_model_derivatives

  abstract interface
     function dlik_func_i(self, para, T) result(dl)
       use, intrinsic :: iso_fortran_env, only: wp => real64
       import :: fortress_abstract_bayesian_model_derivatives
       implicit none
       class(fortress_abstract_bayesian_model_derivatives), intent(inout) :: self
       real(wp), intent(in) :: para(self%npara)
       integer, intent(in), optional :: T
       real(wp) :: dl(self%npara)
     end function dlik_func_i
  end interface

contains

end module fortress_bayesian_model_derivatives_t
