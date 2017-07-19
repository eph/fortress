module test_VAR_t
  use, intrinsic :: iso_fortran_env, only: wp => real64

  use fortress_VAR_t, only: SVAR

  implicit none

  type, public, extends(SVAR) :: model
   contains
     procedure :: para_to_AF
  end type model

  interface model
     module procedure new_model
  end interface model

contains

  type(model) function new_model() result(self)

    character(len=144) :: name, datafile, prior_mu_file, prior_var_file
    integer :: nobs, T, p, constant, nA, nF

    name = 'swz'
    datafile = '/home/eherbst/Dropbox/var_smc_estimation/replication-code/smc_msvar/data.txt'
    prior_mu_file = '/home/eherbst/Dropbox/var_smc_estimation/replication-code/smc_msvar/mu.txt'
    prior_var_file = '/home/eherbst/Dropbox/var_smc_estimation/replication-code/smc_msvar/sigma.txt'

    nobs = 3
    p = 5
    constant = 1
    T = 188
    nA = 6
    nF = 48
    
    call self%construct_SVAR(name, datafile, nobs, T, p, constant, nA, nF, prior_mu_file, prior_VAR_file)

  end function

  subroutine para_to_AF(self, para, A, F)

    class(model), intent(inout) :: self

    real(wp), intent(in) :: para(self%npara)
    real(wp), intent(out) :: A(self%nobs, self%nobs)
    real(wp), intent(out) :: F(self%nobs*self%p+self%constant, self%nobs)


    A(1, 1) = 1.0d0*para(1)
    A(2, 1) = 0.0d0
    A(3, 1) = 0.0d0
    A(1, 2) = 1.0d0*para(2)
    A(2, 2) = 1.0d0*para(3)
    A(3, 2) = 0.0d0
    A(1, 3) = 1.0d0*para(4)
    A(2, 3) = 1.0d0*para(5)
    A(3, 3) = 1.0d0*para(6)

    F(1, 1) = 1.0d0*para(7)
    F(2, 1) = 1.0d0*para(8)
    F(3, 1) = 1.0d0*para(9)
    F(4, 1) = 1.0d0*para(10)
    F(5, 1) = 1.0d0*para(11)
    F(6, 1) = 1.0d0*para(12)
    F(7, 1) = 1.0d0*para(13)
    F(8, 1) = 1.0d0*para(14)
    F(9, 1) = 1.0d0*para(15)
    F(10, 1) = 1.0d0*para(16)
    F(11, 1) = 1.0d0*para(17)
    F(12, 1) = 1.0d0*para(18)
    F(13, 1) = 1.0d0*para(19)
    F(14, 1) = 1.0d0*para(20)
    F(15, 1) = 1.0d0*para(21)
    F(16, 1) = 1.0d0*para(22)
    F(1, 2) = 1.0d0*para(23)
    F(2, 2) = 1.0d0*para(24)
    F(3, 2) = 1.0d0*para(25)
    F(4, 2) = 1.0d0*para(26)
    F(5, 2) = 1.0d0*para(27)
    F(6, 2) = 1.0d0*para(28)
    F(7, 2) = 1.0d0*para(29)
    F(8, 2) = 1.0d0*para(30)
    F(9, 2) = 1.0d0*para(31)
    F(10, 2) = 1.0d0*para(32)
    F(11, 2) = 1.0d0*para(33)
    F(12, 2) = 1.0d0*para(34)
    F(13, 2) = 1.0d0*para(35)
    F(14, 2) = 1.0d0*para(36)
    F(15, 2) = 1.0d0*para(37)
    F(16, 2) = 1.0d0*para(38)
    F(1, 3) = 1.0d0*para(39)
    F(2, 3) = 1.0d0*para(40)
    F(3, 3) = 1.0d0*para(41)
    F(4, 3) = 1.0d0*para(42)
    F(5, 3) = 1.0d0*para(43)
    F(6, 3) = 1.0d0*para(44)
    F(7, 3) = 1.0d0*para(45)
    F(8, 3) = 1.0d0*para(46)
    F(9, 3) = 1.0d0*para(47)
    F(10, 3) = 1.0d0*para(48)
    F(11, 3) = 1.0d0*para(49)
    F(12, 3) = 1.0d0*para(50)
    F(13, 3) = 1.0d0*para(51)
    F(14, 3) = 1.0d0*para(52)
    F(15, 3) = 1.0d0*para(53)
    F(16, 3) = 1.0d0*para(54)

  end subroutine para_to_AF

end module test_VAR_t
