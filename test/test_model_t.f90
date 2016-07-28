module test_model_t
  use, intrinsic :: iso_fortran_env, only: wp => real64

  use fortress, only : fortress_lgss_model

  implicit none

  type, public, extends(fortress_lgss_model) :: model
   contains
     procedure :: system_matrices
  end type model


  interface model
     module procedure new_model
  end interface model
  
  
contains

  type(model) function new_model() result(self)

    character(len=144) :: name, datafile
    integer :: nobs, T, ns, npara, neps
    
    name = 'ss'
    datafile = '/home/eherbst/Dropbox/code/fortress/test/test_data.txt'
    
    nobs = 1
    T = 80
    
    ns = 2
    npara = 2
    neps = 1

    call self%construct_model(name, datafile, npara, nobs, T, ns, neps)

    self%p0 = [0.2_wp, 0.5_wp]

  end function new_model

  subroutine system_matrices(self, para, error)

    class(model), intent(inout) :: self
    real(wp), intent(in) :: para(self%npara)

    integer, intent(out) :: error 
    
    real(wp) :: thet1, thet2
    real(wp) :: phi1, phi2, phi3

    thet1 = para(1)
    thet2 = para(2)

    phi1 = thet1**2
    phi2 = 1.0_wp - thet1**2
    phi3 = phi2 - thet1*thet2

    self%TT = 0.0_wp ; self%RR = 0.0_wp ; self%QQ = 0.0_wp

    self%TT(1,:) = [phi1, 0.0_wp]
    self%TT(2,:) = [phi3, phi2]

    self%RR(1,1) = 1.0_wp
    self%QQ(1,1) = 1.0_wp ** 2

    self%DD = 0.0_wp ; self%ZZ = 0.0_wp ; self%HH = 0.0_wp

    self%ZZ(1,:) = [1.0_wp, 1.0_wp]


    self%HH = 0.1_wp

    
    error = 0
  end subroutine system_matrices


end module test_model_t
