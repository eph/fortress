module test_gensys

  use, intrinsic :: iso_fortran_env, only: wp => real64

  use fruit
  implicit none

  include 'mpif.h'

contains

  subroutine test_gensys_sw

    use gensys, only : do_gensys
    use fortress_util, only : read_array_from_file

    integer :: neta, ns, neps

    real(wp), allocatable  :: GAM0(:, :), GAM1(:, :), C(:), PSI(:, :), PPI(:, :), CC(:), TT(:,:), RR(:,:)
    
    ! gensys 
    real(wp) :: fmat, fwt, ywt, gev, loose, DIV
    integer :: eu(2)


    neta = 12
    ns = 53
    neps = 7




    allocate(GAM0(ns, ns), GAM1(ns, ns), C(ns), PSI(ns, neps), PPI(ns, neta), &
         CC(ns), TT(ns,ns), RR(ns,ns))

    call read_array_from_file('test/sw/GAM0.txt', GAM0)
    call read_array_from_file('test/sw/GAM1.txt', GAM1)
    call read_array_from_file('test/sw/PSI.txt', PSI)
    call read_array_from_file('test/sw/PPI.txt', PPI)
    C = 0.0_wp

    call do_gensys(TT, CC, RR, fmat, fwt, ywt, gev, eu, loose, &
         GAM0, GAM1, C, PSI, PPI, DIV)

    call assert_equals(1, eu(1))
    call assert_equals(1, eu(2))

  end subroutine test_gensys_sw

end module test_gensys
