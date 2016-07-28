module fortress_linalg
  use, intrinsic :: iso_fortran_env, only: wp => real64

  implicit none
  
  real(wp), parameter :: ONE = 1.0_wp

contains

  real(wp) function determinant(matrix, r) result(det)
    ! Computes the determinant of symmetric square matrix, matrix (rank r).
    integer, intent(in) :: r
    real(wp), intent(in) :: matrix(r,r)

    integer :: info, i, piv(r)
    real(wp) :: matrix_copy(r, r)

    call dcopy(r*r, matrix, 1, matrix_copy, 1)
    det = 0.0_wp

    call dpotrf('u', r, matrix_copy, r, info)

    if (info .ne. 0) then
       !write(*,'(a,i4)') 'In determinant(), dgetrf returned error code ', info
       det = -10000.0_wp
       return
    end if

    det = ONE

    do i = 1, r

       if (.true.) then !(piv(i) .ne. i) then
          det = det * matrix_copy(i, i) * matrix_copy(i,i)
       else
          det = det * matrix_copy(i, i)
       end if

    end do
  end function determinant

  subroutine dlyap(TT, RQR, P0, ns, info)
    ! Computes the solution to the discrete Lyapunov equation,
    !      P0 = TT*P0*TT' + RQR
    ! where (inputs) TT, RQR and (output) P0 are ns x ns (real) matrices. 
    !--------------------------------------------------------------------------------
    integer, intent(in) :: ns
    real(wp), intent(in) :: TT(ns,ns), RQR(ns,ns)

    integer, intent(out) :: info
    real(wp), intent(out) :: P0(ns,ns)

    ! for slicot
    real(wp) :: scale, U(ns,ns), UH(ns, ns), rcond, ferr, wr(ns), wi(ns), dwork(14*ns*ns*ns), sepd
    integer :: iwork(ns*ns), ldwork

    integer :: t

    UH = TT
    P0 = -1.0_wp*RQR

    !call sb03md('D','X', 'N', 'T', ns, UH, ns, U, ns, P0, ns, &
    !     scale, sepd, ferr, wr, wi, iwork, dwork, 14*ns*ns*ns, info)
    
    !if (ferr > 0.000001_wp) call dlyap_symm(TT, RQR, P0, ns, info)
    if (info .ne. 0) then
       print*,'SB03MD failed. (info = ', info, ')'
       P0 = 0.0_wp
       info = 1
       do t = 1,ns
          P0(t,t)=1.0_wp
       end do

       return
    else
!       P0 = 0.5_wp*P0 + 0.5_wp*transpose(P0)
       info = 0
    end if

  end subroutine dlyap



end module fortress_linalg
