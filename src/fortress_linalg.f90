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

  subroutine inverse(X,info)

    real(wp), intent(inout) :: X(:,:)
    integer, intent(out) :: info

    integer :: r, c, i, ipiv(size(X,1))
    real(wp) :: work(3*size(X,1))

    r = size(X,1)
    c = size(X,1)

    info = -1
    if (r/=c) return

    call dgetrf(r,r,X,r,ipiv,info)
    call dgetri(r,X,r,ipiv,work,3*r,info)


  end subroutine inverse


  subroutine cholesky(X,info) 

    real(wp), intent(inout) :: X(:,:)
    integer, intent(out) :: info

    integer :: r, c, i

    r = size(X,1)
    c = size(X,1)

    info = -1
    if (r/=c) return

    call dpotrf('l',r,X,r,info)
    do i = 1,r-1
       X(i,i+1:r) = 0.0_wp
    end do
    
  end subroutine

  subroutine Kronecker(ma, na, mb, nb, mc, nc, alpha, beta, A, B, C)

    integer, intent(in) :: ma, na, mb, nb, mc, nc
    real(wp), intent(in) :: alpha, beta
    real(wp), dimension(ma,na), intent(in) :: A
    real(wp), dimension(mb,nb), intent(in) :: B
    real(wp), dimension(mc,nc), intent(inout) :: C

    integer :: i, j, i1, j1  

    do j = 1, na
       do i = 1, ma
          i1 = (i-1)*mb + 1
          j1 = (j-1)*nb + 1 
          C(i1:i1+mb-1,j1:j1+nb-1) = alpha*C(i1:i1+mb-1,j1:j1+nb-1) + beta*A(i,j)*B 
       end do
    end do

  end subroutine Kronecker

#:if defined('SLICOT')
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
#:else
! from elmar
SUBROUTINE DLYAP(A, QQ, Sigma, nx, status)

    ! doubling, calling DSYMM and DGEMM

    ! Sigma = A * Sigma * A' + B * B'

    ! output Sigma is symmetric

    IMPLICIT NONE

    integer, intent(in) :: nx
    
    integer, intent(out) :: status

    real(wp), intent(in) :: QQ(nx,nx), A(nx,nx)
    real(wp), intent(out) :: Sigma(nx,nx)
    
    

    INTEGER, PARAMETER :: maxiter = 100

    DOUBLE PRECISION, PARAMETER :: tol = 1.0d-8

 
    INTEGER :: iter, i
    LOGICAL :: converged

    DOUBLE PRECISION, DIMENSION(Nx,Nx) :: AA, AAA, AASigma, Sigma0


 
    Sigma0 = QQ
    ! Sigma0 = B B'

    ! Sigma0 = 0.0d0

    ! call DSYRK('U','N',Nx,Nw,1.0d0,B,Nx,0.0d0,Sigma0,Nx)

    ! ! fill up lower triangular -- necessary for DGEMM below

    ! FORALL (i=2:Nx) Sigma0(i,1:i-1) = Sigma0(1:i-1,i)

 

    converged = .false.

    iter = 0

 

    AA = A
    DO

 

       iter = iter + 1

 

       ! call sandwichplus(Sigma, AA, Nx, Sigma0, Nx)

       ! MANUAL SANDWICHPLUS: Sigma = AA * Sigma0 * AA' + Sigma

       call DSYMM('R','U',Nx,Nx,1.0d0,Sigma0,Nx,AA,Nx,0.0d0,AASigma,Nx)

       Sigma = Sigma0 ! this line requires Sigma0 to

       call DGEMM('N','T',Nx,Nx,Nx,1.0d0,AASigma,Nx,AA,Nx,1.0d0,Sigma,Nx)

 

       ! balance for symmetry

       Sigma = 0.5d0 * (Sigma + transpose(Sigma))

 

       IF (abs(maxval(Sigma - Sigma0)) < tol) converged = .true.      

 

       ! print *, iter, abs(maxval(Sigma - Sigma0)), tol

       ! Sigma = (Sigma + transpose(Sigma)) / dble(2)

 

       IF (converged .OR. (iter > maxiter)) EXIT

 

       ! AAA = AA * AA

       call DGEMM('N','N',Nx,Nx,Nx,1.0d0,AA,Nx,AA,Nx,0.0d0,AAA,Nx)

       AA     = AAA

       Sigma0 = Sigma

 

    END DO

 

    IF (converged) THEN

       status = 0

    ELSE

       status = -1

    END IF

 

 

  END SUBROUTINE DLYAP
#:endif


end module fortress_linalg
