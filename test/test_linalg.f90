module test_linalg
  use, intrinsic :: iso_fortran_env, only: wp => real64

  use fruit
  implicit none

contains

!  subroutine test_matpow
!
!    use fortress_linalg, only: matpow => matpow_fortran
!
!    ! create a 3x3 matrix and fill it with random numbers
!    real64, dimension(3,3) :: A
!
!    ! fill the matrix with integers 1-9
!    A(1,1) = 1.0_wp
!    A(1,2) = 2.0_wp
!    A(1,3) = 3.0_wp
!    A(2,1) = 4.0_wp
!    A(2,2) = 5.0_wp
!    A(2,3) = 6.0_wp
!    A(3,1) = 7.0_wp
!    A(3,2) = 8.0_wp
!    A(3,3) = 9.0_wp
!
!    ! multiply the matrix by itself twice manually using matmul
!    ! and compare the results with the Fortran routine
!    ! (this is a simple test, but it is useful to see if the
!    ! Fortran routine is working correctly)
!    real64, dimension(3,3) :: B
!    real64, dimension(3,3) :: C
!
!    B = matmul(A,A)
!    B = matmul(B,A)
!
!    call matpower(A,3,C)
!
!    ! compare the results using assert_equals
!    call assert_equals(maxval(abs(B-C)),1.0d-6)
!
!end subroutine test_matpow

  subroutine test_dlyap

    use fortress_linalg, only: dlyap

    real(wp) :: Pt(1,1), TT(1,1), QQ(1,1)

    real(wp) :: Pt2d(2,2), TT2d(2,2), QQ2d(2,2)
    integer :: info

    TT = 0.9_wp
    QQ = 1.0_wp

    call dlyap(TT,QQ,Pt,1,info)
    
    call assert_equals(0, info)
    call assert_equals(1.0_wp / (1.0_wp - 0.9_wp**2), Pt(1,1), 0.000001_wp)

    TT = 1.1_wp
    call dlyap(TT,QQ,Pt,1,info)
    call assert_equals(-1, info)

    TT2d(1,:) = [0.9_wp, 0.01_wp]
    TT2d(2,:) = [0.3_wp, 0.10_wp]

    QQ2d(1,:) = [1.0_wp, 0.5_wp]
    QQ2d(2,:) = [0.5_wp, 2.0_wp]

    call dlyap(TT2d,QQ2d,Pt2d,2,info)

    call assert_equals(0, info)
    call assert_equals(5.47135633_wp, Pt2d(1,1), 0.000001_wp)
    call assert_equals(2.18292845_wp, Pt2d(1,2), 0.000001_wp)
    call assert_equals(2.64989674_wp, Pt2d(2,2), 0.000001_wp)
  end subroutine test_dlyap
  
  subroutine test_determinant

    use fortress_linalg, only: determinant

    integer :: n 
    real(wp) :: res, X(2,2)
    
    n = 2
    X(1,:) = [3.0_wp, 0.5_wp]
    X(2,:) = [0.5_wp, 1.0_wp]

    res = determinant(X, n)

    call assert_equals(2.75_wp, res, 0.000001_wp)
    
  end subroutine test_determinant


  subroutine test_cholesky

    use fortress_linalg, only: cholesky

    integer :: n, info
    real(wp) :: res, X(2,2)
    
    n = 2
    X(1,:) = [1.0_wp, 0.3_wp]
    X(2,:) = [0.3_wp, 2.0_wp]

    call cholesky(X, info)
    
    call assert_equals(1.0_wp, X(1,1),0.000001_wp)
    call assert_equals(0.0_wp, X(1,2),0.000001_wp)
    call assert_equals(0.3_wp, X(2,1),0.000001_wp)
    call assert_equals(1.382075_wp, X(2,2), 0.000001_wp)

  end subroutine test_cholesky

end module test_linalg
