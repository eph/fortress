module test_linalg
  use, intrinsic :: iso_fortran_env, only: wp => real64

  use fruit
  implicit none

contains

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
