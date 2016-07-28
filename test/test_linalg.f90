module test_linalg
  use, intrinsic :: iso_fortran_env, only: wp => real64

  use fruit
  implicit none

contains

  subroutine test_dlyap

    use fortress_linalg, only: dlyap

    real(wp) :: Pt(1,1), TT(1,1), QQ(1,1)

    integer :: info

    TT = 0.9_wp
    QQ = 1.0_wp

    call dlyap(TT,QQ,Pt,1,info)
    
    call assert_equals(0, info)
    call assert_equals(1.0_wp / (1.0_wp - 0.9_wp**2), Pt(1,1), 0.000001_wp)

    TT = 1.1_wp
    call dlyap(TT,QQ,Pt,1,info)
    call assert_equals(-1, info)

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

end module test_linalg
