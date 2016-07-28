module test_linalg
  use, intrinsic :: iso_fortran_env, only: wp => real64

  use fruit
  implicit none

contains

  subroutine test_dlyap

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
