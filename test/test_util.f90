module test_util
  use, intrinsic :: iso_fortran_env, only: wp => real64

  use fruit
  implicit none

contains

  subroutine test_read

    use fortress_util, only: read_array_from_file, write_array_to_file

    character(len=:), allocatable :: array_file, matrix_file
    
    real(wp) :: array_dp(3), array_dp_short(2)
    integer :: array_int(3)
    logical :: array_logical(3)

    real(wp) :: matrix_dp(80,1), matrix_dp_t(1,80)

    array_file = 'test/readfile.txt'

    call read_array_from_file(array_file, array_dp)
    call assert_equals(3.0_wp, array_dp(1))
    call assert_equals(5.0_wp, array_dp(2))
    call assert_equals(1.0_wp, array_dp(3))

    call read_array_from_file(array_file, array_int)
    call assert_equals(3, array_int(1))
    call assert_equals(5, array_int(2))
    call assert_equals(1, array_int(3))
   
    ! call read_array_from_file(array_file, array_logical)
    ! call assert_equals(.true., array_logical(1))
    ! call assert_equals(.true., array_logical(2))
    ! call assert_equals(.true., array_logical(3))
    
    matrix_file = 'test/test_data.txt'
    call read_array_from_file(matrix_file, matrix_dp)
    call assert_equals(0.90512501165_wp, matrix_dp(9,1), 0.000000001_wp)

    call read_array_from_file(matrix_file, matrix_dp_t, transpose=.true.)
    call assert_equals(0.90512501165_wp, matrix_dp_t(1,9), 0.000000001_wp)

    call read_array_from_file('test/test_data.txt', matrix_dp)
    call assert_equals(0.90512501165_wp, matrix_dp(9,1), 0.000000001_wp)
    
    call write_array_to_file('test/write_test.txt', matrix_dp)
    call read_array_from_file('test/test_data.txt', matrix_dp)
    call assert_equals(0.90512501165_wp, matrix_dp(9,1), 0.000000001_wp)

  end subroutine test_read
  
end module test_util
