module fruit
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none

  private

  integer :: total_cases = 0
  integer :: failed_cases = 0
  integer :: failed_assertions = 0
  character(len=256) :: current_case = ''

  public :: init_fruit
  public :: run_test_case
  public :: fruit_summary
  public :: fruit_finalize
  public :: assert_equals

  interface assert_equals
    module procedure assert_equals_real
    module procedure assert_equals_int
    module procedure assert_equals_logical
    module procedure assert_equals_char
  end interface assert_equals

  abstract interface
    subroutine fruit_test_proc()
    end subroutine fruit_test_proc
  end interface

contains

  subroutine init_fruit()
    total_cases = 0
    failed_cases = 0
    failed_assertions = 0
    current_case = ''
  end subroutine init_fruit

  subroutine run_test_case(proc, name)
    procedure(fruit_test_proc) :: proc
    character(len=*), intent(in) :: name
    integer :: before_failures

    current_case = name
    before_failures = failed_assertions

    call report_message('case start')
    call proc()

    total_cases = total_cases + 1
    if (failed_assertions > before_failures) then
      failed_cases = failed_cases + 1
      call report_message('case failed')
    else
      call report_message('case passed')
    end if

    current_case = ''
  end subroutine run_test_case

  subroutine fruit_summary()
    if (failed_assertions == 0) then
      write(*,'(A,I0,A)') 'FRUIT: all ', total_cases, ' assertions passed.'
    else
      write(*,'(A,I0,A,I0,A)') 'FRUIT: ', failed_assertions, ' failed assertions across ', &
        total_cases, ' test cases.'
    end if
  end subroutine fruit_summary

  subroutine fruit_finalize()
    ! No-op placeholder for compatibility with the classic FRUIT API.
  end subroutine fruit_finalize

  subroutine assert_equals_real(expected, actual, tol)
    real(real64), intent(in) :: expected
    real(real64), intent(in) :: actual
    real(real64), intent(in), optional :: tol
    real(real64) :: tolerance

    tolerance = 0.0_real64
    if (present(tol)) tolerance = tol

    if (present(tol)) then
      if (abs(actual - expected) > tolerance) then
        call report_real_failure(expected, actual, tolerance)
      end if
    else
      if (actual /= expected) then
        call report_real_failure(expected, actual, tolerance)
      end if
    end if
  end subroutine assert_equals_real

  subroutine assert_equals_int(expected, actual)
    integer, intent(in) :: expected
    integer, intent(in) :: actual

    if (actual /= expected) then
      call report_int_failure(expected, actual)
    end if
  end subroutine assert_equals_int

  subroutine assert_equals_logical(expected, actual)
    logical, intent(in) :: expected
    logical, intent(in) :: actual

    if (actual .neqv. expected) then
      call report_logical_failure(expected, actual)
    end if
  end subroutine assert_equals_logical

  subroutine assert_equals_char(expected, actual)
    character(len=*), intent(in) :: expected
    character(len=*), intent(in) :: actual

    if (trim(actual) /= trim(expected)) then
      call report_char_failure(expected, actual)
    end if
  end subroutine assert_equals_char

  subroutine report_real_failure(expected, actual, tol)
    real(real64), intent(in) :: expected
    real(real64), intent(in) :: actual
    real(real64), intent(in) :: tol
    character(len=256) :: buffer

    write(buffer,'(A,ES16.8E3,A,ES16.8E3,A,ES16.8E3,A)') &
      'expected=', expected, ' actual=', actual, ' tol=', tol, ' (real)'
    call register_failure(buffer)
  end subroutine report_real_failure

  subroutine report_int_failure(expected, actual)
    integer, intent(in) :: expected
    integer, intent(in) :: actual
    character(len=256) :: buffer

    write(buffer,'(A,I0,A,I0,A)') 'expected=', expected, ' actual=', actual, ' (int)'
    call register_failure(buffer)
  end subroutine report_int_failure

  subroutine report_logical_failure(expected, actual)
    logical, intent(in) :: expected
    logical, intent(in) :: actual
    character(len=256) :: buffer

    write(buffer,'(A,L1,A,L1,A)') 'expected=', expected, ' actual=', actual, ' (logical)'
    call register_failure(buffer)
  end subroutine report_logical_failure

  subroutine report_char_failure(expected, actual)
    character(len=*), intent(in) :: expected
    character(len=*), intent(in) :: actual
    character(len=256) :: buffer

    write(buffer,'(A,A,A,A,A)') 'expected="', trim(expected), '" actual="', trim(actual), '" (char)'
    call register_failure(buffer)
  end subroutine report_char_failure

  subroutine register_failure(message)
    character(len=*), intent(in) :: message

    failed_assertions = failed_assertions + 1
    if (len_trim(current_case) > 0) then
      write(*,'(A,A)') '[FRUIT]['//trim(current_case)//'] ', trim(message)
    else
      write(*,'(A,A)') '[FRUIT] ', trim(message)
    end if
  end subroutine register_failure

  subroutine report_message(summary)
    character(len=*), intent(in) :: summary

    if (len_trim(current_case) > 0) then
      write(*,'(A,A)') '[FRUIT]['//trim(current_case)//'] ', trim(summary)
    else
      write(*,'(A,A)') '[FRUIT] ', trim(summary)
    end if
  end subroutine report_message

end module fruit
