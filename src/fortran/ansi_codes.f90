module ansi_codes
  use iso_fortran_env, only: wp => real64
  implicit none

  character(len=*), parameter :: ESC = achar(27)
  character(len=*), parameter :: CSI = ESC // '['

  ! Clear screen
  character(len=*), parameter :: CLS = CSI // '2J'
  character(len=*), parameter :: CLS_AFTER_CURSOR = CSI // '0J'
  character(len=*), parameter :: CLS_BEFORE_CURSOR = CSI // '1J'

  ! Clear line
  character(len=*), parameter :: CLR_LINE = CSI // '2K'
  character(len=*), parameter :: CLR_LINE_AFTER_CURSOR = CSI // '0K'
  character(len=*), parameter :: CLR_LINE_BEFORE_CURSOR = CSI // '1K'

  ! Cursor movement
  character(len=*), parameter :: CURSOR_UP = CSI // '1A'
  character(len=*), parameter :: CURSOR_DOWN = CSI // '1B'
  character(len=*), parameter :: CURSOR_FORWARD = CSI // '1C'
  character(len=*), parameter :: CURSOR_BACK = CSI // '1D'
  character(len=*), parameter :: CURSOR_HOME = CSI // 'H'

  ! Colors
  character(len=*), parameter :: BLACK = CSI // '30m'
  character(len=*), parameter :: RED = CSI // '31m'
  character(len=*), parameter :: GREEN = CSI // '32m'
  character(len=*), parameter :: YELLOW = CSI // '33m'
  character(len=*), parameter :: BLUE = CSI // '34m'
  character(len=*), parameter :: MAGENTA = CSI // '35m'
  character(len=*), parameter :: CYAN = CSI // '36m'
  character(len=*), parameter :: WHITE = CSI // '37m'
  character(len=*), parameter :: RESET = CSI // '0m'

  ! Styles
  character(len=*), parameter :: BOLD = CSI // '1m'
  character(len=*), parameter :: UNDERLINE = CSI // '4m'
  character(len=*), parameter :: BLINK = CSI // '5m'
  character(len=*), parameter :: REVERSE = CSI // '7m'
  character(len=*), parameter :: HIDDEN = CSI // '8m'

contains

  subroutine ansi_print(message, clear_line)
    character(len=*), intent(in) :: message
    logical, intent(in), optional :: clear_line

    if (present(clear_line) .and. clear_line) then
      write(*, '(A, A)', advance='no') CLR_LINE, message
    else
      write(*, '(A)', advance='no') message
    end if
  end subroutine ansi_print

  function progress_bar(value, max_value, width) result(bar)
    real(wp), intent(in) :: value
    real(wp), intent(in) :: max_value
    integer, intent(in) :: width
    character(len=width) :: bar

    integer :: i, progress
    real(wp) :: percentage

    percentage = value / max_value
    progress = int(width * percentage)

    bar = ''
    do i = 1, width
      if (i <= progress) then
        bar(i:i) = '='
      else
        bar(i:i) = ' '
      end if
    end do
  end function progress_bar

end module ansi_codes
