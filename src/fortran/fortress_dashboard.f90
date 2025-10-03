module fortress_dashboard
  use, intrinsic :: iso_fortran_env, only: stdout => output_unit
  implicit none
  private
  public :: update_dashboard, finish_dashboard, hide_cursor, show_cursor
  character(len=*), parameter :: ESC = achar(27)
  character(len=*), parameter :: CSI = ESC//'['
contains
  subroutine hide_cursor()
    write(stdout,'(a)', advance='no') CSI//'?25l'
    call flush(stdout)
  end subroutine hide_cursor

  subroutine show_cursor()
    write(stdout,'(a)', advance='no') CSI//'?25h'
    call flush(stdout)
  end subroutine show_cursor

  subroutine update_dashboard(msg)
    character(len=*), intent(in) :: msg
    ! \x1b[2K clears the entire line; \x1b[G moves cursor to column 1
    write(stdout,'(a)', advance='no') CSI//'2K'//CSI//'G'//msg
    call flush(stdout)
  end subroutine update_dashboard

  subroutine finish_dashboard()
    write(stdout,*)    ! drop to the next line when finished
  end subroutine finish_dashboard
end module fortress_dashboard
