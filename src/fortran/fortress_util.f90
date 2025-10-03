!> Contains file utilities for fortress. 
!!
module fortress_util
  use, intrinsic :: iso_fortran_env, only: wp => real64

  implicit none


  interface read_array_from_file
#:for dtype in ['int','logical','dreal']
     module procedure read_array_from_file_${dtype}$_1d
     module procedure read_array_from_file_${dtype}$_2d
#:endfor
  end interface read_array_from_file
  
  interface write_array_to_file
#:for dtype in ['int','logical','dreal']
     module procedure write_array_to_file_${dtype}$_1d
     module procedure write_array_to_file_${dtype}$_2d
#:endfor
  end interface write_array_to_file
contains

  !> Create a directory using system commands
  !!
  !! @param[in]  dir_name  Directory name to create
  !! @param[out] retcode   Return code (0=success, 1=already exists, 2=creation failed)
  subroutine mkdir(dir_name, retcode)
    use, intrinsic :: iso_fortran_env, only: error_unit

    character(len=*), intent(in) :: dir_name
    integer, intent(out) :: retcode
    logical :: direxists
    integer :: cmdstat
    character(len=256) :: cmdmsg
    character(len=:), allocatable :: command

    ! Check if directory already exists
    inquire(file=trim(dir_name)//'/.', exist=direxists)
    if (direxists) then
       write(error_unit, '(a)') 'mkdir: directory already exists: '//trim(dir_name)
       retcode = 1
       return
    endif

    ! Build command - use mkdir -p for robustness
    command = 'mkdir -p "'//trim(dir_name)//'"'

    ! Execute command using Fortran 2008 intrinsic
    call execute_command_line(command, wait=.true., exitstat=retcode, cmdstat=cmdstat, cmdmsg=cmdmsg)

    ! Check for execution errors
    if (cmdstat /= 0) then
       write(error_unit, '(a,i0)') 'mkdir: command execution failed with status ', cmdstat
       if (len_trim(cmdmsg) > 0) write(error_unit, '(a)') 'mkdir: '//trim(cmdmsg)
       retcode = 2
       return
    endif

    ! Check mkdir exit status
    if (retcode /= 0) then
       write(error_unit, '(a,i0,a)') 'mkdir: failed to create directory (exit status ', retcode, '): '//trim(dir_name)
       retcode = 2
    endif

  end subroutine mkdir

  ! check if a file exists
   subroutine file_exists(file_name, retcode)
   
      character(len=:), allocatable, intent(in) :: file_name
      integer, intent(out) :: retcode 
      logical :: fileexists
   
      inquire(file=file_name, exist=fileexists)
      if (fileexists) then
          retcode = 1
      else
          retcode = 0
      endif
   
   end subroutine file_exists



#:for dtype, dtype_dec in [('int','integer'),('logical','logical'), ('dreal','real(wp)')]
  subroutine read_array_from_file_${dtype}$_1d(filename, array)

    !character(len=:), allocatable, intent(in) :: filename
    character(len=*), intent(in) :: filename
    ${dtype_dec}$, dimension(:), intent(inout) :: array

    integer :: i, n

    n = size(array)

    open(143,file=filename,action='read')
    do i = 1, n
       read(143,*) array(i)
    end do
    close(143)

  end subroutine read_array_from_file_${dtype}$_1d


  subroutine read_array_from_file_${dtype}$_2d(filename, array, transpose)

    !    character(len=:), allocatable, intent(in) :: filename
    character(len=*), intent(in) :: filename
    ${dtype_dec}$, dimension(:,:), intent(inout) :: array
    logical, intent(in), optional :: transpose

    logical :: transpose_t
    integer :: i, a, b, rows

    transpose_t = .false.

    if (present(transpose)) transpose_t = transpose

    a = size(array,1)
    b = size(array,2)

    if (transpose_t) then
       rows = b
    else
       rows = a
    end if

    open(143,file=filename,action='read')
    do i = 1, rows
       if (transpose_t) then
          read(143,*) array(:,i)
       else
          read(143,*) array(i,:)
       end if
    end do
    close(143)

  end subroutine read_array_from_file_${dtype}$_2d
#:endfor


#:for dtype, dtype_dec in [('int','integer'),('logical','logical'), ('dreal','real(wp)')]
  subroutine write_array_to_file_${dtype}$_1d(filename, array)

    character(len=*), intent(in) :: filename
    ${dtype_dec}$, dimension(:), intent(inout) :: array

    integer :: i, n

    n = size(array)

    open(143,file=filename,action='write')
    do i = 1, n
       write(143,*) array(i)
    end do
    close(143)

  end subroutine write_array_to_file_${dtype}$_1d


  subroutine write_array_to_file_${dtype}$_2d(filename, array, transpose)

    character(len=*), intent(in) :: filename
    ${dtype_dec}$, dimension(:,:), intent(inout) :: array
    logical, intent(in), optional :: transpose

    logical :: transpose_t
    integer :: i, a, b, rows

    transpose_t = .false.

    if (present(transpose)) transpose_t = transpose

    a = size(array,1)
    b = size(array,2)

    if (transpose_t) then
       rows = b
    else
       rows = a
    end if

    open(143,file=filename,action='write')
    do i = 1, rows
       if (transpose_t) then
          write(143,'(1000f16.8)') array(:,i)
       else
          write(143,'(1000f16.8)') array(i,:)
       end if
    end do
    close(143)

  end subroutine write_array_to_file_${dtype}$_2d
#:endfor



end module fortress_util

