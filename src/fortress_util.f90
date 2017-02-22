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

  !TODO: rewrite this
  subroutine mkdir(dir_name, retcode)

    character(len=:), allocatable, intent(in) :: dir_name
    integer, intent(out) :: retcode 
    logical :: direxists

    inquire(file=dir_name, exist=direxists)
    if (direxists) then
       print *, 'directory already exists'
       retcode = 1
    else
       call system('mkdir '//dir_name)
       retcode = 0
    endif

  end subroutine mkdir




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
          write(143,'(1000f)') array(:,i)
       else
          write(143,'(1000f)') array(i,:)
       end if
    end do
    close(143)

  end subroutine write_array_to_file_${dtype}$_2d
#:endfor



end module fortress_util

