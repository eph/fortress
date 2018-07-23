module test_json
  use, intrinsic :: iso_fortran_env, only: wp => real64

  use fruit
  implicit none


contains

  subroutine test_json1


    use json_module, only: json_core, json_value
    type(json_core) :: json
    type(json_value),pointer :: p, inp

    character(len=4) :: my_file
    call json%initialize()

    ! ! initialize the structure:
    call json%create_object(p,'')

    ! ! add an "inputs" object to the structure:
    call json%create_object(inp,'inputs')
    call json%add(p, inp) !add it to the root

    ! add some data to inputs:
    call json%add(inp, 't0', 0.1_wp)
    call json%add(inp, 'tf', 1.1_wp)
    call json%add(inp, 'x0', 9999.0000d0)
    call json%add(inp, 'integer_scalar', 787)
    call json%add(inp, 'integer_array', [2,4,99])
    call json%add(inp, 'names', ['aaa','bbb','ccc'])
    call json%add(inp, 'logical_scalar', .true.)
    call json%add(inp, 'logical_vector', [.true., .false., .true.])
    !nullify(inp)  !don't need this anymore

    ! write the file:
    my_file = 'fdfa'
    open(1,file=my_file,action='write')
    call json%print(p,'fadf')
    close(1)
    ! !cleanup:
    call json%destroy(p)

    
    ! call assert_equals(0.0_wp, 0.0_wp, 0.001_wp)

  end subroutine test_json1

end module test_json
