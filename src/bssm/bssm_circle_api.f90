module bssm_circle_api
  use, intrinsic :: iso_c_binding, only: c_int64_t, c_int, c_double
  use, intrinsic :: iso_fortran_env, only: wp => real64
  use test_model_circle_t, only: model, new_model
  implicit none

  type :: handle_entry
     type(model), allocatable :: m
     logical :: used = .false.
  end type handle_entry

  type(handle_entry), allocatable, save :: pool(:)

contains

  integer(c_int64_t) function bssm_circle_create() bind(C, name="bssm_circle_create")
    integer :: i
    if (.not. allocated(pool)) then
      allocate(pool(1))
    end if
    do i=1, size(pool)
      if (.not. pool(i)%used) then
        allocate(pool(i)%m)
        pool(i)%m = new_model()
        pool(i)%used = .true.
        bssm_circle_create = int(i, c_int64_t)
        return
      end if
    end do
    ! grow pool
    call grow_pool()
    i = size(pool)
    allocate(pool(i)%m)
    pool(i)%m = new_model()
    pool(i)%used = .true.
    bssm_circle_create = int(i, c_int64_t)
  end function bssm_circle_create

  subroutine grow_pool()
    type(handle_entry), allocatable :: tmp(:)
    integer :: n
    if (.not. allocated(pool)) then
      allocate(pool(4))
      return
    end if
    n = size(pool)
    allocate(tmp(n+max(4,n)))
    tmp(:)%used = .false.
    tmp(1:n) = pool
    call move_alloc(tmp, pool)
  end subroutine grow_pool

  subroutine bssm_circle_free(h) bind(C, name="bssm_circle_free")
    integer(c_int64_t), value :: h
    integer :: i
    i = int(h)
    if (allocated(pool) .and. i>=1 .and. i<=size(pool)) then
      if (pool(i)%used) then
        if (allocated(pool(i)%m)) deallocate(pool(i)%m)
        pool(i)%used = .false.
      end if
    end if
  end subroutine bssm_circle_free

  subroutine bssm_circle_prior_logpdf(h, params, n, out) bind(C, name="bssm_circle_prior_logpdf")
    integer(c_int64_t), value :: h
    integer(c_int),    value :: n
    real(c_double), intent(in) :: params(*)
    real(c_double), intent(out) :: out
    integer :: i
    real(wp), allocatable :: p(:)
    i = int(h)
    if (.not. allocated(pool) .or. .not. pool(i)%used) then
      out = -1.0d300
      return
    end if
    allocate(p(n))
    p = real(params(1:n), wp)
    out = pool(i)%m%prior%logpdf(p)
  end subroutine bssm_circle_prior_logpdf

  subroutine bssm_circle_lik_logpdf(h, params, n, t, out) bind(C, name="bssm_circle_lik_logpdf")
    integer(c_int64_t), value :: h
    integer(c_int),    value :: n
    integer(c_int),    value :: t
    real(c_double), intent(in) :: params(*)
    real(c_double), intent(out) :: out
    integer :: i
    real(wp), allocatable :: p(:)
    i = int(h)
    if (.not. allocated(pool) .or. .not. pool(i)%used) then
      out = -1.0d300
      return
    end if
    allocate(p(n))
    p = real(params(1:n), wp)
    out = pool(i)%m%lik(p, T=int(t, kind=kind(t)))
  end subroutine bssm_circle_lik_logpdf

end module bssm_circle_api

