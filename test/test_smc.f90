module test_smc
  use, intrinsic :: iso_fortran_env, only: wp => real64

  use fruit
  implicit none

  include 'mpif.h'

contains

  subroutine test_smc_blocks

    use fortress_smc_t, only: generate_random_blocks

    integer, allocatable :: ind(:), break_points(:)

    integer :: npara, nblocks

    npara = 1
    nblocks = 1

    allocate(ind(1), break_points(2))
    ind = [1]
    call generate_random_blocks(npara, nblocks, ind, break_points)
    call assert_equals(0, break_points(1))
    call assert_equals(1, break_points(2))
    call assert_equals(1, ind(1))
    print*,ind
    deallocate(ind, break_points)

    allocate(ind(2), break_points(3))
    ind = [1,2]
    npara = 2 ; nblocks=2
    call generate_random_blocks(npara, nblocks, ind, break_points, rvs=[0.01_wp,0.9_wp])
    call assert_equals(0, break_points(1))
    call assert_equals(1, break_points(2))
    call assert_equals(2, break_points(3))
    call assert_equals(1, ind(1))

    deallocate(ind, break_points)

    npara=4 ; nblocks=2
    allocate(ind(npara), break_points(nblocks+1))
    ind = [1,2,3,4]
    call generate_random_blocks(npara, nblocks, ind, break_points, rvs=[0.01_wp,0.6_wp, 0.8_wp, 0.9_wp])
    call assert_equals(0, break_points(1))
    call assert_equals(2, break_points(2))
    call assert_equals(4, break_points(3))

    call assert_equals(1, ind(1))
    call assert_equals(2, ind(2))
    call assert_equals(3, ind(3))
    call assert_equals(4, ind(4))
    deallocate(ind, break_points)

    npara=4 ; nblocks=3
    allocate(ind(npara), break_points(nblocks+1))
    ind = [1,2,3,4]
    call generate_random_blocks(npara, nblocks, ind, break_points, rvs=[0.01_wp,0.6_wp, 0.8_wp, 0.9_wp])
    call assert_equals(0, break_points(1))
    call assert_equals(1, break_points(2))
    call assert_equals(2, break_points(3))
    call assert_equals(4, break_points(4))

    call assert_equals(1, ind(1))
    call assert_equals(2, ind(2))
    call assert_equals(3, ind(3))
    call assert_equals(4, ind(4))
    deallocate(ind, break_points)


  end subroutine test_smc_blocks

  subroutine test_smc_initialize

    use fortress, only: fortress_smc
    use test_model_t, only : model 

    type(fortress_smc) :: test_smc
    type(model) :: test_model

    integer :: mpierror, rank, nproc
    
    logical :: is_initialized

    call mpi_initialized(is_initialized, mpierror)
    if (.not. is_initialized) call mpi_init(mpierror)
    call mpi_comm_size(MPI_COMM_WORLD, nproc, mpierror)
    call mpi_comm_rank(MPI_COMM_WORLD, rank, mpierror)

    test_model = model()
    test_smc = fortress_smc(test_model, nproc)

    call assert_equals('ss', test_smc%model%name)
    call assert_equals(1, test_smc%model%nobs)
    call assert_equals(80, test_smc%model%T)
    call assert_equals(2, test_smc%model%npara)

    !call test_smc%estimate(rank)

  end subroutine test_smc_initialize


  subroutine test_smc_endog_ess

    use fortress_smc_particles_t, only: fortress_smc_particles
    use fortress_smc_t, only: ess_gap 

    type(fortress_smc_particles) :: test_particles

    real(wp) :: gap

    test_particles = fortress_smc_particles(npart=4)

    gap = ess_gap(1.0_wp, 0.0_wp, test_particles, 2.0_wp)

    call assert_equals(-1.0_wp, gap, 0.00001_wp)

  end subroutine test_smc_endog_ess


  subroutine test_smc_endog_schedule

    use fortress, only: fortress_smc
    use fortress_smc_t, only: tempering_schedule
    use test_model_t, only : model 

    type(fortress_smc) :: test_smc
    type(model) :: test_model

    integer :: mpierror, rank, nproc, i, j
    
    logical :: is_initialized

    call mpi_initialized(is_initialized, mpierror)
    if (.not. is_initialized) call mpi_init(mpierror)
    call mpi_comm_size(MPI_COMM_WORLD, nproc, mpierror)
    call mpi_comm_rank(MPI_COMM_WORLD, rank, mpierror)

    test_model = model()
    test_smc = fortress_smc(test_model, nproc)
    deallocate(test_smc%temp%T_schedule, test_smc%temp%phi_schedule, &
         test_smc%temp%Z_estimates, test_smc%temp%ESS_estimates)

    test_smc%temp = tempering_schedule(nstages=5000, lambda=2.0_wp, max_T=80)
    test_smc%endog_tempering = .true.

    !call test_smc%estimate(rank)

    test_smc%temp = tempering_schedule(nstages=5000, lambda=2.0_wp, max_T=80)
    test_smc%temp%phi_max = 100000000.0_wp

    !call test_smc%estimate(rank)

  end subroutine test_smc_endog_schedule
  
  subroutine test_smc_schedule

    use fortress, only: fortress_smc
    use fortress_smc_t, only: tempering_schedule
    use test_model_t, only : model 

    type(fortress_smc) :: test_smc
    type(model) :: test_model

    integer :: mpierror, rank, nproc, i, j
    
    logical :: is_initialized

    call mpi_initialized(is_initialized, mpierror)
    if (.not. is_initialized) call mpi_init(mpierror)
    call mpi_comm_size(MPI_COMM_WORLD, nproc, mpierror)
    call mpi_comm_rank(MPI_COMM_WORLD, rank, mpierror)

    test_model = model()
    test_smc = fortress_smc(test_model, nproc)
    deallocate(test_smc%temp%T_schedule, test_smc%temp%phi_schedule, &
         test_smc%temp%Z_estimates, test_smc%temp%ESS_estimates)
    test_smc%temp = tempering_schedule(nstages=81, lambda=0.0_wp, max_T=80)
    do i = 1,81
       test_smc%temp%T_schedule(i) = i-1
    end do

    test_smc%temp%phi_schedule = 1.0_wp
    
    call test_smc%estimate(rank)

    test_smc = fortress_smc(test_model, nproc)
    deallocate(test_smc%temp%T_schedule, test_smc%temp%phi_schedule, &
         test_smc%temp%Z_estimates, test_smc%temp%ESS_estimates)
    test_smc%temp = tempering_schedule(nstages=162, lambda=0.0_wp, max_T=80)
    j = 0
    do i = 1,162,2
       test_smc%temp%T_schedule(i) = j
       test_smc%temp%T_schedule(i+1) = j
       
       test_smc%temp%phi_schedule(i) = 0.5_wp
       test_smc%temp%phi_schedule(i+1) = 1.0_wp
       
       j = j + 1
    end do
    test_smc%temp%phi_schedule(1) = 1.0_wp
    !test_smc%temp%phi_schedule = 1.0_wp
    
    call test_smc%estimate(rank)


  end subroutine test_smc_schedule
  
end module test_smc
