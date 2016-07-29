module test_particle_filter
  use, intrinsic :: iso_fortran_env, only: wp => real64

  use fruit
  implicit none

  include 'mpif.h'
  
contains

  subroutine test_particle_filter_1

    use test_model_t, only : model
    use fortress_particle_filter_t, only: ParallelParticleFilter

    type(model) :: test_model
    type(ParallelParticleFilter) :: ppf

    integer :: rank, nproc, npart
    real(wp) :: lik0 

    npart = 100000
    test_model = model()
    ppf = ParallelParticleFilter(test_model, npart=npart, seed=1848, nproc=1, rank=0)

    call assert_equals(.false., ppf%USE_MPI)
    call assert_equals(npart, ppf%npart)
    call assert_equals(0, ppf%nforeignpart)

    call mpi_init()

    ppf%USE_MPI = .true.    
    lik0 = ppf%lik(test_model%p0,0,1)

    ! verify this 
    call assert_equals(-85.82_wp, lik0, 0.05_wp)


  end subroutine test_particle_filter_1
end module test_particle_filter
