program smc_driver
  use iso_fortran_env, only: wp => real64

  use fortress, only: fortress_smc
  use test_model_circle_t, only: model

  implicit none
  include 'mpif.h'

  type(fortress_smc) :: smc
  type(model) :: smc_model

  integer :: mpierror, rank, nproc, i

  call mpi_init(mpierror)
  call mpi_comm_size(MPI_COMM_WORLD, nproc, mpierror)
  call mpi_comm_rank(MPI_COMM_WORLD, rank, mpierror)

  smc_model = model()
  print*,smc_model%lik(smc_model%p0)
  print*,smc_model%dlik(smc_model%p0)

  smc = fortress_smc(smc_model, nproc)
  call smc%estimate(rank)
  call mpi_finalize(mpierror)
end program smc_driver
