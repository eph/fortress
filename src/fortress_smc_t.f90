module fortress_smc_t
  use, intrinsic :: iso_fortran_env, only: wp => real64

  use flap, only : command_line_interface
  use fortress_bayesian_model_t, only: fortress_abstract_bayesian_model
  !use fortress_particles, only: particle_system
  use fortress_info

  !use mpi

  implicit none

  include 'mpif.h'
  integer, parameter :: DEFAULT_NPART = 4800
  integer, parameter :: DEFAULT_NPHI = 500
  integer, parameter :: DEFAULT_NINTMH = 1
  integer, parameter :: DEFAULT_TRIAL = 0
  real(wp), parameter :: DEFAULT_LAMBDA = 2.0_wp

  integer, parameter :: DEFAULT_WRITE_EVERY = 10
  logical, parameter :: DEFAULT_SAVE_HYPER = .false.
  logical, parameter :: DEFAULT_FIXED_HYPER = .false.

  real(wp), parameter :: DEFAULT_INITIAL_C = 0.4_wp
  real(wp), parameter :: DEFAULT_MIX_ALPHA = 0.0d0
  logical, parameter :: DEFAULT_MCMC_MIX = .false.

  integer, parameter :: DEFAULT_SEED = 1848

  character(len=2), parameter :: DEFAULT_OUTPUT_DIR = './'
  type fortress_smc

     type(command_line_interface) :: cli
     class(fortress_abstract_bayesian_model), allocatable :: model

     integer :: npart, nphi, nintmh, trial
     real(wp) :: lambda

     integer :: write_every
     logical :: save_hyper, fixed_hyper

     real(wp) :: initial_c, mix_alpha
     logical :: mcmc_mix, conditional_covariance

     integer :: seed

     character(len=200) :: output_dir 

   contains
     procedure :: estimate
  end type fortress_smc
     



  interface fortress_smc
     module procedure new_smc
  end interface fortress_smc

contains

  type(fortress_smc) function new_smc(model_p) result(smc)

    class(fortress_abstract_bayesian_model) :: model_p

    integer :: err
    allocate(smc%model, source=model_p)

    smc%cli = initialize_cli()
    call smc%cli%parse()

    ! parse the command line
    call smc%cli%get(switch='-n',val=smc%npart,error=err); if (err/=0) stop 1
    call smc%cli%get(switch='-p',val=smc%nphi,error=err); if (err/=0) stop 1
    call smc%cli%get(switch='-b',val=smc%lambda,error=err); if (err/=0) stop 1
    call smc%cli%get(switch='-m',val=smc%nintmh,error=err); if (err/=0) stop 1
    call smc%cli%get(switch='-i',val=smc%trial,error=err); if (err/=0) stop 1
    call smc%cli%get(switch='-n',val=smc%npart,error=err); if (err/=0) stop 1
    call smc%cli%get(switch='-we',val=smc%write_every,error=err); if (err/=0) stop 1
    call smc%cli%get(switch='-sh',val=smc%save_hyper,error=err); if (err/=0) stop 1
    call smc%cli%get(switch='-fh',val=smc%fixed_hyper,error=err); if (err/=0) stop 1
    call smc%cli%get(switch='-od',val=smc%output_dir,error=err); if (err /=0) stop 1
    print*,'Initializing SMC for model: ', smc%model%name
    print*, smc%model%lik(smc%model%p0)
    print*, smc%npart, smc%nphi

  end function new_smc

  subroutine estimate(self)

    class(fortress_smc), intent(inout) :: self

    character(len=:), allocatable :: estimation_name

    integer :: mpierror, rank, nproc

    real(wp), allocatable :: arate(:,:), loglh(:), prio(:), ESS(:), zi(:), loglhold(:)
    real(wp), allocatable :: incwt(:), wtsq(:)
    integer, allocatable :: acptsim(:,:), resamp(:)

    real(wp), allocatable :: nodepara(:,:), nodeloglh(:), nodewt(:), nodeprio(:)
    real(wp), allocatable :: nodeeps(:,:), nodeu(:), nodevar(:,:)
    integer, allocatable :: nodeacpt(:,:)
    
    estimation_name = self%output_dir//'smc_'//self%model%name

    

    !call write_smc_settings(estimation_name//'.json')

    !call mpi_init(mpierror)
    !call mpi_comm_size(MPI_COMM_WORLD, nproc, mpierror)
    !call mpi_comm_rank(MPI_COMM_WORLD, rank, mpierror)
    ! 
    !call mpi_finalize(mpierror)


  end subroutine estimate


  !subroutine write_smc_settings(

  type(command_line_interface) function initialize_cli() result(cli)

    integer :: error

    call cli%init(progname='SMC Estimation via FORTRESS', &
         version=FORTRESS_VERSION, &
         authors=FORTRESS_AUTHOR)

    call cli%add(switch='--npart', switch_ab='-n', help='Number of particles (integer)', &
         required=.false.,act='store',def='4800', error=error)
    call cli%add(switch='--nphi', switch_ab='-p', help='Number of stages (integer)', &
         required=.false.,act='store',def='500', error=error)
    call cli%add(switch='--bend', switch_ab='-b', help='Lambda coefficient (float)', &
         required=.false.,act='store',def='2.0', error=error)
    call cli%add(switch='--nintmh', switch_ab='-m', help='Number of MH steps (integer)', &
         required=.false.,act='store',def='1', error=error)
    call cli%add(switch='--trial', switch_ab='-i', help='Trial number (integer)', &
         required=.false.,act='store',def='0', error=error)
    call cli%add(switch='--write-every', switch_ab='-we', &
         help='Write every nth stage (integer).  [Always writes first and last]', &
         required=.false.,act='store',def='10',error=error)
    call cli%add(switch='--save-hyper', switch_ab='-sh', &
         help='Write every nth stage (integer).  [Always writes first and last]', &
         required=.false.,act='store_true',def='.false.',error=error)
    call cli%add(switch='--fixed-hyper', switch_ab='-fh', &
         help='Write every nth stage (integer).  [Always writes first and last]', &
         required=.false.,act='store_true',def='.false.',error=error)
    call cli%add(switch='--initial-scale', switch_ab='-c', &
         help='Initial scaling in the random walk metropolis-hastings algorithm.', &
         required=.false.,act='store',def='0.4',error=error)
    call cli%add(switch='--seed', switch_ab='-s', &
         help='Seed for random number generation.', &
         required=.false.,act='store',def='1848',error=error)
    call cli%add(switch='--no-mix', switch_ab='-u', &
         required=.false.,act='store_true',def='.true.',error=error)
    call cli%add(switch='--output_dir', switch_ab='-od', &
         required=.false.,act='store',def='./',error=error, help='The base output directory.')
  end function initialize_cli


end module fortress_smc_t


  
