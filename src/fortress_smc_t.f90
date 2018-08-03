module fortress_smc_t
  use, intrinsic :: iso_fortran_env, only: wp => real64, &
       stdout => output_unit, stderr => error_unit
  !use, intrinsic :: ieee_arithmetic, only: isnan => ieee_is_nan 

  use flap, only : command_line_interface
  use fortress_bayesian_model_t, only: fortress_abstract_bayesian_model
  use fortress_smc_particles_t, only: fortress_smc_particles
  use fortress_random_t, only: fortress_random
  use fortress_linalg, only: cholesky, inverse
  use fortress_util, only: read_array_from_file

  use json_module, only: json_core, json_value
  use fortress_info

  !use mpi

  implicit none

  include 'mpif.h'


  character(len=2), parameter :: DEFAULT_OUTPUT_DIR = './'

  type :: tempering_schedule

     integer :: nstages

     real(wp), allocatable :: phi_schedule(:), Z_estimates(:), ESS_estimates(:)
     integer, allocatable :: T_schedule(:)

     real(wp) :: phi_max

     contains
       procedure :: write_json

  end type tempering_schedule

  interface tempering_schedule
     module procedure new_tempering_schedule
  end interface tempering_schedule


  type fortress_smc
     type(command_line_interface) :: cli

     class(fortress_abstract_bayesian_model), allocatable :: model

     integer :: npart, nphi, nintmh, trial, nblocks, npriorextra

     integer :: start_geweke
     real(wp) :: lambda
     integer :: write_every
     logical :: save_hyper, fixed_hyper
     real(wp) :: initial_c, mix_alpha
     logical :: mcmc_mix, conditional_covariance
     integer :: seed
     character(len=200) :: output_file, init_file

     integer :: nproc
     integer :: ngap


     type(fortress_random) :: rng
     type(tempering_schedule) :: temp

     logical :: endog_tempering
     real(wp) :: resample_tol

   contains
     procedure :: estimate
     procedure :: draw_from_prior
     procedure :: evaluate_time_t_lik
  end type fortress_smc

  interface fortress_smc
     module procedure new_smc
  end interface fortress_smc

contains

  type(tempering_schedule) function new_tempering_schedule(nstages, lambda, max_T) result(self)

    integer, intent(in) :: nstages, max_T
    real(wp), intent(in) :: lambda

    integer :: i


    self%nstages = nstages
    self%phi_max = 1.0_wp

    allocate(self%phi_schedule(nstages), self%T_schedule(nstages), &
         self%Z_estimates(nstages), self%ESS_estimates(nstages))
    do i = 1, self%nstages
       self%phi_schedule(i) = self%phi_max * (i-1.0_wp)/(self%nstages - 1.0_wp)
    end do

    self%phi_schedule = self%phi_schedule**lambda
    self%phi_schedule(nstages) = self%phi_max
    self%Z_estimates(:) = 1.0_wp
    self%ESS_estimates(1) = 0.0_wp !npart
    self%T_schedule = max_T

  end function new_tempering_schedule

  subroutine write_json(self, json_node)

    class(tempering_schedule) :: self

    type(json_value), pointer, intent(inout) :: json_node
    !type(json_value), pointer :: json_temp
    type(json_core) :: json

    call json%add(json_node, 'tempering_phi', self%phi_schedule)
    call json%add(json_node, 'tempering_T', self%T_schedule)
    call json%add(json_node, 'Z_estimates', self%Z_estimates)
    call json%add(json_node, 'ESS_estimates', self%ESS_estimates)

  end subroutine write_json

  type(fortress_smc) function new_smc(model_p, nproc) result(smc)

    class(fortress_abstract_bayesian_model) :: model_p

    integer :: nproc

    integer :: rank
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
    call smc%cli%get(switch='-od',val=smc%output_file,error=err); if (err /=0) stop 1
    call smc%cli%get(switch='-s',val=smc%seed,error=err); if (err /=0) stop 1
    call smc%cli%get(switch='-nb',val=smc%nblocks,error=err); if (err /=0) stop 1
    call smc%cli%get(switch='-cc',val=smc%conditional_covariance,error=err); if (err /=0) stop 1
    call smc%cli%get(switch='-rt',val=smc%resample_tol,error=err); if (err/=0) stop 1
    call smc%cli%get(switch='-et',val=smc%endog_tempering,error=err); if (err/=0) stop 1
    call smc%cli%get(switch='-pf',val=smc%init_file,error=err); if(err/=0) stop 1
    call smc%cli%get(switch='-pe',val=smc%npriorextra,error=err); if(err/=0) stop 1
    
    print*,'Initializing SMC for model: ', smc%model%name
    print*, smc%npart, smc%nphi

    if (.not.(mod(smc%npart, nproc) == 0)) then
       write(stderr, '(a)') 'ERROR: npart must be divisible by nsim!'
       stop
    end if

    smc%nproc = nproc
    smc%ngap= int(smc%npart/smc%nproc)

    smc%temp = tempering_schedule(smc%nphi, smc%lambda, smc%model%T)
    smc%nphi = smc%temp%nstages

  end function new_smc

  subroutine estimate(self, rank)

    class(fortress_smc), intent(inout) :: self

    integer, intent(in) :: rank

    character(len=:), allocatable :: estimation_name

    integer :: mpierror, nproc

    integer :: i, j, k, m, info

    type(json_core) :: json
    type(json_value), pointer :: json_p, json_ip
    integer, allocatable :: resamp(:)

    real(wp) :: uu(self%temp%nstages, 1)
    real(wp) :: u(self%npart*self%nintmh*self%nblocks,1)
    integer :: resample_ind(self%npart)
    real(wp) :: eps(self%model%npara*self%nintmh,self%npart)

    type(fortress_smc_particles) :: parasim, nodepara
    integer :: acptsim(self%model%npara, self%npart)

    real(wp) :: nodeeps(self%model%npara * self%nintmh,self%ngap)
    real(wp) :: nodeu(self%ngap*self%nblocks * self%nintmh, 1)
    integer :: nodeacpt(self%model%npara,self%ngap)
    real(wp) :: nodeloglh(self%ngap), nodepriorsim(self%ngap)

    logical :: add_new_observation
    real(wp) :: phi, phi_old
    integer :: previous_T, current_T

    real(wp) :: mean(self%model%npara), variance(self%model%npara, self%model%npara), chol_var(self%model%npara, self%model%npara)

    real(wp) :: p0(self%model%npara), loglh0, loglhold0, prior0
    real(wp) :: likdiff, likdiffold, prdiff, alp
    integer :: nodeuind

    character(len=3) :: chart

    integer :: npara_old, bj, jj, bsize, break_points(self%nblocks+1), ind(self%model%npara), restsize
    integer, allocatable :: b_ind(:), rest_ind(:)
    real(wp) :: variance_copy(self%model%npara, self%model%npara)
    real(wp), allocatable :: block_variance_chol(:,:), temp(:,:), other_var(:,:)

    real(wp) :: scale, ahat

    real(wp) :: ess_gap1, phi0

    scale = 0.4_wp
    !if (self%endog_tempering) self%resample_tol = 0.0_wp

    self%rng = fortress_random(seed=self%seed+rank)

    ! random numbers for resampling
    if (rank == 0) uu = self%rng%uniform_rvs(self%temp%nstages, 1)

    ! draw from the prior
    parasim = fortress_smc_particles(nvars=self%model%npara, npart=self%npart)
    if (rank == 0) then

       if (self%init_file == 'none') then
          parasim%particles = self%model%prior%rvs(self%npart, rng=self%rng)
       else
          call read_array_from_file(self%init_file, parasim%particles)
          print*,'Reading particles from', self%init_file
       end if
       call json%create_object(json_p,'')
       call json%add(json_p, 'model_name', self%model%name)
       call json%add(json_p, 'nproc', self%nproc)

    end if
    nodepara = fortress_smc_particles(nvars=self%model%npara, npart=self%ngap)

    call scatter_particles(parasim, nodepara)
    call parasim%mean_and_variance(mean, variance)
    call self%draw_from_prior(nodepara, self%temp%T_schedule(1))

    call gather_particles(parasim, nodepara)

    if (rank == 0) then

       call json%create_object(json_ip,'prior')
       call json%add(json_p, json_ip)
       call parasim%write_json(json_ip)

       nullify(json_ip)

    end if

    i = 2
    do while (i <= self%temp%nstages)
    !do i = 2, self%temp%nstages
       previous_T = self%temp%T_schedule(i-1)
       current_T = self%temp%T_schedule(i)

       phi = self%temp%phi_schedule(i)
       phi_old = self%temp%phi_schedule(i-1)

       add_new_observation = abs(previous_T - current_T) > 0

       if (add_new_observation) then
          parasim%loglhold = parasim%loglh
          call scatter_particles(parasim, nodepara)
          call self%evaluate_time_t_lik(nodepara, current_T)
          call gather_particles(parasim, nodepara)
          if (phi_old == 1.0_wp) then
             phi_old = 0.0_wp
          else
             write(stderr, '(a)') 'ERROR: tempering schedule misspecified.'
          end if

       end if

       if (rank==0) then

          if (self%endog_tempering .eqv. .true.) then

             phi0 = self%temp%phi_max
             ess_gap1 = ess_gap(phi0, phi_old, parasim, self%resample_tol)

             if (ess_gap1 > 0.0_wp) then
                phi = self%temp%phi_max
                if (current_T == self%model%T) self%temp%nstages = i
             else
                do while (isnan(ess_gap1))
                   phi0 = max(phi0 / 2.0_wp, phi_old+0.01_wp)
                   if (phi0 < phi_old) stop
                   ess_gap1 = ess_gap(phi0, phi_old, parasim, self%resample_tol)
                end do
                phi = bisection(phi_old, phi0, 0.0001_wp, phi_old, parasim, self%resample_tol)
             end if
             ess_gap1 = ess_gap(phi, phi_old, parasim, self%resample_tol)

             self%temp%phi_schedule(i) = phi

          end if

          !------------------------------------------------------------
          ! Correction
          !------------------------------------------------------------
          do j = 1, parasim%npart
             parasim%weights(j) = parasim%weights(j) * exp( &
                  (phi - phi_old)*(parasim%loglh(j)-parasim%loglhold(j)) )
          end do
             
          call parasim%normalize_weights(self%temp%Z_estimates(i))
          self%temp%ESS_estimates(i) = parasim%ESS()
          print*,'iteration ', i,' of ', self%temp%nstages

          !------------------------------------------------------------
          ! Selection
          !------------------------------------------------------------
          !print*,self%temp%ESS_estimates(i), ahat, scale, phi, current_T
          if (self%temp%ESS_estimates(i) < self%npart * self%resample_tol) then

             call parasim%systematic_resampling(uu(i,1), resample_ind)

             parasim%prior = parasim%prior(resample_ind)
             parasim%loglh = parasim%loglh(resample_ind)
             parasim%loglhold = parasim%loglhold(resample_ind)

          end if

          eps = self%rng%norm_rvs(self%model%npara*self%nintmh, self%npart)
          u = self%rng%uniform_rvs(self%nblocks*self%nintmh*self%npart, 1)
          call parasim%mean_and_variance(mean, variance)
          do bj = 1, self%model%npara
             print*,bj,mean(bj), sqrt(variance(bj,bj))
          end do
          ! blocking
          ind = [ (bj, bj = 1,self%model%npara )]
          call generate_random_blocks(self%model%npara, self%nblocks, ind, break_points)

       end if

       call mpi_scatter(eps, self%ngap*self%nintmh*self%model%npara, MPI_DOUBLE_PRECISION, nodeeps, &
            self%ngap*self%nintmh*self%model%npara, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror)
       call mpi_scatter(u, self%ngap*self%nintmh*self%nblocks, MPI_DOUBLE_PRECISION, nodeu, &
            self%ngap*self%nintmh*self%nblocks, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror)
       call mpi_bcast(variance, self%model%npara*self%model%npara, MPI_DOUBLE_PRECISION, &
            0, MPI_COMM_WORLD, mpierror)

       call mpi_bcast(ind, self%model%npara, MPI_INTEGER, &
            0, MPI_COMM_WORLD, mpierror)
       call mpi_bcast(break_points, self%nblocks+1, MPI_INTEGER, &
            0, MPI_COMM_WORLD, mpierror)
       call scatter_particles(parasim, nodepara)


       !------------------------------------------------------------
       ! mutation
       !------------------------------------------------------------
       chol_var = 0.0_wp
       do bj = 1, self%nblocks
          bsize = break_points(bj+1) - break_points(bj)
          allocate(b_ind(bsize), block_variance_chol(bsize,bsize))
          b_ind = ind(break_points(bj)+1:break_points(bj+1))


          block_variance_chol = variance(b_ind, b_ind)
          if (self%conditional_covariance .and. (self%nblocks>1)) then
             restsize = self%model%npara-bsize
             allocate(rest_ind(restsize))
             rest_ind = pack([ (jj, jj = 1,self%model%npara )], &
                  [( (any(b_ind(:) == jj)), jj = 1,self%model%npara)].eqv..false.)

             allocate(other_var(restsize,restsize), &
                  temp(bsize,restsize))

             other_var = variance(rest_ind, rest_ind)
             call inverse(other_var, info)

             call dgemm('n','n',bsize,restsize,restsize,1.0_wp, variance(b_ind, rest_ind), bsize, &
                   other_var, restsize, 0.0_wp, temp, bsize)

             call dgemm('n','n',bsize,bsize, restsize, -1.0_wp, temp, bsize, &
                  variance(rest_ind, b_ind), restsize, 1.0_wp, block_variance_chol, bsize)

             deallocate(other_var, temp, rest_ind)
          end if

          call cholesky(block_variance_chol, info)

          chol_var(b_ind, b_ind) = block_variance_chol
          deallocate(b_ind, block_variance_chol)
       end do


       nodeuind = 1
       nodeacpt = 0
       do j = 1, self%ngap

          do m = 1, self%nintmh

             npara_old = 0
             do bj = 1, self%nblocks
                bsize = break_points(bj+1) - break_points(bj)
                allocate(b_ind(bsize), block_variance_chol(bsize,bsize))
                b_ind = ind(break_points(bj)+1:break_points(bj+1))

                block_variance_chol = chol_var(b_ind, b_ind)
                ! block_variance_chol = variance(b_ind, b_ind)
                ! if (self%conditional_covariance .and. (self%nblocks>1)) then

                !    rest_ind = pack([ (jj, jj = 1,self%model%npara )], &
                !         [( (any(b_ind(:) == jj)), jj = 1,self%model%npara)].eqv..false.)
                !    restsize = self%model%npara-bsize
                !    allocate(other_var(restsize,restsize), &
                !         temp(bsize,restsize))

                !    other_var = variance(rest_ind, rest_ind)
                !    call inverse(other_var, info)

                !    call dgemm('n','n',bsize,restsize,restsize,1.0_wp, variance(b_ind, rest_ind), bsize, &
                !         other_var, restsize, 0.0_wp, temp, bsize)
                !    call dgemm('n','n',bsize,bsize,restsize, -1.0_wp, temp, bsize, &
                !         variance(rest_ind, b_ind), restsize, 1.0_wp, block_variance_chol, bsize)

                ! end if

                ! call cholesky(block_variance_chol, info)

                p0 = nodepara%particles(:,j)
                p0(b_ind) = p0(b_ind) + scale*matmul(block_variance_chol, nodeeps(b_ind,j))

                loglh0 = self%model%lik(p0, T=current_T)
                if ( nodepara%loglhold(j) /= 0.0_wp) then
                   loglhold0 = self%model%lik(p0, T=previous_T)
                else
                   loglhold0 = 0.0_wp
                end if
                prior0 = self%model%prior%logpdf(p0)

                likdiff = loglh0 - nodepara%loglh(j)
                likdiffold = loglhold0 - nodepara%loglhold(j)
                prdiff = prior0 - nodepara%prior(j)
                alp = exp( phi*(likdiff-likdiffold) + likdiffold + prdiff)

                if (nodeu(nodeuind,1) < alp) then
                   nodeacpt(b_ind,j) = nodeacpt(b_ind,j) + 1
                   nodepara%particles(:,j) = p0
                   nodepara%loglh(j) = loglh0
                   nodepara%loglhold(j) = loglhold0
                   nodepara%prior(j) = prior0
                end if

                deallocate(b_ind, block_variance_chol)

                npara_old = npara_old + bsize
                nodeuind = nodeuind + 1
             end do

          end do
       end do
       previous_T = self%temp%T_schedule(i)

       call mpi_gather(nodeacpt, self%ngap*self%model%npara, MPI_INTEGER, &
            acptsim, self%ngap*self%model%npara, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierror )
       call gather_particles(parasim, nodepara)

       if (rank == 0.0_wp) then

          ahat = sum(acptsim) / (1.0_wp * self%npart * self%nintmh * self%model%npara)
          scale = scale * (0.95_wp + 0.10*exp(16.0_wp*(ahat - 0.25_wp)) &
               / (1.0_wp + exp(16.0_wp*(ahat - 0.25_wp))))
          print*,'acceptance rate', ahat, scale
          print*,sum(log(self%temp%Z_estimates(1:i)))
          print*,sum(parasim%loglh)/self%npart, sum(parasim%prior)/self%npart
          if (phi == 1.0_wp) then
             write(chart, '(I3.3)') current_T
             call json%create_object(json_ip,'posterior.'//trim(chart))
             call json%add(json_p, json_ip)
             call parasim%write_json(json_ip)
             nullify(json_ip)
          end if
       end if

       call mpi_bcast(scale, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror)

       i = i + 1
    end do

    if (rank==0) then
       print*,'Writing Output'
       call self%temp%write_json(json_p)
       call json%print(json_p, self%output_file)
       call json%destroy(json_p)
    end if
  end subroutine estimate


  subroutine scatter_particles(parasim, nodepara)

    type(fortress_smc_particles), intent(inout) :: parasim, nodepara

    integer :: ngap, npara, mpierror

    ngap = nodepara%npart
    npara = nodepara%nvars

    call mpi_scatter(parasim%particles, ngap*npara, MPI_DOUBLE_PRECISION, &
         nodepara%particles, ngap*npara, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror)
    call mpi_scatter(parasim%weights, ngap, MPI_DOUBLE_PRECISION, nodepara%weights, &
         ngap, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror)
    call mpi_scatter(parasim%loglh, ngap, MPI_DOUBLE_PRECISION, nodepara%loglh, ngap, &
         MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror)
    call mpi_scatter(parasim%prior, ngap, MPI_DOUBLE_PRECISION, nodepara%prior, ngap, &
         MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror)
    call mpi_barrier(MPI_COMM_WORLD, mpierror)

  end subroutine scatter_particles

  subroutine gather_particles(parasim, nodepara)
    type(fortress_smc_particles), intent(inout) :: parasim, nodepara

    integer :: ngap, npara, mpierror

    ngap = nodepara%npart
    npara = nodepara%nvars

    call mpi_gather(nodepara%weights, ngap, MPI_DOUBLE_PRECISION, parasim%weights, &
        ngap, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror )
    call mpi_gather(nodepara%loglh, ngap, MPI_DOUBLE_PRECISION, parasim%loglh, &
         ngap, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror )
    call mpi_gather(nodepara%prior, ngap, MPI_DOUBLE_PRECISION, parasim%prior, &
         ngap, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror )
    call mpi_gather(nodepara%particles, ngap*npara, MPI_DOUBLE_PRECISION, parasim%particles, &
         ngap*npara, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror )

    call mpi_barrier(MPI_COMM_WORLD, mpierror)
   end subroutine gather_particles



  type(command_line_interface) function initialize_cli() result(cli)

    integer :: error

    call cli%init(progname='SMC Estimation via FORTRESS', &
         version=FORTRESS_VERSION, authors=FORTRESS_AUTHOR)
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
    call cli%add(switch='--nblocks', switch_ab='-nb', &
         help='Number of blocks', &
         required=.false.,act='store',def='1',error=error)
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
    call cli%add(switch='--output-file', switch_ab='-od', &
         required=.false.,act='store',def='output.json',error=error, help='The name of the output file')
    call cli%add(switch='--conditional-covariance', switch_ab='-cc', &
         required=.false.,act='store_true',def='.false.',error=error, &
         help='Conditional covariance.')
    call cli%add(switch='--endog-tempering', switch_ab='-et', &
         required=.false.,act='store_true',def='.false.',error=error, &
         help='Endogenous tempering.')
    call cli%add(switch='--resample-tol', switch_ab='-rt', &
         required=.false.,act='store',def='0.5',error=error, &
         help='Tolerance for resampling.')
    call cli%add(switch='--initial-particles', switch_ab='-pf', &
         required=.false.,act='store',def='none',error=error, &
         help='file with draws from prior')
    call cli%add(switch='--npriorextra', switch_ab='-pe', &
         help='Number of extra draws from prior (integer)', &
         required=.false.,act='store',def='1000', error=error)

  end function initialize_cli


  subroutine evaluate_time_t_lik(self, smc_particles, likT)

    class(fortress_smc) :: self
    type(fortress_smc_particles), intent(inout) :: smc_particles
    integer, intent(in) :: likT

    integer :: i,j, neval
    real(wp) :: lik0

    real(wp) :: paraextra(smc_particles%nvars, 1000)

    j = 1
    neval = smc_particles%npart

    do i = 1, neval

       lik0 = self%model%lik(smc_particles%particles(:,i), T=likT)


       if (isnan(lik0)) then
          lik0 = -10.0_wp**11
       endif

       smc_particles%loglh(i) = lik0
       ! if (mod(i, 100) == 0) then
       !    write(stdout,'(a$)') '+'
       ! endif

    end do

  end subroutine evaluate_time_t_lik


  subroutine draw_from_prior(self, smc_particles, likT)

    class(fortress_smc) :: self
    type(fortress_smc_particles), intent(inout) :: smc_particles
    integer, intent(in) :: likT

    integer :: i,j, neval
    real(wp) :: lik0, p0(self%model%npara)
    real(wp) :: paraextra(smc_particles%nvars, self%npriorextra)

    j = 1

    print*,'runfdsfas2', smc_particles%nvars
    print*, shape(paraextra)
    paraextra = self%model%prior%rvs(self%npriorextra)
    !print*,'fdsfsad'


    neval = smc_particles%npart
    do i = 1, neval

       p0 = smc_particles%particles(:,i)
       lik0 = self%model%lik(p0, T=likT)
       if (isnan(lik0)) lik0 = -10.0_wp**11


       if ((lik0 < -10.0_wp**9) .and. (self%init_file /= 'none')) then
          print*,'Error with initial particle files'
          stop
       end if

       do while (lik0 < -10.0_wp**9)
          p0 = paraextra(:,j)
          lik0 = self%model%lik(p0, T=likT)

          if (isnan(lik0)) lik0 = -10000000000000.0_wp

          j = j + 1
          if (j > self%npriorextra) then
             write(stderr, *) 'Prior is too far from likelihood ... '
             stop
          end if

       end do

       smc_particles%loglh(i) = lik0
       smc_particles%prior(i) = self%model%prior%logpdf(smc_particles%particles(:,i))
       if (mod(i, 100) == 0) then
          write(stdout,'(a$)') '+'
       endif

    end do

  end subroutine draw_from_prior


  subroutine rperm3(N, p, u)

    integer, intent(in) :: N
    integer, dimension(N), intent(out) :: p
    real(wp), intent(in) :: u(N)
    integer :: j, k


    p = 0

    do j=1,N

       !call random_number(u)
       k = floor(j*u(j)) + 1

       p(j) = p(k)
       p(k) = j

    end do

  end subroutine rperm3

  subroutine generate_random_blocks(npara, nblocks, indices, break_points, rvs)

    integer, intent(in) :: npara, nblocks

    integer, intent(inout) :: indices(npara)
    integer, intent(out) ::  break_points(nblocks+1)

    real(wp), intent(in), optional :: rvs(npara)

    real(wp) :: u(npara)
    integer :: ind2(npara)
    integer :: i, gap

    if (.not. (present(rvs))) then
       do i = 1,npara
          call random_number(u(i))
       end do
    else
       u = rvs
    end if

    call rperm3(npara, ind2, u)

    indices(ind2) = indices

    gap = int(1.0_wp*npara / nblocks)

    do i = 2, nblocks
       break_points(i) = (i-1)*gap

    end do

    break_points(1) = 0
    break_points(nblocks+1) = npara


  end subroutine generate_random_blocks


  real(wp) function ess_gap(phi, phi_old, parasim, r) result(f)

    real(wp), intent(in) :: phi, phi_old, r
    type(fortress_smc_particles) :: parasim

    !type(fortress_smc_particles) :: nw

    real(wp) :: ESS, temp, a1, a2
    real(wp) :: new_weight(parasim%npart), new_weight2(parasim%npart)
    real(wp) ::  max1, max2
    !nw = fortress_smc_particles(npart=parasim%npart)
    !nw%weights = exp( (phi - phi_old) * parasim%loglh) !* nw%weights
    new_weight = (phi-phi_old)*parasim%loglh
    new_weight2 = 2.0_wp*(phi-phi_old)*parasim%loglh

    max1 = maxval(new_weight)
    max2 = maxval(new_weight2)

    a1 = sum(exp(new_weight2-max2)) * parasim%npart
    a2 = sum(exp(new_weight-max1))**2

    f = exp(max2-2.0_wp*max1)*a2/a1  - r

  end function ess_gap


    double precision function bisection(lb1, ub1, tol, phi_old, parasim, rstar)

      double precision, intent(in) :: lb1, ub1, tol, phi_old, rstar
      type(fortress_smc_particles) :: parasim

      logical :: bisection_converged
      integer :: bisection_loops

      double precision :: x1, essx, lb, ub

      lb = lb1
      ub = ub1
      x1 = (lb+ub)/2
      essx =  ess_gap(x1, phi_old, parasim, rstar)
      bisection_converged = abs(essx) < tol
      bisection_loops = 1
      do while (.not. bisection_converged)

         if (essx < 0.0) then
            ub = x1
            x1 = (x1 + lb) / 2.0d0
         else
            lb = x1
            x1 = (x1 + ub) / 2.0d0
         endif

         essx =  ess_gap(x1, phi_old, parasim, rstar)

         bisection_converged = abs(essx) < tol
         bisection_loops = bisection_loops + 1



      end do
      !print*,bisection_loops
      bisection = x1

    end function bisection


end module fortress_smc_t
