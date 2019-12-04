module fortress_smc_t
  use, intrinsic :: iso_fortran_env, only: wp => real64, &
       stdout => output_unit, stderr => error_unit
  !use, intrinsic :: ieee_arithmetic, only: isnan => ieee_is_nan 

  use flap, only : command_line_interface
  use fortress_bayesian_model_t, only: fortress_abstract_bayesian_model, fortress_lgss_model
  use fortress_prior_t, only: logmvnormpdf
  use fortress_smc_particles_t, only: fortress_smc_particles
  use fortress_random_t, only: fortress_random
  use fortress_linalg, only: cholesky, inverse
  use fortress_util, only: read_array_from_file, write_array_to_file
  use fortress_constants, only: BAD_LOG_LIKELIHOOD
  use fortress_smc_utilities, only: generate_random_blocks, ess_gap, bisection
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
     integer :: write_every, T_write_thresh
     logical :: save_hyper, fixed_hyper
     real(wp) :: initial_c, mix_alpha
     logical :: mcmc_mix, conditional_covariance
     integer :: seed
     character(len=200) :: output_file, init_file

     character(len=200) :: mutation_type 
     integer :: L 

     integer :: nproc
     integer :: ngap

     type(fortress_random) :: rng
     type(tempering_schedule) :: temp

     logical :: endog_tempering
     logical :: resample_every_period
     
     real(wp) :: resample_tol
     real(wp) :: bisection_thresh
     real(wp) :: target_acpt = 0.25_wp
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

    call json%add(json_node, 'tempering_phi', self%phi_schedule(1:self%nstages))
    call json%add(json_node, 'tempering_T', self%T_schedule(1:self%nstages))
    call json%add(json_node, 'Z_estimates', self%Z_estimates(1:self%nstages))
    call json%add(json_node, 'ESS_estimates', self%ESS_estimates(1:self%nstages))

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
    call smc%cli%get(switch='-c',val=smc%initial_c,error=err); if (err/=0) stop 1
    call smc%cli%get(switch='-n',val=smc%npart,error=err); if (err/=0) stop 1
    call smc%cli%get(switch='-we',val=smc%write_every,error=err); if (err/=0) stop 1
    call smc%cli%get(switch='-sh',val=smc%save_hyper,error=err); if (err/=0) stop 1
    call smc%cli%get(switch='-fh',val=smc%fixed_hyper,error=err); if (err/=0) stop 1
    call smc%cli%get(switch='-od',val=smc%output_file,error=err); if (err /=0) stop 1
    call smc%cli%get(switch='-s',val=smc%seed,error=err); if (err /=0) stop 1
    call smc%cli%get(switch='-nb',val=smc%nblocks,error=err); if (err /=0) stop 1
    call smc%cli%get(switch='-cc',val=smc%conditional_covariance,error=err); if (err /=0) stop 1
    call smc%cli%get(switch='-rep',val=smc%resample_every_period,error=err); if (err /=0) stop 1
    call smc%cli%get(switch='-rt',val=smc%resample_tol,error=err); if (err/=0) stop 1
    call smc%cli%get(switch='-et',val=smc%endog_tempering,error=err); if (err/=0) stop 1
    call smc%cli%get(switch='-L',val=smc%L,error=err); if (err/=0) stop 1
    call smc%cli%get(switch='-bt',val=smc%bisection_thresh,error=err); if (err/=0) stop 1
    call smc%cli%get(switch='-mt',val=smc%mutation_type,error=err); if (err/=0) stop 1
    call smc%cli%get(switch='-pf',val=smc%init_file,error=err); if(err/=0) stop 1
    call smc%cli%get(switch='-pe',val=smc%npriorextra,error=err); if(err/=0) stop 1
    call smc%cli%get(switch='-twt',val=smc%T_write_thresh,error=err); if(err/=0) stop 1
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

    if (smc%T_write_thresh == -1) smc%T_write_thresh = smc%model%T

    if (smc%endog_tempering) then
       print*,'Using endogenous tempering (on full sample only).'
       smc%temp = tempering_schedule(nstages=100000, lambda=0.0_wp, max_T=32)
       select type(mod => smc%model )
       class is (fortress_lgss_model)
          smc%temp%T_schedule = mod%T
       class default
          smc%temp%T_schedule = mod%T
       end select
       smc%temp%phi_schedule = 0.0_wp
    end if

    if (smc%mutation_type(1:3)=="HMC") then
       smc%target_acpt = 0.65_wp
    endif

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

    real(wp) :: nodeeps(self%model%npara ,self%ngap*self%nintmh)
    real(wp) :: nodeu(self%ngap*self%nblocks * self%nintmh, 1)
    integer :: nodeacpt(self%model%npara,self%ngap)
    real(wp) :: nodeloglh(self%ngap), nodepriorsim(self%ngap)

    logical :: add_new_observation
    real(wp) :: phi, phi_old
    integer :: previous_T, current_T

    real(wp) :: mean(self%model%npara), variance(self%model%npara, self%model%npara), chol_var(self%model%npara, self%model%npara)

    real(wp) :: p0(self%model%npara), loglh0, loglhold0, prior0
    real(wp) :: likdiff, likdiffold, prdiff, alp
    integer :: nodeuind, nodemind

    real(wp) :: loglhvec(self%model%T)

    character(len=3) :: chart

    integer :: npara_old, bj, jj, bsize, break_points(self%nblocks+1), ind(self%model%npara), restsize
    integer, allocatable :: b_ind(:), rest_ind(:)
    real(wp) :: variance_copy(self%model%npara, self%model%npara)
    real(wp), allocatable :: block_variance_chol(:,:), temp(:,:), other_var(:,:)

    real(wp) :: scale, ahat

    real(wp) :: ess_gap1, phi0

    real(wp) :: momemtum(self%model%npara), zeros_vec(self%model%npara), momemtum_dens_diff, scale_rand, momemtum0(self%model%npara)
    integer :: i_hmc

    real(wp) :: old_particles(self%model%npara, self%npart), mu_old, mu_new, sig_old, sig_new, average_corr(self%model%npara)

    real(wp) :: incwt(self%npart), maxincwt, Zt

    zeros_vec = 0.0_wp

    scale = self%initial_c
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
                   phi0 = max(phi0 , phi_old+0.01_wp)
                   if (phi0 < phi_old) stop
                   ess_gap1 = ess_gap(phi0, phi_old, parasim, self%resample_tol)
                   print*,phi0,ess_gap1
                end do
                phi = bisection(phi_old, phi0, self%bisection_thresh, phi_old, parasim, self%resample_tol)
             end if
             ess_gap1 = ess_gap(phi, phi_old, parasim, self%resample_tol)
             self%temp%phi_schedule(i) = phi

          end if

          !------------------------------------------------------------
          ! Correction
          !------------------------------------------------------------
          ! print*,phi - phi_old

          ! do j = 1, parasim%npart
          !    print*,parasim%weights(j),parasim%loglh(j),parasim%loglhold(j)
          !    parasim%weights(j) = parasim%weights(j) * exp( &
          !         (phi - phi_old)*(parasim%loglh(j)-parasim%loglhold(j)) )

          !    if (isnan(parasim%weights(j))) then
          !       print*,'Particle ', j, 'is nan'
          !       print *,parasim%particles(:,j)
          !       parasim%weights(j) = 0.0_wp
          !    end if
          ! end do
          ! print*,parasim%weights
          ! call parasim%normalize_weights(self%temp%Z_estimates(i))
          incwt = (phi - phi_old)*(parasim%loglh-parasim%loglhold)
          maxincwt = maxval(incwt)
          Zt = sum(exp(incwt - maxincwt)*parasim%weights)
          parasim%weights = exp(incwt - maxincwt)*parasim%weights / Zt
          self%temp%Z_estimates = log(Zt) + maxincwt

          self%temp%ESS_estimates(i) = parasim%ESS()
          if (isnan(self%temp%ESS_estimates(i) )) stop
          print*,'============================================================='
          write(*,'(A,I8,A,I9)') 'Iteration', i,' of ', self%temp%nstages
          write(*,'(A,I4,A,F8.5,A)') 'Current  (T,phi): (', current_T, ',', phi, ')'
          write(*,'(A,I4,A,F8.5,A)') 'Previous (T,phi): (', previous_T, ',', phi_old, ')'
          write(*,'(A,F8.2, A, I8.0)',advance='no')        'Current ESS     :', self%temp%ESS_estimates(i), ' /', self%npart
          !------------------------------------------------------------
          ! Selection
          !------------------------------------------------------------
          !print*,self%temp%ESS_estimates(i), ahat, scale, phi, current_T
          if ((self%temp%ESS_estimates(i) < self%npart * 0.5_wp) .or. &
              (self%resample_every_period .eqv. .true.))then 
             !.or. (self%endog_tempering)) then
             write(*,'(A)') ' [resampling]'
             call parasim%systematic_resampling(uu(i,1), resample_ind)

             parasim%prior = parasim%prior(resample_ind)
             parasim%loglh = parasim%loglh(resample_ind)
             parasim%loglhold = parasim%loglhold(resample_ind)

          else
             print*, ''
          end if

          eps = self%rng%norm_rvs(self%model%npara, self%npart*self%nintmh)
          u = self%rng%uniform_rvs(self%nblocks*self%nintmh*self%npart, 1)
          call parasim%mean_and_variance(mean, variance)
          write(*,'(A)') '             mu    std'
          write(*,'(A)') '           -----  -----'
          do bj = 1, self%model%npara
             write(*, '(A,I3,A,F7.3,F7.3)') 'para[',bj,']',mean(bj), sqrt(variance(bj,bj))
          end do
          ! blocking
          ind = [ (bj, bj = 1,self%model%npara )]
          call generate_random_blocks(self%model%npara, self%nblocks, ind, break_points)
          old_particles = parasim%particles

       end if

       ! for endogenous tempering wiht multiprocessing
       call mpi_barrier(MPI_COMM_WORLD, mpierror)
       call mpi_bcast(self%temp%nstages, 1, MPI_DOUBLE_PRECISION, &
            0, MPI_COMM_WORLD, mpierror)
       call mpi_bcast(phi, 1, MPI_DOUBLE_PRECISION, &
            0, MPI_COMM_WORLD, mpierror)
       self%temp%phi_schedule(i) = phi

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
       !call generate_subblock_cholesky(self%npara, self%nblocks, break_points, variance, chol_var)
       
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

       if (self%mutation_type(1:3) == "HMC") then
 
          if (self%mutation_type == 'HMC-diagM') then 
             chol_var = 0.0_wp
             do j = 1,self%model%npara
                chol_var(j,j) = variance(j,j) 
             end do
             variance = chol_var
          elseif (self%mutation_type == 'HMC-M') then
             chol_var = variance
          elseif (self%mutation_type == 'HMC-I') then
             chol_var = 0.0_wp
             variance = 0.0_wp
             do j = 1,self%model%npara
                chol_var(j,j) = 1.0_wp
                variance(j,j) = 1.0_wp
             end do
          else
             chol_var = 0.0_wp
             do j = 1,self%model%npara
                chol_var(j,j) = variance(j,j) 
             end do
             variance = chol_var
          end if

          call inverse(chol_var, info)
          call cholesky(chol_var, info)

       end if
          


       nodeuind = 1
       nodemind = 1
       nodeacpt = 0

       do j = 1, self%ngap
          
          do m = 1, self%nintmh

             npara_old = 0
             do bj = 1, self%nblocks
                bsize = break_points(bj+1) - break_points(bj)
                allocate(b_ind(bsize), block_variance_chol(bsize,bsize))
                b_ind = ind(break_points(bj)+1:break_points(bj+1))

                block_variance_chol = chol_var(b_ind, b_ind) 

                p0 = nodepara%particles(:,j)

                if (self%mutation_type == "RWMH") then
                   p0(b_ind) = p0(b_ind) + scale*matmul(block_variance_chol, nodeeps(b_ind,nodemind))
                 elseif (self%mutation_type(1:3) == "HMC") then 
                    call dgemv('N', self%model%npara, self%model%npara, 1.0_wp, chol_var, self%model%npara, nodeeps(b_ind,nodemind), 1, 0.0_wp, momemtum0, 1)

                    momemtum = momemtum0 + scale*(phi*self%model%dlik(p0) + self%model%prior%dlogpdf(p0))/2.0_wp
                    do i_hmc = 1, self%L
                       call dgemv('N', self%model%npara, self%model%npara, scale, variance, self%model%npara, momemtum, 1, 1.0_wp, p0, 1)
                       if (i_hmc < self%L) momemtum = momemtum + scale*(phi*self%model%dlik(p0) + self%model%prior%dlogpdf(p0))
                       
                    end do
                    momemtum = momemtum + scale*(phi*self%model%dlik(p0) + self%model%prior%dlogpdf(p0))/2.0_wp
                    momemtum = -momemtum
                    momemtum_dens_diff = logmvnormpdf(momemtum, zeros_vec, chol_var) &
                         - logmvnormpdf(momemtum0, zeros_vec, chol_var)
                end if


                ! if ( self%endog_tempering .eqv. .true.) then !thennodepara%loglhold(j) /= 0.0_wp) then
                !    select type(mod => self%model )
                !    class is (fortress_lgss_model)
                !       loglhvec(1:current_T) = mod%lik_filter_vec(p0, T=current_T)
                !       loglh0 = sum(loglhvec(1:current_T))
                !       loglhold0 = sum(loglhvec(1:maxval([current_T-1,0])))

                !       if ((isnan(loglh0))) loglh0 = BAD_LOG_LIKELIHOOD
                !       if ((isnan(loglhold0))) loglhold0 = BAD_LOG_LIKELIHOOD
                !    class default
                !       loglh0 = mod%lik(p0, T=current_T)
                !       loglhold0 = 0.0_wp !mod%lik(p0, T=maxval([current_T-1,0]))
                !    end select
                ! else
                   loglh0 = self%model%lik(p0, T=current_T)
                   loglhold0 = 0.0_wp
                !end if

                prior0 = self%model%prior%logpdf(p0)

                likdiff = loglh0 - nodepara%loglh(j)


                likdiffold = loglhold0 - nodepara%loglhold(j)
                prdiff = prior0 - nodepara%prior(j)
                alp = exp( phi*(likdiff-likdiffold) + likdiffold + prdiff)

                if (self%mutation_type(1:3) == "HMC")  alp = exp( phi*(likdiff-likdiffold) + likdiffold + prdiff + momemtum_dens_diff)
                if (.not.(self%model%inbounds(p0))) alp = 0.0_wp
                if (isnan(alp)) alp = 0.0_wp

                if (nodeu(nodeuind,1) < alp) then
                   !print*,loglh0,loglhold0,p0,alp
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
             nodemind = nodemind + 1

          end do


       end do


       call mpi_gather(nodeacpt, self%ngap*self%model%npara, MPI_INTEGER, &
            acptsim, self%ngap*self%model%npara, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierror )
       call gather_particles(parasim, nodepara)

       if (rank == 0) then
          ahat = sum(acptsim) / (1.0_wp * self%npart * self%nintmh * self%model%npara)
          !scale = scale * (0.95_wp + 0.10*exp(16.0_wp*(ahat - self%target_acpt)) &
          !     / (1.0_wp + exp(16.0_wp*(ahat - self%target_acpt))))
          scale = scale * (0.80_wp + 0.40*exp(16.0_wp*(ahat - self%target_acpt)) &
               / (1.0_wp + exp(16.0_wp*(ahat - self%target_acpt))))
          print*,(0.80_wp + 0.40*exp(16.0_wp*(ahat - self%target_acpt)) &
               / (1.0_wp + exp(16.0_wp*(ahat - self%target_acpt))))
          write(*,'(A,F4.2,A,F8.4,A)') 'MCMC average acceptance: ', ahat, ' [c = ', scale, ']'
          write(*,'(A,F 10.2)') 'Log MDD estimate:', sum(self%temp%Z_estimates(1:i))
          write(*,'(A,F 10.2)') 'Avg. Log Lik    :', sum(parasim%loglh * parasim%weights)
          write(*,'(A,F 10.2)') 'Avg. Log Prior  :', sum(parasim%prior * parasim%weights)

          if ((phi == 1.0_wp) .and. (current_T >= self%T_write_thresh) )then
             write(chart, '(I3.3)') current_T
             call json%create_object(json_ip,'posterior.'//trim(chart))
             call json%add(json_p, json_ip)
             call parasim%write_json(json_ip)
             nullify(json_ip)
          end if

          do j = 1, self%model%npara
          
             mu_old = sum(old_particles(j,:))/self%npart
             sig_old = sum( (old_particles(j,:)-mu_old)**2)/self%npart
          
             mu_new = sum(parasim%particles(j,:))/self%npart
             sig_new = sum( (parasim%particles(j,:)-mu_new)**2)/self%npart
          
             average_corr(j) = (sum(old_particles(j,:) * parasim%particles(j,:))/self%npart - mu_old*mu_new) / sqrt(sig_new * sig_old)
          end do
          
          print*,minval(average_corr), maxval(average_corr)
          
          ! if ((maxval(average_corr) > 0.75_wp)) then
          !    self%L = floor(self%L * 1.25_wp)
          !    !scale = scale * 1.1_wp
          ! elseif (maxval(average_corr) < 0.25_wp) then 
          !    self%L = floor(self%L * 0.75_wp)+1
          ! end if
          ! self%L = min(self%L, 100)
          ! print*,self%L

       end if

       if ((self%endog_tempering) .and. (phi == 1.0_wp)) then
          self%temp%T_schedule(i+1:self%temp%nstages) = current_T + 1
       end if
       call mpi_bcast(scale, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror)
       !call mpi_bcast(self%L, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierror)
       call mpi_barrier(MPI_COMM_WORLD, mpierror)

       i = i + 1
    end do

    if (rank==0) then
       print*,''
       print*,'Finished SMC sampler! Writing output to '//self%output_file
       call self%temp%write_json(json_p)
       call json%print(json_p, self%output_file)
       call json%destroy(json_p)
    end if

    call mpi_barrier(MPI_COMM_WORLD, mpierror)
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
    call cli%add(switch='--mutation-type', switch_ab='-mt', help='Mutation Type', &
         required=.false.,act='store',def='RWMH', error=error)
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
    call cli%add(switch='--resample-every-period', switch_ab='-rep', &
         required=.false.,act='store_true',def='.false.',error=error, &
         help='Resample every period')
    call cli%add(switch='--endog-tempering', switch_ab='-et', &
         required=.false.,act='store_true',def='.false.',error=error, &
         help='Endogenous tempering.')
    call cli%add(switch='--resample-tol', switch_ab='-rt', &
         required=.false.,act='store',def='0.5',error=error, &
         help='Tolerance for resampling.')
    call cli%add(switch='--bisection-thresh', switch_ab='-bt', &
         required=.false.,act='store',def='0.001',error=error, &
         help='Tolerance for resampling.')
    call cli%add(switch='--initial-particles', switch_ab='-pf', &
         required=.false.,act='store',def='none',error=error, &
         help='file with draws from prior')
    call cli%add(switch='--npriorextra', switch_ab='-pe', &
         help='Number of extra draws from prior (integer)', &
         required=.false.,act='store',def='1000', error=error)
    call cli%add(switch='--HMCL', switch_ab='-L', &
         help='Steps for hmc', &
         required=.false.,act='store',def='10', error=error)
    call cli%add(switch='--T-write-thresh', switch_ab='-twt', &
         help='T to start writing posterior at', &
         required=.false.,act='store',def='-1', error=error)



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


    paraextra = self%model%prior%rvs(self%npriorextra, rng=self%rng)


    neval = smc_particles%npart
    do i = 1, neval

       p0 = smc_particles%particles(:,i)
       lik0 = self%model%lik(p0, T=likT)
       if (isnan(lik0)) lik0 = BAD_LOG_LIKELIHOOD


       if ((lik0 <= BAD_LOG_LIKELIHOOD ) .and. (self%init_file /= 'none')) then
          print*,'Error with initial particle files'
          stop
       end if

       do while (lik0 <= BAD_LOG_LIKELIHOOD)
          p0 = paraextra(:,j)
          lik0 = self%model%lik(p0, T=likT)

          if (isnan(lik0)) lik0 = BAD_LOG_LIKELIHOOD

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




end module fortress_smc_t
