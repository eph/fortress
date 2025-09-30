module fortress_particle_filter_t
  use, intrinsic :: iso_fortran_env, only: wp => real64
#:if defined('ENABLE_MPI')
  use mpi_f08, only: MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_MIN, MPI_SUM, MPI_Request, MPI_Status, &
       MPI_Allgather, MPI_Allreduce, MPI_Barrier, MPI_Irecv, MPI_Isend, MPI_Wait
#:endif
#:if defined('ENABLE_OPENMP')
  use omp_lib
#:endif

  use fortress_bayesian_model_t, only: fortress_ss_model
  use fortress_particles_t, only: fortress_particles
  use fortress_random_t, only: fortress_random

  implicit none

  type ParallelParticleFilter

     integer :: npart = 1000
     integer :: nproc = 1
     integer :: nlocalpart, nforeignpart, naltpart

     integer :: initialization_T = 1

     logical, allocatable :: adjusted_proposal(:)
     double precision, allocatable :: adjusted_proposal_mu(:,:), adjusted_proposal_std(:)

     class(fortress_ss_model), allocatable :: m
     type(fortress_random) :: rng

     logical :: USE_MPI = .true.

   contains

     procedure :: lik 

  end type ParallelParticleFilter

  interface ParallelParticleFilter
     module procedure new_ParallelParticleFilter
  end interface ParallelParticleFilter

contains

  type(ParallelParticleFilter) function new_ParallelParticleFilter & 
       (m, npart, seed, nproc, rank) result(ppf)

    class(fortress_ss_model) :: m

    integer, optional, intent(in) :: npart, nproc, seed, rank

    integer :: rng_seed, rank_opt
    
    allocate(ppf%m, source=m)

    if (present(npart)) ppf%npart = npart
    if (present(nproc)) ppf%nproc = nproc

    !if (present(nlocalpart)) ppf%nlocalpart = nlocalpart
    !if (present(nforeignpart)) ppf%nforeignpart = nforeignpart
    !if (present(naltpart)) ppf%naltpart = naltpart

#:if defined('ENABLE_MPI')
    ppf%nlocalpart = ppf%npart / max(ppf%nproc, 1)

    if (ppf%nproc > 1) then
       ppf%nforeignpart = ppf%nlocalpart / 1.5d0
       ppf%naltpart = ppf%nlocalpart / 2.0d0
       ppf%USE_MPI = .true.
    else
       ppf%nforeignpart = 0.0_wp
       ppf%naltpart = 0.0_wp
       ppf%USE_MPI = .false.
    end if
#:else
    ppf%nproc = 1
    ppf%nlocalpart = ppf%npart
    ppf%nforeignpart = 0.0_wp
    ppf%naltpart = 0.0_wp
    ppf%USE_MPI = .false.
#:endif
    
    if (rank==0) then
       print *, 'Initializing a parallel particle filter '
       print *, 'nlocalpart = ', ppf%nlocalpart, 'nproc = ', ppf%nproc
       print *, 'nforeignpart = ', ppf%nforeignpart, 'naltpart = ', ppf%naltpart
    end if

    if (present(seed)) then
       rng_seed = seed
    else
       rng_seed = 0
    end if

    ppf%rng = fortress_random(seed=rng_seed)

    allocate(ppf%adjusted_proposal(ppf%m%T), ppf%adjusted_proposal_mu(ppf%m%T, ppf%m%neps), &
         ppf%adjusted_proposal_std(ppf%m%T))

    ppf%adjusted_proposal = .false.
    ppf%adjusted_proposal_mu = 0.0d0
    ppf%adjusted_proposal_std = 1.0d0

  end function new_ParallelParticleFilter


  double precision function lik(ppf, para, rank, nproc, save_states) 
    class(ParallelParticleFilter) :: ppf

    double precision, intent(in) :: para(ppf%m%npara)

    integer, intent(in) :: rank
    integer, intent(in) :: nproc
    logical, optional, intent(in) :: save_states

    logical :: converged 


    type(fortress_particles) :: old_local, old_copy, old_foreign, new_local, alt_foreign, alt_copy
    double precision :: incwt, py, py_across_procs, relative_weight, py_vector(nproc)

    double precision :: ess, min_ess, randu(ppf%m%T,1), effective_procs
    integer :: foreignpartstart, ranked_procs(nproc), partner
    logical :: remaining_procs(nproc)
    double precision :: endogsteady(ppf%m%ns), shocks(ppf%m%neps, ppf%nlocalpart)
    double precision :: minwt, maxwt

    integer :: i, j, k, t

    double precision :: is_sampling_std, is_sampling_mu(ppf%m%neps), kap

    integer :: explode_count

    integer :: regime(2)
    integer :: left, right
    integer :: mpierror
#:if defined('ENABLE_MPI')
    type(MPI_Status) :: mpi_recv_status1, mpi_recv_status2, mpi_recv_status3, mpi_recv_status4
    type(MPI_Status) :: mpi_send_status1, mpi_send_status2, mpi_send_status3, mpi_send_status4

    type(MPI_Request) :: mpi_req_recv1
    type(MPI_Request) :: mpi_req_recv2
    type(MPI_Request) :: mpi_req_recv3
    type(MPI_Request) :: mpi_req_recv4

    type(MPI_Request) :: mpi_req_send1
    type(MPI_Request) :: mpi_req_send2
    type(MPI_Request) :: mpi_req_send3
    type(MPI_Request) :: mpi_req_send4
#:endif

    character(len=5) :: rb, char_rank, char_t
    character(len=200) :: out_str
    logical :: save_states_opt

    save_states_opt = .false.
    if (present(save_states)) save_states_opt = save_states

#:if defined('ENABLE_MPI')
    if (ppf%USE_MPI) then
       converged = ppf%m%solve(para, nproc, rank)
    else
       converged = ppf%m%solve(para)
    end if
#:else
    converged = ppf%m%solve(para)
#:endif
       
    if (rank==0) print*,'entering particle filter'
    lik = 0.0d0


#:if defined('ENABLE_MPI')
    if (ppf%USE_MPI) call MPI_Barrier(MPI_COMM_WORLD, mpierror)
#:endif

    if (converged .eqv. .false.) then
       lik = -1000000000000.0d0
       return 
    end if

    randu = ppf%rng%uniform_rvs(ppf%m%T, 1)


    right = modulo(rank+1, nproc)
    left = modulo(rank-1, nproc) 
    foreignpartstart = ppf%nlocalpart - ppf%nforeignpart

    ! initialize part
    old_local = fortress_particles(nvars=ppf%m%ns, npart=ppf%nlocalpart)
    new_local = fortress_particles(nvars=ppf%m%ns, npart=ppf%nlocalpart)

    old_copy = fortress_particles(nvars=ppf%m%ns, npart=ppf%nforeignpart)
    old_foreign = fortress_particles(nvars=ppf%m%ns, npart=ppf%nforeignpart)

    alt_copy = fortress_particles(nvars=ppf%m%ns, npart=ppf%naltpart)
    alt_foreign = fortress_particles(nvars=ppf%m%ns, npart=ppf%naltpart)

    endogsteady = ppf%m%steadystate(para)
    old_local%particles = 0.0d0

#:if defined('ENABLE_OPENMP')
    !$omp parallel do if (.not. ppf%USE_MPI) default(shared) private(j)
#:endif
    do j = 1, ppf%m%ns
       old_local%particles(j,:) = endogsteady(j)
    end do
#:if defined('ENABLE_OPENMP')
    !$omp end parallel do
#:endif


    !if (rank==0) call ppf%m%describe_params()
    !if (rank==0) print*, old_local%particles(:,1)
#:if defined('ENABLE_MPI')
    if (ppf%USE_MPI) call MPI_Barrier(MPI_COMM_WORLD, mpierror)
#:endif

    !if (rank==0) print*,'entering first initializaztion'
    ! first initialization
    do t = 1, ppf%initialization_T
       shocks = ppf%rng%norm_rvs(ppf%m%neps, ppf%nlocalpart)

       !new_local%particles(:,1) = ppf%m%policy_function(old_local%particles(:,ppf%nlocalpart), shocks(:,1), [1,1])
       j = 2
       explode_count = 1
       do while (j <= ppf%nlocalpart)

          !do j = 2, ppf%nlocalpart
          old_local%particles(:,j) = ppf%m%policy_function(old_local%particles(:,j-1), shocks(:,j))
          ! check for explosions
          if ( isnan(sum(old_local%particles(:,j))) ) then
             ! write(rb,'(I3.3)') rank
             ! open(rank+27,file='explode_'//trim(adjustl(rb))//'.txt',action='write')
             ! write(rank+27,'(100f)'),old_local%particles(:,j)
             ! write(rank+27,'(100f)') old_local%particles(:,j-1)
             ! write(rank+27,'(100f)') shocks(:,j)
             ! close(rank+27)
             ! stop
             j = max(j-3, 2)
             shocks = ppf%rng%norm_rvs(ppf%m%neps, ppf%nlocalpart)
             explode_count = explode_count + 1

             if (explode_count > 1) then
                print*,'particle',j,'explodes on rank',rank

                old_local%particles(:,j-1) = 0.0d0
                do k = 1, ppf%m%ns
                   old_local%particles(k,j-1) = endogsteady(k)
                end do


             end if
          else
             j = j + 1
             explode_count = 1
          end if

       end do
    end do


    ! second initialization
#:if defined('ENABLE_MPI')
    if (ppf%USE_MPI) call MPI_Barrier(MPI_COMM_WORLD, mpierror)
#:endif
    !if (rank==0) print*,'entering second initialization'

    do t = 1, ppf%initialization_T
       shocks = ppf%rng%norm_rvs(ppf%m%neps, ppf%nlocalpart)

       do j = 1, ppf%nlocalpart
          new_local%particles(:,j) = ppf%m%policy_function(old_local%particles(:,j), shocks(:,j))
       end do
       old_local = new_local
    end do

    do j = 1, ppf%nlocalpart
       do k = 1, ppf%m%ns
          if (isnan(old_local%particles(k,j))) then
             old_local%particles(:,j) = 0.0d0
             old_local%weights(j) = 0.0d0
          end if
       end do
    end do

    if (sum(old_local%weights) == 0.0d0) then
       print*,'divergence on rank ', rank
       lik = -1000000000.0
       stop
    end if
    old_copy = old_local 

#:if defined('ENABLE_MPI')
    if (ppf%USE_MPI) call MPI_Barrier(MPI_COMM_WORLD, mpierror)
#:endif

    old_copy%weights = 1.0d0 / (ppf%nlocalpart * nproc)
    old_local%weights = 1.0d0 / (ppf%nlocalpart * nproc)




    ! ! the particle filter
    !if (rank==0) print*,'entering loop'
    do t = 1, ppf%m%T


#:if defined('ENABLE_MPI')
       if (ppf%USE_MPI) then
          call MPI_Irecv(old_foreign%particles, ppf%m%ns*ppf%nforeignpart, MPI_DOUBLE_PRECISION, &
               left, 1, MPI_COMM_WORLD, mpi_req_recv1, mpierror)
          call MPI_Irecv(old_foreign%weights, ppf%nforeignpart, MPI_DOUBLE_PRECISION, &
               left, 1, MPI_COMM_WORLD, mpi_req_recv2, mpierror)
          call MPI_Isend(old_copy%particles(:,1:ppf%nforeignpart), ppf%m%ns*ppf%nforeignpart, MPI_DOUBLE_PRECISION, &
               right, 1, MPI_COMM_WORLD, mpi_req_send1, mpierror)
          call MPI_Isend(old_copy%weights(1:ppf%nforeignpart), ppf%nforeignpart, MPI_DOUBLE_PRECISION, &
               right, 1, MPI_COMM_WORLD, mpi_req_send2, mpierror)
       end if
#:endif
       ! draw shocks and regimes
       shocks = ppf%rng%norm_rvs(ppf%m%neps, ppf%nlocalpart)

       if (ppf%adjusted_proposal(t)) then
          is_sampling_std = ppf%adjusted_proposal_std(t)
          is_sampling_mu = ppf%adjusted_proposal_mu(t,:)

          ! should be speed up
          shocks = shocks * is_sampling_std
          do j = 1, ppf%m%neps
             shocks(j,:) = shocks(j,:) + is_sampling_mu(j)
          end do
       end if

       do j = 1, ppf%nlocalpart - ppf%nforeignpart
          new_local%particles(:,j) = ppf%m%policy_function(old_local%particles(:,ppf%nforeignpart+j), &
               shocks(:,j))
       end do

       !calculate unnormalized weights 
       do j = 1, ppf%nlocalpart - ppf%nforeignpart
          incwt = ppf%m%pdfy(t, new_local%particles(:,j), old_local%particles(:,ppf%nforeignpart+j), para)
          new_local%weights(j) = old_local%weights(ppf%nforeignpart+j) * incwt 
       end do

#:if defined('ENABLE_MPI')
       if (ppf%USE_MPI) call MPI_Wait(mpi_req_recv1, mpi_recv_status1, mpierror)
#:endif

       do j = 1, ppf%nforeignpart
          new_local%particles(:,foreignpartstart+j) = ppf%m%policy_function(old_foreign%particles(:,j), &
               shocks(:,foreignpartstart+j))
       end do

#:if defined('ENABLE_MPI')
       if (ppf%USE_MPI) call MPI_Wait(mpi_req_recv2, mpi_recv_status2, mpierror)
#:endif

       do j = 1, ppf%nforeignpart
          incwt = ppf%m%pdfy(t, new_local%particles(:, foreignpartstart+j), &
               old_foreign%particles(:,j), para)
          !if (incwt==0) print*,'***',j,rank,incwt,new_local%particles(:,foreignpartstart+j)
          new_local%weights(foreignpartstart+j) = old_foreign%weights(j) * incwt
       end do

#:if defined('ENABLE_MPI')
       if (ppf%USE_MPI) call MPI_Wait(mpi_req_send1, mpi_send_status1, mpierror)
       if (ppf%USE_MPI) call MPI_Wait(mpi_req_send2, mpi_send_status2, mpierror)
#:endif

       ! could be speed up
       if (ppf%adjusted_proposal(t)) then

#:if defined('ENABLE_OPENMP')
       !$omp parallel do if (.not. ppf%USE_MPI) default(shared) private(j)
#:endif
       do j = 1, ppf%nlocalpart
          kap = is_sampling_std ** ppf%m%neps &
               * exp(-0.5d0 * dot_product(shocks(:,j),shocks(:,j))) &
               / exp(-0.5d0 * dot_product(shocks(:,j)-is_sampling_mu, shocks(:,j)-is_sampling_mu) / is_sampling_std ** 2)
          new_local%weights(j) = new_local%weights(j) * kap
       end do
#:if defined('ENABLE_OPENMP')
       !$omp end parallel do
#:endif
       end if


       minwt = minval(new_local%weights)
       maxwt = maxval(new_local%weights)

       do j = 1,ppf%nlocalpart
          if (isnan(new_local%weights(j))) new_local%weights(j) = 0.000000000001d0
       end do

       ! check the particles
       call new_local%normalize_weights(py)
       ess = new_local%ess()
#:if defined('ENABLE_MPI')
       if (ppf%USE_MPI) then
          call MPI_Allreduce(ess, min_ess, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, mpierror)
       else
          min_ess = ess
       end if
#:else
       min_ess = ess
#:endif
       !new_local%weights = new_local%weights * relative_weight 

       ! check lnpy

       ! if (sumproby ==0.0d0) then
       !    write(*,*) 'divergence of particle filter on rank', rank, 'at time ', t
       ! end if

       ! if (sumproby == nlocalparticles*minweight_param) then
       !    write(*,*) 'divergence of particle filter on rank', rank, 'avoided at time ', t
       ! end if

       ! if (isnan(sumproby)==.TRUE.) then
       !    write(*,*) 'prob(y) is nan on rank', rank
       ! end if
       ! totsumproby = 0.0d0
       py_across_procs = 0.0d0
       if (ppf%USE_MPI) then
#:if defined('ENABLE_MPI')
          call MPI_Allreduce(py, py_across_procs, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
#:else
          py_across_procs = py
#:endif
       else
          py_across_procs = py
       end if
       relative_weight = py / py_across_procs

       lik = lik + log(py_across_procs)

       if (min_ess < ppf%nlocalpart) then
          !if (rank==0) print*,'min_ess',min_ess!, !ppf%nlocalpart, t
          call new_local%systematic_resampling(randu(t,1))
          new_local%weights = relative_weight / ppf%nlocalpart
       end if
       old_local = new_local


#:if defined('ENABLE_MPI')
       if (ppf%USE_MPI) then
          py = sum(old_local%weights)
          py_vector = 0.0d0
          call MPI_Allgather(py, 1, MPI_DOUBLE_PRECISION, py_vector, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpierror)
          effective_procs = 1.0d0 / sum(py_vector**2)

          !if (rank==0) print*,'before passing', t, effective_procs 
          if (effective_procs < (nproc/2.0d0)) then 
             ! sort processors by relative weight, partner 
             remaining_procs = .true.
             do i = 1, nproc
                ranked_procs(i) = minloc(py_vector, dim=1, mask=remaining_procs) - 1
                remaining_procs(ranked_procs(i)+1) = .false.
             end do
             partner = nproc - minloc(ranked_procs, dim=1, mask=ranked_procs==rank) + 1
             partner = ranked_procs(partner)
             !print*,'left', left, 'partner', partner

             ! pass the last nalt particlse
             call MPI_Barrier(MPI_COMM_WORLD, mpierror)

             alt_copy%particles = old_local%particles(:,ppf%nlocalpart-ppf%naltpart+1:ppf%nlocalpart)
             alt_copy%weights = old_local%weights(ppf%nlocalpart-ppf%naltpart+1:ppf%nlocalpart)

             call MPI_Irecv(alt_foreign%particles, ppf%m%ns*ppf%naltpart, MPI_DOUBLE_PRECISION, &
                  partner, 2, MPI_COMM_WORLD, mpi_req_recv3, mpierror)
             call MPI_Irecv(alt_foreign%weights, ppf%naltpart, MPI_DOUBLE_PRECISION, partner, &
                  2, MPI_COMM_WORLD, mpi_req_recv4, mpierror)

             call MPI_Isend(alt_copy%particles, ppf%m%ns*ppf%naltpart, MPI_DOUBLE_PRECISION, &
                  partner, 2, MPI_COMM_WORLD, mpi_req_send3, mpierror)
             call MPI_Isend(alt_copy%weights, ppf%naltpart, MPI_DOUBLE_PRECISION, partner, 2, &
                  MPI_COMM_WORLD, mpi_req_send4, mpierror)

             call MPI_Wait(mpi_req_send3, mpi_send_status3, mpierror)
             call MPI_Wait(mpi_req_send4, mpi_send_status4, mpierror)
             call MPI_Wait(mpi_req_recv3, mpi_recv_status3, mpierror)
             call MPI_Wait(mpi_req_recv4, mpi_recv_status4, mpierror)

             old_local%particles(:, ppf%nlocalpart-ppf%naltpart+1:ppf%nlocalpart) = alt_foreign%particles
             old_local%weights(ppf%nlocalpart-ppf%naltpart+1:ppf%nlocalpart) = alt_foreign%weights

            py = sum(old_local%weights)
            py_vector = 0.0d0
            call MPI_Allgather(py, 1, MPI_DOUBLE_PRECISION, py_vector, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpierror)
            effective_procs = 1.0d0 / sum(py_vector**2)
            !if (rank==0) print*, 'after passing', t, effective_procs 

         end if

         old_copy%particles = old_local%particles(:,1:ppf%nforeignpart)
         old_copy%weights = old_local%weights(1:ppf%nforeignpart)

       end if
#:endif
       if (.not. ppf%USE_MPI) effective_procs = 1.0d0
       if (save_states_opt) then
          write(char_rank,'(I3.3)') rank
          write(char_t, '(I3.3)') t
          out_str = 'test-output/t_'//trim(adjustl(char_t))//'_rank_'//trim(adjustl(char_rank))//'_'
          call old_copy%write(out_str)
       end if

#:if defined('ENABLE_MPI')
       call MPI_Barrier(MPI_COMM_WORLD, mpierror)
#:endif
       !print*,rank,py_vector
    end do

    ! deallocate
    call old_local%free()
    call old_copy%free()
    call old_foreign%free()
    call new_local%free()
    call alt_foreign%free()
    call alt_copy%free()

  end function lik
end module fortress_particle_filter_t
