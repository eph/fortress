!*******************************************************************************
!>author: Ed Herbst
!
!# Description
!  Abstract class for Bayesian models in Fortran.
!
module fortress_bayesian_model_t

  use, intrinsic :: iso_fortran_env, only: wp => real64, stderr => error_unit
  use fortress_prior_t, only: fortress_abstract_prior, model_prior => prior
  use fortress_util, only: read_array_from_file
  use filter, only: kalman_filter, chand_recursion, kalman_filter_missing, REALLY_NEG


  implicit none
 
  real(wp), parameter :: M_PI = 3.14159265358979323846d0


  type,public,abstract :: fortress_abstract_bayesian_model

     character(len=:), allocatable :: name
     character(len=:), allocatable :: datafile

     logical :: MISSING_DATA = .false.
     integer :: npara
     integer :: nobs, T

     class(fortress_abstract_prior), allocatable :: prior
     real(wp), allocatable :: yy(:,:)
     real(wp), allocatable :: p0(:)

   contains

     procedure :: read_data
     procedure(lik_func), deferred :: lik
     procedure :: dlik
     procedure :: inbounds
     !final :: finalize
     procedure :: construct_abstract_bayesian_model
     generic :: construct_model => construct_abstract_bayesian_model

  end type fortress_abstract_bayesian_model

  abstract interface
     function lik_func(self, para, T) result(l)
       use, intrinsic :: iso_fortran_env, only: wp => real64
       import :: fortress_abstract_bayesian_model
       implicit none
       class(fortress_abstract_bayesian_model), intent(inout) :: self
       real(wp), intent(in) :: para(self%npara)
       integer, intent(in), optional :: T
       real(wp) :: l
       
     end function lik_func
  end interface


  type,public,extends(fortress_abstract_bayesian_model), abstract :: fortress_ss_model

     integer :: ns, neps
     integer :: t0 = 0

     real(wp), allocatable :: HH(:,:)

   contains

     procedure(policy_function_i), public, deferred :: policy_function
     procedure(pdfy_i), public, deferred :: pdfy
     procedure(logpdfy_kernel_i), public, deferred :: logpdfy_kernel
     procedure :: steadystate

     procedure(solve_serial_i), public, deferred :: solve_serial
     procedure(solve_parallel_i), public, deferred :: solve_parallel
     generic :: solve => solve_serial, solve_parallel

  end type fortress_ss_model

  abstract interface
     function policy_function_i(self, states_old, shocks_new) result(states_new)
       use, intrinsic :: iso_fortran_env, only: wp => real64
       import :: fortress_ss_model
       implicit none
       class(fortress_ss_model), intent(inout) :: self
       real(wp), intent(in) :: states_old(self%ns), shocks_new(self%neps)
       real(wp) :: states_new(self%ns)
     end function policy_function_i

     real(wp) function pdfy_i(m, t, states_new, states_old, para) result(pdf)
       use, intrinsic :: iso_fortran_env, only: wp => real64
       import
       implicit none
       class(fortress_ss_model), intent(inout) :: m
       integer, intent(in) :: t
       real(wp), intent(in) :: states_new(m%ns), states_old(m%ns)
       real(wp), intent(in) :: para(m%npara)
     end function pdfy_i

     real(wp) function logpdfy_kernel_i(m, t, states_new, states_old, para) result(logpdf)
       use, intrinsic :: iso_fortran_env, only: wp => real64
       import
       implicit none
       class(fortress_ss_model), intent(inout) :: m
       integer, intent(in) :: t
       real(wp), intent(in) :: states_new(m%ns), states_old(m%ns)
       real(wp), intent(in) :: para(m%npara)
       double precision :: z
       integer :: i
     end function logpdfy_kernel_i

     logical function solve_serial_i(self, para) result(converged)
       use, intrinsic :: iso_fortran_env, only: wp => real64
       import
       implicit none
       class(fortress_ss_model), intent(inout) :: self
       real(wp), intent(in) :: para(self%npara)

     end function solve_serial_i

     logical function solve_parallel_i(self, para, nproc, rank) result(converged)
       use, intrinsic :: iso_fortran_env, only: wp => real64
       import
       implicit none
       class(fortress_ss_model), intent(inout) :: self
       real(wp), intent(in) :: para(self%npara)
       integer, intent(in) :: nproc, rank

     end function solve_parallel_i


  end interface

  type,public,extends(fortress_ss_model), abstract :: fortress_lgss_model

     logical :: USE_CR = .true.
     logical :: NONSTATIONARY = .false.

     character(len=:), allocatable :: priorfile

     real(wp), allocatable :: TT(:,:), RR(:,:), QQ(:,:)
     real(wp), allocatable :: DD(:), ZZ(:,:)

   contains

     procedure :: lik => lik_filter
     procedure :: lik_filter_vec
     procedure :: policy_function => policy_function_lgss
     procedure :: pdfy => pdfy_lgss
     procedure :: logpdfy_kernel => logpdfy_kernel_lgss
     procedure :: construct_lgss_model
     procedure :: construct_lgss_model_noprior
     procedure :: construct_lgss_model_noprior_nodata
     procedure(system_matrices_i), public, deferred :: system_matrices
     generic :: construct_model => construct_lgss_model

     procedure :: solve_serial => solve_serial_lgss
     procedure :: solve_parallel => solve_parallel_lgss

  end type fortress_lgss_model


  abstract interface
     subroutine system_matrices_i(self, para, error)
       use, intrinsic :: iso_fortran_env, only: wp => real64
       import 
       implicit none
       class(fortress_lgss_model), intent(inout) :: self
       real(wp), intent(in) :: para(self%npara)
       integer, intent(out) :: error
     end subroutine system_matrices_i

  end interface

  !type,public,extends(fortress_lgss_model) :: fortress_lgss_model_with_sv
  !   real(wp), allocatable :: time_varying_QQ(:,:,:)
  !
  ! contains
  !   procedure :: lik_conditional_on_sv
  !end type fortress_lgss_model_cr


contains 

    
  !*******************************************************************************
  !>
  !  This routine reads the data file (after allocating yy attribute.)
  !  Note that yy is stored in a column-major fashion, while the textfile read in
  !  is stored rowwise.
  !
  subroutine read_data(self)

    class(fortress_abstract_bayesian_model), intent(inout) :: self

    integer :: i

    allocate(self%yy(self%nobs,self%T))

    call read_array_from_file(self%datafile, self%yy, transpose=.true.)

    do i = 1, self%T
       if ( any(isnan(self%yy(:,i))) ) self%MISSING_DATA = .true.
    end do

  end subroutine read_data
  !*******************************************************************************  


  subroutine construct_abstract_bayesian_model(self, name, datafile, npara, nobs, T)

    class(fortress_abstract_bayesian_model), intent(inout) :: self

    character(len=144), intent(in) :: name, datafile
    integer, intent(in) :: nobs, T, npara

    self%name = name
    self%datafile = datafile

    self%nobs = nobs
    self%T = T
    self%npara = npara
    allocate(self%p0(npara))

    call self%read_data()

  end subroutine construct_abstract_bayesian_model

  function steadystate(self, para) result(ss)

    class(fortress_ss_model) :: self
    real(wp), intent(in) :: para(self%npara)
    real(wp) :: ss(self%ns)

    ss = 0.0_wp

  end function steadystate

  subroutine construct_lgss_model(self, name, datafile, priorfile, npara, nobs, T, ns, neps)

    class(fortress_lgss_model), intent(inout) :: self

    character(len=144), intent(in) :: name, datafile, priorfile
    integer, intent(in) :: nobs, T, ns, neps, npara

    allocate(self%prior, source=model_prior(priorfile))
    call self%construct_model(name, datafile, npara, nobs, T)
    
    self%ns = ns
    self%neps = neps

    allocate(self%TT(self%ns, self%ns), self%RR(self%ns, self%neps), self%QQ(self%neps, self%neps), &
         self%ZZ(self%nobs, self%ns), self%DD(self%nobs), self%HH(self%nobs, self%nobs))

  end subroutine construct_lgss_model


  subroutine construct_lgss_model_noprior(self, name, datafile, npara, nobs, T, ns, neps)

    class(fortress_lgss_model), intent(inout) :: self

    character(len=144), intent(in) :: name, datafile
    integer, intent(in) :: nobs, T, ns, neps, npara

    !allocate(self%prior, source=model_prior(priorfile))
    call self%construct_model(name, datafile, npara, nobs, T)
    
    self%ns = ns
    self%neps = neps

    allocate(self%TT(self%ns, self%ns), self%RR(self%ns, self%neps), self%QQ(self%neps, self%neps), &
         self%ZZ(self%nobs, self%ns), self%DD(self%nobs), self%HH(self%nobs, self%nobs))

  end subroutine construct_lgss_model_noprior


  subroutine construct_lgss_model_noprior_nodata(self, name, npara, nobs, T, ns, neps)

    class(fortress_lgss_model), intent(inout) :: self

    character(len=144), intent(in) :: name
    integer, intent(in) :: nobs, T, ns, neps, npara

    ! Initialize model parameters without reading prior or data files
    self%name = name
    self%datafile = ''  ! No datafile needed

    self%nobs = nobs
    self%T = T
    self%npara = npara
    self%ns = ns
    self%neps = neps

    allocate(self%p0(npara))

    ! Note: self%yy and self%prior are allocated/initialized by the caller
    ! This allows for hardcoded data and custom prior types

    allocate(self%TT(self%ns, self%ns), self%RR(self%ns, self%neps), self%QQ(self%neps, self%neps), &
         self%ZZ(self%nobs, self%ns), self%DD(self%nobs), self%HH(self%nobs, self%nobs))

  end subroutine construct_lgss_model_noprior_nodata


  !*******************************************************************************
  !>
  !
  function lik_filter_vec(self, para, T) result(l)

    class(fortress_lgss_model), intent(inout) :: self
    real(wp), intent(in) :: para(self%npara)
    integer, intent(in) :: T
    real(wp) :: l(T)
    integer :: error


    l = 0.0_wp
    if (T==0) return
    l(1) = REALLY_NEG

    if ( (T < 0) .or. (T>self%T)) then
       write(stderr,'(a,i0,a,i0)') 'ERROR: Invalid time period T=', T, ' (must be between 0 and ', self%T, ')'
       stop 1
    end if

    if (.not. self%inbounds(para)) return


    call self%system_matrices(para, error)
    if (error /= 0) return


    associate(yy => self%yy(:,1:T), TT => self%TT, RR => self%RR, QQ => self%QQ, &
         DD => self%DD, ZZ => self%ZZ, HH => self%HH, ny => self%nobs, &
         ns => self%ns, neps => self%neps, t0 => self%t0)

    if (self%MISSING_DATA) then
       l = kalman_filter_missing(yy,TT,RR,QQ,DD,ZZ,HH,ny,T,neps,ns,t0)
    elseif (self%USE_CR) then
       l = chand_recursion(yy,TT,RR,QQ,DD,ZZ,HH,ny,T,neps,ns,t0)
    else
       l = kalman_filter(yy,TT,RR,QQ,DD,ZZ,HH,ny,T,neps,ns,t0)
    end if
  end associate
end function lik_filter_vec

real(wp) function lik_filter(self, para, T) result(l)

  class(fortress_lgss_model), intent(inout) :: self
  real(wp), intent(in) :: para(self%npara)
  integer, intent(in), optional :: T
  integer :: use_T
  
  use_T = self%T
  if (present(T)) use_T = T

  l = sum(self%lik_filter_vec(para, T=use_T))

  if ((isnan(l))) l = REALLY_NEG


end function

function dlik(self, para, T) result(dl)

  class(fortress_abstract_bayesian_model), intent(inout) :: self
  real(wp), intent(in) :: para(self%npara)
  integer, intent(in), optional :: T

  real(wp) :: dl(self%npara)

  write(stderr,'(a)') 'ERROR: dlik() function not implemented for this model type'
  stop 1


end function dlik



double precision function pdfy_lgss(m, t, states_new, states_old, para) result(pdf)
    class(fortress_lgss_model), intent(inout) :: m 

    integer, intent(in) :: t

    double precision, intent(in) :: states_new(m%ns), states_old(m%ns)
    double precision, intent(in) :: para(m%npara)

    double precision :: z
    integer :: i

    pdf = 1.0d0
    do i = 1, m%nobs
       z = (m%yy(i,t) - m%DD(i) - dot_product(m%ZZ(i,:), states_new)) / sqrt(m%HH(i,i))
       pdf = pdf / sqrt(m%HH(i,i)) * exp(-0.5d0 * z**2)
    end do

    pdf = pdf * (1.0d0 / sqrt(2.0d0*M_PI))**m%nobs
  end function pdfy_lgss


  real(wp) function logpdfy_kernel_lgss(m, t, states_new, states_old, para) result(logpdf)
    class(fortress_lgss_model), intent(inout) :: m 

    integer, intent(in) :: t

    real(wp), intent(in) :: states_new(m%ns), states_old(m%ns)
    real(wp), intent(in) :: para(m%npara)

    double precision :: z
    integer :: i

    logpdf = 0.0d0
    do i = 1, m%nobs
       z = (m%yy(i,t) - m%DD(i) - dot_product(m%ZZ(i,:), states_new)) / sqrt(m%HH(i,i))
       logpdf =  logpdf - 0.5d0 * z**2
    end do

  end function logpdfy_kernel_lgss

  logical function solve_serial_lgss(self, para) result(converged)

    class(fortress_lgss_model), intent(inout) :: self
    real(wp), intent(in) :: para(self%npara)
    integer :: error
    
    call self%system_matrices(para, error)
    converged = error == 0

  end function solve_serial_lgss

  logical function solve_parallel_lgss(self, para, nproc, rank) result(converged)

    class(fortress_lgss_model), intent(inout) :: self
    real(wp), intent(in) :: para(self%npara)
    integer, intent(in) :: nproc, rank
    integer :: error
    
    call self%system_matrices(para, error)
    converged = error == 0

  end function solve_parallel_lgss



 function policy_function_lgss(self, states_old, shocks_new) result(states_curr)
    class(fortress_lgss_model), intent(inout) :: self

    real(wp), intent(in) :: states_old(self%ns), shocks_new(self%neps)
    real(wp) :: states_curr(self%ns)

    real(wp) :: tmp(self%neps)

    integer :: i

    states_curr = 0.0d0
    do i = 1, self%neps
       tmp(i) = sqrt(self%QQ(i,i))*shocks_new(i)
    end do


    call dgemv('n', self%ns, self%neps, 1.0d0, self%RR, self%ns, tmp, 1, 0.0d0, states_curr, 1)
    call dgemv('n', self%ns, self%ns, 1.0d0, self%TT, self%ns, states_old, 1, 1.0d0, states_curr, 1)

  end function policy_function_lgss


!*******************************************************************************
!>  whether the function is inbounds
!
logical function inbounds(self, para)

  class(fortress_abstract_bayesian_model), intent(in) :: self
  real(wp) :: para(self%npara)

  integer :: j
  
  inbounds = .true.


  associate(prior => self%prior)
  select type(prior)
  class is (model_prior)
 
     do j = 1, prior%npara
        if ( (para(j) < prior%plower(j)) .or. &
             (para(j) > prior%pupper(j))) then
           inbounds = .false.
        end if
     end do


  class default
     inbounds = .true.
  end select
  end associate
end function inbounds







!*******************************************************************************
!>
!  Finalizer for [[fortress_abstract_bayesian_model]].
! pure subroutine finalize(self)
!   type(fortress_abstract_bayesian_model), intent(inout) :: self

!   if (allocated(self%yy)) deallocate(self%yy)

  ! end subroutine finalize
!*******************************************************************************


!*******************************************************************************  
end module fortress_bayesian_model_t
!*******************************************************************************
