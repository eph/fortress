!*******************************************************************************
!>
! author: Ed Herbst
!
!# Description
!  Class for priors.

module fortress_prior_t

  use, intrinsic :: iso_fortran_env, only: wp => real64
  use, intrinsic :: iso_fortran_env, only: error_unit 
  use fortress_random_t, only : fortress_random
  use logbeta, only : betaln, gamln

  implicit none


  type, abstract :: fortress_abstract_prior
     integer :: npara
   contains
     procedure(rvs_interface), deferred :: rvs
     procedure(logpdf_interface), deferred :: logpdf
     procedure :: dlogpdf
  end type fortress_abstract_prior

  interface fortress_abstract_prior

     function rvs_interface(self,nsim,seed,rng) result(parasim)
       use, intrinsic :: iso_fortran_env, only: wp => real64
       import
       class(fortress_abstract_prior), intent(inout) :: self
       integer, intent(in) :: nsim
       integer, optional :: seed 
       type(fortress_random), optional, intent(inout) :: rng

       real(wp) :: parasim(self%npara,nsim)
     end function rvs_interface

     function logpdf_interface(self, para) result(lpdf)
       use, intrinsic :: iso_fortran_env, only: wp => real64
       import fortress_abstract_prior

       class(fortress_abstract_prior), intent(inout) :: self
       real(wp), intent(in) :: para(self%npara)
       real(wp) :: lpdf
     end function logpdf_interface

  end interface fortress_abstract_prior


  integer, parameter :: PARA_BETA = 1
  integer, parameter :: PARA_GAMMA = 2
  integer, parameter :: PARA_NORMAL = 3
  integer, parameter :: PARA_INVGAMMA = 4
  integer, parameter :: PARA_UNIFORM = 5
  integer, parameter :: PARA_FIXED = 6
  integer, parameter :: PARA_INVGAMMA1 = 7

  integer, parameter :: PRIOR_FILE_UNIT = 26

  real(wp), parameter :: M_PI = 3.14159265358979d0

  type, extends(fortress_abstract_prior) :: prior

     integer :: seed = 1848

     integer, allocatable :: ptype(:), pfix(:)
     real(wp), allocatable :: pmean(:), pstdd(:), pval(:)
     real(wp), allocatable :: plower(:), pupper(:)

     type(fortress_random) :: rn

   contains
     procedure :: logpdf
     procedure :: rvs
     procedure :: inbounds
     procedure :: dlogpdf => dlogpdf_frank
  end type prior

  interface prior
     module procedure new_prior
  end interface prior

contains

  type(prior) function new_prior(frankfile) result(pr)

    integer :: nlines, unit, io, i

    character(len=*), intent(in) :: frankfile


    open(PRIOR_FILE_UNIT, file=frankfile)
    nlines = 0

    do
       read(PRIOR_FILE_UNIT,*,iostat=io)
       if (io/=0) exit
       nlines = nlines + 1
    end do

    rewind(PRIOR_FILE_UNIT)

    pr%npara = nlines
    allocate(pr%ptype(pr%npara), pr%pmean(pr%npara), pr%pstdd(pr%npara), &
         pr%pfix(pr%npara), pr%pval(pr%npara), &
         pr%plower(pr%npara), pr%pupper(pr%npara))

    do i = 1, pr%npara
       read(PRIOR_FILE_UNIT,*) pr%ptype(i), pr%pmean(i), pr%pstdd(i), pr%pfix(i), pr%pval(i)
       select case( pr%ptype(i) )
       case( PARA_BETA )
          pr%plower(i) = 0.000001d0
          pr%pupper(i) = 0.999999d0
       case( PARA_GAMMA )
          pr%plower(i) = 0.000001d0
          pr%pupper(i) = 50.00000d0
       case( PARA_NORMAL )
          pr%plower(i) = -9999999.0d0
          pr%pupper(i) =  9999999.0d0
       case( PARA_INVGAMMA )
          pr%plower(i) = 0.0000001d0
          pr%pupper(i) = 50.000000d0
       case( PARA_UNIFORM )
          pr%plower(i) = pr%pmean(i)
          pr%pupper(i) = pr%pstdd(i)
       case( PARA_FIXED )
          pr%plower(i) = pr%pval(i)
          pr%pupper(i) = pr%pval(i)
       case( PARA_INVGAMMA1 )
          pr%plower(i) = 0.0000001d0
          pr%pupper(i) = 50.000000d0

       case default
          print*,'ERROR: in prior(.), parameter', i, 'has a misspecified prior'
          stop
       end select
    end do

    close(PRIOR_FILE_UNIT)
  end function new_prior

  function dlogpdf(self, para) result (dl)
    class(fortress_abstract_prior), intent(inout) :: self
    real(wp), intent(in) :: para(self%npara)

    real(wp) :: dl(self%npara)

    dl = 0.0_wp
  end function dlogpdf


  logical function inbounds(self, para)
    class(prior), intent(inout) :: self
    real(wp), intent(in) :: para(self%npara)

    integer :: i

    inbounds = .true.

    do i=1,self%npara

       if ( (para(i) < self%plower(i)) .or. (para(i) > self%pupper(i)) ) then
          inbounds = .false.
          return
       end if
    end do
  end function inbounds

  real(wp) function logpdf(self, para) result(logprior)

    class(prior), intent(inout) :: self
    real(wp), intent(in) :: para(self%npara)

    real(wp) :: a, b
    integer :: i

    logprior = 0.0d0

    if (.not. self%inbounds(para)) then
       logprior = -1000000000000.0_wp
       return
    end if

    associate(pmean => self%pmean, pstdd => self%pstdd )
      do i = 1, self%npara

         select case ( self%ptype(i) )
         case( PARA_BETA )
            a = (1-pmean(i))*pmean(i)**2/pstdd(i)**2 - pmean(i)
            b = a*(1/pmean(i) - 1)
            logprior = logprior + logbetapdf(para(i),a,b)
         case ( PARA_GAMMA )
            b = pstdd(i)**2/pmean(i) !theta
            a = pmean(i)/b           !k
            logprior = logprior + loggampdf(para(i),a,b)
         case ( PARA_NORMAL )
            a = pmean(i)
            b = pstdd(i)
            logprior = logprior + lognorpdf(para(i),a,b)
         case ( PARA_INVGAMMA )
            a = pmean(i)
            b = pstdd(i)
            logprior = logprior + logigpdf(para(i),a,b)
          case ( PARA_UNIFORM )
             a = pmean(i)
             b = pstdd(i)
             logprior = logprior + log(1.0d0/(b - a))
          !case ( PARA_INVGAMMA1)
          !    a = pmean(i)
          !    b = pstdd(i)
          !    logprior = logprior + logig1pdf(para(i), a, b)
         end select
      end do
    end associate

  end function logpdf

  function rvs(self, nsim, seed, rng) result(parasim)

    class(prior), intent(inout) :: self

    integer, intent(in) :: nsim
    integer, optional :: seed 
    type(fortress_random), optional, intent(inout) :: rng

    integer :: prior_seed
    type(fortress_random) :: prior_rvs

    real(wp) :: parasim(self%npara, nsim)
    real(wp) :: a, b, temp(nsim,1)

    integer :: i 

    prior_seed = 1848
    if (present(seed)) prior_seed = seed

    if (present(rng)) then
       prior_rvs = rng
    else
       prior_rvs = fortress_random(seed=seed)
    end if

    do i = 1, self%npara
       if (self%pfix(i) == 1) then
          parasim(i,:) = self%pfix(i)
       else

          select case (self%ptype(i))
          case( PARA_BETA )
             a = (1-self%pmean(i))*self%pmean(i)**2/self%pstdd(i)**2 - self%pmean(i)
             b = a*(1/self%pmean(i) - 1)
             temp= prior_rvs%beta_rvs(nsim, 1, a=a, b=b)
             parasim(i,:) = temp(:,1)
         case ( PARA_GAMMA )
            b = self%pstdd(i)**2/self%pmean(i) !theta
            a = self%pmean(i)/b                !k
            temp = prior_rvs%gamma_rvs(nsim, 1, theta=a, k=b)
            parasim(i,:) = temp(:,1)
         case ( PARA_NORMAL )
            a = self%pmean(i)
            b = self%pstdd(i)
            temp = prior_rvs%norm_rvs(nsim, 1, mu=a, sig=b)
            parasim(i,:) = temp(:,1)
         case ( PARA_INVGAMMA )
            a = self%pmean(i)
            b = self%pstdd(i)
            temp = prior_rvs%inv_gamma_rvs(nsim, 1, a, b)
            parasim(i,:) = temp(:,1)
         case ( PARA_UNIFORM )
            a = self%pmean(i)
            b = self%pstdd(i)
            temp = prior_rvs%uniform_rvs(nsim, 1, lb=a, ub=b)
            parasim(i,:)  = temp(:,1)
         end select
      end if
   end do
    

  end function rvs

  function dlogpdf_frank(self, para) result(dlogprior)

    class(prior), intent(inout) :: self
    real(wp), intent(in) :: para(self%npara)

    real(wp) :: dlogprior(self%npara)
    real(wp) :: a, b
    integer :: i

    dlogprior = 0.0d0

    if (.not. self%inbounds(para)) then
       dlogprior = 0.0_wp
       return
    end if

    associate(pmean => self%pmean, pstdd => self%pstdd )
      do i = 1, self%npara

         select case ( self%ptype(i) )
         case( PARA_BETA )
            a = (1-pmean(i))*pmean(i)**2/pstdd(i)**2 - pmean(i)
            b = a*(1/pmean(i) - 1)
            dlogprior(i) = dlogbetapdf(para(i),a,b)
         case ( PARA_GAMMA )
            b = pstdd(i)**2/pmean(i) !theta
            a = pmean(i)/b           !k
            dlogprior(i) = dloggampdf(para(i),a,b)
         case ( PARA_NORMAL )
            a = pmean(i)
            b = pstdd(i)
            dlogprior(i) = dlognorpdf(para(i),a,b)
         case ( PARA_INVGAMMA )
            a = pmean(i)
            b = pstdd(i)
            dlogprior(i) = dlogigpdf(para(i),a,b)
          case ( PARA_UNIFORM )
             a = pmean(i)
             b = pstdd(i)
             dlogprior(i) = 0.0_wp
         end select
      end do
    end associate

  end function dlogpdf_frank


  real(wp) function logbetapdf(x, a, b) result(retval)

    real(wp), intent(in) :: x, a, b

    retval = -betaln(a,b) + (a-1.0d0)*log(x) + (b-1.0d0)*log(1.0-x)

  end function logbetapdf

  real(wp) function loggampdf(x, a, b) result(retval)

    real(wp), intent(in) :: x, a, b

    retval = -gamln(a) -a*log(b) + (a-1.0d0)*log(x) - x/b

  end function loggampdf

  real(wp) function lognorpdf(x, a, b) result(retval)

    real(wp), intent(in) :: x, a, b

    retval = -0.5d0*log(2.0d0*M_PI) - log(b) - 0.5d0*(x-a)**2/b**2

  end function lognorpdf

  real(wp) function logigpdf(x,a,b) result(retval)

    real(wp), intent(in) :: x, a, b

    retval = log(2.0d0) - gamln(b/2.0d0) + b/2.0d0*log(b*a**2/2.0d0) &
         -(b+1.0d0)/2.0d0*log(x**2) - b*a**2.0/(2.0d0*x**2)

  end function logigpdf

  function dlogbetapdf(x,a,b);
    real(wp), intent(in) :: x, a, b
    real(wp) :: dlogbetapdf
    dlogbetapdf = (a-1.0_wp)/x - (b-1.0_wp)/(1.0_wp-x);
  endfunction dlogbetapdf

  function dloggampdf(x,a,b);
    real(wp), intent(in) :: x, a, b
    real(wp) :: dloggampdf
    dloggampdf = (a-1.0_wp)/x - 1.0_wp/b;
  endfunction dloggampdf

  function dlognorpdf(x,a,b);
    real(wp), intent(in) :: x, a, b
    real(wp) :: dlognorpdf
    dlognorpdf =  -(x-a)/b**2;
  endfunction dlognorpdf

  function dlogigpdf(x,a,b);
    real(wp), intent(in) :: x, a, b
    real(wp) :: dlogigpdf
    dlogigpdf = -(b+1.0_wp)/x + b*a**2/x**3;
  endfunction dlogigpdf

    function logmvnormpdf(x, mu, chol_sigma) result(logq)
    !! Computes the log of the n-dimensional multivariate normal pdf at x
    !! with mean mu and variance = chol_sigma*chol_sigma'.
    !!
    real(wp), intent(in) :: x(:), mu(:), chol_sigma(:,:)
    real(wp), external :: ddot

    real(wp) :: logq
    real(wp) :: a(size(x, 1)), det_sigma

    integer :: n, i
    n = size(x, 1)

    if (n /= size(mu, 1)) then
       print*, 'mvnormal pdf, size error'
       stop
    endif

    det_sigma = 1.0_wp

    do i = 1, n

       det_sigma = det_sigma*chol_sigma(i,i)

    end do
    det_sigma = det_sigma**2
    a = x - mu
    call dtrsv('l','n', 'n', n, chol_sigma, n, a, 1)

    logq = -n*0.5_wp*log(2.0_wp*3.14159_wp) - 0.5_wp*log(det_sigma) - 0.5*ddot(n, a, 1, a, 1)

  end function logmvnormpdf


end module fortress_prior_t
