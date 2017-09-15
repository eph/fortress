module tvt_t
  use, intrinsic :: iso_fortran_env, only: wp => real64

  use fortress_bayesian_model_t, only: fortress_abstract_bayesian_model
  use fortress_prior_t, only: fortress_abstract_prior
  use fortress_random_t, only: fortress_random
  use fortress_util!, only: read_array_from_file, write_array_to_file
  use fortress_linalg, only: cholesky, determinant, Kronecker
  use logbeta, only: gamln

  implicit none

  real(wp), parameter :: M_PI = 3.14159265358979323846d0

  type, public, extends(fortress_abstract_prior) :: MinnesotaPrior
     integer  :: p, constant, ny
     integer  :: nA, nF

     real(wp), allocatable :: hyper_phistar(:), hyper_Omega_inv(:,:), hyper_iw_Psi(:,:)
     integer :: hyper_iw_nu


   contains
     procedure :: rvs => rvs_mn
     procedure :: logpdf => logpdf_mn
     procedure :: para_to_sigma_phi
     procedure :: construct_prior_hyper

     final :: clean_MinnesotaPrior
  end type MinnesotaPrior

  interface MinnesotaPrior
      module procedure new_mn_prior
  end interface MinnesotaPrior


contains

  type(MinnesotaPrior) function new_mn_prior(nobs, p, cons, lam1, lam2, lam3, lam4, lam5, lamxx, tau, ybar, sbar) result(pr)

    integer, intent(in) :: nobs, p, cons
    real(wp), intent(in) :: lam1, lam2, lam3, lam4, lam5, lamxx, tau, ybar(nobs), sbar(nobs)

    pr%p = p
    pr%constant = cons
    pr%nA = nobs*(nobs+1)/2
    pr%nF = nobs*(cons + nobs*p)
    pr%ny = nobs
    pr%npara = pr%nA + pr%nF

    allocate(pr%hyper_phistar(pr%nF), pr%hyper_Omega_inv(pr%nF/pr%ny,pr%nF/pr%ny), &
         pr%hyper_iw_Psi(pr%ny,pr%ny))
    call pr%construct_prior_hyper(lam1, lam2, lam3, lam4, lam5, tau, lamxx, ybar, sbar, &
         pr%hyper_phistar, pr%hyper_Omega_inv, pr%hyper_iw_Psi, pr%hyper_iw_nu)

  end function new_mn_prior

  subroutine clean_MinnesotaPrior(self)
    type(MinnesotaPrior) :: self
    if (allocated(self%hyper_phistar)) deallocate(self%hyper_phistar)
    if (allocated(self%hyper_Omega_inv)) deallocate(self%hyper_Omega_inv)
    if (allocated(self%hyper_iw_Psi)) deallocate(self%hyper_iw_Psi)
  end subroutine clean_MinnesotaPrior

  function mvnormal_pdf(x, mu, chol_sigma) result(logq)
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

  end function mvnormal_pdf

  function logiwishpdf(nu,S,X,n)

    integer, intent(in) :: nu, n
    real(wp), intent(in) :: S(n,n), X(n,n)
    real(wp) :: logiwishpdf
    integer :: i, ipiv(n), info
    real(wp) :: kap, iX(n,n), tr, work(n), ret(n,n),x0

    iX = X
    call dgetrf(n, n, iX, n, ipiv, info)
    call dgetri(n, iX, n, ipiv, work, n, info)
    call dgemm('n','n',n,n,n,1.0_wp,iX,n,S,n,0.0_wp,ret,n)

    tr = 0.0_wp
    kap = 0.0_wp
    x0 = (nu*1.0_wp)/2.0_wp

    do i = 1,n
       tr = tr + ret(i,i)
       kap = kap + gamln(x0+(1-i)*0.5_wp)!gamln(0.5_wp*(nu*1.0_wp-i*1.0_wp-1.0_wp))
    end do

    kap = -kap - 0.5_wp*nu*n*log(2.0_wp) - 0.25_wp*n*(n-1)*log(M_PI)


    logiwishpdf = kap + 0.5_wp*nu*log(determinant(S,n)) &
         - 0.5_wp*(nu+n+1)*log(determinant(X,n)) -0.5_wp*tr

  end function logiwishpdf



  subroutine construct_prior_hyper(self, lam1,lam2,lam3,lam4,lam5,tau,lamxx,ybar,sbar, &
       hyper_phistar, hyper_Omega_inv, hyper_iw_Psi, hyper_iw_nu)

    class(MinnesotaPrior), intent(inout) :: self
    real(wp), intent(in) :: lam1,lam2,lam3,lam4,lam5,lamxx,tau, ybar(self%ny), sbar(self%ny)

    integer :: dumr, dumc, subt, disp, diag, l, j
    real(wp) :: Omega(self%nF/self%ny,self%nF/self%ny), mvn_mu_mat(self%nF/self%ny,self%ny)
    real(wp), allocatable :: dumx(:,:), dumy(:,:), dumxp_dumy(:,:)

    real(wp), intent(out) :: hyper_phistar(self%nF), hyper_Omega_inv(self%nF/self%ny,self%nF/self%ny), hyper_iw_Psi(self%ny,self%ny)
    integer, intent(out):: hyper_iw_nu

    integer :: ipiv(self%nF/self%ny), info, ny, p, constant, nF
    real(wp) :: work(self%nF)
    real(wp):: use_lamxx

    use_lamxx = 0.0_wp

    if (lamxx > 0.0_wp) use_lamxx = (1.0_wp)

    ny = self%ny
    constant = self%constant
    p = self%p
    nF = self%nF

    dumr = ny*2 + int(lam3)*ny + ny*(p-1) + constant + int(use_lamxx)
    dumc = ny*p+constant
    allocate(dumy(dumr,ny), dumx(dumr,dumc), dumxp_dumy(dumc,ny))


    dumx = 0.0_wp
    dumy = 0.0_wp

    ! tightness on prior of own lag coefficients
    do j = 1,ny
       dumy(j,j) = lam1 * tau * sbar(j)
       dumx(j,j) = lam1 * sbar(j)
    end do

    ! scaling coefficient for higher order lags
    subt = ny
    do l = 1,p-1
       do j = 1,ny
          disp = subt + j
          diag = ny * l + j
          dumx(disp, diag) = lam1 * sbar(j) * ((l+1)*1.0_wp)**lam2!(2.0_wp ** lam2)
       end do
       subt = subt + ny
    end do

    ! prior for the covariance matrix
    do l = 1,int(lam3)
       do j = 1,ny
          dumy(subt+j,j) = sbar(j)
       end do
       subt = subt + ny
    end do

    ! co-persistent dummy observations
    subt = subt+1;

    dumy(subt,:) = lam5 * ybar
    do j = 1,p
       dumx(subt,(j-1)*ny+1:j*ny) = lam5 * ybar
    end do

    if (constant==1) dumx(subt,ny*p+constant) = lam5
       


    ! sum of coefficients dummies
    do j = 1,ny
       disp = subt+j
       dumy(disp,j) = lam4 * ybar(j)
       do l = 1,p
          diag = (l-1)*ny + j
          dumx(disp,diag) = lam4 * ybar(j)
       end do
    end do

    subt = subt + ny
    ! mark's thing
    if (use_lamxx>0.0_wp) then
       dumx(dumr,dumc) = lam1/lamxx
    end if

    ! Omega^-1 = (dumx'dumx)^-1
    Omega = 0.0_wp
    call dgemm('t','n',dumc,dumc,dumr,1.0_wp,dumx,dumr,dumx,dumr,0.0_wp,Omega,dumc)
    hyper_Omega_inv = Omega
    call dgetrf(dumc,dumc,hyper_Omega_inv,dumc,ipiv,info)
    call dgetri(dumc,hyper_Omega_inv,dumc,ipiv,work,nF,info)

    ! mvn_mu = (dumx'*dumx)^-1 * dumx ' dumy
    dumxp_dumy = 0.0_wp; mvn_mu_mat = 0.0_wp
    call dgemm('t','n',dumc,ny,dumr,1.0_wp,dumx,dumr,dumy,dumr,0.0_wp,dumxp_dumy,dumc)

    call dgemm('n','n',dumc,ny,dumc,1.0_wp,hyper_Omega_inv,dumc,dumxp_dumy,dumc,0.0_wp,mvn_mu_mat,dumc)
    do j = 1,ny
       hyper_phistar(dumc*(j-1)+1:dumc*j) = mvn_mu_mat(:,j)
    end do

    ! Psi = dumy'dumy - mvn_mu ' omega * mvn_mu
    hyper_iw_Psi = 0.0_wp
    call dgemm('t','n',ny,ny,dumr,1.0_wp,dumy,dumr,dumy,dumr,0.0_wp,hyper_iw_Psi(1:ny,:),ny)
    call dgemm('t','n',ny,ny,dumc,-1.0_wp,hyper_phistar,dumc,dumxp_dumy,dumc,1.0_wp,hyper_iw_Psi(1:ny,:),ny)

    hyper_iw_nu = dumr - dumc

    ! do j = 1, dumr
    !    write(*,'(100f8.3)') dumx(j,:)
    ! end do

    ! print *, ''

    ! do j = 1, dumr
    !    write(*,'(100f8.3)') real(j), dumy(j,:)
    ! end do


    deallocate(dumx,dumy,dumxp_dumy)
  end subroutine construct_prior_hyper


  subroutine para_to_sigma_phi(self, para, sigma, phi)

    class(MinnesotaPrior), intent(inout) :: self
    real(wp), intent(in) :: para(self%npara)

    real(wp), intent(out) :: sigma(self%ny,self%ny), phi(self%nF/self%ny, self%ny)

    integer :: k, ind0

    ind0 = 1
    do k = 1, self%ny
       sigma(1:k,k) = para(ind0:ind0+k-1)
       ind0 = ind0 + k
       sigma(k,1:k) = sigma(1:k,k)
    end do


    do k = 1,self%ny
       phi(:,k) = para(self%nA + (k-1)*self%nF/self%ny+1:self%nA+k*self%nF/self%ny)
    end do


  end subroutine para_to_sigma_phi


  function rvs_mn(self, nsim, seed, rng) result(parasim)

    class(MinnesotaPrior), intent(inout) :: self

    integer, intent(in) :: nsim
    integer, optional :: seed
    type(fortress_random), optional :: rng

    type(fortress_random) :: use_rng
    real(wp) :: parasim(self%npara, nsim)

    real(wp) :: iW(self%ny, self%ny)
    integer :: j, ind0, ny, k, info, nF

    real(wp) :: Fvec(self%nF), dev_F(self%nF, 1), mvn_covar(self%nF, self%nF)

    real(wp) :: F(self%nF/self%ny, self%ny)

    use_rng = fortress_random()
    if (present(rng)) use_rng = rng

    ny = self%ny
    nF = self%nF

    do j = 1, nsim

       iW = use_rng%iw_rvs(self%hyper_iw_Psi, self%hyper_iw_nu, self%ny)

       ind0 = 1
       do k=1,ny
          parasim(ind0:ind0+k-1,j) = iW(1:k,k)
          ind0 = ind0 + k
       end do

       mvn_covar = 0.0_wp
       call Kronecker(ny,ny,nF/ny,nF/ny,nF,nF,0.0_wp,1.0_wp,iW,self%hyper_Omega_inv,mvn_covar)
       info = 0
       call cholesky(mvn_covar, info)

       dev_F = use_rng%norm_rvs(self%nF, 1)
       Fvec = self%hyper_phistar

       call dgemv('n',nF,nF,1.0_wp,mvn_covar,nF,dev_F(:,1),1,1.0_wp,Fvec,1)
       parasim(self%nA+1:self%npara,j) = Fvec

       ! do k = 1,ny
       !    !F(:,k) = Fvec((k-1)*nF/ny+1:k*nF/ny)
       ! end do

    end do
  end function rvs_mn

  real(wp) function logpdf_mn(self, para) result(lpdf)

    class(MinnesotaPrior), intent(inout) :: self
    real(wp), intent(in) :: para(self%npara)

    real(wp) :: sigma(self%ny,self%ny), phi(self%nF/self%ny, self%ny)
    real(wp) :: mvn_covar(self%nF, self%nF)

    integer :: ny, nF, info

    ny = self%ny
    nF = self%nF

    lpdf = 0.0_wp


    call self%para_to_sigma_phi(para, sigma, phi)

    lpdf = logiwishpdf(self%hyper_iw_nu, self%hyper_iw_Psi, sigma, self%ny)

    mvn_covar = 0.0_wp
    call Kronecker(ny,ny,nF/ny,nF/ny,nF,nF,0.0_wp,1.0_wp,sigma,self%hyper_Omega_inv,mvn_covar)
    call cholesky(mvn_covar, info)
    
    lpdf = lpdf + mvnormal_pdf(para(self%nA+1:self%npara), self%hyper_phistar, mvn_covar)


  end function logpdf_mn


end module tvt_t

module fortress_VAR_t

  use, intrinsic :: iso_fortran_env, only: wp => real64

  use fortress_bayesian_model_t, only: fortress_abstract_bayesian_model
  use fortress_prior_t, only: fortress_abstract_prior
  use fortress_random_t, only: fortress_random
  use fortress_util!, only: read_array_from_file, write_array_to_file
  use fortress_linalg, only: cholesky, determinant, Kronecker
  use logbeta, only: gamln

  implicit none

  real(wp), parameter :: M_PI = 3.14159265358979323846d0

  type, public, extends(fortress_abstract_prior) :: SimsZhaSVARPrior
     real(wp), allocatable :: mu(:), var(:,:), chol_var(:,:)
   contains
     procedure rvs
     procedure logpdf

  end type SimsZhaSVARPrior

  interface SimsZhaSVARPrior
     module procedure new_prior
  end interface SimsZhaSVARPrior


  type, public, extends(fortress_abstract_bayesian_model), abstract :: SVAR
     integer :: p, constant
     integer :: nA, nF

     integer :: likT
     real(wp), allocatable :: YY_lik(:,:), XX_lik(:,:)


   contains
     procedure :: construct_SVAR
     generic :: construct_model => construct_SVAR!, construct_SVAR_mn
     procedure(para_to_AF_i), public, deferred :: para_to_AF
     procedure :: lik

  end type SVAR

  abstract interface
     subroutine para_to_AF_i(self, para, A, F)
       use, intrinsic :: iso_fortran_env, only: wp => real64
       import

       implicit none

       class(SVAR), intent(inout) :: self

       real(wp), intent(in) :: para(self%npara)
       real(wp), intent(out) :: A(self%nobs, self%nobs)
       real(wp), intent(out) :: F(self%nobs**2*self%p+self%nobs*self%nobs, self%nobs)


     end subroutine para_to_AF_i
  end interface

contains

 !  subroutine construct_SVAR_mn(self, name, datafile, nobs, T, p, cons, &
 !     lam1, lam2, lam3, lam4, lam5, lamxx, tau, ybar, sbar)

 !    class(SVAR), intent(inout) :: self

 !    character(len=144), intent(in) :: name, datafile
 !    integer, intent(in) :: nobs, T, p, cons, nA, nF

 !    real(wp), intent(in) :: lam1, lam2, lam3, lam4, lam5, lamxx, tau, ybar, sbar

 !    integer :: i,j
 !    self%p = p
 !    self%constant = cons
 !    self%nA = nobs*(nobs+1)/2
 !    self%nF = nobs*(cons + nobs*p)

 !    allocate(self%prior, &
 !    source=MinnesotaPrior(nobs,p,cons, lam1, lam2, lam3, lam4, lam5, lamxx, tau, ybar, sbar))

 !    call self%construct_model(name, datafile, self%nA+self%nF, nobs, T)

 !        self%likT = self%T - self%p

 !    allocate(self%YY_lik(self%likT, self%nobs), &
 !         self%XX_lik(self%likT, self%nobs*self%p+self%constant))

 !    self%XX_lik = 1.0_wp
 !    do i = 1, self%likT
 !       self%YY_lik(i,:) = self%YY(:,self%p+i)
 !       do j = 1, self%p
 !          self%XX_lik(i,(j-1)*self%nobs+1:j*self%nobs) = self%YY(:,i-j+self%P)
 !       end do
 !    end do

 !    self%T = self%likT

 ! end subroutine construct_SVAR_mn



  subroutine construct_SVAR(self, name, datafile, nobs, T, p, cons, nA, nF, &
       prior_mu_file, prior_VAR_file)

    class(SVAR), intent(inout) :: self

    character(len=144), intent(in) :: name, datafile, prior_mu_file, prior_VAR_file
    integer, intent(in) :: nobs, T, p, cons, nA, nF

    integer :: i,j
    self%p = p
    self%constant = cons
    self%nA = nA
    self%nF = nF

    allocate(self%prior, source=SimsZhaSVARPrior(nA+nF, prior_mu_file, prior_var_file))
    call self%construct_model(name, datafile, nA+nF, nobs, T)

        self%likT = self%T - self%p

    allocate(self%YY_lik(self%likT, self%nobs), &
         self%XX_lik(self%likT, self%nobs*self%p+self%constant))


    self%XX_lik = 1.0_wp
    do i = 1, self%likT
       self%YY_lik(i,:) = self%YY(:,self%p+i)
       do j = 1, self%p
          self%XX_lik(i,(j-1)*self%nobs+1:j*self%nobs) = self%YY(:,i-j+self%P)
       end do
    end do

    self%T = self%likT

end subroutine construct_SVAR

  real(wp) function lik(self, para, T) result(l)

    class(SVAR), intent(inout) :: self
    real(wp), intent(in) :: para(self%npara)

    integer, intent(in), optional :: T

    integer :: use_T, ny, nF,  p, constant, TT, info, i

    real(wp) :: likvec(self%likT), A0(self%nobs, self%nobs), F(self%nF/self%nobs, self%nobs), AXiAXip(self%nobs,self%nobs)
    real(wp) :: resid(self%likT, self%nobs)

    use_T = self%likT
    ny = self%nobs
    nF = self%nF
    p = self%p
    TT= self%likT
    constant = self%constant
    if (present(T)) use_T = T

    call self%para_to_AF(para, A0, F)

    likvec = - ny / 2.0_wp*log(2.0_wp*M_PI)
    call dgemm('n','n',TT,ny,ny,1.0_wp,self%YY_lik,TT,A0,ny,0.0_wp,resid,TT)
    call dgemm('n','n',TT,ny,ny*p+constant,-1.0_wp,self%XX_lik,TT,F,ny*p+constant,1.0_wp,resid,TT)
    call dgemm('n','t',ny,ny,ny,1.0_wp,A0,ny,A0,ny,0.0_wp,AXiAXip,ny)
    likvec = likvec + 0.5_wp*log(determinant(AXiAXip,ny))
    do i = 1,ny
       likvec = likvec - 0.5_wp * resid(:,i)**2
    end do

    l = sum(likvec(1:use_T))

  end function lik


  type(SimsZhaSVARPrior) function new_prior(npara, mufile, varfile) result(pr)

    integer, intent(in) :: npara
    character(len=144), intent(in) :: mufile, varfile
    integer :: info


    pr%npara = npara

    allocate(pr%mu(pr%npara), pr%var(pr%npara,pr%npara), pr%chol_var(pr%npara,pr%npara))

    call read_array_from_file(mufile, pr%mu)
    call read_array_from_file(varfile, pr%var)

    pr%chol_var = pr%var
    call cholesky(pr%chol_var, info)

    if (info/=0) print*,'In SimsZhaSVARPrior, variance not pd'

  end function new_prior

  function rvs(self, nsim, seed, rng) result(parasim)

    class(SimsZhaSVARPrior), intent(inout) :: self
    integer, intent(in) :: nsim
    integer, optional :: seed

    type(fortress_random), optional :: rng

    type(fortress_random) :: use_rng
    real(wp) :: parasim(self%npara, nsim), parasimT(nsim, self%npara)

    use_rng = fortress_random()
    if (present(rng)) use_rng = rng

    parasimT = use_rng%mv_norm_rvs(nsim, self%npara, self%mu, self%chol_var)
    parasim = transpose(parasimT)

  end function rvs

  real(wp) function logpdf(self, para) result(lpdf)

    class(SimsZhaSVARPrior), intent(inout) :: self
    real(wp), intent(in) :: para(self%npara)

    lpdf = mvnormal_pdf(para, self%mu, self%chol_var)

  end function logpdf


  function mvnormal_pdf(x, mu, chol_sigma) result(logq)
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

  end function mvnormal_pdf

  function logiwishpdf(nu,S,X,n)

    integer, intent(in) :: nu, n
    real(wp), intent(in) :: S(n,n), X(n,n)
    real(wp) :: logiwishpdf
    integer :: i, ipiv(n), info
    real(wp) :: kap, iX(n,n), tr, work(n), ret(n,n),x0

    iX = X
    call dgetrf(n, n, iX, n, ipiv, info)
    call dgetri(n, iX, n, ipiv, work, n, info)
    call dgemm('n','n',n,n,n,1.0_wp,iX,n,S,n,0.0_wp,ret,n)

    tr = 0.0_wp
    kap = 0.0_wp
    x0 = (nu*1.0_wp)/2.0_wp

    do i = 1,n
       tr = tr + ret(i,i)
       kap = kap + gamln(x0+(1-i)*0.5_wp)!gamln(0.5_wp*(nu*1.0_wp-i*1.0_wp-1.0_wp))
    end do

    kap = -kap - 0.5_wp*nu*n*log(2.0_wp) - 0.25_wp*n*(n-1)*log(M_PI)


    logiwishpdf = kap + 0.5_wp*nu*log(determinant(S,n)) &
         - 0.5_wp*(nu+n+1)*log(determinant(X,n)) -0.5_wp*tr

  end function logiwishpdf


end module fortress_VAR_t
