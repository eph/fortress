#:setvar FC defined('IFORT')
#:if FC>0
     include 'mkl_vsl.fi'
#:endif
module fortress_random_t
  use, intrinsic :: iso_fortran_env, only: wp => real64

  use fortress_linalg, only: cholesky

#:if FC == 0
  use randlib, only : rand_uniform, rand_normal, rand_gamma, rand_beta

#:else
  use mkl_vsl
  !use mkl_vsl_type
#:endif

  implicit none


  type fortress_random

     integer :: seed = 1848        ! year of the revolutions

     double precision :: normal_mean = 0.0d0
     double precision :: normal_std  = 1.0d0

     double precision :: uniform_lb = 0.0d0
     double precision :: uniform_ub = 1.0d0


#:if FC>0
     integer :: brng = vsl_brng_mt19937

     integer :: methodu = VSL_RNG_METHOD_UNIFORM_STD
     integer :: methodn = VSL_RNG_METHOD_GAUSSIAN_BOXMULLER
     integer :: methodg = VSL_RNG_METHOD_GAMMA_GNORM
     integer :: methodb = VSL_RNG_METHOD_BETA_CJA
     type(vsl_stream_state) :: stream
#:endif

   contains

     procedure :: norm_rvs
     procedure :: uniform_rvs
     procedure :: beta_rvs
     procedure :: gamma_rvs
     procedure :: inv_gamma_rvs
     procedure :: iw_rvs
     procedure :: mv_norm_rvs

  end type fortress_random

  interface fortress_random
     module procedure new_fortress_random
  end interface fortress_random

contains

  type(fortress_random) function new_fortress_random(seed) result(rn)

    integer, intent(in), optional :: seed

    integer :: errcode, n, i
    integer, allocatable :: nseed(:)

    if (present(seed)) rn%seed = seed

#:if FC > 0
    errcode = vslnewstream(rn%stream, rn%brng, rn%seed)
#:else
    call random_seed(size=n)
    allocate(nseed(n))
    nseed = rn%seed + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put=nseed)
    deallocate(nseed)
#:endif

  end function new_fortress_random

  function norm_rvs(rn, dim_a, dim_b, mu, sig) result(rvs)
    class(fortress_random) :: rn

    integer, intent(in) :: dim_a, dim_b
    double precision :: rvs(dim_a, dim_b)

    double precision, intent(in), optional :: mu, sig
    double precision :: rvs_mu, rvs_sig

    integer :: i, j, errcode

    rvs_mu = rn%normal_mean
    rvs_sig = rn%normal_std

    if (present(mu)) rvs_mu = mu
    if (present(sig)) rvs_sig = sig

#:if FC==0
    do j = 1, dim_b
       do i = 1, dim_a
          rvs(i,j) = rand_normal(rvs_mu, rvs_sig)
       end do
    end do
#:elif FC>0
       !print*,rn%methodn, rn%methodg, dim_a, dim_b, rvs_mu, rvs_sig!, rn%stream
       errcode = vdrnggaussian( rn%methodn, rn%stream, dim_b*dim_a, rvs, rvs_mu, rvs_sig)
#:endif
  end function norm_rvs

  function mv_norm_rvs(rn, dim_a, dim_b, mu, L, use_cholesky) result(rvs)

    class(fortress_random), intent(inout) :: rn
    integer, intent(in) :: dim_a, dim_b

    double precision, intent(in), optional :: mu(dim_b), L(dim_b, dim_b)
    double precision :: rvs(dim_a, dim_b)

    logical, intent(in), optional :: use_cholesky
    logical :: do_cholesky

    double precision :: eps(dim_a, dim_b), L_copy(dim_b, dim_b)
    integer :: i , info

    do i = 1, dim_b
       rvs(:,i) = mu(i)
    end do
    eps = rn%norm_rvs(dim_a, dim_b)

    do_cholesky = .true.
    if (present(use_cholesky)) do_cholesky = use_cholesky

    if (do_cholesky .eqv. .false.) then
       L_copy = L
       call cholesky(L_copy, info)
       call dgemm('n', 't', dim_a, dim_b, dim_b, 1.0_wp, eps, dim_a, L_copy, dim_b, 1.0_wp, rvs, dim_a)
    else
       call dgemm('n', 't', dim_a, dim_b, dim_b, 1.0_wp, eps, dim_a, L, dim_b, 1.0_wp, rvs, dim_a)
    end if




  end function mv_norm_rvs


  function uniform_rvs(rn, dim_a, dim_b, lb, ub) result(rvs)
    class(fortress_random) :: rn
    integer, intent(in) :: dim_a, dim_b

    double precision :: rvs(dim_a, dim_b)

    double precision, intent(in), optional :: lb, ub
    double precision :: rvs_lb, rvs_ub

    integer :: errcode, i, j

    rvs_lb = rn%uniform_lb
    rvs_ub = rn%uniform_ub

    if (present(lb)) rvs_lb = lb
    if (present(ub)) rvs_ub = ub
#:if FC==0
    do j = 1, dim_b
       do i = 1, dim_a
          rvs(i,j) = rand_uniform(rvs_lb, rvs_ub)
       end do

    end do
#:elif FC > 0
    errcode = vdrnguniform( rn%methodn, rn%stream, dim_a*dim_b, rvs, rvs_lb, rvs_ub)
#:endif
  end function uniform_rvs

  function gamma_rvs(rn, dim_a, dim_b, theta, k) result(rvs)
    class(fortress_random) :: rn
    integer, intent(in) :: dim_a, dim_b

    double precision :: rvs(dim_a, dim_b)

    double precision, intent(in) :: theta, k
    double precision :: rvs_a, rvs_b

    integer :: errcode, i, j

#:if FC==0
    do j = 1, dim_b
       do i = 1, dim_a
          rvs(i,j) = rand_gamma(k,1.0_wp/theta)
       end do
    end do
#:elif FC > 0
    rvs = 0.0_wp
    errcode = vdrnggamma(rn%methodg, rn%stream, dim_a*dim_b, rvs, k, 0.0_wp, theta)
#:endif


  end function gamma_rvs


  function inv_gamma_rvs(rn, dim_a, dim_b, a, b) result(rvs)

    class(fortress_random) :: rn
    integer, intent(in) :: dim_a, dim_b

    double precision :: rvs(dim_a, dim_b)

    double precision, intent(in) :: a, b
    double precision :: rvs_a, rvs_b

    integer :: int_b, i

    real(wp), allocatable :: rand_norm(:,:)

    allocate(rand_norm(int(b), dim_a))

    do i = 1, dim_b
       rand_norm = rn%norm_rvs(int(b), dim_a)
       rvs(:,i) = sqrt(b*a**2 / sum(rand_norm**2, 1))
    end do
    deallocate(rand_norm)
  end function inv_gamma_rvs


  function beta_rvs(rn, dim_a, dim_b, a, b) result(rvs)

    class(fortress_random) :: rn
    integer, intent(in) :: dim_a, dim_b

    double precision :: rvs(dim_a, dim_b)

    real(wp), intent(in), optional :: a, b
    double precision :: rvs_a, rvs_b

    integer :: errcode, i, j

    !rvs_a = rn%uniform_a
    !rvs_b = rn%uniform_b


    if (present(a)) rvs_a = a
    if (present(b)) rvs_b = b
#:if FC==0
    do j = 1, dim_b
       do i = 1, dim_a
          rvs(i,j) = rand_beta(rvs_a, rvs_b)
       end do
    end do
#:elif FC > 0
    rvs = 0.0_wp
    errcode = vdrngbeta(rn%methodb, rn%stream, dim_a*dim_b, rvs, rvs_a, rvs_b, 0.0_wp, 1.0_wp)
#:endif


  end function beta_rvs

  function iw_rvs(rn, S, nu, n) result(iW)

    class(fortress_random) :: rn

    integer, intent(in) :: nu, n
    real(wp), intent(in) :: S(n,n)



    real(wp) :: iW(n, n), chol_S(n,n)

    real(wp) :: ny

    real(wp), allocatable :: dev_iw(:,:), W(:,:)

    real(wp) :: work(3*n)
    integer(wp) :: ipiv(n), info

    integer :: jj


    chol_S = S

    ! C = chol(inv(S),'lower')
    call dgetrf(n,n,chol_S,n,ipiv,info)
    call dgetri(n,chol_S,n,ipiv,work,3*n,info)
    call dpotrf('l',n,chol_S,n,info)

    do jj = 1,n-1
       chol_S(jj+1:n,jj) = 0.0_wp
    end do

    allocate(dev_iw(n,nu), W(n,nu))

    dev_iw = rn%norm_rvs(n, nu)

    call dgemm('n','n',n,nu,n,1.0_wp,chol_S,n,dev_iw,n,0.0_wp,W,n)
    call dgemm('n','t',n,n,nu,1.0_wp,W,n,W,n,0.0_wp,iW,n)

    call dgetrf(n,n,iW,n,ipiv,info)
    call dgetri(n,iW,n,ipiv,work,3*n,info)


    deallocate(dev_iw, W)

  end function iw_rvs


end module fortress_random_t
