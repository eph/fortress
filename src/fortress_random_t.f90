module fortress_random_t
  use, intrinsic :: iso_fortran_env, only: wp => real64



#:setvar FC defined('IFORT')

  #:if FC == 0
  use randlib, only : rand_uniform, rand_normal, rand_gamma, rand_beta
  #:endif

  implicit none

  
  type fortress_random

     integer :: seed = 1848        ! year of the revolutions

     double precision :: normal_mean = 0.0d0
     double precision :: normal_std  = 1.0d0

     double precision :: uniform_lb = 0.0d0
     double precision :: uniform_ub = 1.0d0


#:if FC>0
     integer :: methodu = VSL_METHOD_DUNIFORM_STD
     integer :: methodn = VSL_METHOD_DGAUSSIAN_BOXMULLER
     type(vsl_stream_state) :: stream 
#:endif

   contains

     procedure :: norm_rvs
     procedure :: uniform_rvs
     procedure :: beta_rvs
     procedure :: gamma_rvs
     procedure :: inv_gamma_rvs
     procedure :: iw_rvs
     
  end type fortress_random

  interface fortress_random
     module procedure new_fortress_random
  end interface fortress_random

contains

  type(fortress_random) function new_fortress_random(seed) result(rn)

    integer, intent(in), optional :: seed

    integer :: errcode

    if (present(seed)) rn%seed = seed

#:if FC > 0
    errcode = vslnewstream(rn%stream, rn%brng, rn%seed)
#:endif

  end function new_fortress_random

  function norm_rvs(rn, dim_a, dim_b, mu, sig) result(rvs)
    class(fortress_random) :: rn

    integer, intent(in) :: dim_a, dim_b
    double precision :: rvs(dim_a, dim_b)

    double precision, intent(in), optional :: mu, sig
    double precision :: rvs_mu, rvs_sig

    integer :: i, j

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
    errcode = vdrnggaussian( rn%methodn, rn%stream, dim_a*dim_b, rvs, rvs_mu, rvs_sig)
#:endif
  end function norm_rvs

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

    integer :: i, j


#:if FC==0
    do j = 1, dim_b
       do i = 1, dim_a
          rvs(i,j) = rand_gamma(theta,k)
       end do
    end do
#:elif FC > 0
    print*,'not implemented'
    stop
    !errcode = vdrnguniform( rn%methodn, rn%stream, dim_a*dim_b, rvs, rvs_lb, rvs_ub)
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
    print*,'not implemented'
    stop
    !errcode = vdrnguniform( rn%methodn, rn%stream, dim_a*dim_b, rvs, rvs_lb, rvs_ub)
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
