module fortress_random_t

#:setvar FC defined('IFORT')

  #:if FC == 0
  use randlib, only : rand_uniform, rand_normal
  #:endif

  
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
      
      integer :: errcode 

      rvs_lb = rn%uniform_lb
      rvs_ub = rn%uniform_ub

      if (present(lb)) rvs_lb = lb
      if (present(ub)) rvs_ub = ub

#:if FC > 0
      errcode = vdrnguniform( rn%methodn, rn%stream, dim_a*dim_b, rvs, rvs_lb, rvs_ub)
#:endif
    end function uniform_rvs

  ! function beta_rvs(self, dim_a, dim_b, mean, std) result(ans)

  ! end function beta_rvs

end module fortress_random_t
