module fortress_smc_utilities
  use, intrinsic :: iso_fortran_env, only: wp => real64
  use fortress_smc_particles_t, only: fortress_smc_particles


contains
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

    if (nblocks == 1) then
       break_points = [0,npara]
       return
    endif
    

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
    new_weight = (phi-phi_old)*(parasim%loglh - parasim%loglhold) 
    new_weight2 = 2.0_wp*(phi-phi_old)*(parasim%loglh- parasim%loglhold) 

    max1 = maxval(new_weight)
    max2 = maxval(new_weight2)

    a1 = sum(exp(new_weight2-max2)*parasim%weights) 
    a2 = sum(exp(new_weight-max1)*parasim%weights)**2 

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



end module fortress_smc_utilities
