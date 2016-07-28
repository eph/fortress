! file: filter.f90
!
!  The file contains the filtering module with three functions
!  - kalman_filter(y, TT, RR, QQ, DD, ZZ, HH, ny, nobs, neps, ns, t0)
!  - kalman_filter_missing(y, TT, RR, QQ, DD, ZZ, HH, ny, nobs, neps, ns, t0)
!  - chand_recursion(y, TT, RR, QQ, DD, ZZ, HH, ny, nobs, neps, ns, t0)
!
!  Both of which evaluate the log likelihood of a LGSS
!
!    s_t = TT*s_{t-1} + RR*eps_t,  eps_t ~ iidN(0, QQ)
!    y_t = DD + ZZ*s_t + eta_t, eta_t ~ iidN(0, HH)
! 
!  where 
!  YY  - ny x nobs
!  TT  - ns x ns
!  RR  - ns x neps
!  QQ  - neps x neps
!  DD  - ny
!  ZZ  - ny x ns
!  HH  - ny x ny
! 
!  t0 is the number of initial observations to condition on.
!
!  The state is initialized as the stationary distribution.
!  kalman_filter_missing can handle missing data entered as `NaN.'
!
!  returns:
!  loglh -- scalar containing the log likelihood.
!
!  Helper routines:
!   - determinant(): computes determinant of matrix.
!   - dlyap(): solves discrete Lyapunov equation via SLICOT.
! 
!  In addition, the module contains the subroutine `reduce_system', which 
!   removes redundant states from the LGSS.
!
!  Gradients can be computed via:
!  - chand_recursion_derivative
!  
! Implemented by Ed Herbst <edward.p.herbst@frb.gov>.
!-------------------------------------------------------------------------------
module filter

  use, intrinsic :: iso_fortran_env, only: wp => real64
  use fortress_linalg, only: determinant, dlyap
  use omp_lib

  implicit none 

  !integer, parameter :: wp = kind(1.0d0)

  real(wp), parameter :: ONE = 1.0_wp, ZERO = 0.0_wp, NEG_ONE = -1.0_wp, M_PI = 3.141592653589793_wp
  real(wp), parameter :: really_small = 1e-10

  real(wp), parameter :: REALLY_NEG = -100000000.0_wp

contains



  function kalman_filter(y, TT, RR, QQ, DD, ZZ, HH, ny, nobs, neps, ns, t0) result(loglh) 
    ! Evaluating the likelihood of LGSS via Kalman Filter.
    !--------------------------------------------------------------------------------
    integer, intent(in) :: ny, nobs, neps, ns, t0
    real(wp), intent(in) :: y(:,:), TT(:,:), RR(:,:), QQ(:,:), DD(:), ZZ(:,:), HH(:,:)
    real(wp) :: loglh

    real(wp) :: At(ns), Pt(ns,ns), RQR(ns,ns), Kt(ns,ny), QQRRp(neps,ns)
    real(wp) :: yhat(ny), nut(ny), Ft(ny,ny), iFt(ny,ny), detFt
    integer :: t, info

    real(wp) :: ZZP0(ny,ns), iFtnut(ny,ny), gain(ns), C(ns,ns), KtiFt(ns,ny), TTPt(ns,ns)

    ! BLAS functions
    real(wp) :: ddot

    ! initialization 
    At = 0.0_wp
    call dgemm('n','t', neps, ns, neps, 1.0_wp, QQ, neps, RR, ns, 0.0_wp, QQRRp, neps)
    call dgemm('n','n', ns, ns, neps, 1.0_wp, RR, ns, QQRRp, neps, 0.0_wp, RQR, ns)

    call dlyap(TT, RQR, Pt, ns, info)

    ! Pt = TT*Pt*TT' + RQR
    call dgemm('n','n', ns, ns, ns, 1.0_wp, TT, ns, Pt, ns, 0.0_wp, TTPt, ns)
    Pt = RQR
    call dgemm('n','t', ns, ns, ns, 1.0_wp, TTPt, ns, TT, ns, 1.0_wp, Pt, ns)

    loglh = 0.0_wp
    do t = 1, nobs

       ! yhat = ZZ*At + DD
       call dcopy(ny, DD, 1, yhat, 1)
       call dgemv('n', ny, ns, ONE, ZZ, ny, At, 1, ONE, yhat, 1)

       ! nut = yt - yhat
       nut = y(:, t) - yhat

       ! Ft = ZZ*Pt*ZZ' + HH
       call dcopy(ny*ny, HH, 1, Ft, 1)
       call dsymm('r', 'l', ny, ns, ONE, Pt, ns, ZZ, ny, ZERO, ZZP0, ny)
       call dgemm('n', 't', ny, ny, ns, ONE, ZZP0, ny, ZZ, ny, ONE, Ft, ny)

       ! iFt = inv(Ft)
       call dcopy(ny*ny, Ft, 1, iFt, 1)
       call dpotrf('u', ny, iFt, ny, info)
       call dpotri('u', ny, iFt, ny, info)

       ! det(Ft)
       detFt = determinant(Ft, ny);

       call dsymv('u', ny, ONE, iFt, ny, nut, 1, ZERO, iFtnut, 1)

       if (t > t0) then
          loglh = loglh - 0.5_wp*ny*log(2*M_PI) - 0.5_wp*log(detFt) & 
               - 0.5_wp*ddot(ny, nut, 1, iFtnut, 1)
       endif

       ! Kt = TT*Pt*ZZ'
       call dgemm('n','t', ns, ny, ns, ONE, TT, ns, ZZP0, ny, ZERO, Kt, ns)

       ! At = TT*At + Kt*iFt*nut'
       call dgemv('n', ns, ny, ONE, Kt, ns, iFtnut, 1, ZERO, gain, 1)
       call dgemv('n', ns, ns, ONE, TT, ns, At, 1, ONE, gain, 1)
       call dcopy(ns, gain, 1, At, 1)

       ! Pt = TT*Pt*TT' + RQR - Kt*iFt*Kt'
       call dgemm('n','n', ns, ns, ns, ONE, TT, ns, Pt, ns, ZERO, C, ns)
       call dsymm('r', 'u', ns, ns, ONE, Pt, ns, TT, ns, ZERO, C, ns)
       call dcopy(ns*ns, RQR, 1, Pt, 1)
       call dgemm('n', 't', ns, ns, ns, ONE, C, ns, TT, ns, ONE, Pt, ns)

       call dsymm('r', 'u', ns, ny, ONE, iFt, ny, Kt, ns, ZERO, KtiFt, ns)
       call dgemm('n', 't', ns, ns, ny, NEG_ONE, KtiFt, ns, Kt, ns, ONE, Pt, ns)

    end do


  end function kalman_filter




  function kalman_filter_missing(y, TT, RR, QQ, DD, ZZ, HH, ny, nobs, neps, ns, t0) result(loglh)
    ! Evaluating the log likelihood of LGSS via the Kalman Filter.
    !  Can handle missing data entered as `NaN.'
    !--------------------------------------------------------------------------------
    integer, intent(in) :: ny, nobs, neps, ns, t0
    real(wp), intent(in) :: y(:,:), TT(:,:), RR(:,:), QQ(:,:), DD(:), ZZ(:,:), HH(:,:)
    real(wp) :: loglh

    real(wp) :: At(ns), Pt(ns,ns), RQR(ns,ns), QQRRp(neps,ns)
    real(wp) :: detFt
    integer :: t,i,j, info

    integer :: nmiss, ngood
    integer, allocatable :: oind(:)
    real(wp), allocatable :: iFtnut(:),ZZPt(:,:), Ft(:,:), iFt(:,:), Kt(:,:), KtiFt(:,:), nut(:),yhat(:)
    real(wp) :: ZZP0(ny,ns), gain(ns), C(ns,ns),  TTPt(ns,ns), HHmiss(ny,ny)

    ! BLAS functions
    real(wp) :: ddot

    ! initialization 
    At = 0.0_wp
    call dgemm('n','t', neps, ns, neps, 1.0_wp, QQ, neps, RR, ns, 0.0_wp, QQRRp, neps)
    call dgemm('n','n', ns, ns, neps, 1.0_wp, RR, ns, QQRRp, neps, 0.0_wp, RQR, ns)

    call dlyap(TT, RQR, Pt, ns, info)

    ! Pt = TT*Pt*TT' + RQR

    call dgemm('n','n', ns, ns, ns, 1.0_wp, TT, ns, Pt, ns, 0.0_wp, TTPt, ns)
    Pt = RQR
    call dgemm('n','t', ns, ns, ns, 1.0_wp, TTPt, ns, TT, ns, 1.0_wp, Pt, ns)


    loglh = 0.0_wp

    do t = 1, nobs
       HHmiss = 0.0_wp

       nmiss = count(isnan(y(:,t)))

       ngood = ny - nmiss
       allocate(oind(ngood),iFtnut(ngood),ZZPt(ngood, ns),Ft(ngood, ngood), &
            iFt(ngood,ngood),Kt(ns,ngood), KtiFt(ns,ngood),nut(ngood),yhat(ngood))
       j = 1;
       do i = 1,ny
          if (isnan(y(i,t))) then 
          else
             oind(j) = i
             j = j + 1;
          end if
       end do

       if (ngood > 0) then

          ! yhat = ZZ*At + DD
          call dcopy(ngood, DD(oind), 1, yhat, 1)
          call dgemv('n', ngood, ns, ONE, ZZ(oind,:), ngood, At, 1, ONE, yhat, 1)

          ! nut = yt - yhat
          nut = y(oind, t) - yhat

          ! Ft = ZZ*Pt*ZZ' + HH
          call dcopy(ngood*ngood, HH(oind, oind), 1, Ft, 1)

          !call dsymm('r', 'l', ngood, ns, ONE, Pt, ns, ZZ(oind,:), ngood, ZERO, ZZPt, ngood)
          call dgemm('n','n', ngood, ns, ns, ONE, ZZ(oind, :), ngood, Pt, ns, ZERO, ZZPt, ngood)
          call dgemm('n', 't', ngood, ngood, ns, ONE, ZZPt, ngood, ZZ(oind,:), ngood, ONE, Ft, ngood)

          ! iFt = inv(Ft)
          call dcopy(ngood*ngood, Ft, 1, iFt, 1)
          call dpotrf('u', ngood, iFt, ngood, info)
          call dpotri('u', ngood, iFt, ngood, info)


          ! det(Ft)
          detFt = determinant(Ft, ngood);

          call dsymv('u', ngood, ONE, iFt, ngood, nut, 1, ZERO, iFtnut, 1)

          if (t > t0) then
             loglh = loglh - 0.5_wp*ngood*log(2*M_PI) - 0.5_wp*log(detFt) & 
                  - 0.5_wp*ddot(ngood, nut, 1, iFtnut, 1)
          endif


          ! Kt = TT*Pt*ZZ'
          call dgemm('n','t', ns, ngood, ns, ONE, TT, ns, ZZPt, ngood, ZERO, Kt, ns)

          ! At = TT*At + Kt*iFt*nut'
          call dgemv('n', ns, ngood, ONE, Kt, ns, iFtnut, 1, ZERO, gain, 1)
          call dgemv('n', ns, ns, ONE, TT, ns, At, 1, ONE, gain, 1)
          call dcopy(ns, gain, 1, At, 1)

          ! Pt = TT*Pt*TT' + RQR - Kt*iFt*Kt'
          call dgemm('n','n', ns, ns, ns, ONE, TT, ns, Pt, ns, ZERO, C, ns)
          call dsymm('r', 'u', ns, ns, ONE, Pt, ns, TT, ns, ZERO, C, ns)
          call dcopy(ns*ns, RQR, 1, Pt, 1)
          call dgemm('n', 't', ns, ns, ns, ONE, C, ns, TT, ns, ONE, Pt, ns)

          call dsymm('r', 'u', ns, ngood, ONE, iFt, ngood, Kt, ns, ZERO, KtiFt, ns)
          call dgemm('n', 't', ns, ns, ngood, NEG_ONE, KtiFt, ns, Kt, ns, ONE, Pt, ns)

       else
          ! At = TT *At-1
          call dgemv('n', ns, ns, ONE, TT, ns, At, 1, ZERO, gain, 1)
          call dcopy(ns, gain, 1, At, 1)


          call dgemm('n','n', ns, ns, ns, ONE, TT, ns, Pt, ns, ZERO, C, ns)
          call dsymm('r', 'u', ns, ns, ONE, Pt, ns, TT, ns, ZERO, C, ns)
          call dcopy(ns*ns, RQR, 1, Pt, 1)
          call dgemm('n', 't', ns, ns, ns, ONE, C, ns, TT, ns, ONE, Pt, ns)

       end if



       deallocate(oind)

       deallocate(Ft,iFt,iFtnut,Kt,KtiFt,ZZPt,nut,yhat)
    end do


  end function kalman_filter_missing


  subroutine kalman_filter_missing_with_states(y, TT, RR, QQ, DD, ZZ, HH, ny, nobs, neps, ns, t0, &
    loglh, smooth_states, smooth_vars, smooth_shocks) 
    ! Evaluating the log likelihood of LGSS via the Kalman Filter.
    !  Can handle missing data entered as `NaN.'
    !--------------------------------------------------------------------------------
    integer, intent(in) :: ny, nobs, neps, ns, t0
    real(wp), intent(in) :: y(:,:), TT(:,:), RR(:,:), QQ(:,:), DD(:), ZZ(:,:), HH(:,:)
    real(wp),intent(out) :: loglh, smooth_states(nobs,ns), smooth_shocks(nobs,neps), smooth_vars(ns,ns,nobs)

    real(wp) :: At(ns), Pt(ns,ns), RQR(ns,ns), QQRRp(neps,ns)
    real(wp) :: detFt
    integer :: t,i,j, info

    integer :: nmiss, ngood
    integer, allocatable :: oind(:)
    real(wp), allocatable :: iFtnut(:),ZZPt(:,:), Ft(:,:), iFt(:,:), Kt(:,:), KtiFt(:,:), nut(:),yhat(:), L(:,:)
    real(wp) :: ZZP0(ny,ns), gain(ns), C(ns,ns),  TTPt(ns,ns), HHmiss(ny,ny), beta_(ns,ns,nobs), eye(ns,ns)

    real(wp) :: fcst_error(ny,nobs), L_mat(ny,ns,nobs), P_mat(ns,ns,nobs), cc(ns),cc_old(ns), psi_(ns,ns)
    ! BLAS functions
    real(wp) :: ddot

    ! initialization 
    At = 0.0_wp
    call dgemm('n','t', neps, ns, neps, 1.0_wp, QQ, neps, RR, ns, 0.0_wp, QQRRp, neps)
    call dgemm('n','n', ns, ns, neps, 1.0_wp, RR, ns, QQRRp, neps, 0.0_wp, RQR, ns)

    call dlyap(TT, RQR, Pt, ns, info)

    ! Pt = TT*Pt*TT' + RQR

    call dgemm('n','n', ns, ns, ns, 1.0_wp, TT, ns, Pt, ns, 0.0_wp, TTPt, ns)
    Pt = RQR
    call dgemm('n','t', ns, ns, ns, 1.0_wp, TTPt, ns, TT, ns, 1.0_wp, Pt, ns)

    eye = 0.0_wp
    do t = 1,ns
       eye(t,t) = 1.0_wp
    end do
    fcst_error = 0.0_wp  
    
    loglh = 0.0_wp
    P_mat = 0.0_wp
    L_mat = 0.0_wp
    do t = 1, nobs
       HHmiss = 0.0_wp

       nmiss = count(isnan(y(:,t)))

       ngood = ny - nmiss
       allocate(oind(ngood),iFtnut(ngood),ZZPt(ngood, ns),Ft(ngood, ngood), &
            iFt(ngood,ngood),Kt(ns,ngood), KtiFt(ns,ngood),nut(ngood),yhat(ngood),&
            L(ngood,ns))

       j = 1;
       do i = 1,ny
          if (isnan(y(i,t))) then 
          else
             oind(j) = i
             j = j + 1;
          end if
       end do

       smooth_states(t,:) = At
       P_mat(:,:,t) = Pt

       if (ngood > 0) then 

          ! yhat = ZZ*At + DD
          call dcopy(ngood, DD(oind), 1, yhat, 1)
          call dgemv('n', ngood, ns, ONE, ZZ(oind,:), ngood, At, 1, ONE, yhat, 1)

          ! nut = yt - yhat
          nut = y(oind, t) - yhat
          fcst_error(1:ngood,t) = nut

          ! Ft = ZZ*Pt*ZZ' + HH
          call dcopy(ngood*ngood, HH(oind, oind), 1, Ft, 1)

          !call dsymm('r', 'l', ngood, ns, ONE, Pt, ns, ZZ(oind,:), ngood, ZERO, ZZPt, ngood)
          call dgemm('n','n', ngood, ns, ns, ONE, ZZ(oind, :), ngood, Pt, ns, ZERO, ZZPt, ngood)
          call dgemm('n', 't', ngood, ngood, ns, ONE, ZZPt, ngood, ZZ(oind,:), ngood, ONE, Ft, ngood)

          ! iFt = inv(Ft)
          call dcopy(ngood*ngood, Ft, 1, iFt, 1)
          call dpotrf('u', ngood, iFt, ngood, info)
          call dpotri('u', ngood, iFt, ngood, info)


          ! det(Ft)
          detFt = determinant(Ft, ngood);

          call dsymv('u', ngood, ONE, iFt, ngood, nut, 1, ZERO, iFtnut, 1)

          if (t > t0) then
             loglh = loglh - 0.5_wp*ngood*log(2*M_PI) - 0.5_wp*log(detFt) & 
                  - 0.5_wp*ddot(ngood, nut, 1, iFtnut, 1)
          endif

          ! Kt = TT*Pt*ZZ'
          call dgemm('n','t', ns, ngood, ns, ONE, TT, ns, ZZPt, ngood, ZERO, Kt, ns)

          ! At = TT*At + Kt*iFt*nut'
          call dgemv('n', ns, ngood, ONE, Kt, ns, iFtnut, 1, ZERO, gain, 1)
          call dgemv('n', ns, ns, ONE, TT, ns, At, 1, ONE, gain, 1)
          call dcopy(ns, gain, 1, At, 1)

          ! Pt = TT*Pt*TT' + RQR - Kt*iFt*Kt'
          call dgemm('n','n', ns, ns, ns, ONE, TT, ns, Pt, ns, ZERO, C, ns)
          call dsymm('r', 'u', ns, ns, ONE, Pt, ns, TT, ns, ZERO, C, ns)
          call dcopy(ns*ns, RQR, 1, Pt, 1)
          call dgemm('n', 't', ns, ns, ns, ONE, C, ns, TT, ns, ONE, Pt, ns)

          call dsymm('r', 'u', ns, ngood, ONE, iFt, ngood, Kt, ns, ZERO, KtiFt, ns)
          call dgemm('n', 't', ns, ns, ngood, NEG_ONE, KtiFt, ns, Kt, ns, ONE, Pt, ns)

          ! from Hess's code
          !call dgemm('t','n',ns,ngood,ngood,ONE,ZZ(oind,:),ngood,iFt,ngood,0.0_wp,L,ns)
          call dsymm('l','u',ngood,ns,ONE,iFt,ngood,ZZ(oind,:),ngood,0.0_wp,L,ngood)
          L_mat(1:ngood,:,t) = L

          beta_(:,:,t) = eye
          !call dcopy(ns*ns,eye, 1, beta_(t,:,:),1)
          call dgemm('t','n',ns,ns,ngood,-1.0_wp,L,ngood,ZZPt,ngood,ONE,beta_(:,:,t),ns)

          ! if (t == nobs-1) then
          !    open(2,file='beta_pre.txt',action='write')
          !    do j = 1,ns
          !       write(2,'(200f32.8)') beta_(j,:,t)
          !    end do
          !    close(2)

          !    open(1,file='L2.txt',action='write')
          !    do j = 1,ngood
          !       write(1,'(200f32.8)') L(j,:)
          !    end do
          !    close(1)

          !    open(1,file='ZZPt.txt',action='write')
          !    do j = 1,ngood
          !       write(1,'(200f32.8)') ZZPt(j,:)
          !    end do
          !    close(1)

          ! end if

          if (t > 1) then
             call dgemm('n','t',ns,ns,ns,1.0_wp,beta_(:,:,t-1),ns,TT,ns,ZERO,TTPt,ns)
             !call dcopy(ns*ns, TTPt, 1, beta_(:,:,t-1), 1)
             beta_(:,:,t-1) = TTPt
          end if


    else
       ! At = TT *At-1
       call dgemv('n', ns, ns, ONE, TT, ns, At, 1, ZERO, gain, 1)
       call dcopy(ns, gain, 1, At, 1)


       call dgemm('n','n', ns, ns, ns, ONE, TT, ns, Pt, ns, ZERO, C, ns)
       call dsymm('r', 'u', ns, ns, ONE, Pt, ns, TT, ns, ZERO, C, ns)
       call dcopy(ns*ns, RQR, 1, Pt, 1)
       call dgemm('n', 't', ns, ns, ns, ONE, C, ns, TT, ns, ONE, Pt, ns)

    end if


       ! if (t == nobs) then
       !    open(2,file='beta_post.txt',action='write')
       !    do j = 1,ns
       !       write(2,'(200f32.8)') beta_(:,j,t-1)
       !    end do
       !    close(2)
       ! end if

       deallocate(oind)

       deallocate(Ft,iFt,iFtnut,Kt,KtiFt,ZZPt,nut,yhat,L)
    end do


    ! very rudimentary smoother
    cc = 0.0_wp
    cc_old = 0.0_wp

    psi_ = 0.0_wp

    do t=nobs,1,-1

       nmiss = count(isnan(y(:,t)))

       ngood = ny - nmiss
       allocate(oind(ngood))
       j = 1;
       do i = 1,ny
          if (isnan(y(i,t))) then 
          else
             oind(j) = i
             j = j + 1;
          end if
       end do

       ! psi = beta*psi_*beta' + L'ZZ;
       call dgemm('n','t',ns,ns,ns,ONE,psi_,ns,beta_(:,:,t),ns,ZERO,TTPt,ns)
       call dgemm('n','n',ns,ns,ns,ONE,beta_(:,:,t),ns,TTPt,ns,ZERO,psi_,ns)
       call dgemm('t','n',ns,ns,ngood,ONE,L_mat(1:ngood,:,t),ngood,ZZ(oind,:),ngood,ONE,psi_,ns)

       ! c = L_t*eta_t + beta_t*c
       ! a_tT = a_t + P_t*c
       cc_old = cc
       call dgemv('n',ns,ns,ONE,beta_(:,:,t),ns,cc_old,1,ZERO,cc,1)
       call dgemv('t',ngood,ns,ONE,L_mat(1:ngood,:,t),ngood,fcst_error(1:ngood,t),1,ONE,cc,1)
       call dgemv('n',ns,ns,ONE,P_mat(:,:,t),ns,cc,1,ONE,smooth_states(t,:),1)

       call dgemv('n',neps,ns,ONE,QQRRp,neps,cc,1,ONE,smooth_shocks(t,:),1)

       smooth_vars(:,:,t) = P_mat(:,:,t)
       call dgemm('n','n',ns,ns,ns,NEG_ONE,psi_,ns,P_mat(:,:,t),ns,ZERO,TTPt,ns)
       call dgemm('n','n',ns,ns,ns,ONE,P_mat(:,:,t),ns,psi_,ns,ONE,smooth_vars,ns)
       

       deallocate(oind)

    end do


  end subroutine kalman_filter_missing_with_states



  function chand_recursion(y, TT, RR, QQ, DD, ZZ, HH, ny, nobs, neps, ns, t0) result(loglh)
    ! Evaluating the log likelihood of a LGSS via the Chandrasekhar Recursions.
    !--------------------------------------------------------------------------------
    integer, intent(in) :: ny, nobs, neps, ns, t0
    real(wp), intent(in) :: y(:,:), TT(:,:), RR(:,:), QQ(:,:), DD(:), ZZ(:,:), HH(:,:)
    real(wp) :: loglh

    real(wp) :: Kt(ns, ny), St(ns, ny), Mt(ny, ny), yhat(ny), nut(ny), At(ns), RQR(ns,ns)
    real(wp) :: Ft(ny, ny), iFt(ny, ny), Ft1(ny, ny), iFt1(ny, ny), MSpZp(ny, ny), P0(ns,ns)
    real(wp) :: detFt, Ft2(2,2), Kt1(ns,ny)
    integer :: t

    ! for blas / lapack
    real(wp) :: ZZP0(ny, ns), TTP0(ns,ns), gain(ns), iFtnut(ny), ZZSt(ny,ny), ZZStMt(ny,ny)
    real(wp) :: KtFt(ns, ny), TTSt(ns,ny), QQRRp(neps,ns)
    real(wp) :: ddot
    integer :: info

    ! for slicot
    real(wp) :: scale, U(ns,ns), UH(ns, ns), rcond, ferr, wr(ns), wi(ns), dwork(5*ns*ns), sepd
    integer :: iwork(ns*ns), ldwork
    double precision :: kk

    real(wp) :: tol
    logical :: converged

    loglh = 0.0_wp

    ! initialization 
    At = 0.0_wp

    call dgemm('n','t', neps, ns, neps, 1.0_wp, QQ, neps, RR, ns, 0.0_wp, QQRRp, neps)
    call dgemm('n','n', ns, ns, neps, 1.0_wp, RR, ns, QQRRp, neps, 0.0_wp, RQR, ns)

    call dlyap(TT, RQR, P0, ns, info)

    ! Pt = TT*Pt*TT' + RQR
    !call dgemm('n','n', ns, ns, ns, 1.0_wp, TT, ns, P0, ns, 0.0_wp, TTP0, ns)
    call dsymm('r', 'l', ns, ns, 1.0_wp, P0, ns, TT, ns, 0.0_wp, TTP0, ns)
    P0 = RQR
    call dgemm('n','t', ns, ns, ns, 1.0_wp, TTP0, ns, TT, ns, 1.0_wp, P0, ns)

    ! Ft = ZZ*P0*ZZ' + HH
    call dcopy(ny*ny, HH, 1, Ft, 1)
    call dsymm('r', 'l', ny, ns, ONE, P0, ns, ZZ, ny, ZERO, ZZP0, ny)
    call dgemm('n', 't', ny, ny, ns, ONE, ZZP0, ny, ZZ, ny, ONE, Ft, ny)

    ! iFt = inv(Ft)
    !call dcopy(ny*ny, Ft, 1, iFt, 1)
    iFt = 0.0_wp
    call dlacpy('u', ny, ny, Ft, ny, iFt, ny)
    call dpotrf('u', ny, iFt, ny, info)
    call dpotri('u', ny, iFt, ny, info)


    !St = TT*P0*ZZ.transpose();
    call dsymm('r', 'u', ns, ns, ONE, P0, ns, TT, ns, ZERO, TTP0, ns)
    call dgemm('n', 't', ns, ny, ns, ONE, TTP0, ns, ZZ, ny, ZERO, St, ns)


    ! Mt = -iFt
    Mt = 0.0_wp                 ! this initialization is necessary!
    call daxpy(ny*ny, NEG_ONE, iFt, 1, Mt, 1)

    ! Kt = St*iFt
    call dsymm('r', 'u', ns, ny, ONE, iFt, ny, St, ns, ZERO, Kt, ns)
    !    call dgemm('n','n',ns, ny, ny, ONE, St, ns, iFt, ny, ZERO, Kt, ns)

    tol = 1e-5
    converged = .false.


    do t = 1, nobs

       ! yhat = ZZ*At + DD
       call dcopy(ny, DD, 1, yhat, 1)
       call dgemv('n', ny, ns, ONE, ZZ, ny, At, 1, ONE, yhat, 1)

       ! nut = yt - yhat
       nut = y(:, t) - yhat


       ! det
       detFt = determinant(Ft, ny)


       call dsymv('u', ny, ONE, iFt, ny, nut, 1, ZERO, iFtnut, 1)
       if (t > t0) then
          loglh = loglh - 0.5_wp*ny*log(2*M_PI) - 0.5_wp*log(detFt) & 
               - 0.5_wp*ddot(ny, nut, 1, iFtnut, 1)

       endif

       ! At = TT*At + Kt*nut
       call dgemv('n', ns, ny, ONE, Kt, ns, nut, 1, ZERO, gain, 1)
       call dgemv('n', ns, ns, ONE, TT, ns, At, 1, ONE, gain, 1)
       call dcopy(ns, gain, 1, At, 1)

       if (converged .eqv. .false.) then 
          ! Ft1 = Ft + ZZ*St*Mt*St'ZZ';
          call dcopy(ny*ny, Ft, 1, Ft1, 1)
          call dgemm('n', 'n', ny, ny, ns, ONE, ZZ, ny, St, ns, ZERO, ZZSt, ny)
          call dsymm('r', 'u', ny, ny, ONE, Mt, ny, ZZSt, ny, ZERO, ZZStMt, ny)
          call dgemm('n', 't', ny, ny, ny, ONE, ZZStMt, ny, ZZSt, ny, ONE, Ft1, ny)
          !Ft1 = 0.5_wp*Ft1 + 0.5_wp*transpose(Ft1)

          ! iFt1 = inv(Ft1)
          call dcopy(ny*ny, Ft1, 1, iFt1, 1)
          call dpotrf('u', ny, iFt1, ny, info)
          call dpotri('u', ny, iFt1, ny, info)

          ! Kt = (Kt*Ft + TT*St*MSpZp)*iFt1;
          call dsymm('r', 'u', ns, ny, ONE, Ft, ny, Kt, ns, ZERO, KtFt, ns)
          call dgemm('n', 'n', ns, ny, ns, ONE, TT, ns, St, ns, ZERO, TTSt, ns)
          call dgemm('n', 't', ns, ny, ny, ONE, TTSt, ns, ZZStMt, ny, ONE, KtFt, ns)
          call dsymm('r', 'u', ns, ny, ONE, iFt1, ny, KtFt, ns, ZERO, Kt1, ns)


          ! St = (TT - Kt*ZZ)*St = TT*St - Kt*ZZ*St
          call dgemm('n','n', ns, ny, ny, NEG_ONE, Kt1, ns, ZZSt, ny, ONE, TTSt, ns)
          call dcopy(ny*ns, TTSt, 1, St, 1)


          ! Mt = Mt + MSpZp*iFt*MSpZp.transpose();
          call dsymm('l','u', ny, ny, ONE, iFt, ny, ZZStMt, ny, ZERO, Ft, ny)
          call dgemm('t','n', ny, ny, ny, ONE, ZZStMt, ny, Ft, ny, ONE, Mt, ny)

          if (maxval(abs(Kt1 - Kt)) < tol) then
             converged = .true.
          endif

          ! Ft = Ft1, iFt = iFt1
          call dcopy(ny*ny, Ft1, 1, Ft, 1)
          call dcopy(ny*ny, iFt1, 1, iFt, 1)

          call dcopy(ny*ns, Kt1, 1, Kt, 1)

       end if
    end do



  end function chand_recursion










  !----------------------------------------------------------------------------------------------------
  !
  ! State space system reduction.
  !
  !----------------------------------------------------------------------------------------------------
  
  subroutine reduce_system(TT,RR,ZZ,TTred,RRred,ZZred,ny,ns,neps,ns1,check_dim)

    integer, intent(in) :: ny, ns, neps
    integer, intent(out) :: ns1

    real(wp), intent(in) :: TT(ns,ns), RR(ns,neps), ZZ(ny,ns)
    real(wp), intent(out) :: TTred(:,:), RRred(:,:), ZZred(:,:)

    logical, intent(in) :: check_dim 

    ! for Schur Decomposition
    integer :: nmin, info, nwork
    real(wp) :: U(ns,ns), Q(ns,ns)
    logical :: bwork(ns)
    real(wp) :: wr(ns), wi(ns), work(3*ns)

    real(wp) :: RRtemp(ns,neps)
    integer :: i, j

    real(wp), allocatable :: temp(:,:), ZZQ(:,:)
    integer :: nconteps, conteps(neps)


    U = TT

    call dgees('v','s', selctg, ns, U, ns, nmin, wr, wi, Q, ns, work, 3*ns, bwork, info)

    call dgemm('t','n',ns,neps,ns,1.0_wp, Q, ns, RR, ns, 0.0_wp, RRtemp, ns)

    nconteps = count(sum(abs(RRtemp(nmin+1:ns,:)),dim=1)>really_small)

    ns1 = nmin + nconteps
    if (check_dim .eqv. .true.) return

    j = 1
    do i = 1,neps
       if (sum(abs(RRtemp(nmin+1:ns,i))) > really_small) then
          conteps(j) = i
          j = j + 1
       end if
    end do

    TTred = 0.0_wp
    RRred = 0.0_wp
    ZZred = 0.0_wp

    if (nconteps == 0) then

       TTred = U(1:ns1,1:ns1)
       RRred = RRtemp(1:ns,:);
       call dgemm('n','n',ny, ns1, ns, 1.0_wp, ZZ, ny, Q(:,1:ns1), ns, 0.0_wp, ZZred, ny)

    else

       allocate(temp(ns,ns1), ZZQ(ny,ns))

       TTred(1:nmin,1:nmin) = U(1:nmin,1:nmin)
       call dgemm('n','n', nmin, nconteps, ns-nmin, 1.0_wp, U(1:nmin,nmin+1:ns), &
            nmin, RRtemp(nmin+1:ns,conteps(1:nconteps)), ns-nmin, 0.0_wp, &
            TTred(1:nmin,nmin+1:ns1), nmin)

       RRred(1:nmin,:) = RRtemp(1:nmin,:)
       j = 1;
       do i = nmin+1,ns1
          RRred(i,conteps(j)) = 1.0_wp
          j = j + 1
       end do

       temp = 0.0_wp
       do i = 1,nmin
          temp(i,i) = 1.0_wp
       end do
       temp(nmin+1:ns,nmin+1:ns1) = RRtemp(nmin+1:ns,conteps(1:nconteps))

       call dgemm('n','n', ny, ns, ns, 1.0_wp, ZZ, ny, Q, ns, 0.0_wp, ZZQ, ny)
       call dgemm('n','n', ny, ns1, ns, 1.0_wp, ZZQ, ny, temp, ns, 0.0_wp, ZZred, ny)

       deallocate(temp,ZZQ)

    end if

  end subroutine reduce_system

  logical function selctg(wi, wr)

    real(wp) :: wi, wr

    if (sqrt(wi**2 + wr**2) > 1e-15) then 
       selctg = .true.
    else
       selctg = .false.
    end if

  end function selctg




end module filter
