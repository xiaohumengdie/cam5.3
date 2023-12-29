#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef _COLLAPSE_AND_ALIGN
#define OMP_COLLAPSE_SIMD  $OMP SIMD COLLAPSE(2)
#define DIR_VECTOR_ALIGNED DIR$ VECTOR ALIGNED
#endif

!#define _DBG_ print *,"File:",__FILE__," at ",__LINE__
!#define _DBG_ !DBG
!
!
module prim_advance_mod
  use edgetype_mod, only : EdgeBuffer_t
  use kinds, only : real_kind, iulog
  use perf_mod, only: t_startf, t_stopf, t_barrierf, t_adj_detailf ! _EXTERNAL
  use parallel_mod, only : abortmp, parallel_t
  implicit none
  private
  save
  public :: prim_advance_exp, prim_advance_init, preq_robert3,&
       applyCAMforcing_dynamics, applyCAMforcing, smooth_phis, overwrite_SEdensity

  type (EdgeBuffer_t) :: edge1
  type (EdgeBuffer_t) :: edge2
  type (EdgeBuffer_t) :: edge3p1

  real (kind=real_kind) :: initialized_for_dt   = 0

  real (kind=real_kind), allocatable :: ur_weights(:)

contains

  subroutine prim_advance_init(par,elem,integration)
    use edge_mod, only : initEdgeBuffer
    use element_mod, only : element_t
    use dimensions_mod, only : nlev
    use control_mod, only : qsplit,rsplit
    type(parallel_t) :: par
    type(element_t), intent(inout), target :: elem(:)
    character(len=*)    , intent(in) :: integration
    integer :: i

    if (rsplit==0) then
       call initEdgeBuffer(par,edge3p1,elem,3*nlev+1)
    else
       ! need extra buffer space for dp3d
       call initEdgeBuffer(par,edge3p1,elem,4*nlev+1)
    endif

    if(integration == 'semi_imp') then
       call initEdgeBuffer(par,edge1,elem,nlev)
       call initEdgeBuffer(par,edge2,elem,2*nlev)
    end if

    ! compute averaging weights for RK+LF (tstep_type=1) timestepping:
    allocate(ur_weights(qsplit))
    ur_weights(:)=0.0d0

    if(mod(qsplit,2).NE.0)then
       ur_weights(1)=1.0d0/qsplit
       do i=3,qsplit,2
         ur_weights(i)=2.0d0/qsplit
       enddo
    else
       do i=2,qsplit,2
         ur_weights(i)=2.0d0/qsplit
       enddo
    endif

  end subroutine prim_advance_init


  subroutine prim_advance_exp(elem, deriv, hvcoord, hybrid,&
       dt, tl,  nets, nete, compute_diagnostics)
    use bndry_mod, only : bndry_exchangev
    use control_mod, only : prescribed_wind, qsplit, tstep_type, rsplit, qsplit, moisture, integration
    use derivative_mod, only : derivative_t, vorticity, divergence, gradient, gradient_wk
    use dimensions_mod, only : np, nlev, nlevp, nvar
!    use prim_state_mod, only : prim_printstate
    use edgetype_mod, only : EdgeBuffer_t
    use edge_mod, only : edgevpack, edgevunpack, initEdgeBuffer
    use element_mod, only : element_t
    use hybvcoord_mod, only : hvcoord_t
    use hybrid_mod, only : hybrid_t
    use reduction_mod, only : reductionbuffer_ordered_1d_t
    use time_mod, only : TimeLevel_t,  timelevel_qdp

#ifndef CAM
    use asp_tests, only : asp_advection_vertical
#endif

    implicit none

    type (element_t), intent(inout), target   :: elem(:)
    type (derivative_t), intent(in)   :: deriv
    type (hvcoord_t)                  :: hvcoord

    type (hybrid_t)    , intent(in):: hybrid

    real (kind=real_kind), intent(in) :: dt
    type (TimeLevel_t)   , intent(in) :: tl
    integer              , intent(in) :: nets
    integer              , intent(in) :: nete
    logical, intent(in)               :: compute_diagnostics

    ! =================
    ! Local
    ! =================
    real (kind=real_kind) ::  dt2, time, dt_vis, eta_ave_w
    real (kind=real_kind) ::  eta_dot_dpdn(np,np,nlevp)
    real (kind=real_kind) ::  dp(np,np)
    integer :: ie,nm1,n0,np1,nstep,method,qsplit_stage,k, qn0
    integer :: n,i,j,lx,lenx


!JMD    call t_barrierf('sync_prim_advance_exp', hybrid%par%comm)
    call t_adj_detailf(+1)
    call t_startf('prim_advance_exp')
    nm1   = tl%nm1
    n0    = tl%n0
    np1   = tl%np1
    nstep = tl%nstep

    ! timelevel to use for accessing Qdp() to compute virtual temperature
    qn0 = -1    ! -1 = disabled (assume dry dynamics)
    if ( moisture /= "dry") then
       call TimeLevel_Qdp(tl, qsplit, qn0)  ! compute current Qdp() timelevel
    endif


! integration = "explicit"
!
!   tstep_type=0  pure leapfrog except for very first timestep   CFL=1
!                    typically requires qsplit=4 or 5
!   tstep_type=1  RK2 followed by qsplit-1 leapfrog steps        CFL=close to qsplit
!                    typically requires qsplit=4 or 5
!   tstep_type=2  RK2-SSP 3 stage (as used by tracers)           CFL=.58
!                    optimal in terms of SSP CFL, but not        CFLSSP=2
!                    optimal in terms of CFL
!                    typically requires qsplit=3
!                    but if windspeed > 340m/s, could use this
!                    with qsplit=1
!   tstep_type=3  classic RK3                                    CFL=1.73 (sqrt(3))
!
!   tstep_type=4  Kinnmark&Gray RK4 4 stage                      CFL=sqrt(8)=2.8
!                 should we replace by standard RK4 (CFL=sqrt(8))?
!                 (K&G 1st order method has CFL=3)
!   tstep_type=5  Kinnmark&Gray RK3 5 stage 3rd order            CFL=3.87  (sqrt(15))
!                 From Paul Ullrich.  3rd order for nonlinear terms also
!                 K&G method is only 3rd order for linear
!                 optimal: for windspeeds ~120m/s,gravity: 340m/2
!                 run with qsplit=1
!                 (K&G 2nd order method has CFL=4. tiny CFL improvement not worth 2nd order)
!
! integration = "full_imp"
!
!   tstep_type=1  Backward Euler or BDF2 implicit dynamics
!

! default weights for computing mean dynamics fluxes
    eta_ave_w = 1d0/qsplit

    if(tstep_type==0)then
       method=0                ! pure leapfrog
       if (nstep==0) method=1  ! but use RK2 on first step
    else if (tstep_type==1) then
       method=0                           ! LF
       qsplit_stage = mod(nstep,qsplit)
       if (qsplit_stage==0) method=1      ! RK2 on first of qsplit steps
       ! RK2 + LF scheme has tricky weights:
       eta_ave_w=ur_weights(qsplit_stage+1)
    else
       method = tstep_type                ! other RK variants
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    ! fix dynamical variables, skip dynamics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    if (1==prescribed_wind) then
       time=tl%nstep*dt
       do ie=nets,nete
#ifdef CAM
          ! read in omega met data and store it in derived%omega_prescribed
          ! --- compute eta_dot_dpdn from omega_prescribed ...
          eta_dot_dpdn(:,:,nlev+1) = 0.0d0

          do k = nlev,2,-1
             eta_dot_dpdn(:,:,k) = eta_dot_dpdn(:,:,k+1) - 2.d0*elem(ie)%derived%omega_prescribed(:,:,k)
          enddo

          eta_dot_dpdn(:,:,1) = 0.0d0
#else
          ! assume most fields are constant in time
          elem(ie)%state%ps_v(:,:,np1) = elem(ie)%state%ps_v(:,:,n0)
          elem(ie)%state%lnps(:,:,np1) = elem(ie)%state%lnps(:,:,n0)
          elem(ie)%state%T(:,:,:,np1)   = elem(ie)%state%T(:,:,:,n0)
          elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,n0)
          elem(ie)%derived%div = 0
	  !get eta_dot_dpdn at n0, ps_v is not used or calculated in this routine
          call asp_advection_vertical(time,hvcoord,elem(ie)%state%ps_v(:,:,n0),&
              eta_dot_dpdn)
#endif
          ! accumulate mean fluxes for advection
          elem(ie)%derived%eta_dot_dpdn(:,:,:) = &
               elem(ie)%derived%eta_dot_dpdn(:,:,:) + eta_dot_dpdn(:,:,:)*eta_ave_w

          ! subcycling code uses a mean flux to advect tracers
#if (defined COLUMN_OPENMP)
          !$omp parallel do private(k,dp)
#endif
          do k=1,nlev
             dp(:,:) =&
                  ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                  ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,tl%n0)

             elem(ie)%derived%vn0(:,:,1,k)=elem(ie)%derived%vn0(:,:,1,k)+&
                  eta_ave_w*elem(ie)%state%v(:,:,1,k,n0)*dp(:,:)
             elem(ie)%derived%vn0(:,:,2,k)=elem(ie)%derived%vn0(:,:,2,k)+&
                  eta_ave_w*elem(ie)%state%v(:,:,2,k,n0)*dp(:,:)
          enddo
       end do
       call t_stopf('prim_advance_exp')
       call t_adj_detailf(-1)
       return
    endif


    ! ==================================
    ! Take timestep
    ! ==================================
    dt_vis = dt
    if (method==0) then
       ! regular LF step
       dt2 = 2*dt
       call compute_and_apply_rhs(np1,nm1,n0,qn0,dt2,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,eta_ave_w)
       dt_vis = dt2  ! dt to use for time-split dissipation
    else if (method==1) then
       ! RK2
       ! forward euler to u(dt/2) = u(0) + (dt/2) RHS(0)  (store in u(np1))
       call compute_and_apply_rhs(np1,n0,n0,qn0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,0d0)
       ! leapfrog:  u(dt) = u(0) + dt RHS(dt/2)     (store in u(np1))
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,eta_ave_w)
    else if (method==2) then
       ! RK2-SSP 3 stage.  matches tracer scheme. optimal SSP CFL, but
       ! not optimal for regular CFL
       ! u1 = u0 + dt/2 RHS(u0)
       call compute_and_apply_rhs(np1,n0,n0,qn0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,eta_ave_w/3)
       ! u2 = u1 + dt/2 RHS(u1)
       call compute_and_apply_rhs(np1,np1,np1,qn0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,eta_ave_w/3)
       ! u3 = u2 + dt/2 RHS(u2)
       call compute_and_apply_rhs(np1,np1,np1,qn0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,eta_ave_w/3)
       ! unew = u/3 +2*u3/3  = u + 1/3 (RHS(u) + RHS(u1) + RHS(u2))
       do ie=nets,nete
          elem(ie)%state%v(:,:,:,:,np1)= elem(ie)%state%v(:,:,:,:,n0)/3 &
               + 2*elem(ie)%state%v(:,:,:,:,np1)/3
          elem(ie)%state%T(:,:,:,np1)= elem(ie)%state%T(:,:,:,n0)/3 &
               + 2*elem(ie)%state%T(:,:,:,np1)/3
          if (rsplit==0) then
             elem(ie)%state%ps_v(:,:,np1)= elem(ie)%state%ps_v(:,:,n0)/3 &
                  + 2*elem(ie)%state%ps_v(:,:,np1)/3
          else
             elem(ie)%state%dp3d(:,:,:,np1)= elem(ie)%state%dp3d(:,:,:,n0)/3 &
                  + 2*elem(ie)%state%dp3d(:,:,:,np1)/3
          endif
       enddo
    else if (method==3) then
       ! classic RK3  CFL=sqrt(3)
       ! u1 = u0 + dt/3 RHS(u0)
       call compute_and_apply_rhs(np1,n0,n0,qn0,dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,0d0)
       ! u2 = u0 + dt/2 RHS(u1)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0)
       ! u3 = u0 + dt RHS(u2)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,eta_ave_w)
    else if (method==4) then
       ! KG 4th order 4 stage:   CFL=sqrt(8)
       ! low storage version of classic RK4
       ! u1 = u0 + dt/4 RHS(u0)
       call compute_and_apply_rhs(np1,n0,n0,qn0,dt/4,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,0d0)
       ! u2 = u0 + dt/3 RHS(u1)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0)
       ! u3 = u0 + dt/2 RHS(u2)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0)
       ! u4 = u0 + dt RHS(u3)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,eta_ave_w)
    else if (method==5) then
#if 0
       ! KG 3nd order 5 stage:   CFL=sqrt( 4^2 -1) = 3.87
       ! but nonlinearly only 2nd order
       ! u1 = u0 + dt/5 RHS(u0)
       call compute_and_apply_rhs(np1,n0,n0,qn0,dt/5,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,0d0)
       ! u2 = u0 + dt/5 RHS(u1)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt/5,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0)
       ! u3 = u0 + dt/3 RHS(u2)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0)
       ! u4 = u0 + dt/2 RHS(u3)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0)
       ! u5 = u0 + dt RHS(u4)
       call compute_and_apply_rhs(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,eta_ave_w)
#else
       ! Ullrich 3nd order 5 stage:   CFL=sqrt( 4^2 -1) = 3.87
       ! u1 = u0 + dt/5 RHS(u0)  (save u1 in timelevel nm1)

       !
       ! phl: rhs: t=t
       !
       call compute_and_apply_rhs(nm1,n0,n0,qn0,dt/5,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,eta_ave_w/4)
       ! u2 = u0 + dt/5 RHS(u1)

       !
       ! phl: rhs: t=t+dt/5
       !
       call compute_and_apply_rhs(np1,n0,nm1,qn0,dt/5,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0)
       ! u3 = u0 + dt/3 RHS(u2)
       !
       ! phl: rhs: t=t+2*dt/5
       !

       call compute_and_apply_rhs(np1,n0,np1,qn0,dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0)
       ! u4 = u0 + 2dt/3 RHS(u3)

       !
       ! phl: rhs: t=t+2*dt/5+dt/3
       !

       call compute_and_apply_rhs(np1,n0,np1,qn0,2*dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0)

       ! compute (5*u1/4 - u0/4) in timelevel nm1:
       do ie=nets,nete
          elem(ie)%state%v(:,:,:,:,nm1)= (5*elem(ie)%state%v(:,:,:,:,nm1) &
               - elem(ie)%state%v(:,:,:,:,n0) ) /4
          elem(ie)%state%T(:,:,:,nm1)= (5*elem(ie)%state%T(:,:,:,nm1) &
               - elem(ie)%state%T(:,:,:,n0) )/4
          if (rsplit==0) then
             elem(ie)%state%ps_v(:,:,nm1)= ( 5*elem(ie)%state%ps_v(:,:,nm1) &
                  - elem(ie)%state%ps_v(:,:,n0) )/4
          else
             elem(ie)%state%dp3d(:,:,:,nm1)= (5*elem(ie)%state%dp3d(:,:,:,nm1) &
                  - elem(ie)%state%dp3d(:,:,:,n0) )/4
          endif
       enddo
       ! u5 = (5*u1/4 - u0/4) + 3dt/4 RHS(u4)
       !
       ! phl: rhs: t=t+2*dt/5+dt/3+3*dt/4         -wrong RK times ...
       !

       call compute_and_apply_rhs(np1,nm1,np1,qn0,3*dt/4,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,3*eta_ave_w/4)
       ! final method is the same as:
       ! u5 = u0 +  dt/4 RHS(u0)) + 3dt/4 RHS(u4)
#endif

    else if ((method==11).or.(method==12)) then
       ! Fully implicit JFNK method (vertically langragian not active yet)
       if (rsplit > 0) then
       call abortmp('ERROR: full_imp integration not yet coded for vert lagrangian adv option')
       end if
!      if (hybrid%masterthread) print*, "fully implicit integration is still under development"

    else
       call abortmp('ERROR: bad choice of tstep_type')
    endif

!    call prim_printstate(elem,tl,hybrid,hvcoord,nets,nete)

    ! ==============================================
    ! Time-split Horizontal diffusion: nu.del^2 or nu.del^4
    ! U(*) = U(t+1)  + dt2 * HYPER_DIFF_TERM(t+1)
    ! ==============================================
#ifdef ENERGY_DIAGNOSTICS
    if (compute_diagnostics) then
       do ie = nets,nete
          elem(ie)%accum%DIFF(:,:,:,:)=elem(ie)%state%v(:,:,:,:,np1)
          elem(ie)%accum%DIFFT(:,:,:)=elem(ie)%state%T(:,:,:,np1)
       enddo
    endif
#endif

    ! note:time step computes u(t+1)= u(t*) + RHS.
    ! for consistency, dt_vis = t-1 - t*, so this is timestep method dependent
    if (tstep_type==0) then
       ! leapfrog special case
       call advance_hypervis_lf(edge3p1,elem,hvcoord,hybrid,deriv,nm1,n0,np1,nets,nete,dt_vis)
    else
       if (rsplit==0) then
          ! forward-in-time, maybe hypervis applied to PS
          call advance_hypervis(edge3p1,elem,hvcoord,hybrid,deriv,np1,nets,nete,dt_vis,eta_ave_w)
       else
          ! forward-in-time, hypervis applied to dp3d
          call advance_hypervis_dp(edge3p1,elem,hvcoord,hybrid,deriv,np1,nets,nete,dt_vis,eta_ave_w)
       endif
    endif

#ifdef ENERGY_DIAGNOSTICS
    if (compute_diagnostics) then
       do ie = nets,nete
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
         do k=1,nlev  !  Loop index added (AAM)
          elem(ie)%accum%DIFF(:,:,:,k)=( elem(ie)%state%v(:,:,:,k,np1) -&
               elem(ie)%accum%DIFF(:,:,:,k) ) / dt_vis
          elem(ie)%accum%DIFFT(:,:,k)=( elem(ie)%state%T(:,:,k,np1) -&
               elem(ie)%accum%DIFFT(:,:,k) ) / dt_vis
         enddo
       enddo
    endif
#endif

    call t_stopf('prim_advance_exp')
    call t_adj_detailf(-1)
    end subroutine prim_advance_exp

  subroutine preq_robert3(nm1,n0,np1,elem,hvcoord,nets,nete)
  use dimensions_mod, only : np, nlev, qsize
  use hybvcoord_mod, only : hvcoord_t
  use element_mod, only : element_t
  use time_mod, only: smooth
  use control_mod, only : integration

  implicit none
  integer              , intent(in) :: nm1,n0,np1,nets,nete
  type (hvcoord_t), intent(in)      :: hvcoord
  type (element_t)     , intent(inout) :: elem(:)


  integer :: i,j,k,ie,q
  real (kind=real_kind) :: dp
  logical :: filter_ps = .false.
  if (integration == "explicit") filter_ps = .true.

  call t_adj_detailf(+1)
  call t_startf('preq_robert')
  do ie=nets,nete
     if (filter_ps) then
        elem(ie)%state%ps_v(:,:,n0) = elem(ie)%state%ps_v(:,:,n0) + smooth*(elem(ie)%state%ps_v(:,:,nm1) &
             - 2.0D0*elem(ie)%state%ps_v(:,:,n0)   + elem(ie)%state%ps_v(:,:,np1))
        elem(ie)%state%lnps(:,:,n0) = LOG(elem(ie)%state%ps_v(:,:,n0))
     else
        elem(ie)%state%lnps(:,:,n0) = elem(ie)%state%lnps(:,:,n0) + smooth*(elem(ie)%state%lnps(:,:,nm1) &
             - 2.0D0*elem(ie)%state%lnps(:,:,n0)   + elem(ie)%state%lnps(:,:,np1))
        elem(ie)%state%ps_v(:,:,n0) = EXP(elem(ie)%state%lnps(:,:,n0))
     endif

     elem(ie)%state%T(:,:,:,n0) = elem(ie)%state%T(:,:,:,n0) + smooth*(elem(ie)%state%T(:,:,:,nm1) &
          - 2.0D0*elem(ie)%state%T(:,:,:,n0)   + elem(ie)%state%T(:,:,:,np1))
     elem(ie)%state%v(:,:,:,:,n0) = elem(ie)%state%v(:,:,:,:,n0) + smooth*(elem(ie)%state%v(:,:,:,:,nm1) &
          - 2.0D0*elem(ie)%state%v(:,:,:,:,n0) + elem(ie)%state%v(:,:,:,:,np1))

  end do
  call t_stopf('preq_robert')
  call t_adj_detailf(-1)

  end subroutine preq_robert3




  subroutine applyCAMforcing(elem,fvm,hvcoord,np1,np1_qdp,dt_q,nets,nete)
  use dimensions_mod, only : np, nc, nlev, qsize, ntrac
  use element_mod, only : element_t
  use hybvcoord_mod, only : hvcoord_t
  use control_mod, only : moisture, tracer_grid_type
  use control_mod, only : TRACER_GRIDTYPE_GLL, TRACER_GRIDTYPE_FVM
  use physical_constants, only: Cp
  use fvm_control_volume_mod, only : fvm_struct
  implicit none
  type (element_t)     , intent(inout) :: elem(:)
  type(fvm_struct)     , intent(inout) :: fvm(:)
  real (kind=real_kind), intent(in) :: dt_q
  type (hvcoord_t), intent(in)      :: hvcoord
  integer,  intent(in) :: np1,nets,nete,np1_qdp

  ! local
  integer :: i,j,k,ie,q
  real (kind=real_kind) :: v1,dp
  real (kind=real_kind) :: beta(np,np),E0(np,np),ED(np,np),dp0m1(np,np),dpsum(np,np)
  logical :: wet

  wet = (moisture /= "dry")

  do ie=nets,nete
     ! apply forcing to Qdp
     elem(ie)%derived%FQps(:,:,1)=0
#if (defined COLUMN_OPENMP)
!$omp parallel do private(q,k,i,j,v1)
#endif
     if (tracer_grid_type == TRACER_GRIDTYPE_GLL) then
        do q=1,qsize
           do k=1,nlev
              do j=1,np
                 do i=1,np
                    v1 = dt_q*elem(ie)%derived%FQ(i,j,k,q,1)
                    !if (elem(ie)%state%Qdp(i,j,k,q,np1) + v1 < 0 .and. v1<0) then
                    if (elem(ie)%state%Qdp(i,j,k,q,np1_qdp) + v1 < 0 .and. v1<0) then
                       !if (elem(ie)%state%Qdp(i,j,k,q,np1) < 0 ) then
                       if (elem(ie)%state%Qdp(i,j,k,q,np1_qdp) < 0 ) then
                          v1=0  ! Q already negative, dont make it more so
                       else
                          !v1 = -elem(ie)%state%Qdp(i,j,k,q,np1)
                          v1 = -elem(ie)%state%Qdp(i,j,k,q,np1_qdp)
                       endif
                    endif
                    !elem(ie)%state%Qdp(i,j,k,q,np1) = elem(ie)%state%Qdp(i,j,k,q,np1)+v1
                    elem(ie)%state%Qdp(i,j,k,q,np1_qdp) = elem(ie)%state%Qdp(i,j,k,q,np1_qdp)+v1
                    if (q==1) then
                       elem(ie)%derived%FQps(i,j,1)=elem(ie)%derived%FQps(i,j,1)+v1/dt_q
                    endif
                 enddo
              enddo
           enddo
        enddo
     else if (tracer_grid_type == TRACER_GRIDTYPE_FVM) then
        ! Repeat for the fvm tracers
        do q = 1, ntrac
           do k = 1, nlev
              do j = 1, nc
                 do i = 1, nc
                    v1 = fvm(ie)%fc(i,j,k,q)
                    if (fvm(ie)%c(i,j,k,q,np1_qdp) + v1 < 0 .and. v1<0) then
                       if (fvm(ie)%c(i,j,k,q,np1_qdp) < 0 ) then
                          v1 = 0  ! C already negative, dont make it more so
                       else
                          v1 = -fvm(ie)%c(i,j,k,q,np1_qdp)
                       end if
                    end if
                    fvm(ie)%c(i,j,k,q,np1_qdp) = fvm(ie)%c(i,j,k,q,np1_qdp) + v1
!                    if (q == 1) then
!!XXgoldyXX: Should update the pressure forcing here??!!??
!                    elem(ie)%derived%FQps(i,j,1)=elem(ie)%derived%FQps(i,j,1)+v1/dt_q
!                  end if
                 end do
              end do
           end do
        end do
     else
        call abortmp('applyCAMforcing: ERROR: unsupported value for tracer_grid_type')
     end if

     if (wet .and. qsize>0) then
        ! to conserve dry mass in the precese of Q1 forcing:
        elem(ie)%state%ps_v(:,:,np1) = elem(ie)%state%ps_v(:,:,np1) + &
             dt_q*elem(ie)%derived%FQps(:,:,1)
     endif

#if 0
     ! disabled - energy fixers will be moving into CAM physics
     ! energy fixer for FQps term
     ! dp1 = dp0 + d(FQps)
     ! dp0-dp1 = -d(FQps)
     ! E0-E1 = sum( dp0*ED) - sum( dp1*ED) = sum( dp0-dp1) * ED )
     ! compute E0-E1
     E0=0
     do k=1,nlev
        ED(:,:) = ( 0.5d0* &
             (elem(ie)%state%v(:,:,1,k,np1)**2 + elem(ie)%state%v(:,:,2,k,np1)**2)&
             + cp*elem(ie)%state%T(:,:,k,np1)  &
             + elem(ie)%state%phis(:,:) )

        dp0m1(:,:) = -dt_q*( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%derived%FQps(:,:,1)

        E0(:,:) = E0(:,:) + dp0m1(:,:)*ED(:,:)
     enddo
     ! energy fixer:
     ! Tnew = T + beta
     ! cp*dp*beta  = E0-E1   beta = (E0-E1)/(cp*sum(dp))

     dpsum(:,:) = ( hvcoord%hyai(nlev+1) - hvcoord%hyai(1) )*hvcoord%ps0 + &
          ( hvcoord%hybi(nlev+1) - hvcoord%hybi(1) )*elem(ie)%state%ps_v(:,:,np1)

     beta(:,:)=E0(:,:)/(dpsum(:,:)*cp)
     do k=1,nlev
        elem(ie)%state%T(:,:,k,np1)=elem(ie)%state%T(:,:,k,np1)+beta(:,:)
     enddo
#endif




     ! Qdp(np1) and ps_v(np1) were updated by forcing - update Q(np1)
#if (defined COLUMN_OPENMP)
!$omp parallel do private(q,k,i,j,dp)
#endif
     do q=1,qsize
        do k=1,nlev
           do j=1,np
              do i=1,np
                 dp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                      ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,np1)
                 elem(ie)%state%Q(i,j,k,q) = elem(ie)%state%Qdp(i,j,k,q,np1_qdp)/dp
              enddo
           enddo
        enddo
     enddo

     elem(ie)%state%T(:,:,:,np1) = elem(ie)%state%T(:,:,:,np1) + &
          dt_q*elem(ie)%derived%FT(:,:,:,1)
     elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,np1) + &
          dt_q*elem(ie)%derived%FM(:,:,:,:,1)
  enddo
  end subroutine applyCAMforcing



  subroutine applyCAMforcing_dynamics(elem,hvcoord,np1,dt_q,nets,nete)
  use dimensions_mod, only : np, nlev, qsize
  use element_mod, only : element_t
  use hybvcoord_mod, only : hvcoord_t

  implicit none
  type (element_t)     , intent(inout) :: elem(:)
  real (kind=real_kind), intent(in) :: dt_q
  type (hvcoord_t), intent(in)      :: hvcoord
  integer,  intent(in) :: np1,nets,nete

  ! local
  integer :: i,j,k,ie,q
  real (kind=real_kind) :: v1,dp
  logical :: wet

  do ie=nets,nete
     elem(ie)%state%T(:,:,:,np1) = elem(ie)%state%T(:,:,:,np1) + &
          dt_q*elem(ie)%derived%FT(:,:,:,1)
     elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,np1) + &
          dt_q*elem(ie)%derived%FM(:,:,:,:,1)
  enddo
  end subroutine applyCAMforcing_dynamics



  subroutine advance_hypervis(edge3,elem,hvcoord,hybrid,deriv,nt,nets,nete,dt2,eta_ave_w)
  !
  !  take one timestep of:
  !          u(:,:,:,np) = u(:,:,:,np) +  dt2*nu*laplacian**order ( u )
  !          T(:,:,:,np) = T(:,:,:,np) +  dt2*nu_s*laplacian**order ( T )
  !
  !
  !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
  !
  !
  use dimensions_mod, only : np, np, nlev
  use control_mod, only : nu, nu_div, nu_s, hypervis_order, hypervis_subcycle, nu_p, nu_top, psurf_vis
  use hybrid_mod, only : hybrid_t
  use hybvcoord_mod, only : hvcoord_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t, laplace_sphere_wk, vlaplace_sphere_wk
  use edgetype_mod, only : EdgeBuffer_t
  use edge_mod, only : edgevpack, edgevunpack
  use bndry_mod, only : bndry_exchangev
  use viscosity_mod, only : biharmonic_wk
  use physical_constants, only: Cp
!  use time_mod, only : TimeLevel_t
  implicit none

  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (EdgeBuffer_t)  , intent(inout) :: edge3
  type (derivative_t)  , intent(in) :: deriv
  type (hvcoord_t), intent(in)      :: hvcoord
!  type (TimeLevel_t)   , intent(in) :: tl

  real (kind=real_kind) :: dt2
  integer :: nets,nete

  ! local
  real (kind=real_kind) :: eta_ave_w  ! weighting for mean flux terms
  real (kind=real_kind) :: nu_scale, dpdn,dpdn0, nu_scale_top
  integer :: k,kptr,i,j,ie,ic,nt
  real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)      :: vtens
  real (kind=real_kind), dimension(np,np,nlev,nets:nete)        :: ptens
  real (kind=real_kind), dimension(np,np,nets:nete) :: pstens
  real (kind=real_kind), dimension(np,np,nlev) :: p
  real (kind=real_kind), dimension(np,np) :: dptemp1,dptemp2


! NOTE: PGI compiler bug: when using spheremp, rspheremp and ps as pointers to elem(ie)% members,
  !       data is incorrect (offset by a few numbers actually)
  !       removed for now.
  !       real (kind=real_kind), dimension(:,:), pointer :: spheremp,rspheremp
  !       real (kind=real_kind), dimension(:,:,:), pointer   :: ps

  real (kind=real_kind), dimension(np,np) :: lap_p
  real (kind=real_kind), dimension(np,np,2) :: lap_v
  real (kind=real_kind) :: v1,v2,dt,heating,utens_tmp,vtens_tmp,ptens_tmp


  if (nu_s == 0 .and. nu == 0 .and. nu_p==0 ) return;
!JMD  call t_barrierf('sync_advance_hypervis', hybrid%par%comm)
  call t_adj_detailf(+1)
  call t_startf('advance_hypervis')


  dt=dt2/hypervis_subcycle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  regular viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hypervis_order == 1) then
     if (nu_p>0) call abortmp( 'ERROR: hypervis_order == 1 not coded for nu_p>0')
     do ic=1,hypervis_subcycle
        do ie=nets,nete

#if (defined COLUMN_OPENMP)
! Not sure about deriv here
!$omp parallel do private(k,lap_p,lap_v,deriv,i,j)
#endif
           do k=1,nlev
              ! lap_p=laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              call laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),lap_p,var_coef=.false.)
              ! lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              call vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),lap_v,var_coef=.false.)
              ! advance in time.  (note: DSS commutes with time stepping, so we
              ! can time advance and then DSS.  this has the advantage of
              ! not letting any discontinuties accumulate in p,v via roundoff
              do j=1,np
                 do i=1,np
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt)*elem(ie)%spheremp(i,j)  +  dt*nu_s*lap_p(i,j)
                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%state%v(i,j,1,k,nt)*elem(ie)%spheremp(i,j) + dt*nu*lap_v(i,j,1)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%state%v(i,j,2,k,nt)*elem(ie)%spheremp(i,j) + dt*nu*lap_v(i,j,2)
                 enddo
              enddo
           enddo

           kptr=0
           call edgeVpack(edge3, elem(ie)%state%T(:,:,:,nt),nlev,kptr,ie)
           kptr=nlev
           call edgeVpack(edge3,elem(ie)%state%v(:,:,:,:,nt),2*nlev,kptr,ie)
        enddo

        call bndry_exchangeV(hybrid,edge3)

        do ie=nets,nete

           kptr=0
           call edgeVunpack(edge3, elem(ie)%state%T(:,:,:,nt), nlev, kptr, ie)
           kptr=nlev
           call edgeVunpack(edge3, elem(ie)%state%v(:,:,:,:,nt), 2*nlev, kptr, ie)

           ! apply inverse mass matrix
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j)
#endif
           do k=1,nlev
              do j=1,np
                 do i=1,np
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%T(i,j,k,nt)
                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,1,k,nt)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,2,k,nt)
                 enddo
              enddo
           enddo
        enddo
     enddo  ! subcycle
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  hyper viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! nu_p=0:
!   scale T dissipaton by dp  (conserve IE, dissipate T^2)
! nu_p>0
!   dont scale:  T equation IE dissipation matches (to truncation error)
!                IE dissipation from continuity equation
!                (1 deg: to about 0.1 W/m^2)
!
  if (hypervis_order == 2) then
     do ic=1,hypervis_subcycle
        call biharmonic_wk(elem,pstens,ptens,vtens,deriv,edge3,hybrid,nt,nets,nete)
        do ie=nets,nete

           ! comptue mean flux
           if (nu_p>0) then
#if 0
              elem(ie)%derived%psdiss_ave(:,:)=&
                   elem(ie)%derived%psdiss_ave(:,:)+eta_ave_w*elem(ie)%state%ps_v(:,:,nt)/hypervis_subcycle
              elem(ie)%derived%psdiss_biharmonic(:,:)=&
                   elem(ie)%derived%psdiss_biharmonic(:,:)+eta_ave_w*pstens(:,:,ie)/hypervis_subcycle
#else
              do k=1,nlev
                 dptemp1(:,:) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                      ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,nt)
                 elem(ie)%derived%dpdiss_ave(:,:,k)=elem(ie)%derived%dpdiss_ave(:,:,k)+eta_ave_w*dptemp1(:,:)/hypervis_subcycle

                 dptemp2(:,:) = (hvcoord%hybi(k+1)-hvcoord%hybi(k))*pstens(:,:,ie)
                 elem(ie)%derived%dpdiss_biharmonic(:,:,k)=&
                      elem(ie)%derived%dpdiss_biharmonic(:,:,k)+eta_ave_w*dptemp2(:,:)/hypervis_subcycle
              enddo
#endif
           endif
           nu_scale=1
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,lap_p,lap_v,nu_scale_top,dpdn,dpdn0,nu_scale,utens_tmp,vtens_tmp,ptens_tmp)
#endif
           do k=1,nlev
              ! advace in time.
              ! note: DSS commutes with time stepping, so we can time advance and then DSS.
              ! note: weak operators alreayd have mass matrix "included"

              ! add regular diffusion in top 3 layers:
              if (nu_top>0 .and. k<=3) then
                 ! lap_p=laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                 call laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),lap_p,var_coef=.false.)
                 !lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                 call vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),lap_v,var_coef=.false.)
              endif
              nu_scale_top = 1
              if (k==1) nu_scale_top=4
              if (k==2) nu_scale_top=2

              do j=1,np
                 do i=1,np
                    if (nu_p==0) then
                       ! normalize so as to conserve IE
                       ! scale by 1/rho (normalized to be O(1))
                       ! dp/dn = O(ps0)*O(delta_eta) = O(ps0)/O(nlev)
                       dpdn = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                            ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,nt)
                       dpdn0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                            ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
                       nu_scale = dpdn0/dpdn
                    endif

                    ! biharmonic terms need a negative sign:
                    if (nu_top>0 .and. k<=3) then
                       utens_tmp=(-nu*vtens(i,j,1,k,ie) + nu_scale_top*nu_top*lap_v(i,j,1))
                       vtens_tmp=(-nu*vtens(i,j,2,k,ie) + nu_scale_top*nu_top*lap_v(i,j,2))
                       ptens_tmp=nu_scale*(-nu_s*ptens(i,j,k,ie) + nu_scale_top*nu_top*lap_p(i,j) )
                    else
                       utens_tmp=-nu*vtens(i,j,1,k,ie)
                       vtens_tmp=-nu*vtens(i,j,2,k,ie)
                       ptens_tmp=-nu_scale*nu_s*ptens(i,j,k,ie)
                    endif

                    ptens(i,j,k,ie) = ptens_tmp
                    vtens(i,j,1,k,ie)=utens_tmp
                    vtens(i,j,2,k,ie)=vtens_tmp
                 enddo
              enddo
           enddo

           pstens(:,:,ie)  =  -nu_p*pstens(:,:,ie)
           kptr=0
           call edgeVpack(edge3, ptens(:,:,:,ie),nlev,kptr, ie)
           kptr=nlev
           call edgeVpack(edge3,vtens(:,:,:,:,ie),2*nlev,kptr,ie)
           kptr=3*nlev
           call edgeVpack(edge3,pstens(:,:,ie),1,kptr,ie)
        enddo


        call bndry_exchangeV(hybrid,edge3)

        do ie=nets,nete

           kptr=0
           call edgeVunpack(edge3, ptens(:,:,:,ie), nlev, kptr, ie)
           kptr=nlev
           call edgeVunpack(edge3, vtens(:,:,:,:,ie), 2*nlev, kptr, ie)


           ! apply inverse mass matrix, accumulate tendencies
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
           do k=1,nlev
              vtens(:,:,1,k,ie)=dt*vtens(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)
              vtens(:,:,2,k,ie)=dt*vtens(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)
              ptens(:,:,k,ie)=dt*ptens(:,:,k,ie)*elem(ie)%rspheremp(:,:)
           enddo

           ! apply hypervis to u -> u+utens:
           ! E0 = dpdn * .5*u dot u + dpdn * T  + dpdn*PHIS
           ! E1 = dpdn * .5*(u+utens) dot (u+utens) + dpdn * (T-X) + dpdn*PHIS
           ! E1-E0:   dpdn (u dot utens) + dpdn .5 utens dot utens   - dpdn X
           !      X = (u dot utens) + .5 utens dot utens
           !  alt:  (u+utens) dot utens
           do k=1,nlev
              do j=1,np
                 do i=1,np
                    ! update v first (gives better results than updating v after heating)
                    elem(ie)%state%v(i,j,:,k,nt)=elem(ie)%state%v(i,j,:,k,nt) + &
                         vtens(i,j,:,k,ie)

                    v1=elem(ie)%state%v(i,j,1,k,nt)
                    v2=elem(ie)%state%v(i,j,2,k,nt)
                    heating = (vtens(i,j,1,k,ie)*v1  + vtens(i,j,2,k,ie)*v2 )
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt) &
                         +ptens(i,j,k,ie)-heating/cp

                 enddo
              enddo
           enddo

           if (nu_p>0) then
              kptr=3*nlev
              call edgeVunpack(edge3, pstens(:,:,ie), 1, kptr, ie)
              pstens(:,:,ie)=dt*pstens(:,:,ie)*elem(ie)%rspheremp(:,:)
              elem(ie)%state%ps_v(:,:,nt)=elem(ie)%state%ps_v(:,:,nt) + pstens(:,:,ie)
           endif

        enddo
     enddo
  endif

  call t_stopf('advance_hypervis')
  call t_adj_detailf(-1)

  end subroutine advance_hypervis




  subroutine advance_hypervis_dp(edge3,elem,hvcoord,hybrid,deriv,nt,nets,nete,dt2,eta_ave_w)
  !
  !  take one timestep of:
  !          u(:,:,:,np) = u(:,:,:,np) +  dt2*nu*laplacian**order ( u )
  !          T(:,:,:,np) = T(:,:,:,np) +  dt2*nu_s*laplacian**order ( T )
  !
  !
  !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
  !
  !
  use dimensions_mod, only : np, np, nlev
  use control_mod, only : nu, nu_div, nu_s, hypervis_order, hypervis_subcycle, nu_p, nu_top, psurf_vis
  use hybrid_mod, only : hybrid_t
  use hybvcoord_mod, only : hvcoord_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t, laplace_sphere_wk, vlaplace_sphere_wk
  use edgetype_mod, only : EdgeBuffer_t
  use edge_mod, only : edgevpack, edgevunpack
  use bndry_mod, only : bndry_exchangev
  use viscosity_mod, only : biharmonic_wk_dp3d
  use physical_constants, only: Cp
!  use time_mod, only : TimeLevel_t
  implicit none

  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (EdgeBuffer_t)  , intent(inout) :: edge3
  type (derivative_t)  , intent(in) :: deriv
  type (hvcoord_t), intent(in)      :: hvcoord
!  type (TimeLevel_t)   , intent(in) :: tl

  real (kind=real_kind) :: dt2
  integer :: nets,nete

  ! local
  real (kind=real_kind) :: eta_ave_w  ! weighting for mean flux terms
  real (kind=real_kind) :: dpdn,dpdn0, nu_scale_top
  integer :: k,kptr,i,j,ie,ic,nt
  real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)      :: vtens
  real (kind=real_kind), dimension(np,np,nlev,nets:nete)        :: ttens
  real (kind=real_kind), dimension(np,np,nlev,nets:nete)        :: dptens
  real (kind=real_kind), dimension(np,np,nlev) :: p
  real (kind=real_kind), dimension(np,np) :: dptemp1,dptemp2


! NOTE: PGI compiler bug: when using spheremp, rspheremp and ps as pointers to elem(ie)% members,
  !       data is incorrect (offset by a few numbers actually)
  !       removed for now.
  !       real (kind=real_kind), dimension(:,:), pointer :: spheremp,rspheremp
  !       real (kind=real_kind), dimension(:,:,:), pointer   :: ps

  real (kind=real_kind), dimension(np,np) :: lap_t,lap_dp
  real (kind=real_kind), dimension(np,np,2) :: lap_v
  real (kind=real_kind) :: v1,v2,dt,heating,utens_tmp,vtens_tmp,ttens_tmp,dptens_tmp


  if (nu_s == 0 .and. nu == 0 .and. nu_p==0 ) return;
!JMD  call t_barrierf('sync_advance_hypervis', hybrid%par%comm)
  call t_adj_detailf(+1) 
  call t_startf('advance_hypervis_dp')


  dt=dt2/hypervis_subcycle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  regular viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hypervis_order == 1) then
     if (nu_p>0) call abortmp( 'ERROR: hypervis_order == 1 not coded for nu_p>0')
     do ic=1,hypervis_subcycle
        do ie=nets,nete

#if (defined COLUMN_OPENMP)
! Not sure about deriv here
!$omp parallel do private(k,lap_t,lap_v,deriv,i,j)
#endif
           do k=1,nlev
              ! lap_t=laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              call laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),lap_t,var_coef=.false.)
              !lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              call vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),lap_v,var_coef=.false.)
              ! advace in time.  (note: DSS commutes with time stepping, so we
              ! can time advance and then DSS.  this has the advantage of
              ! not letting any discontinuties accumulate in p,v via roundoff
              do j=1,np
                 do i=1,np
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt)*elem(ie)%spheremp(i,j)  +  dt*nu_s*lap_t(i,j)
                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%state%v(i,j,1,k,nt)*elem(ie)%spheremp(i,j) + dt*nu*lap_v(i,j,1)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%state%v(i,j,2,k,nt)*elem(ie)%spheremp(i,j) + dt*nu*lap_v(i,j,2)
                 enddo
              enddo
           enddo

           kptr=0
           call edgeVpack(edge3, elem(ie)%state%T(:,:,:,nt),nlev,kptr,ie)
           kptr=nlev
           call edgeVpack(edge3,elem(ie)%state%v(:,:,:,:,nt),2*nlev,kptr,ie)
        enddo

        call bndry_exchangeV(hybrid,edge3)

        do ie=nets,nete

           kptr=0
           call edgeVunpack(edge3, elem(ie)%state%T(:,:,:,nt), nlev, kptr, ie)
           kptr=nlev
           call edgeVunpack(edge3, elem(ie)%state%v(:,:,:,:,nt), 2*nlev, kptr, ie)

           ! apply inverse mass matrix
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j)
#endif
           do k=1,nlev
              do j=1,np
                 do i=1,np
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%T(i,j,k,nt)
                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,1,k,nt)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,2,k,nt)
                 enddo
              enddo
           enddo
        enddo
     enddo  ! subcycle
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  hyper viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! nu_p=0:
!   scale T dissipaton by dp  (conserve IE, dissipate T^2)
! nu_p>0
!   dont scale:  T equation IE dissipation matches (to truncation error)
!                IE dissipation from continuity equation
!                (1 deg: to about 0.1 W/m^2)
!
  if (hypervis_order == 2) then
     do ic=1,hypervis_subcycle
        call biharmonic_wk_dp3d(elem,dptens,ttens,vtens,deriv,edge3,hybrid,nt,nets,nete)
        do ie=nets,nete

           ! comptue mean flux
           if (nu_p>0) then
              do k=1,nlev
                !OMP_COLLAPSE_SIMD
                !DIR_VECTOR_ALIGNED
                do j=1,np
                do i=1,np
              elem(ie)%derived%dpdiss_ave(i,j,k)=elem(ie)%derived%dpdiss_ave(i,j,k)+&
                   eta_ave_w*elem(ie)%state%dp3d(i,j,k,nt)/hypervis_subcycle
              elem(ie)%derived%dpdiss_biharmonic(i,j,k)=elem(ie)%derived%dpdiss_biharmonic(i,j,k)+&
                   eta_ave_w*dptens(i,j,k,ie)/hypervis_subcycle
                enddo
              enddo
              enddo
           endif
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,lap_t,lap_dp,lap_v,nu_scale_top,dpdn,dpdn0,utens_tmp,vtens_tmp,ttens_tmp,dptens_tmp)
#endif
           do k=1,nlev
              ! advace in time.
              ! note: DSS commutes with time stepping, so we can time advance and then DSS.
              ! note: weak operators alreayd have mass matrix "included"

              ! add regular diffusion in top 3 layers:
              if (nu_top>0 .and. k<=3) then
                 ! lap_t=laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                 call laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),lap_t,var_coef=.false.)
                 ! lap_dp=laplace_sphere_wk(elem(ie)%state%dp3d(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                 call laplace_sphere_wk(elem(ie)%state%dp3d(:,:,k,nt),deriv,elem(ie),lap_dp, var_coef=.false.)
                 ! lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                 call vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),lap_v, var_coef=.false.)
              endif
              nu_scale_top = 1
              if (k==1) nu_scale_top=4
              if (k==2) nu_scale_top=2

              ! biharmonic terms need a negative sign:
              if (nu_top>0 .and. k<=3) then
                 !OMP_COLLAPSE_SIMD
                 !DIR_VECTOR_ALIGNED
                 do j=1,np
                    do i=1,np
                       ttens(i,j,k,ie)=(-nu_s*ttens(i,j,k,ie) + nu_scale_top*nu_top*lap_t(i,j) )
                       dptens(i,j,k,ie)=(-nu_p*dptens(i,j,k,ie) + nu_scale_top*nu_top*lap_dp(i,j) )
                       vtens(i,j,1,k,ie)=(-nu*vtens(i,j,1,k,ie) + nu_scale_top*nu_top*lap_v(i,j,1))
                       vtens(i,j,2,k,ie)=(-nu*vtens(i,j,2,k,ie) + nu_scale_top*nu_top*lap_v(i,j,2))
                    enddo
                 enddo
              else
                 !OMP_COLLAPSE_SIMD
                 !DIR_VECTOR_ALIGNED
                 do j=1,np
                    do i=1,np
                       ttens(i,j,k,ie)=-nu_s*ttens(i,j,k,ie)
                       dptens(i,j,k,ie)=-nu_p*dptens(i,j,k,ie)
                       vtens(i,j,1,k,ie)=-nu*vtens(i,j,1,k,ie)
                       vtens(i,j,2,k,ie)=-nu*vtens(i,j,2,k,ie)
                    enddo
                 enddo
              endif
           enddo


           kptr=0
           call edgeVpack(edge3, ttens(:,:,:,ie),nlev,kptr,ie)
           kptr=nlev
           call edgeVpack(edge3,vtens(:,:,:,:,ie),2*nlev,kptr,ie)
           kptr=3*nlev
           call edgeVpack(edge3,dptens(:,:,:,ie),nlev,kptr, ie)
        enddo


        call bndry_exchangeV(hybrid,edge3)

        do ie=nets,nete

           kptr=0
           call edgeVunpack(edge3, ttens(:,:,:,ie), nlev, kptr, ie)
           kptr=nlev
           call edgeVunpack(edge3, vtens(:,:,:,:,ie), 2*nlev, kptr, ie)
           kptr=3*nlev
           call edgeVunpack(edge3, dptens(:,:,:,ie), nlev, kptr, ie)


           ! apply inverse mass matrix, accumulate tendencies
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
           do k=1,nlev
              vtens(:,:,1,k,ie)=dt*vtens(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)
              vtens(:,:,2,k,ie)=dt*vtens(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)
              ttens(:,:,k,ie)=dt*ttens(:,:,k,ie)*elem(ie)%rspheremp(:,:)
              dptens(:,:,k,ie)=dt*dptens(:,:,k,ie)*elem(ie)%rspheremp(:,:)
           enddo

           ! apply hypervis to u -> u+utens:
           ! E0 = dpdn * .5*u dot u + dpdn * T  + dpdn*PHIS
           ! E1 = dpdn * .5*(u+utens) dot (u+utens) + dpdn * (T-X) + dpdn*PHIS
           ! E1-E0:   dpdn (u dot utens) + dpdn .5 utens dot utens   - dpdn X
           !      X = (u dot utens) + .5 utens dot utens
           !  alt:  (u+utens) dot utens
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,v1,v2,heating)
#endif
           do k=1,nlev
              do j=1,np
                 do i=1,np
                    ! update v first (gives better results than updating v after heating)
                    elem(ie)%state%v(i,j,:,k,nt)=elem(ie)%state%v(i,j,:,k,nt) + &
                         vtens(i,j,:,k,ie)

                    v1=elem(ie)%state%v(i,j,1,k,nt)
                    v2=elem(ie)%state%v(i,j,2,k,nt)
                    heating = (vtens(i,j,1,k,ie)*v1  + vtens(i,j,2,k,ie)*v2 )
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt) &
                         +ttens(i,j,k,ie)-heating/cp
                    elem(ie)%state%dp3d(i,j,k,nt)=elem(ie)%state%dp3d(i,j,k,nt) + &
                         dptens(i,j,k,ie)
                 enddo
              enddo
           enddo
        enddo
     enddo
  endif

  call t_stopf('advance_hypervis_dp')
  call t_adj_detailf(-1)

  end subroutine advance_hypervis_dp






  subroutine advance_hypervis_lf(edge3,elem,hvcoord,hybrid,deriv,nm1,n0,nt,nets,nete,dt2)
  !
  !  take one timestep of:
  !          u(:,:,:,np) = u(:,:,:,np) +  dt2*nu*laplacian**order ( u )
  !          T(:,:,:,np) = T(:,:,:,np) +  dt2*nu_s*laplacian**order ( T )
  !
  !
  !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
  !
  !
  use dimensions_mod, only : np, np, nlev
  use control_mod, only : nu, nu_div, nu_s, hypervis_order, hypervis_subcycle, nu_p, nu_top, psurf_vis
  use hybrid_mod, only : hybrid_t
  use hybvcoord_mod, only : hvcoord_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t, laplace_sphere_wk, vlaplace_sphere_wk
  use edgetype_mod, only : EdgeBuffer_t
  use edge_mod, only : edgevpack, edgevunpack
  use bndry_mod, only : bndry_exchangev
  use viscosity_mod, only : biharmonic_wk
  use physical_constants, only: Cp
!  use time_mod, only : TimeLevel_t
  implicit none

  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (EdgeBuffer_t)  , intent(inout) :: edge3
  type (derivative_t)  , intent(in) :: deriv
  type (hvcoord_t), intent(in)      :: hvcoord
!  type (TimeLevel_t)   , intent(in) :: tl

  real (kind=real_kind) :: dt2
  integer :: nets,nete

  ! local
  real (kind=real_kind) :: nu_scale, dpdn,dpdn0, nu_scale_top
  integer :: k,kptr,i,j,ie,ic,n0,nt,nm1
  real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)      :: vtens
  real (kind=real_kind), dimension(np,np,nlev,nets:nete)        :: ptens
  real (kind=real_kind), dimension(np,np,nets:nete) :: pstens
  real (kind=real_kind), dimension(np,np,nlev) :: p
  real (kind=real_kind), dimension(np,np) :: dXdp


! NOTE: PGI compiler bug: when using spheremp, rspheremp and ps as pointers to elem(ie)% members,
  !       data is incorrect (offset by a few numbers actually)
  !       removed for now.
  !       real (kind=real_kind), dimension(:,:), pointer :: spheremp,rspheremp
  !       real (kind=real_kind), dimension(:,:,:), pointer   :: ps

  real (kind=real_kind), dimension(np,np) :: lap_p
  real (kind=real_kind), dimension(np,np,2) :: lap_v
  real (kind=real_kind) :: v1,v2,dt,heating,utens_tmp,vtens_tmp,ptens_tmp


  if (nu_s == 0 .and. nu == 0 .and. nu_p==0 ) return;
!JMD  call t_barrierf('sync_advance_hypervis_lf', hybrid%par%comm)
  call t_adj_detailf(+1)
  call t_startf('advance_hypervis_lf')

! for non-leapfrog,nt=n0=nmt
!
!  nm1 = tl%nm1   ! heating term uses U,V at average of nt and nm1 levels
!  n0 = tl%n0     ! timelevel used for ps scaling.  use n0 for leapfrog.
!  nt = tl%np1    ! apply viscosity to this timelevel  (np1)


  dt=dt2/hypervis_subcycle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  regular viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hypervis_order == 1) then
     if (nu_p>0) stop 'ERROR: hypervis_order == 1 not coded for nu_p>0'
     do ic=1,hypervis_subcycle
        do ie=nets,nete

#if (defined COLUMN_OPENMP)
! Not sure about deriv here
!$omp parallel do private(k,lap_p,lap_v,deriv,i,j)
#endif
           do k=1,nlev
              !lap_p=laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              call laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),lap_p,var_coef=.false.)
              !lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              call vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),lap_v,var_coef=.false.)
              ! advace in time.  (note: DSS commutes with time stepping, so we
              ! can time advance and then DSS.  this has the advantage of
              ! not letting any discontinuties accumulate in p,v via roundoff
              do j=1,np
                 do i=1,np
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt)*elem(ie)%spheremp(i,j)  +  dt*nu_s*lap_p(i,j)
                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%state%v(i,j,1,k,nt)*elem(ie)%spheremp(i,j) + dt*nu*lap_v(i,j,1)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%state%v(i,j,2,k,nt)*elem(ie)%spheremp(i,j) + dt*nu*lap_v(i,j,2)
                 enddo
              enddo
           enddo

           kptr=0
           call edgeVpack(edge3, elem(ie)%state%T(:,:,:,nt),nlev,kptr,ie)
           kptr=nlev
           call edgeVpack(edge3,elem(ie)%state%v(:,:,:,:,nt),2*nlev,kptr,ie)
        enddo

        call bndry_exchangeV(hybrid,edge3)

        do ie=nets,nete

           kptr=0
           call edgeVunpack(edge3, elem(ie)%state%T(:,:,:,nt), nlev, kptr, ie)
           kptr=nlev
           call edgeVunpack(edge3, elem(ie)%state%v(:,:,:,:,nt), 2*nlev, kptr, ie)

           ! apply inverse mass matrix
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j)
#endif
           do k=1,nlev
              do j=1,np
                 do i=1,np
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%T(i,j,k,nt)
                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,1,k,nt)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,2,k,nt)
                 enddo
              enddo
           enddo
        enddo
     enddo  ! subcycle
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  hyper viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hypervis_order == 2) then
     do ic=1,hypervis_subcycle
        call biharmonic_wk(elem,pstens,ptens,vtens,deriv,edge3,hybrid,nt,nets,nete)
        do ie=nets,nete

           nu_scale=1
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,lap_p,lap_v,nu_scale_top,dpdn,dpdn0,nu_scale,utens_tmp,vtens_tmp,ptens_tmp)
#endif
           do k=1,nlev
              ! advace in time.
              ! note: DSS commutes with time stepping, so we can time advance and then DSS.
              ! note: weak operators alreayd have mass matrix "included"

              ! add regular diffusion in top 3 layers:
              if (nu_top>0 .and. k<=3) then
                 ! lap_p=laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                 call laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),lap_p, var_coef=.false.)
                 ! lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                 call vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),lap_v, var_coef=.false.)
              endif
              nu_scale_top = 1
              if (k==1) nu_scale_top=4
              if (k==2) nu_scale_top=2

              do j=1,np
                 do i=1,np
                    if (psurf_vis==0) then
                       ! normalize so as to conserve IE  (not needed when using p-surface viscosity)
                       ! scale velosity by 1/rho (normalized to be O(1))
                       ! dp/dn = O(ps0)*O(delta_eta) = O(ps0)/O(nlev)
                       dpdn = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                            ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,n0)  ! nt ?
                       dpdn0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                            ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
                       nu_scale = dpdn0/dpdn
                    endif

                    ! biharmonic terms need a negative sign:
                    if (nu_top>0 .and. k<=3) then
                       utens_tmp=(-nu*vtens(i,j,1,k,ie) + nu_scale_top*nu_top*lap_v(i,j,1))
                       vtens_tmp=(-nu*vtens(i,j,2,k,ie) + nu_scale_top*nu_top*lap_v(i,j,2))
                       ptens_tmp=nu_scale*(-nu_s*ptens(i,j,k,ie) + nu_scale_top*nu_top*lap_p(i,j) )
                    else
                       utens_tmp=-nu*vtens(i,j,1,k,ie)
                       vtens_tmp=-nu*vtens(i,j,2,k,ie)
                       ptens_tmp=-nu_scale*nu_s*ptens(i,j,k,ie)
                    endif

                    ptens(i,j,k,ie) = ptens_tmp
                    vtens(i,j,1,k,ie)=utens_tmp
                    vtens(i,j,2,k,ie)=vtens_tmp
                 enddo
              enddo
           enddo

           pstens(:,:,ie)  =  -nu_p*pstens(:,:,ie)
           kptr=0
           call edgeVpack(edge3, ptens(:,:,:,ie),nlev,kptr,ie)
           kptr=nlev
           call edgeVpack(edge3,vtens(:,:,:,:,ie),2*nlev,kptr,ie)
           kptr=3*nlev
           call edgeVpack(edge3,pstens(:,:,ie),1,kptr,ie)
        enddo


        call bndry_exchangeV(hybrid,edge3)

        do ie=nets,nete

           kptr=0
           call edgeVunpack(edge3, ptens(:,:,:,ie), nlev, kptr, ie)
           kptr=nlev
           call edgeVunpack(edge3, vtens(:,:,:,:,ie), 2*nlev, kptr, ie)
           kptr=3*nlev
           call edgeVunpack(edge3, pstens(:,:,ie), 1, kptr, ie)

           if (psurf_vis == 1 ) then
              ! apply p-surface correction
              do k=1,nlev
                 p(:,:,k)   = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%ps_v(:,:,nt)
              enddo
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,dXdp)
#endif
              do k=1,nlev
                 if (k.eq.1) then
                    ! no correction needed
                 else if (k.eq.nlev) then
                    ! one-sided difference
                    dXdp = (elem(ie)%state%T(:,:,k,nt) - elem(ie)%state%T(:,:,k-1,nt)) / &
                        (p(:,:,k)-p(:,:,k-1))
                    ptens(:,:,k,ie) = ptens(:,:,k,ie) - dXdp(:,:)*hvcoord%hybm(k)*pstens(:,:,ie)
                 else
                    dXdp = (elem(ie)%state%T(:,:,k+1,nt) - elem(ie)%state%T(:,:,k-1,nt)) / &
                         (p(:,:,k+1)-p(:,:,k-1))
                    ptens(:,:,k,ie) = ptens(:,:,k,ie) - dXdp(:,:)*hvcoord%hybm(k)*pstens(:,:,ie)
                 endif
              enddo
           endif


           ! apply inverse mass matrix, accumulate tendencies
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,v1,v2,heating)
#endif
           do k=1,nlev
              do j=1,np
                 do i=1,np

                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%state%v(i,j,1,k,nt) + &
                         dt*elem(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%state%v(i,j,2,k,nt) +  &
                         dt*elem(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)

                    ! better E conservation if we use v after adding in vtens:
                    v1=.5*(elem(ie)%state%v(i,j,1,k,nt)+elem(ie)%state%v(i,j,1,k,nm1))
                    v2=.5*(elem(ie)%state%v(i,j,2,k,nt)+elem(ie)%state%v(i,j,2,k,nm1))
                    heating = (vtens(i,j,1,k,ie)*v1  + vtens(i,j,2,k,ie)*v2 )

                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt)     + &
                         dt*elem(ie)%rspheremp(i,j)*(cp*ptens(i,j,k,ie) - heating)/cp

                 enddo
              enddo
           enddo
           elem(ie)%state%ps_v(:,:,nt)=elem(ie)%state%ps_v(:,:,nt) + dt*elem(ie)%rspheremp(:,:)*pstens(:,:,ie)
        enddo
     enddo
  endif

  call t_stopf('advance_hypervis_lf')
  call t_adj_detailf(-1)

  end subroutine advance_hypervis_lf


  !
  ! phl notes: output is stored in first argument. Advances from 2nd argument using tendencies evaluated at 3rd rgument: 
  ! phl: for offline winds use time at 3rd argument (same as rhs currently)
  !
  subroutine compute_and_apply_rhs(np1,nm1,n0,qn0,dt2,elem,hvcoord,hybrid,&
       deriv,nets,nete,compute_diagnostics,eta_ave_w)
  ! ===================================
  ! compute the RHS, accumulate into u(np1) and apply DSS
  !
  !           u(np1) = u(nm1) + dt2*DSS[ RHS(u(n0)) ]
  !
  ! This subroutine is normally called to compute a leapfrog timestep
  ! but by adjusting np1,nm1,n0 and dt2, many other timesteps can be
  ! accomodated.  For example, setting nm1=np1=n0 this routine will
  ! take a forward euler step, overwriting the input with the output.
  !
  !    qn0 = timelevel used to access Qdp() in order to compute virtual Temperature
  !          qn0=-1 for the dry case
  !
  ! if  dt2<0, then the DSS'd RHS is returned in timelevel np1
  !
  ! Combining the RHS and DSS pack operation in one routine
  ! allows us to fuse these two loops for more cache reuse
  !
  ! Combining the dt advance and DSS unpack operation in one routine
  ! allows us to fuse these two loops for more cache reuse
  !
  ! note: for prescribed velocity case, velocity will be computed at
  ! "real_time", which should be the time of timelevel n0.
  !
  !
  ! ===================================
  use kinds, only : real_kind
  use dimensions_mod, only : np, np, nlev
  use hybrid_mod, only : hybrid_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t, divergence_sphere, gradient_sphere, vorticity_sphere
  use edge_mod, only : edgevpack, edgevunpack
  use bndry_mod, only : bndry_exchangev
  use control_mod, only : moisture, qsplit, use_cpstar, rsplit
  use hybvcoord_mod, only : hvcoord_t

  use physical_constants, only : cp, cpwater_vapor, Rgas, kappa
  use physics_mod, only : virtual_specific_heat, virtual_temperature
  use prim_si_mod, only : preq_vertadv, preq_omega_ps, preq_hydrostatic


  implicit none
  integer, intent(in) :: np1,nm1,n0,qn0,nets,nete
  real*8, intent(in) :: dt2
  logical, intent(in)  :: compute_diagnostics

  type (hvcoord_t)     , intent(in) :: hvcoord
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv
  real (kind=real_kind) :: eta_ave_w  ! weighting for eta_dot_dpdn mean flux

  ! local
  real (kind=real_kind), pointer, dimension(:,:)      :: ps         ! surface pressure for current tiime level
  real (kind=real_kind), pointer, dimension(:,:,:)   :: phi

  real (kind=real_kind), dimension(np,np,nlev)   :: omega_p
  real (kind=real_kind), dimension(np,np,nlev)   :: T_v
  real (kind=real_kind), dimension(np,np,nlev)   :: divdp
  real (kind=real_kind), dimension(np,np,nlev+1)   :: eta_dot_dpdn  ! half level vertical velocity on p-grid
  real (kind=real_kind), dimension(np,np)      :: sdot_sum   ! temporary field
  real (kind=real_kind), dimension(np,np,2)    :: vtemp     ! generic gradient storage
  real (kind=real_kind), dimension(np,np)      :: vgrad_T    ! v.grad(T)
  real (kind=real_kind), dimension(np,np)      :: Ephi       ! kinetic energy + PHI term
  real (kind=real_kind), dimension(np,np,2)      :: grad_ps    ! lat-lon coord version
  real (kind=real_kind), dimension(np,np,2,nlev) :: grad_p
  real (kind=real_kind), dimension(np,np,nlev)   :: vort       ! vorticity
  real (kind=real_kind), dimension(np,np,nlev)   :: p          ! pressure
  real (kind=real_kind), dimension(np,np,nlev)   :: dp         ! delta pressure
  real (kind=real_kind), dimension(np,np,nlev)   :: rdp        ! inverse of delta pressure
  real (kind=real_kind), dimension(np,np,nlev)   :: T_vadv     ! temperature vertical advection
  real (kind=real_kind), dimension(np,np,nlev)   :: vgrad_p    ! v.grad(p)
  real (kind=real_kind), dimension(np,np,nlev+1) :: ph               ! half level pressures on p-grid
  real (kind=real_kind), dimension(np,np,2,nlev) :: v_vadv   ! velocity vertical advection
  real (kind=real_kind) ::  kappa_star(np,np,nlev)
  real (kind=real_kind) ::  vtens1(np,np,nlev)
  real (kind=real_kind) ::  vtens2(np,np,nlev)
  real (kind=real_kind) ::  ttens(np,np,nlev)

  real (kind=real_kind) ::  cp2,cp_ratio,E,de,Qt,v1,v2
  real (kind=real_kind) ::  glnps1,glnps2,gpterm
  integer :: i,j,k,kptr,ie

!JMD  call t_barrierf('sync_compute_and_apply_rhs', hybrid%par%comm)
  call t_adj_detailf(+1)
  call t_startf('compute_and_apply_rhs')
  do ie=nets,nete
     !ps => elem(ie)%state%ps_v(:,:,n0)
     phi => elem(ie)%derived%phi(:,:,:)

     ! ==================================================
     ! compute pressure (p) on half levels from ps
     ! using the hybrid coordinates relationship, i.e.
     ! e.g. equation (3.a.92) of the CCM-2 description,
     ! (NCAR/TN-382+STR), June 1993, p. 24.
     ! ==================================================
     ! vertically eulerian only needs grad(ps)
     if (rsplit==0) &
          call gradient_sphere(elem(ie)%state%ps_v(:,:,n0),deriv,elem(ie)%Dinv,grad_ps)


     ! ============================
     ! compute p and delta p
     ! ============================
!#if (defined COLUMN_OPENMP)
!!$omp parallel do private(k,i,j)
!#endif
!     do k=1,nlev+1
!       ph(:,:,k)   = hvcoord%hyai(k)*hvcoord%ps0 + hvcoord%hybi(k)*elem(ie)%state%ps_v(:,:,n0)
!     end do

#if (defined COLUMN_OPENMP_BROKEN)
!  Note that the following does not work with OpenMP threading turned on.
!  The line 
!              p(:,:,k)=p(:,:,k-1) + dp(:,:,k-1)/2 + dp(:,:,k)/2
!  is sequential in that it depends upon the previous p(:,:,k-1) and therefore
!  gives a race condition.
!  !$omp parallel do private(k,i,j,v1,v2,vtemp)
#endif
     do k=1,nlev
        if (rsplit==0) then
           dp(:,:,k) = (hvcoord%hyai(k+1)*hvcoord%ps0 + hvcoord%hybi(k+1)*elem(ie)%state%ps_v(:,:,n0)) &
                - (hvcoord%hyai(k)*hvcoord%ps0 + hvcoord%hybi(k)*elem(ie)%state%ps_v(:,:,n0))
           p(:,:,k)   = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%ps_v(:,:,n0)
           grad_p(:,:,:,k) = hvcoord%hybm(k)*grad_ps(:,:,:)
        else
           ! vertically lagrangian code: we advect dp3d instead of ps_v
           ! we also need grad(p) at all levels (not just grad(ps))
           !p(k)= hyam(k)*ps0 + hybm(k)*ps
           !    = .5*(hyai(k+1)+hyai(k))*ps0 + .5*(hybi(k+1)+hybi(k))*ps
           !    = .5*(ph(k+1) + ph(k) )  = ph(k) + dp(k)/2
           !
           ! p(k+1)-p(k) = ph(k+1)-ph(k) + (dp(k+1)-dp(k))/2
           !             = dp(k) + (dp(k+1)-dp(k))/2 = (dp(k+1)+dp(k))/2
           dp(:,:,k) = elem(ie)%state%dp3d(:,:,k,n0)
           if (k==1) then
              p(:,:,k)=hvcoord%hyai(k)*hvcoord%ps0 + dp(:,:,k)/2
           else
              p(:,:,k)=p(:,:,k-1) + dp(:,:,k-1)/2 + dp(:,:,k)/2
           endif
           ! grad_p(:,:,:,k) = gradient_sphere(p(:,:,k),deriv,elem(ie)%Dinv)
           call gradient_sphere(p(:,:,k),deriv,elem(ie)%Dinv,grad_p(:,:,:,k))
        endif

        rdp(:,:,k) = 1.0D0/dp(:,:,k)

        ! ============================
        ! compute vgrad_lnps
        ! ============================
        do j=1,np
           do i=1,np
              v1 = elem(ie)%state%v(i,j,1,k,n0)
              v2 = elem(ie)%state%v(i,j,2,k,n0)
!              vgrad_p(i,j,k) = &
!                   hvcoord%hybm(k)*(v1*grad_ps(i,j,1) + v2*grad_ps(i,j,2))
              vgrad_p(i,j,k) = (v1*grad_p(i,j,1,k) + v2*grad_p(i,j,2,k))
              vtemp(i,j,1) = v1*dp(i,j,k)
              vtemp(i,j,2) = v2*dp(i,j,k)
           end do
        end do



        ! ================================
        ! Accumulate mean Vel_rho flux in vn0
        ! ================================
        elem(ie)%derived%vn0(:,:,:,k)=elem(ie)%derived%vn0(:,:,:,k)+eta_ave_w*vtemp(:,:,:)


        ! =========================================
        !
        ! Compute relative vorticity and divergence
        !
        ! =========================================
        ! divdp(:,:,k)=divergence_sphere(vtemp,deriv,elem(ie))
        call divergence_sphere(vtemp,deriv,elem(ie),divdp(:,:,k))
        ! vort(:,:,k)=vorticity_sphere(elem(ie)%state%v(:,:,:,k,n0),deriv,elem(ie))
        call vorticity_sphere(elem(ie)%state%v(:,:,:,k,n0),deriv,elem(ie),vort(:,:,k))

     enddo

     ! compute T_v for timelevel n0
     !if ( moisture /= "dry") then
     if (qn0 == -1 ) then
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j)
#endif
        do k=1,nlev
           do j=1,np
              do i=1,np
                 T_v(i,j,k) = elem(ie)%state%T(i,j,k,n0)
                 kappa_star(i,j,k) = kappa
              end do
           end do
        end do
     else
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,Qt)
#endif

        do k=1,nlev
           do j=1,np
              do i=1,np
                 ! Qt = elem(ie)%state%Q(i,j,k,1)
                 Qt = elem(ie)%state%Qdp(i,j,k,1,qn0)/dp(i,j,k)
!!XXgoldyXX
!Qt=0._real_kind
!!XXgoldyXX
                 T_v(i,j,k) = Virtual_Temperature(elem(ie)%state%T(i,j,k,n0),Qt)
                 if (use_cpstar==1) then
                    kappa_star(i,j,k) =  Rgas/Virtual_Specific_Heat(Qt)
                 else
                    kappa_star(i,j,k) = kappa
                 endif
              end do
           end do
        end do
     end if


     ! ====================================================
     ! Compute Hydrostatic equation, modeld after CCM-3
     ! ====================================================
     !call geopotential_t(p,dp,T_v,Rgas,phi)
     call preq_hydrostatic(phi,elem(ie)%state%phis,T_v,p,dp)

     ! ====================================================
     ! Compute omega_p according to CCM-3
     ! ====================================================
     call preq_omega_ps(omega_p,hvcoord,p,vgrad_p,divdp)


     ! ==================================================
     ! zero partial sum for accumulating sum
     !    (div(v_k) + v_k.grad(lnps))*dsigma_k = div( v dp )
     ! used by eta_dot_dpdn and lnps tendency
     ! ==================================================
     sdot_sum=0


     ! ==================================================
     ! Compute eta_dot_dpdn
     ! save sdot_sum as this is the -RHS of ps_v equation
     ! ==================================================
     if (rsplit>0) then
        ! VERTICALLY LAGRANGIAN:   no vertical motion
        eta_dot_dpdn=0
        T_vadv=0
        v_vadv=0
     else
        do k=1,nlev
           ! ==================================================
           ! add this term to PS equation so we exactly conserve dry mass
           ! ==================================================
           sdot_sum(:,:) = sdot_sum(:,:) + divdp(:,:,k)
           eta_dot_dpdn(:,:,k+1) = sdot_sum(:,:)
        end do


        ! ===========================================================
        ! at this point, eta_dot_dpdn contains integral_etatop^eta[ divdp ]
        ! compute at interfaces:
        !    eta_dot_dpdn = -dp/dt - integral_etatop^eta[ divdp ]
        ! for reference: at mid layers we have:
        !    omega = v grad p  - integral_etatop^eta[ divdp ]
        ! ===========================================================
#if (defined COLUMN_OPENMP)
        !$omp parallel do private(k)
#endif
        do k=1,nlev-1
           eta_dot_dpdn(:,:,k+1) = hvcoord%hybi(k+1)*sdot_sum(:,:) - eta_dot_dpdn(:,:,k+1)
        end do

        eta_dot_dpdn(:,:,1     ) = 0.0D0
        eta_dot_dpdn(:,:,nlev+1) = 0.0D0

        ! ===========================================================
        ! Compute vertical advection of T and v from eq. CCM2 (3.b.1)
        ! ==============================================
        call preq_vertadv(elem(ie)%state%T(:,:,:,n0),elem(ie)%state%v(:,:,:,:,n0), &
             eta_dot_dpdn,rdp,T_vadv,v_vadv)
     endif


     ! ================================
     ! accumulate mean vertical flux:
     ! ================================
#if (defined COLUMN_OPENMP)
     !$omp parallel do private(k)
#endif
     do k=1,nlev  !  Loop index added (AAM)
        elem(ie)%derived%eta_dot_dpdn(:,:,k) = &
             elem(ie)%derived%eta_dot_dpdn(:,:,k) + eta_ave_w*eta_dot_dpdn(:,:,k)
        elem(ie)%derived%omega_p(:,:,k) = &
             elem(ie)%derived%omega_p(:,:,k) + eta_ave_w*omega_p(:,:,k)
     enddo
     elem(ie)%derived%eta_dot_dpdn(:,:,nlev+1) = &
          elem(ie)%derived%eta_dot_dpdn(:,:,nlev+1) + eta_ave_w*eta_dot_dpdn(:,:,nlev+1)





     ! ==============================================
     ! Compute phi + kinetic energy term: 10*nv*nv Flops
     ! ==============================================
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,v1,v2,E,Ephi,vtemp,vgrad_T,gpterm,glnps1,glnps2)
#endif
     do k=1,nlev
        do j=1,np
           do i=1,np
              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)
              E = 0.5D0*( v1*v1 + v2*v2 )
              Ephi(i,j)=E+phi(i,j,k)+elem(ie)%derived%pecnd(i,j,k)
           end do
        end do
        ! ================================================
        ! compute gradp term (ps/p)*(dp/dps)*T
        ! ================================================
        ! vtemp(:,:,:)   = gradient_sphere(elem(ie)%state%T(:,:,k,n0),deriv,elem(ie)%Dinv)
        call gradient_sphere(elem(ie)%state%T(:,:,k,n0),deriv,elem(ie)%Dinv,vtemp(:,:,:))
        do j=1,np
           do i=1,np
              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)
              vgrad_T(i,j) =  v1*vtemp(i,j,1) + v2*vtemp(i,j,2)
           end do
        end do


        ! vtemp = grad ( E + PHI )
        ! vtemp = gradient_sphere(Ephi(:,:),deriv,elem(ie)%Dinv)
        call gradient_sphere(Ephi(:,:),deriv,elem(ie)%Dinv,vtemp)

        do j=1,np
           do i=1,np
!              gpterm = hvcoord%hybm(k)*T_v(i,j,k)/p(i,j,k)
!              glnps1 = Rgas*gpterm*grad_ps(i,j,1)
!              glnps2 = Rgas*gpterm*grad_ps(i,j,2)
              gpterm = T_v(i,j,k)/p(i,j,k)
              glnps1 = Rgas*gpterm*grad_p(i,j,1,k)
              glnps2 = Rgas*gpterm*grad_p(i,j,2,k)

              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)

              vtens1(i,j,k) =   - v_vadv(i,j,1,k)                           &
                   + v2*(elem(ie)%fcor(i,j) + vort(i,j,k))        &
                   - vtemp(i,j,1) - glnps1
              !
              ! phl: add forcing term to zonal wind u
              !
              vtens2(i,j,k) =   - v_vadv(i,j,2,k)                            &
                   - v1*(elem(ie)%fcor(i,j) + vort(i,j,k))        &
                   - vtemp(i,j,2) - glnps2
              !
              ! phl: add forcing term to meridional wind v
              !
              ttens(i,j,k)  = - T_vadv(i,j,k) - vgrad_T(i,j) + kappa_star(i,j,k)*T_v(i,j,k)*omega_p(i,j,k)
              !
              ! phl: add forcing term to T
              !
           end do
        end do

     end do

#ifdef ENERGY_DIAGNOSTICS
     ! =========================================================
     !
     ! diagnostics
     ! recomputes some gradients that were not saved above
     ! uses:  sdot_sum(), eta_dot_dpdn(), grad_ps()
     ! grad_phi(), dp(), p(), T_vadv(), v_vadv(), divdp()
     ! =========================================================

     ! =========================================================
     ! (AAM) - This section has accumulations over vertical levels.
     !   Be careful if implementing OpenMP
     ! =========================================================

     if (compute_diagnostics) then
        elem(ie)%accum%KEhorz1=0
        elem(ie)%accum%KEhorz2=0
        elem(ie)%accum%IEhorz1=0
        elem(ie)%accum%IEhorz2=0
        elem(ie)%accum%IEhorz1_wet=0
        elem(ie)%accum%IEhorz2_wet=0
        elem(ie)%accum%KEvert1=0
        elem(ie)%accum%KEvert2=0
        elem(ie)%accum%IEvert1=0
        elem(ie)%accum%IEvert2=0
        elem(ie)%accum%IEvert1_wet=0
        elem(ie)%accum%IEvert2_wet=0
        elem(ie)%accum%T1=0
        elem(ie)%accum%T2=0
        elem(ie)%accum%T2_s=0
        elem(ie)%accum%S1=0
        elem(ie)%accum%S1_wet=0
        elem(ie)%accum%S2=0

        do j=1,np
           do i=1,np
              elem(ie)%accum%S2(i,j) = elem(ie)%accum%S2(i,j) - &
                   sdot_sum(i,j)*elem(ie)%state%phis(i,j)
           enddo
        enddo

        do k=1,nlev
           ! vtemp = grad_E(:,:,k)
           do j=1,np
              do i=1,np
                 v1     = elem(ie)%state%v(i,j,1,k,n0)
                 v2     = elem(ie)%state%v(i,j,2,k,n0)
                 Ephi(i,j)=0.5D0*( v1*v1 + v2*v2 )
              enddo
           enddo
           ! vtemp = gradient_sphere(Ephi,deriv,elem(ie)%Dinv)
           call gradient_sphere(Ephi,deriv,elem(ie)%Dinv,vtemp)
           do j=1,np
              do i=1,np
                 ! dp/dn u dot grad(E)
                 v1     = elem(ie)%state%v(i,j,1,k,n0)
                 v2     = elem(ie)%state%v(i,j,2,k,n0)
                 elem(ie)%accum%KEhorz2(i,j) = elem(ie)%accum%KEhorz2(i,j) + &
                      (v1*vtemp(i,j,1)  + v2*vtemp(i,j,2))*dp(i,j,k)
                 ! E div( u dp/dn )
                 elem(ie)%accum%KEhorz1(i,j) = elem(ie)%accum%KEhorz1(i,j) + Ephi(i,j)*divdp(i,j,k)

                 ! Cp T div( u dp/dn)   ! dry horizontal advection component
                 elem(ie)%accum%IEhorz1(i,j) = elem(ie)%accum%IEhorz1(i,j) + Cp*elem(ie)%state%T(i,j,k,n0)*divdp(i,j,k)


              enddo
           enddo


           ! vtemp = grad_phi(:,:,k)
           ! vtemp = gradient_sphere(phi(:,:,k),deriv,elem(ie)%Dinv)
           call gradient_sphere(phi(:,:,k),deriv,elem(ie)%Dinv,vtemp)
           do j=1,np
              do i=1,np
                 v1     = elem(ie)%state%v(i,j,1,k,n0)
                 v2     = elem(ie)%state%v(i,j,2,k,n0)
                 E = 0.5D0*( v1*v1 + v2*v2 )
                 ! NOTE:  Cp_star = Cp + (Cpv-Cp)*q
                 ! advection terms can thus be broken into two components: dry and wet
                 ! dry components cancel exactly
                 ! wet components should cancel exactly
                 !
                 ! some diagnostics
                 ! e = eta_dot_dpdn()
                 de =  eta_dot_dpdn(i,j,k+1)-eta_dot_dpdn(i,j,k)
                 ! Cp T de/dn, integral dn:
                 elem(ie)%accum%IEvert1(i,j)=elem(ie)%accum%IEvert1(i,j) + Cp*elem(ie)%state%T(i,j,k,n0)*de
                 ! E de/dn
                 elem(ie)%accum%KEvert1(i,j)=elem(ie)%accum%KEvert1(i,j) + E*de
                 ! Cp T_vadv dp/dn
                 elem(ie)%accum%IEvert2(i,j)=elem(ie)%accum%IEvert2(i,j) + Cp*T_vadv(i,j,k)*dp(i,j,k)
                 ! dp/dn V dot V_vadv
                 elem(ie)%accum%KEvert2(i,j)=elem(ie)%accum%KEvert2(i,j) + (v1*v_vadv(i,j,1,k) + v2*v_vadv(i,j,2,k)) *dp(i,j,k)

                 ! IEvert1_wet():  (Cpv-Cp) T Qdp_vadv  (Q equation)
                 ! IEvert2_wet():  (Cpv-Cp) Qdp T_vadv   T equation
                 if (use_cpstar==1) then
                 elem(ie)%accum%IEvert2_wet(i,j)=elem(ie)%accum%IEvert2_wet(i,j) +&
                      (Cpwater_vapor-Cp)*elem(ie)%state%Q(i,j,k,1)*T_vadv(i,j,k)*dp(i,j,k)
                 endif

                 gpterm = T_v(i,j,k)/p(i,j,k)
                 elem(ie)%accum%T1(i,j) = elem(ie)%accum%T1(i,j) - &
                      Rgas*gpterm*(grad_p(i,j,1,k)*v1 + grad_p(i,j,2,k)*v2)*dp(i,j,k)

                 elem(ie)%accum%T2(i,j) = elem(ie)%accum%T2(i,j) - &
                      (vtemp(i,j,1)*v1 + vtemp(i,j,2)*v2)*dp(i,j,k)

                 ! S1 = < Cp_star dp/dn , RT omega_p/cp_star >
                 elem(ie)%accum%S1(i,j) = elem(ie)%accum%S1(i,j) + &
                      Rgas*T_v(i,j,k)*omega_p(i,j,k)*dp(i,j,k)

                 ! cp_star = cp + cp2
                 if (use_cpstar==1) then
                 cp2 = (Cpwater_vapor-Cp)*elem(ie)%state%Q(i,j,k,1)
                 cp_ratio = cp2/(cp+cp2)
                 elem(ie)%accum%S1_wet(i,j) = elem(ie)%accum%S1_wet(i,j) + &
                      cp_ratio*(Rgas*T_v(i,j,k)*omega_p(i,j,k)*dp(i,j,k))
                 endif

                 elem(ie)%accum%CONV(i,j,:,k)=-Rgas*gpterm*grad_p(i,j,:,k)-vtemp(i,j,:)
              enddo
           enddo

           ! vtemp(:,:,:) = gradient_sphere(elem(ie)%state%phis(:,:),deriv,elem(ie)%Dinv)
           call gradient_sphere(elem(ie)%state%phis(:,:),deriv,elem(ie)%Dinv,vtemp)
           do j=1,np
              do i=1,np
                 v1     = elem(ie)%state%v(i,j,1,k,n0)
                 v2     = elem(ie)%state%v(i,j,2,k,n0)
                 elem(ie)%accum%T2_s(i,j) = elem(ie)%accum%T2_s(i,j) - &
                      (vtemp(i,j,1)*v1 + vtemp(i,j,2)*v2)*dp(i,j,k)
              enddo
           enddo

           !vtemp(:,:,:)   = gradient_sphere(elem(ie)%state%T(:,:,k,n0),deriv,elem(ie)%Dinv)
           call gradient_sphere(elem(ie)%state%T(:,:,k,n0),deriv,elem(ie)%Dinv,vtemp)
           do j=1,np
              do i=1,np
                 v1     = elem(ie)%state%v(i,j,1,k,n0)
                 v2     = elem(ie)%state%v(i,j,2,k,n0)

                 ! Cp dp/dn u dot gradT
                 elem(ie)%accum%IEhorz2(i,j) = elem(ie)%accum%IEhorz2(i,j) + &
                      Cp*(v1*vtemp(i,j,1) + v2*vtemp(i,j,2))*dp(i,j,k)

                 if (use_cpstar==1) then
                 elem(ie)%accum%IEhorz2_wet(i,j) = elem(ie)%accum%IEhorz2_wet(i,j) + &
                      (Cpwater_vapor-Cp)*elem(ie)%state%Q(i,j,k,1)*&
                      (v1*vtemp(i,j,1) + v2*vtemp(i,j,2))*dp(i,j,k)
                 endif

              enddo
           enddo

        enddo
     endif
#endif
     ! =========================================================
     ! local element timestep, store in np1.
     ! note that we allow np1=n0 or nm1
     ! apply mass matrix
     ! =========================================================
     if (dt2<0) then
        ! calling program just wanted DSS'd RHS, skip time advance
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
        do k=1,nlev
           elem(ie)%state%v(:,:,1,k,np1) = elem(ie)%spheremp(:,:)*vtens1(:,:,k)
           elem(ie)%state%v(:,:,2,k,np1) = elem(ie)%spheremp(:,:)*vtens2(:,:,k)
           elem(ie)%state%T(:,:,k,np1) = elem(ie)%spheremp(:,:)*ttens(:,:,k)
           if (rsplit>0) &
                elem(ie)%state%dp3d(:,:,k,np1) = -elem(ie)%spheremp(:,:)*&
                (divdp(:,:,k) + eta_dot_dpdn(:,:,k+1)-eta_dot_dpdn(:,:,k))
        enddo
        elem(ie)%state%ps_v(:,:,np1) = -elem(ie)%spheremp(:,:)*sdot_sum
     else
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
        do k=1,nlev
           elem(ie)%state%v(:,:,1,k,np1) = elem(ie)%spheremp(:,:)*( elem(ie)%state%v(:,:,1,k,nm1) + dt2*vtens1(:,:,k) )
           elem(ie)%state%v(:,:,2,k,np1) = elem(ie)%spheremp(:,:)*( elem(ie)%state%v(:,:,2,k,nm1) + dt2*vtens2(:,:,k) )
           elem(ie)%state%T(:,:,k,np1) = elem(ie)%spheremp(:,:)*(elem(ie)%state%T(:,:,k,nm1) + dt2*ttens(:,:,k))
           if (rsplit>0) &
                elem(ie)%state%dp3d(:,:,k,np1) = elem(ie)%spheremp(:,:)*&
                (elem(ie)%state%dp3d(:,:,k,nm1)-dt2*&
                (divdp(:,:,k) + eta_dot_dpdn(:,:,k+1)-eta_dot_dpdn(:,:,k)))
        enddo
        elem(ie)%state%ps_v(:,:,np1) = elem(ie)%spheremp(:,:)*( elem(ie)%state%ps_v(:,:,nm1) - dt2*sdot_sum )

     endif


     ! =========================================================
     !
     ! Pack ps(np1), T, and v tendencies into comm buffer
     !
     ! =========================================================
     kptr=0
     call edgeVpack(edge3p1, elem(ie)%state%ps_v(:,:,np1),1,kptr,ie)

     kptr=1
     call edgeVpack(edge3p1, elem(ie)%state%T(:,:,:,np1),nlev,kptr,ie)

     kptr=nlev+1
     call edgeVpack(edge3p1, elem(ie)%state%v(:,:,:,:,np1),2*nlev,kptr,ie)

     if (rsplit>0) then
        kptr=kptr+2*nlev
        call edgeVpack(edge3p1, elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr,ie)
     endif
  end do

  ! =============================================================
    ! Insert communications here: for shared memory, just a single
  ! sync is required
  ! =============================================================
  call bndry_exchangeV(hybrid,edge3p1)
  do ie=nets,nete
     ! ===========================================================
     ! Unpack the edges for vgrad_T and v tendencies...
     ! ===========================================================
     kptr=0
     call edgeVunpack(edge3p1, elem(ie)%state%ps_v(:,:,np1), 1, kptr, ie)

     kptr=1
     call edgeVunpack(edge3p1, elem(ie)%state%T(:,:,:,np1), nlev, kptr, ie) 

     kptr=nlev+1
     call edgeVunpack(edge3p1, elem(ie)%state%v(:,:,:,:,np1), 2*nlev, kptr, ie)

     if (rsplit>0) then
        kptr=kptr+2*nlev
        call edgeVunpack(edge3p1, elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr,ie)
     endif

     ! ====================================================
     ! Scale tendencies by inverse mass matrix
     ! ====================================================

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
     do k=1,nlev
        elem(ie)%state%T(:,:,k,np1)   = elem(ie)%rspheremp(:,:)*elem(ie)%state%T(:,:,k,np1)
        elem(ie)%state%v(:,:,1,k,np1) = elem(ie)%rspheremp(:,:)*elem(ie)%state%v(:,:,1,k,np1)
        elem(ie)%state%v(:,:,2,k,np1) = elem(ie)%rspheremp(:,:)*elem(ie)%state%v(:,:,2,k,np1)
     end do

     if (rsplit>0) then
        ! vertically lagrangian: complete dp3d timestep:
        do k=1,nlev
           elem(ie)%state%dp3d(:,:,k,np1)= elem(ie)%rspheremp(:,:)*elem(ie)%state%dp3d(:,:,k,np1)
        enddo
        ! when debugging: also update ps_v
        !elem(ie)%state%ps_v(:,:,np1) = elem(ie)%rspheremp(:,:)*elem(ie)%state%ps_v(:,:,np1)
     else
        ! vertically eulerian: complete ps_v timestep:
        elem(ie)%state%ps_v(:,:,np1) = elem(ie)%rspheremp(:,:)*elem(ie)%state%ps_v(:,:,np1)
     endif

  end do


  call t_stopf('compute_and_apply_rhs')
  call t_adj_detailf(-1)

  end subroutine compute_and_apply_rhs

  subroutine smooth_phis(phis,elem,hybrid,deriv,nets,nete,minf,numcycle)
  use dimensions_mod, only : np, np, nlev
  use control_mod, only : smooth_phis_nudt, hypervis_scaling
  use hybrid_mod, only : hybrid_t
  use edgetype_mod, only : EdgeBuffer_t
  use edge_mod, only : edgevpack, edgevunpack, edgevunpackmax, edgevunpackmin
  use bndry_mod, only : bndry_exchangev
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t , laplace_sphere_wk
  use time_mod, only : TimeLevel_t
  implicit none

  integer :: nets,nete
  real (kind=real_kind), dimension(np,np,nets:nete), intent(inout)   :: phis
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv
  real (kind=real_kind), intent(in)   :: minf
  integer,               intent(in) :: numcycle

  ! local
  real (kind=real_kind), dimension(np,np,nets:nete) :: pstens
  real (kind=real_kind), dimension(nets:nete) :: pmin,pmax
  real (kind=real_kind) :: mx,mn
  integer :: nt,ie,ic,i,j,order,order_max, iuse
  logical :: use_var_coef


  ! compute local element neighbor min/max
  do ie=nets,nete
     pstens(:,:,ie)=minval(phis(:,:,ie))
     call edgeVpack(edge3p1,pstens(:,:,ie),1,0,ie)
  enddo
  call bndry_exchangeV(hybrid,edge3p1)
  do ie=nets,nete
     call edgeVunpackMin(edge3p1, pstens(:,:,ie), 1, 0, ie)
     pmin(ie)=minval(pstens(:,:,ie))
  enddo
  do ie=nets,nete
     pstens(:,:,ie)=maxval(phis(:,:,ie))
     call edgeVpack(edge3p1,pstens(:,:,ie),1,0,ie)
  enddo
  call bndry_exchangeV(hybrid,edge3p1)
  do ie=nets,nete
     call edgeVunpackMax(edge3p1, pstens(:,:,ie), 1, 0, ie)
     pmax(ie)=maxval(pstens(:,:,ie))
  enddo

  ! order = 1   grad^2 laplacian
  ! order = 2   grad^4 (need to add a negative sign)
  ! order = 3   grad^6
  ! order = 4   grad^8 (need to add a negative sign)
  order_max = 1


  use_var_coef=.true.
  if (hypervis_scaling/=0) then
     ! for tensorHV option, we turn off the tensor except for *last* laplace operator
     use_var_coef=.false.
     if (hypervis_scaling>=3) then
        ! with a 3.2 or 4 scaling, assume hyperviscosity
        order_max = 2
     endif
  endif


  do ic=1,numcycle
     pstens=phis

     do order=1,order_max-1

        do ie=nets,nete
           ! pstens(:,:,ie)=laplace_sphere_wk(pstens(:,:,ie),deriv,elem(ie),var_coef=use_var_coef)
           call laplace_sphere_wk(pstens(:,:,ie),deriv,elem(ie),pstens(:,:,ie), var_coef=use_var_coef)
           call edgeVpack(edge3p1,pstens(:,:,ie),1,0,ie)
        enddo
        call bndry_exchangeV(hybrid,edge3p1)
        do ie=nets,nete
           call edgeVunpack(edge3p1, pstens(:,:,ie), 1, 0, ie)
           pstens(:,:,ie)=pstens(:,:,ie)*elem(ie)%rspheremp(:,:)
        enddo
     enddo
     do ie=nets,nete
        ! pstens(:,:,ie)=laplace_sphere_wk(pstens(:,:,ie),deriv,elem(ie),var_coef=.true.)
        call laplace_sphere_wk(pstens(:,:,ie),deriv,elem(ie),pstens(:,:,ie), var_coef=.true.)
     enddo
     if (mod(order_max,2)==0) pstens=-pstens

     do ie=nets,nete
        !  ps(t+1) = ps(t) + Minv * DSS * M * RHS
        !  ps(t+1) = Minv * DSS * M [ ps(t) +  RHS ]
        ! but output of biharminc_wk is of the form M*RHS.  rewrite as:
        !  ps(t+1) = Minv * DSS * M [ ps(t) +  M*RHS/M ]
        ! so we can apply limiter to ps(t) +  (M*RHS)/M
#if 1
        mn=pmin(ie)
        mx=pmax(ie)
        iuse = numcycle+1  ! always apply min/max limiter
#endif
        phis(:,:,ie)=phis(:,:,ie) + &
           smooth_phis_nudt*pstens(:,:,ie)/elem(ie)%spheremp(:,:)


        ! remove new extrema.  could use conservative reconstruction from advection
        ! but no reason to conserve mean PHI.
        if (ic < iuse) then
        do i=1,np
        do j=1,np
           if (phis(i,j,ie)>mx) phis(i,j,ie)=mx
           if (phis(i,j,ie)<mn) phis(i,j,ie)=mn
        enddo
        enddo
        endif


        ! user specified minimum
        do i=1,np
        do j=1,np
           if (phis(i,j,ie)<minf) phis(i,j,ie)=minf
        enddo
        enddo

        phis(:,:,ie)=phis(:,:,ie)*elem(ie)%spheremp(:,:)
        call edgeVpack(edge3p1,phis(:,:,ie),1,0,ie)
     enddo
     call bndry_exchangeV(hybrid,edge3p1)
     do ie=nets,nete
        call edgeVunpack(edge3p1, phis(:,:,ie), 1, 0, ie)
        phis(:,:,ie)=phis(:,:,ie)*elem(ie)%rspheremp(:,:)
     enddo
  enddo
  end subroutine smooth_phis

  subroutine overwrite_SEdensity(elem, fvm, dt_q, hybrid,nets,nete, np1)
    use fvm_reconstruction_mod, only: reconstruction
    use fvm_filter_mod, only: monotonic_gradient_cart, recons_val_cart
    use dimensions_mod, only : np, nlev, nc,nhe
    use hybrid_mod, only : hybrid_t
    use edgetype_mod, only : EdgeBuffer_t
    use edge_mod, only : edgevpack, edgevunpack, edgevunpackmax, edgevunpackmin
    use bndry_mod, only : bndry_exchangev
    use element_mod, only : element_t
    use derivative_mod, only : derivative_t , laplace_sphere_wk
    use time_mod, only : TimeLevel_t
    use fvm_control_volume_mod, only : fvm_struct


    type (element_t) , intent(inout)        :: elem(:)

    type(fvm_struct), intent(inout)         :: fvm(:)
    type (hybrid_t), intent(in)             :: hybrid  ! distributed parallel structure (shared)

    integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete  ! ending thread element number   (private)
    integer, intent(in)                     :: np1
    integer :: ie, k

    real (kind=real_kind)             :: xp,yp, tmpval, dt_q
    integer                           :: i, j,ix, jy, starti,endi,tmpi

    real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe)      :: recons

    if ((nc .ne. 4) .or. (np .ne. 4)) then
      if(hybrid%masterthread) then
        print *,"You are in OVERWRITE SE AIR DENSITY MODE"
        print *,"This only works for nc=4 and np=4"
        print *,"Write a new search algorithm or pay $10000!"
      endif
      stop
    endif
#if defined(_FVM)
    do ie=nets,nete
      call reconstruction(fvm(ie)%psc, fvm(ie),recons)
      call monotonic_gradient_cart(fvm(ie)%psc, fvm(ie),recons, elem(ie)%desc)
      do j=1,np
        do i=1,np
          xp=tan(elem(ie)%cartp(i,j)%x)
          yp=tan(elem(ie)%cartp(i,j)%y)
          ix=i
          jy=j
          ! Search index along "x"  (bisection method)
!           starti = 1
!           endi = nc+1
!           do
!              if  ((endi-starti) <=  1)  exit
!              tmpi = (endi + starti)/2
!              if (xp  >  fvm%acartx(tmpi)) then
!                 starti = tmpi
!              else
!                 endi = tmpi
!              endif
!           enddo
!           ix = starti
!
!         ! Search index along "y"
!           starti = 1
!           endi = nc+1
!           do
!              if  ((endi-starti) <=  1)  exit
!              tmpi = (endi + starti)/2
!              if (yp  >  fvm%acarty(tmpi)) then
!                 starti = tmpi
!              else
!                 endi = tmpi
!              endif
!           enddo
!           jy = starti

          call recons_val_cart(fvm(ie)%psc, xp,yp,fvm(ie)%spherecentroid,recons,ix,jy,tmpval)
          elem(ie)%state%ps_v(i,j,np1)= elem(ie)%state%ps_v(i,j,np1) +&
               dt_q*(tmpval - elem(ie)%state%ps_v(i,j,np1) )/(7*24*60*60)
        end do
      end do
      elem(ie)%state%ps_v(:,:,np1)=elem(ie)%state%ps_v(:,:,np1)*elem(ie)%spheremp(:,:)
     call edgeVpack(edge3p1,elem(ie)%state%ps_v(:,:,np1),1,0,ie)
  enddo
  call bndry_exchangeV(hybrid,edge3p1)
  do ie=nets,nete
     call edgeVunpack(edge3p1, elem(ie)%state%ps_v(:,:,np1), 1, 0, ie)
     elem(ie)%state%ps_v(:,:,np1)=elem(ie)%state%ps_v(:,:,np1)*elem(ie)%rspheremp(:,:)
  enddo
#endif
  end subroutine overwrite_SEdensity


end module prim_advance_mod

