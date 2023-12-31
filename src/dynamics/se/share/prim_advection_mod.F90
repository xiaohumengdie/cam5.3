#define OVERLAP 1
module prim_advection_mod
!
! two formulations.  both are conservative
! u grad Q formulation:
!
!    d/dt[ Q] +  U grad Q  +  eta_dot dp/dn dQ/dp  = 0
!                            ( eta_dot dQ/dn )
!
!    d/dt[ dp/dn ] = div( dp/dn U ) + d/dn ( eta_dot dp/dn )
!
! total divergence formulation:
!    d/dt[dp/dn Q] +  div( U dp/dn Q ) + d/dn ( eta_dot dp/dn Q ) = 0
!
! for convience, rewrite this as dp Q:  (since dn does not depend on time or the horizonal):
! equation is now:
!    d/dt[dp Q] +  div( U dp Q ) + d( eta_dot_dpdn Q ) = 0
!
!
  use shr_kind_mod,           only: r8=>shr_kind_r8
  use dimensions_mod,         only: nlev, np, qsize, nc
  use physconst,              only: cpair
  use derivative_mod,         only: derivative_t
  use element_mod,            only: element_t

  use hybvcoord_mod,          only: hvcoord_t
  use time_mod,               only: TimeLevel_t, TimeLevel_Qdp
  use control_mod,            only: hypervis_order
  use control_mod,            only: nu_q, nu_p, limiter_option, hypervis_subcycle_q, rsplit
  use edge_mod,               only: edgevpack, edgevunpack, initedgebuffer, initedgesbuffer

  use edgetype_mod,           only: EdgeBuffer_t
  use hybrid_mod,             only: hybrid_t
  use viscosity_mod,          only: biharmonic_wk_scalar, neighbor_minmax, &
                                 neighbor_minmax_start, neighbor_minmax_finish
  use perf_mod,               only: t_startf, t_stopf, t_barrierf
  use cam_abortutils,         only: endrun

  implicit none

  private
  save

  public :: Prim_Advec_Init1, Prim_Advec_Init2
  public :: Prim_Advec_Tracers_remap

  type (EdgeBuffer_t)      :: edgeAdv, edgeAdvp1, edgeAdvQminmax, edgeAdv1,  edgeveloc

  integer,parameter :: DSSeta = 1
  integer,parameter :: DSSomega = 2
  integer,parameter :: DSSdiv_vdp_ave = 3
  integer,parameter :: DSSno_var = -1

  real(kind=r8), allocatable :: qmin(:,:,:), qmax(:,:,:)

  type (derivative_t), public, allocatable   :: deriv(:) ! derivative struct (nthreads)


contains


  subroutine Prim_Advec_Init1(par, elem, n_domains)
    use dimensions_mod, only : nlev, qsize, nelemd
    use parallel_mod,   only : parallel_t
    type(parallel_t)    :: par
    type (element_t)    :: elem(:)
    integer, intent(in) :: n_domains
    !
    ! Shared buffer pointers.
    ! Using "=> null()" in a subroutine is usually bad, because it makes
    ! the variable have an implicit "save", and therefore shared between
    ! threads. But in this case we want shared pointers.
    real(kind=r8), pointer :: buf_ptr(:) => null()
    real(kind=r8), pointer :: receive_ptr(:) => null()












    ! remaining edge buffers can share %buf and %receive with edgeAdvQ3
    ! (This is done through the optional 1D pointer arguments.)
    call initEdgeBuffer(par,edgeAdvp1,elem,qsize*nlev + nlev)

    call initEdgeBuffer(par,edgeAdv,elem,qsize*nlev)

    ! This is a different type of buffer pointer allocation
    ! used for determine the minimum and maximum value from
    ! neighboring  elements
    call initEdgeSBuffer(par,edgeAdvQminmax,elem,qsize*nlev*2)


    call initEdgeBuffer(par,edgeAdv1,elem,nlev)
    call initEdgeBuffer(par,edgeveloc,elem,2*nlev)


    ! Don't actually want these saved, if this is ever called twice.
    nullify(buf_ptr)
    nullify(receive_ptr)

    allocate(deriv(0:n_domains-1))

    ! this static array is shared by all threads, so dimension for all threads (nelemd), not nets:nete:
    allocate (qmin(nlev,qsize,nelemd))
    allocate (qmax(nlev,qsize,nelemd))

  end subroutine Prim_Advec_Init1

  subroutine Prim_Advec_Init2(hybrid)
    use dimensions_mod, only : nc
    use derivative_mod, only : derivinit

    type (hybrid_t), intent(in) :: hybrid

    ! ==================================
    ! Initialize derivative structure
    ! ==================================
    call derivinit(deriv(hybrid%ithr))
  end subroutine Prim_Advec_Init2

  subroutine Prim_Advec_Tracers_remap( elem , deriv , hvcoord , hybrid , dt , tl , nets , nete )
    implicit none
    type (element_t)     , intent(inout) :: elem(:)
    type (derivative_t)  , intent(in   ) :: deriv
    type (hvcoord_t)     , intent(in   ) :: hvcoord
    type (hybrid_t)      , intent(in   ) :: hybrid
    real(kind=r8) , intent(in   ) :: dt
    type (TimeLevel_t)   , intent(inout) :: tl
    integer              , intent(in   ) :: nets
    integer              , intent(in   ) :: nete


    !print *,'prim_Advec_Tracers_remap: qsize: ',qsize
    call Prim_Advec_Tracers_remap_rk2( elem , deriv , hvcoord , hybrid , dt , tl , nets , nete )
  end subroutine Prim_Advec_Tracers_remap


  subroutine euler_step_driver(np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier )


  integer              , intent(in   )         :: np1_qdp, n0_qdp
  real (kind=r8), intent(in   )         :: dt
  type (element_t)     , intent(inout)         :: elem(:)
  type (hvcoord_t)     , intent(in   )         :: hvcoord
  type (hybrid_t)      , intent(in   )         :: hybrid
  type (derivative_t)  , intent(in   )         :: deriv
  integer              , intent(in   )         :: nets
  integer              , intent(in   )         :: nete
  integer              , intent(in   )         :: DSSopt
  integer              , intent(in   )         :: rhs_multiplier

  call euler_step( np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier)

  end subroutine euler_step_driver

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! forward-in-time 2 level vertically lagrangian step
!  this code takes a lagrangian step in the horizontal
! (complete with DSS), and then applies a vertical remap
!
! This routine may use dynamics fields at timelevel np1
! In addition, other fields are required, which have to be
! explicitly saved by the dynamics:  (in elem(ie)%derived struct)
!
! Fields required from dynamics: (in
!    omega_p   it will be DSS'd here, for later use by CAM physics
!              we DSS omega here because it can be done for "free"
!    Consistent mass/tracer-mass advection (used if subcycling turned on)
!       dp()   dp at timelevel n0
!       vn0()  mean flux  < U dp > going from n0 to np1
!
! 3 stage
!    Euler step from t     -> t+.5
!    Euler step from t+.5  -> t+1.0
!    Euler step from t+1.0 -> t+1.5
!    u(t) = u(t)/3 + u(t+2)*2/3
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
  subroutine Prim_Advec_Tracers_remap_rk2( elem , deriv , hvcoord , hybrid , dt , tl , nets , nete )
    use derivative_mod, only : divergence_sphere
    use control_mod   , only : qsplit
    use hybrid_mod    , only : get_loop_ranges

    type (element_t)     , intent(inout) :: elem(:)
    type (derivative_t)  , intent(in   ) :: deriv
    type (hvcoord_t)     , intent(in   ) :: hvcoord
    type (hybrid_t)      , intent(in   ) :: hybrid
    real(kind=r8) , intent(in   ) :: dt
    type (TimeLevel_t)   , intent(inout) :: tl
    integer              , intent(in   ) :: nets
    integer              , intent(in   ) :: nete

    real (kind=r8), dimension(np,np,2     ) :: gradQ
    integer :: k,ie
    integer :: rkstage,rhs_multiplier
    integer :: n0_qdp, np1_qdp
    integer :: kbeg,kend,qbeg,qend

!    call t_barrierf('sync_prim_advec_tracers_remap_k2', hybrid%par%comm)
!    call t_startf('prim_advec_tracers_remap_rk2')
    call TimeLevel_Qdp( tl, qsplit, n0_qdp, np1_qdp) !time levels for qdp are not the same
    rkstage = 3 !   3 stage RKSSP scheme, with optimal SSP CFL

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! RK2 2D advection step
    ! note: stage 3 we take the oppertunity to DSS omega
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! use these for consistent advection (preserve Q=1)
    ! derived%vdp_ave        =  mean horiz. flux:   U*dp
    ! derived%eta_dot_dpdn    =  mean vertical velocity (used for remap)
    ! derived%omega_p         =  advection code will DSS this for the physics, but otherwise
    !                            it is not needed
    ! Also: save a copy of div(U dp) in derived%div(:,:,:,1), which will be DSS'd
    !       and a DSS'ed version stored in derived%div(:,:,:,2)

    call get_loop_ranges(hybrid,kbeg=kbeg,kend=kend,qbeg=qbeg,qend=qend)

    do ie=nets,nete
      do k=kbeg,kend
        ! div( U dp Q),
        gradQ(:,:,1)=elem(ie)%derived%vn0(:,:,1,k)
        gradQ(:,:,2)=elem(ie)%derived%vn0(:,:,2,k)
        ! elem(ie)%derived%divdp(:,:,k) = divergence_sphere(gradQ,deriv,elem(ie))
        call divergence_sphere(gradQ,deriv,elem(ie),elem(ie)%derived%divdp(:,:,k))
        elem(ie)%derived%divdp_proj(:,:,k) = elem(ie)%derived%divdp(:,:,k)
      enddo
    enddo


    !rhs_multiplier is for obtaining dp_tracers at each stage:
    !dp_tracers(stage) = dp - rhs_multiplier*dt*divdp_proj
!    call t_startf('euler_step')

    rhs_multiplier = 0
    call euler_step_driver( np1_qdp, n0_qdp , dt/2, elem, hvcoord, hybrid, deriv, nets, nete, DSSdiv_vdp_ave, rhs_multiplier )

    rhs_multiplier = 1
    call euler_step_driver( np1_qdp, np1_qdp, dt/2, elem, hvcoord, hybrid, deriv, nets, nete, DSSeta         , rhs_multiplier )

    rhs_multiplier = 2
    call euler_step_driver( np1_qdp, np1_qdp, dt/2, elem, hvcoord, hybrid, deriv, nets, nete, DSSomega      , rhs_multiplier )

!    call t_stopf ('euler_step')

    !to finish the 2D advection step, we need to average the t and t+2 results to get a second order estimate for t+1.
    call qdp_time_avg( elem , rkstage , n0_qdp , np1_qdp , hybrid, nets , nete )

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Dissipation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ( limiter_option == 8  ) then
      ! dissipation was applied in RHS.
    else
      call advance_hypervis_scalar(edgeadv,elem,hvcoord,hybrid,deriv,tl%np1,np1_qdp,nets,nete,dt)
    endif

!    call t_stopf('prim_advec_tracers_remap_rk2')

  end subroutine prim_advec_tracers_remap_rk2

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine qdp_time_avg( elem , rkstage , n0_qdp , np1_qdp , hybrid , nets , nete )
    use hybrid_mod, only : hybrid_t, get_loop_ranges
    implicit none
    type(element_t)     , intent(inout) :: elem(:)
    integer             , intent(in   ) :: rkstage , n0_qdp , np1_qdp , nets , nete 
    type(hybrid_t) :: hybrid
    integer :: i,j,ie,q,k
    integer :: kbeg,kend,qbeg,qend
    real(kind=r8) :: rrkstage

    call get_loop_ranges(hybrid,kbeg=kbeg,kend=kend,qbeg=qbeg,qend=qend)

    rrkstage=1.0_r8/real(rkstage,kind=r8)
    do ie=nets,nete
      do q=qbeg,qend
        do k=kbeg,kend
          do j=1,np
          do i=1,np
           elem(ie)%state%Qdp(i,j,k,q,np1_qdp) =               &
               ( elem(ie)%state%Qdp(i,j,k,q,n0_qdp) + &
               (rkstage-1)*elem(ie)%state%Qdp(i,j,k,q,np1_qdp) ) * rrkstage
          enddo
          enddo
        enddo
      enddo
    enddo
  end subroutine qdp_time_avg

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine euler_step( np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier )
  ! ===================================
  ! This routine is the basic foward
  ! euler component used to construct RK SSP methods
  !
  !           u(np1) = u(n0) + dt2*DSS[ RHS(u(n0)) ]
  !
  ! n0 can be the same as np1.
  !
  ! DSSopt = DSSeta or DSSomega:   also DSS eta_dot_dpdn or omega
  !
  ! ===================================
  use dimensions_mod , only : np, nlev
  use hybrid_mod     , only : hybrid_t
  use hybrid_mod     , only : get_loop_ranges
  use element_mod    , only : element_t
  use derivative_mod , only : derivative_t, divergence_sphere, limiter_optim_iter_full
  use edge_mod       , only : edgevpack, edgevunpack
  use bndry_mod      , only : bndry_exchange
  use hybvcoord_mod  , only : hvcoord_t

  integer              , intent(in   )         :: np1_qdp, n0_qdp
  real (kind=r8), intent(in   )         :: dt
  type (element_t)     , intent(inout), target :: elem(:)
  type (hvcoord_t)     , intent(in   )         :: hvcoord
  type (hybrid_t)      , intent(in   )         :: hybrid
  type (derivative_t)  , intent(in   )         :: deriv
  integer              , intent(in   )         :: nets
  integer              , intent(in   )         :: nete
  integer              , intent(in   )         :: DSSopt
  integer              , intent(in   )         :: rhs_multiplier

  ! local
  real(kind=r8), dimension(np,np                       ) :: dpdiss
  real(kind=r8), dimension(np,np,nlev)                   :: dpdissk
  real(kind=r8), dimension(np,np,2                     ) :: gradQ
  real(kind=r8), dimension(np,np,2,nlev                ) :: Vstar
  real(kind=r8), dimension(np,np  ,nlev                ) :: Qtens
  real(kind=r8), dimension(np,np  ,nlev                ) :: dp
  real(kind=r8), dimension(np,np  ,nlev,qsize,nets:nete) :: Qtens_biharmonic
  real(kind=r8), dimension(np,np)                        :: div
  real(kind=r8), pointer, dimension(:,:,:)               :: DSSvar
  real(kind=r8) :: dp0(nlev),qmin_val(nlev),qmax_val(nlev)
  integer :: ie,q,i,j,k, kptr
  integer :: rhs_viss = 0
  integer :: kblk,qblk   ! The per thead size of the vertical and tracers
  integer :: kbeg, kend, qbeg, qend

  call get_loop_ranges(hybrid,kbeg=kbeg,kend=kend,qbeg=qbeg,qend=qend)

  kblk = kend - kbeg + 1   ! calculate size of the block of vertical levels
  qblk = qend - qbeg + 1   ! calculate size of the block of tracers

  do k = kbeg, kend
    dp0(k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
          ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
  enddo

!  call t_startf('euler_step')

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   compute Q min/max values for lim8
  !   compute biharmonic mixing term f
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  rhs_viss = 0
  if ( limiter_option == 8  ) then
    ! when running lim8, we also need to limit the biharmonic, so that term needs
    ! to be included in each euler step.  three possible algorithms here:
    ! 1) most expensive:
    !     compute biharmonic (which also computes qmin/qmax) during all 3 stages
    !     be sure to set rhs_viss=1
    !     cost:  3 biharmonic steps with 3 DSS
    !
    ! 2) cheapest:
    !     compute biharmonic (which also computes qmin/qmax) only on first stage
    !     be sure to set rhs_viss=3
    !     reuse qmin/qmax for all following stages (but update based on local qmin/qmax)
    !     cost:  1 biharmonic steps with 1 DSS
    !     main concern:  viscosity
    !
    ! 3)  compromise:
    !     compute biharmonic (which also computes qmin/qmax) only on last stage
    !     be sure to set rhs_viss=3
    !     compute qmin/qmax directly on first stage
    !     reuse qmin/qmax for 2nd stage stage (but update based on local qmin/qmax)
    !     cost:  1 biharmonic steps, 2 DSS
    !
    !  NOTE  when nu_p=0 (no dissipation applied in dynamics to dp equation), we should
    !        apply dissipation to Q (not Qdp) to preserve Q=1
    !        i.e.  laplace(Qdp) ~  dp0 laplace(Q)
    !        for nu_p=nu_q>0, we need to apply dissipation to Q * diffusion_dp
    !
    ! initialize dp, and compute Q from Qdp (and store Q in Qtens_biharmonic)
    do ie = nets, nete
      ! add hyperviscosity to RHS.  apply to Q at timelevel n0, Qdp(n0)/dp
      do k = kbeg, kend
        do j=1,np
        do i=1,np
          dp(i,j,k) = elem(ie)%derived%dp(i,j,k) - rhs_multiplier*dt*elem(ie)%derived%divdp_proj(i,j,k)
          do q = 1, qsize
            Qtens_biharmonic(i,j,k,q,ie) = elem(ie)%state%Qdp(i,j,k,q,n0_qdp)/dp(i,j,k)
          enddo
        enddo
        enddo
      enddo

      do q = 1, qsize
        do k= 1, nlev
          qmin_val(k) = +1.0e+24
          qmax_val(k) = -1.0e+24
          do j=1,np
            do i=1,np
!             Qtens_biharmonic(i,j,k,q,ie) = elem(ie)%state%Qdp(i,j,k,q,n0_qdp)/dp(i,j,k)
              qmin_val(k) = min(qmin_val(k),Qtens_biharmonic(i,j,k,q,ie))
              qmax_val(k) = max(qmax_val(k),Qtens_biharmonic(i,j,k,q,ie))
            enddo
          enddo

          if ( rhs_multiplier == 1 ) then
              qmin(k,q,ie)=min(qmin(k,q,ie),qmin_val(k))
              qmin(k,q,ie)=max(qmin(k,q,ie),0d0)
              qmax(k,q,ie)=max(qmax(k,q,ie),qmax_val(k))
          else
              qmin(k,q,ie)=max(qmin_val(k),0d0)
              qmax(k,q,ie)=qmax_val(k)
          endif
        enddo
      enddo
    enddo

    ! compute element qmin/qmax
    if ( rhs_multiplier == 0 ) then
      ! update qmin/qmax based on neighbor data for lim8
!      call t_startf('euler_neighbor_minmax1')
      call neighbor_minmax(hybrid,edgeAdvQminmax,nets,nete,qmin(:,:,nets:nete),qmax(:,:,nets:nete))
!      call t_stopf('euler_neighbor_minmax1')
    endif

    ! get niew min/max values, and also compute biharmonic mixing term
    if ( rhs_multiplier == 2 ) then
      rhs_viss = 3
      ! two scalings depending on nu_p:
      ! nu_p=0:    qtens_biharmonic *= dp0                   (apply viscsoity only to q)
      ! nu_p>0):   qtens_biharmonc *= elem()%psdiss_ave      (for consistency, if nu_p=nu_q)
      if ( nu_p > 0 ) then
        do ie = nets, nete
          do k = kbeg, kend
            do j=1,np
            do i=1,np
               dpdiss(i,j) = elem(ie)%derived%dpdiss_ave(i,j,k)
            enddo
            enddo

            do q = qbeg,qend
              ! NOTE: divide by dp0 since we multiply by dp0 below
              do j=1,np
              do i=1,np
                 Qtens_biharmonic(i,j,k,q,ie)=Qtens_biharmonic(i,j,k,q,ie)*dpdiss(i,j)/dp0(k)
              enddo
              enddo
            enddo
          enddo
        enddo
      endif

!   Previous version of biharmonic_wk_scalar_minmax included a min/max
!   calculation into the boundary exchange.  This was causing cache issues.
!   Split the single operation into two separate calls
!      call neighbor_minmax()
!      call biharmonic_wk_scalar()
!
#ifdef OVERLAP
      call neighbor_minmax_start(hybrid,edgeAdvQminmax,nets,nete,qmin(:,:,nets:nete),qmax(:,:,nets:nete))
      call biharmonic_wk_scalar(elem,qtens_biharmonic,deriv,edgeAdv,hybrid,nets,nete)
      do ie = nets, nete
        do q = qbeg, qend
          do k = kbeg, kend
            !OMP_COLLAPSE_SIMD 
            !DIR_VECTOR_ALIGNED
            do j=1,np
            do i=1,np
               ! note: biharmonic_wk() output has mass matrix already applied. Un-apply since we apply again below:
               qtens_biharmonic(i,j,k,q,ie) = &
                     -rhs_viss*dt*nu_q*dp0(k)*Qtens_biharmonic(i,j,k,q,ie) / elem(ie)%spheremp(i,j)
            enddo
            enddo
          enddo
        enddo
      enddo
      call neighbor_minmax_finish(hybrid,edgeAdvQminmax,nets,nete,qmin(:,:,nets:nete),qmax(:,:,nets:nete))
#else
      call t_startf('euler_neighbor_minmax2')
      call neighbor_minmax(hybrid,edgeAdvQminmax,nets,nete,qmin(:,:,nets:nete),qmax(:,:,nets:nete))
      call t_stopf('euler_neighbor_minmax2')
      call biharmonic_wk_scalar(elem,qtens_biharmonic,deriv,edgeAdv,hybrid,nets,nete)

      do ie = nets, nete
        do q = qbeg, qend
          do k = kbeg, kend
            !OMP_COLLAPSE_SIMD 
            !DIR_VECTOR_ALIGNED
            do j=1,np
            do i=1,np
               ! note: biharmonic_wk() output has mass matrix already applied. Un-apply since we apply again below:
               qtens_biharmonic(i,j,k,q,ie) = &
                     -rhs_viss*dt*nu_q*dp0(k)*Qtens_biharmonic(i,j,k,q,ie) / elem(ie)%spheremp(i,j)
            enddo
            enddo
          enddo
        enddo
      enddo
#endif


    endif
  endif  ! compute biharmonic mixing term and qmin/qmax


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   2D Advection step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ie = nets, nete


    ! Compute velocity used to advance Qdp
    do k = kbeg, kend
      ! derived variable divdp_proj() (DSS'd version of divdp) will only be correct on 2nd and 3rd stage
      ! but that's ok because rhs_multiplier=0 on the first stage:
      do j=1,np
      do i=1,np
         dp(i,j,k) = elem(ie)%derived%dp(i,j,k) - rhs_multiplier * dt * elem(ie)%derived%divdp_proj(i,j,k)
         Vstar(i,j,1,k) = elem(ie)%derived%vn0(i,j,1,k) / dp(i,j,k)
         Vstar(i,j,2,k) = elem(ie)%derived%vn0(i,j,2,k) / dp(i,j,k)
      enddo
      enddo
    enddo
    if ( limiter_option == 8) then
        ! Note that the term dpdissk is independent of Q
        do k = kbeg, kend
          ! UN-DSS'ed dp at timelevel n0+1:
          do j=1,np
          do i=1,np
             dpdissk(i,j,k) = dp(i,j,k) - dt * elem(ie)%derived%divdp(i,j,k)
          enddo
          enddo
          if ( nu_p > 0 .and. rhs_viss /= 0 ) then
            ! add contribution from UN-DSS'ed PS dissipation
!            dpdiss(:,:) = ( hvcoord%hybi(k+1) - hvcoord%hybi(k) ) *
!            elem(ie)%derived%psdiss_biharmonic(:,:)


            do j=1,np
            do i=1,np
               dpdiss(i,j) = elem(ie)%derived%dpdiss_biharmonic(i,j,k)
               dpdissk(i,j,k) = dpdissk(i,j,k) - rhs_viss * dt * nu_q * dpdiss(i,j) / elem(ie)%spheremp(i,j)
            enddo
            enddo
          endif
        enddo
    endif  ! limiter == 8


    ! advance Qdp
    do q = qbeg, qend
      do k = kbeg, kend
        ! div( U dp Q),
        do j=1,np
        do i=1,np
           gradQ(i,j,1) = Vstar(i,j,1,k) * elem(ie)%state%Qdp(i,j,k,q,n0_qdp)
           gradQ(i,j,2) = Vstar(i,j,2,k) * elem(ie)%state%Qdp(i,j,k,q,n0_qdp)
        enddo
        enddo


        call divergence_sphere( gradQ , deriv , elem(ie),div )

        do j=1,np
          do i=1,np
            Qtens(i,j,k) = elem(ie)%state%Qdp(i,j,k,q,n0_qdp) - dt * div(i,j)
          enddo
        enddo

        ! optionally add in hyperviscosity computed above:
        if ( rhs_viss /= 0 ) then 
          do j=1,np
          do i=1,np
              Qtens(i,j,k) = Qtens(i,j,k) + Qtens_biharmonic(i,j,k,q,ie)
          enddo
          enddo
        endif
      enddo

      if ( limiter_option == 8) then
        ! apply limiter to Q = Qtens / dp_star
        call limiter_optim_iter_full( Qtens(:,:,:) , elem(ie)%spheremp(:,:) , qmin(:,q,ie) , &
                                                          qmax(:,q,ie) , dpdissk)
      endif


      ! apply mass matrix, overwrite np1 with solution:
      ! dont do this earlier, since we allow np1_qdp == n0_qdp
      ! and we dont want to overwrite n0_qdp until we are done using it
      do k = kbeg, kend
        do j=1,np
        do i=1,np
            elem(ie)%state%Qdp(i,j,k,q,np1_qdp) = elem(ie)%spheremp(i,j) * Qtens(i,j,k)
        enddo
        enddo
      enddo

      if ( limiter_option == 4 ) then
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        ! sign-preserving limiter, applied after mass matrix
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!JMD !$OMP BARRIER
!JMD !$OMP MASTER
        call limiter2d_zero(elem(ie)%state%Qdp(:,:,:,q,np1_qdp))
!JMD !$OMP END MASTER
!JMD !$OMP BARRIER
      endif

      kptr = nlev*(q-1) + kbeg - 1
      call edgeVpack(edgeAdvp1 , elem(ie)%state%Qdp(:,:,kbeg:kend,q,np1_qdp) , kblk , kptr , ie )
    enddo
    ! only perform this operation on thread which owns the first tracer
    if (DSSopt>0) then

        if ( DSSopt == DSSeta         ) DSSvar => elem(ie)%derived%eta_dot_dpdn(:,:,:)
        if ( DSSopt == DSSomega       ) DSSvar => elem(ie)%derived%omega_p(:,:,:)
        if ( DSSopt == DSSdiv_vdp_ave ) DSSvar => elem(ie)%derived%divdp_proj(:,:,:)
        ! also DSS extra field
        do k = kbeg, kend
          do j=1,np
          do i=1,np
             DSSvar(i,j,k) = elem(ie)%spheremp(i,j) * DSSvar(i,j,k)
          enddo
          enddo
        enddo
        kptr = nlev*qsize + kbeg - 1
        call edgeVpack( edgeAdvp1 , DSSvar(:,:,kbeg:kend), kblk, kptr, ie)

    end if
  enddo

  call bndry_exchange( hybrid , edgeAdvp1,location='edgeAdvp1')

  do ie = nets, nete
    ! only perform this operation on thread which owns the first tracer
    if (DSSopt>0) then

        if ( DSSopt == DSSeta         ) DSSvar => elem(ie)%derived%eta_dot_dpdn(:,:,:)
        if ( DSSopt == DSSomega       ) DSSvar => elem(ie)%derived%omega_p(:,:,:)
        if ( DSSopt == DSSdiv_vdp_ave ) DSSvar => elem(ie)%derived%divdp_proj(:,:,:)
        kptr = qsize*nlev + kbeg -1
        call edgeVunpack( edgeAdvp1 , DSSvar(:,:,kbeg:kend) , kblk , kptr , ie )
        do k = kbeg, kend
           !OMP_COLLAPSE_SIMD 
           !DIR_VECTOR_ALIGNED
           do j=1,np
           do i=1,np
               DSSvar(i,j,k) = DSSvar(i,j,k) * elem(ie)%rspheremp(i,j)
           enddo
           enddo
        enddo

    end if
    do q = qbeg, qend
      kptr = nlev*(q-1) + kbeg - 1
      call edgeVunpack( edgeAdvp1 , elem(ie)%state%Qdp(:,:,kbeg:kend,q,np1_qdp) , kblk , kptr , ie )
        do k = kbeg, kend
          !OMP_COLLAPSE_SIMD 
          !DIR_VECTOR_ALIGNED
          do j=1,np
          do i=1,np
             elem(ie)%state%Qdp(i,j,k,q,np1_qdp) = elem(ie)%rspheremp(i,j) * elem(ie)%state%Qdp(i,j,k,q,np1_qdp)
          enddo
          enddo
        enddo
    enddo
  enddo
!  call t_stopf('euler_step')

  end subroutine euler_step



  subroutine limiter2d_zero(Q)
  ! mass conserving zero limiter (2D only).  to be called just before DSS
  !
  ! this routine is called inside a DSS loop, and so Q had already
  ! been multiplied by the mass matrix.  Thus dont include the mass
  ! matrix when computing the mass = integral of Q over the element
  !
  ! ps is only used when advecting Q instead of Qdp
  ! so ps should be at one timelevel behind Q
  implicit none
  real (kind=r8), intent(inout) :: Q(np,np,nlev)

  ! local
!  real (kind=r8) :: dp(np,np)
  real (kind=r8) :: mass,mass_new,ml
  integer i,j,k

  do k = nlev , 1 , -1
    mass = 0
    do j = 1 , np
      do i = 1 , np
        !ml = Q(i,j,k)*dp(i,j)*spheremp(i,j)  ! see above
        ml = Q(i,j,k)
        mass = mass + ml
      enddo
    enddo

    ! negative mass.  so reduce all postive values to zero
    ! then increase negative values as much as possible
    if ( mass < 0 ) Q(:,:,k) = -Q(:,:,k)
    mass_new = 0
    do j = 1 , np
      do i = 1 , np
        if ( Q(i,j,k) < 0 ) then
          Q(i,j,k) = 0
        else
          ml = Q(i,j,k)
          mass_new = mass_new + ml
        endif
      enddo
    enddo

    ! now scale the all positive values to restore mass
    if ( mass_new > 0 ) Q(:,:,k) = Q(:,:,k) * abs(mass) / mass_new
    if ( mass     < 0 ) Q(:,:,k) = -Q(:,:,k)
  enddo
  end subroutine limiter2d_zero

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine advance_hypervis_scalar( edgeAdv , elem , hvcoord , hybrid , deriv , nt , nt_qdp , nets , nete , dt2 )
  !  hyperviscsoity operator for foward-in-time scheme
  !  take one timestep of:
  !          Q(:,:,:,np) = Q(:,:,:,np) +  dt2*nu*laplacian**order ( Q )
  !
  !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
  use dimensions_mod , only : np, nlev
  use hybrid_mod     , only : hybrid_t!, PrintHybrid
  use hybrid_mod     , only : get_loop_ranges
  use element_mod    , only : element_t
  use derivative_mod , only : derivative_t
  use edge_mod       , only : edgevpack, edgevunpack
  use edgetype_mod   , only : EdgeBuffer_t
  use bndry_mod      , only : bndry_exchange

  implicit none
  type (EdgeBuffer_t)  , intent(inout)         :: edgeAdv
  type (element_t)     , intent(inout), target :: elem(:)
  type (hvcoord_t)     , intent(in   )         :: hvcoord
  type (hybrid_t)      , intent(in   )         :: hybrid
  type (derivative_t)  , intent(in   )         :: deriv
  integer              , intent(in   )         :: nt
  integer              , intent(in   )         :: nt_qdp
  integer              , intent(in   )         :: nets
  integer              , intent(in   )         :: nete
  real (kind=r8), intent(in   )         :: dt2

  ! local
  real (kind=r8), dimension(np,np,nlev,qsize,nets:nete) :: Qtens
  real (kind=r8), dimension(np,np,nlev                ) :: dp
!  real (kind=r8), dimension(      nlev,qsize,nets:nete) :: min_neigh
!  real (kind=r8), dimension(      nlev,qsize,nets:nete) :: max_neigh
  integer :: k,kptr,ie,ic,q,i,j
  integer :: kbeg,kend,qbeg,qend

! NOTE: PGI compiler bug: when using spheremp, rspheremp and ps as pointers to elem(ie)% members,
!       data is incorrect (offset by a few numbers actually)
!       removed for now.
!  real (kind=r8), dimension(:,:), pointer :: spheremp,rspheremp
  real (kind=r8) :: dt,dp0
  integer :: kblk,qblk   ! The per thead size of the vertical and tracers

  call get_loop_ranges(hybrid,kbeg=kbeg,kend=kend,qbeg=qbeg,qend=qend)

  if ( nu_q           == 0 ) return
  if ( hypervis_order /= 2 ) return

  kblk = kend - kbeg + 1   ! calculate size of the block of vertical levels
  qblk = qend - qbeg + 1   ! calculate size of the block of tracers

  call t_startf('advance_hypervis_scalar')

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  hyper viscosity
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dt = dt2 / hypervis_subcycle_q

  do ic = 1 , hypervis_subcycle_q
    do ie = nets, nete
      ! Qtens = Q/dp   (apply hyperviscsoity to dp0 * Q, not Qdp)
      do k = kbeg, kend
         ! various options:
         !   1)  biharmonic( Qdp )
         !   2)  dp0 * biharmonic( Qdp/dp )
         !   3)  dpave * biharmonic(Q/dp)
         ! For trace mass / mass consistenciy, we use #2 when nu_p=0
         ! and #e when nu_p>0, where dpave is the mean mass flux from the nu_p
         ! contribution from dynamics.
         dp0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) ) * hvcoord%ps0 + &
              ( hvcoord%hybi(k+1) - hvcoord%hybi(k) ) * hvcoord%ps0
         dp(:,:,k) = elem(ie)%derived%dp(:,:,k) - dt2*elem(ie)%derived%divdp_proj(:,:,k)
         if (nu_p>0) then
            do q = qbeg, qend
               Qtens(:,:,k,q,ie) = elem(ie)%derived%dpdiss_ave(:,:,k)*&
                    elem(ie)%state%Qdp(:,:,k,q,nt_qdp) / dp(:,:,k)
            enddo
         else
            do q = qbeg, qend
               Qtens(:,:,k,q,ie) = dp0*elem(ie)%state%Qdp(:,:,k,q,nt_qdp) / dp(:,:,k)
            enddo
         endif
      enddo
   enddo

    ! compute biharmonic operator. Qtens = input and output
    call biharmonic_wk_scalar( elem , Qtens , deriv , edgeAdv , hybrid ,  nets , nete )

    do ie = nets, nete
      !spheremp     => elem(ie)%spheremp
      do q = qbeg, qend
        do k = kbeg, kend
          dp0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) ) * hvcoord%ps0 + &
                ( hvcoord%hybi(k+1) - hvcoord%hybi(k) ) * hvcoord%ps0
         do j = 1 , np
            do i = 1 , np

              ! advection Qdp.  For mass advection consistency:
              ! DIFF( Qdp) ~   dp0 DIFF (Q)  =  dp0 DIFF ( Qdp/dp )
              elem(ie)%state%Qdp(i,j,k,q,nt_qdp) = elem(ie)%state%Qdp(i,j,k,q,nt_qdp) * elem(ie)%spheremp(i,j) &
                                                                                 - dt * nu_q * Qtens(i,j,k,q,ie)
        enddo
        enddo
        enddo





           ! smooth some of the negativities introduced by diffusion:
           call limiter2d_zero( elem(ie)%state%Qdp(:,:,:,q,nt_qdp) )
      enddo
      do q = qbeg, qend
         kptr = nlev*(q-1) + kbeg - 1
         call edgeVpack( edgeAdv , elem(ie)%state%Qdp(:,:,kbeg:kend,q,nt_qdp) , kblk, kptr, ie )
      enddo
    enddo

    call bndry_exchange( hybrid , edgeAdv,location='advance_hypervis_scalar')

    do ie = nets, nete
      do q = qbeg, qend
         kptr = nlev*(q-1) + kbeg - 1
         call edgeVunpack( edgeAdv , elem(ie)%state%Qdp(:,:,kbeg:kend,q,nt_qdp) , kblk, kptr, ie )
      enddo
      !rspheremp     => elem(ie)%rspheremp
      do q = qbeg, qend
        ! apply inverse mass matrix
        do k = kbeg, kend
          elem(ie)%state%Qdp(:,:,k,q,nt_qdp) = elem(ie)%rspheremp(:,:) * elem(ie)%state%Qdp(:,:,k,q,nt_qdp)
        enddo
      enddo
    enddo

  enddo
  call t_stopf('advance_hypervis_scalar')
  end subroutine advance_hypervis_scalar

end module prim_advection_mod
