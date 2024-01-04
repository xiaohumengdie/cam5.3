#define _DBG_
module prim_driver_mod
  use shr_kind_mod,           only: r8=>shr_kind_r8
  use cam_logfile,            only: iulog
  use cam_abortutils,         only: endrun
  use dimensions_mod,         only: np, nlev, nlevp, nelem, nelemd, nelemdmax, GlobalUniqueCols, qsize, nc,nhc
  use hybrid_mod,             only: hybrid_t
  use derivative_mod,         only: derivative_t
  use quadrature_mod,         only : quadrature_t, test_gauss, test_gausslobatto, gausslobatto
  use derivative_mod, only : derivative_t
  use reduction_mod, only : reductionbuffer_ordered_1d_t, red_min, red_max, &
         red_sum, red_sum_int, red_flops, initreductionbuffer
  use element_mod, only : element_t, timelevels,  allocate_element_desc
  use thread_mod, only : omp_get_num_threads
  implicit none
  private
  public :: prim_init2 , prim_run_subcycle, prim_finalize
  public :: smooth_topo_datasets



contains

!=============================================================================!

  subroutine prim_init2(elem, hybrid, nets, nete, tl, hvcoord)

    use parallel_mod, only : parallel_t, syncmp
    use time_mod, only : timelevel_t, tstep, phys_tscale, timelevel_init, nendstep, smooth, nsplit, TimeLevel_Qdp
    use control_mod, only : runtype, &
         topology,columnpackage, moisture, rsplit, qsplit, rk_stage_user,&
         limiter_option, nu, nu_q, nu_div, tstep_type, hypervis_subcycle, &
         hypervis_subcycle_q
    use control_mod, only : tracer_transport_type
    use control_mod, only : TRACERTRANSPORT_SE_GLL
    use prim_si_ref_mod, only: prim_set_mass
    use thread_mod, only : max_num_threads, omp_get_thread_num
    use derivative_mod, only : derivinit, v2pinit
    use global_norms_mod, only : test_global_integral, print_cfl
    use hybvcoord_mod, only : hvcoord_t
    use prim_advection_mod, only: prim_advec_init2, deriv

    type (element_t), intent(inout) :: elem(:)
    type (hybrid_t), intent(in) :: hybrid

    type (TimeLevel_t), intent(inout)    :: tl              ! time level struct
    type (hvcoord_t), intent(inout)      :: hvcoord         ! hybrid vertical coordinate struct

     integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete  ! ending thread element number   (private)


    ! ==================================
    ! Local variables
    ! ==================================

    real (kind=r8) :: dt              ! "timestep dependent" timestep
!   variables used to calculate CFL
    real (kind=r8) :: dtnu            ! timestep*viscosity parameter
    real (kind=r8) :: dt_dyn_vis      ! viscosity timestep used in dynamics
    real (kind=r8) :: dt_tracer_vis      ! viscosity timestep used in tracers

    real (kind=r8) :: dp


    real (kind=r8) :: ps(np,np)       ! surface pressure

    character(len=80)     :: fname
    character(len=8)      :: njusn
    character(len=4)      :: charnum

    integer :: simday
    integer :: i,j,k,ie,iptr,t,q
    integer :: ierr
    integer :: nfrc
    integer :: n0_qdp

    ! ==========================
    ! begin executable code
    ! ==========================
    if (topology == "cube") then
       call test_global_integral(elem, hybrid,nets,nete)
    end if


    ! compute most restrictive dt*nu for use by variable res viscosity:
    if (tstep_type == 0) then
       ! LF case: no tracers, timestep seen by viscosity is 2*tstep
       dt_tracer_vis = 0
       dt_dyn_vis = 2*tstep
       dtnu = 2.0d0*tstep*max(nu,nu_div)
    else
       ! compute timestep seen by viscosity operator:
       dt_dyn_vis = tstep
       if (qsplit>1 .and. tstep_type == 1) then
          ! tstep_type==1: RK2 followed by LF.  internal LF stages apply viscosity at 2*dt
          dt_dyn_vis = 2*tstep
       endif
       dt_tracer_vis=tstep*qsplit

       ! compute most restrictive condition:
       ! note: dtnu ignores subcycling
       dtnu=max(dt_dyn_vis*max(nu,nu_div), dt_tracer_vis*nu_q)
       ! compute actual viscosity timesteps with subcycling
       dt_tracer_vis = dt_tracer_vis/hypervis_subcycle_q
       dt_dyn_vis = dt_dyn_vis/hypervis_subcycle
    endif


    ! ==================================
    ! Initialize derivative structure
    ! ==================================
    call Prim_Advec_Init2(hybrid)

    if (hybrid%ithr==0) then
       call syncmp(hybrid%par)
    end if

    if (topology /= "cube") then
       call endrun('Error: only cube topology supported for primaitve equations')
    endif

    ! For new runs, and branch runs, convert state variable to (Qdp)
    ! because initial conditon reads in Q, not Qdp
    ! restart runs will read dpQ from restart file
    ! need to check what CAM does on a branch run
    if (runtype==0 .or. runtype==2) then
       do ie=nets,nete
          elem(ie)%derived%omega_p(:,:,:) = 0D0
       end do
       do ie=nets,nete
#if (defined COLUMN_OPENMP)
!$omp parallel do default(shared), private(k, t, q, i, j, dp)
#endif
          do k=1,nlev    !  Loop inversion (AAM)
             do q=1,qsize
                do i=1,np
                   do j=1,np
                      dp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                           ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,tl%n0)
                      
                      elem(ie)%state%Qdp(i,j,k,q,1)=elem(ie)%state%Q(i,j,k,q)*dp
                      elem(ie)%state%Qdp(i,j,k,q,2)=elem(ie)%state%Q(i,j,k,q)*dp
                      
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif

    ! for restart runs, we read in Qdp for exact restart, and rederive Q
    if (runtype==1) then
       call TimeLevel_Qdp( tl, qsplit, n0_qdp)
       do ie=nets,nete
#if (defined COLUMN_OPENMP)
!$omp parallel do default(shared), private(k, t, q, i, j, dp)
#endif
          do k=1,nlev    !  Loop inversion (AAM)
             do t=tl%n0,tl%n0
                do q=1,qsize
                   do i=1,np
                      do j=1,np
                         dp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                              ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,t)
                         elem(ie)%state%Q(i,j,k,q)=elem(ie)%state%Qdp(i,j,k,q, n0_qdp)/dp
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif


    ! timesteps to use for advective stability:  tstep*qsplit and tstep
    call print_cfl(elem,hybrid,nets,nete,dtnu)

    if (hybrid%masterthread) then
       ! CAM has set tstep based on dtime before calling prim_init2(),
       ! so only now does HOMME learn the timstep.  print them out:
       write(iulog,'(a,2f9.2)') "dt_remap: (0=disabled)   ",tstep*qsplit*rsplit

       if (tracer_transport_type == TRACERTRANSPORT_SE_GLL) then
          write(iulog,'(a,2f9.2)') "dt_tracer, per RK stage: ",tstep*qsplit,(tstep*qsplit)/(rk_stage_user-1)
       end if
       write(iulog,'(a,2f9.2)')    "dt_dyn:                  ",tstep
       write(iulog,'(a,2f9.2)')    "dt_dyn (viscosity):      ",dt_dyn_vis
       write(iulog,'(a,2f9.2)')    "dt_tracer (viscosity):   ",dt_tracer_vis


       if (phys_tscale/=0) then
          write(iulog,'(a,2f9.2)') "CAM physics timescale:       ",phys_tscale
       endif
       write(iulog,'(a,2f9.2)') "CAM dtime (dt_phys):         ",tstep*nsplit*qsplit*max(rsplit,1)
    end if


    if (hybrid%masterthread) write(iulog,*) "initial state:"

  end subroutine prim_init2

!=======================================================================================================!

  subroutine prim_run_subcycle(elem, hybrid,nets,nete, dt, tl, hvcoord,nsubstep)
!
!   advance all variables (u,v,T,ps,Q,C) from time t to t + dt_q
!
!   input:
!       tl%nm1   not used
!       tl%n0    data at time t
!       tl%np1   new values at t+dt_q
!
!   then we update timelevel pointers:
!       tl%nm1 = tl%n0
!       tl%n0  = tl%np1
!   so that:
!       tl%nm1   tracers:  t    dynamics:  t+(qsplit-1)*dt
!       tl%n0    time t + dt_q
!
!
    use hybvcoord_mod, only : hvcoord_t
    use time_mod, only : TimeLevel_t, timelevel_update, timelevel_qdp, nsplit
    use control_mod, only: statefreq,&
           energy_fixer, ftype, qsplit, rsplit
    use prim_advance_mod, only : applycamforcing, &
                                 applycamforcing_dynamics
    use reduction_mod, only : parallelmax
    use prim_advection_mod, only : vertical_remap


    type (element_t) , intent(inout)        :: elem(:)

    type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)

    type (hvcoord_t), intent(in)      :: hvcoord         ! hybrid vertical coordinate struct

    integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete  ! ending thread element number   (private)
    real(kind=r8), intent(in)        :: dt  ! "timestep dependent" timestep
    type (TimeLevel_t), intent(inout)       :: tl
    integer, intent(in)                     :: nsubstep  ! nsubstep = 1 .. nsplit
    real(kind=r8) :: st, st1, dp, dt_q, dt_remap
    integer :: ie, t, q,k,i,j,n, n_Q
    integer :: n0_qdp,np1_qdp,r, nstep_end

    real (kind=r8)                          :: maxcflx, maxcfly
    real (kind=r8) :: dp_np1(np,np)

    ! ===================================
    ! Main timestepping loop
    ! ===================================
    dt_q = dt*qsplit
    dt_remap = dt_q
    nstep_end = tl%nstep + qsplit
    if (rsplit>0) then
       dt_remap=dt_q*rsplit   ! rsplit=0 means use eulerian code, not vert. lagrange
       nstep_end = tl%nstep + qsplit*rsplit  ! nstep at end of this routine
    endif

    ! ftype=2  Q was adjusted by physics, but apply u,T forcing here
    ! ftype=1  forcing was applied time-split in CAM coupling layer
    ! ftype=0 means forcing apply here
    ! ftype=-1 do not apply forcing
    call TimeLevel_Qdp( tl, qsplit, n0_qdp)
    if (ftype==0) then
      call ApplyCAMForcing(elem, hvcoord,tl%n0,n0_qdp, dt_remap,nets,nete)
    end if
    if (ftype==2) call ApplyCAMForcing_dynamics(elem, hvcoord,tl%n0,dt_remap,nets,nete)

    ! initialize dp3d from ps
    if (rsplit>0) then
    do ie=nets,nete
       do k=1,nlev
          elem(ie)%state%dp3d(:,:,k,tl%n0)=&
               ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,tl%n0)
       enddo
       ! DEBUGDP step: ps_v should not be used for rsplit>0 code during prim_step
       ! vertical_remap.  so to this for debugging:
       elem(ie)%state%ps_v(:,:,tl%n0)=-9e9
    enddo
    endif


    ! loop over rsplit vertically lagrangian timesteps
    call prim_step(elem, hybrid,nets,nete, dt, tl, hvcoord)
    do r=2,rsplit
       call TimeLevel_update(tl,"leapfrog")
       call prim_step(elem, hybrid,nets,nete, dt, tl, hvcoord)
    enddo
    ! defer final timelevel update until after remap and diagnostics


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !  apply vertical remap
    !  always for tracers
    !  if rsplit>0:  also remap dynamics and compute reference level ps_v
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !compute timelevels for tracers (no longer the same as dynamics)
    call TimeLevel_Qdp( tl, qsplit, n0_qdp, np1_qdp)
    call vertical_remap(hybrid,elem,hvcoord,dt_remap,tl%np1,np1_qdp,nets,nete)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! time step is complete.  update some diagnostic variables:
    ! lnps (we should get rid of this)
    ! Q    (mixing ratio)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do ie=nets,nete
       elem(ie)%state%lnps(:,:,tl%np1)= LOG(elem(ie)%state%ps_v(:,:,tl%np1))
#if (defined COLUMN_OPENMP)
       !$omp parallel do default(shared), private(k,q,dp_np1)
#endif
       do k=1,nlev    !  Loop inversion (AAM)
          !if (k == 1) then
           !write(*,*) "In prim run there are ", omp_get_num_threads(), " in the current team in parallel region"
          !endif
          dp_np1(:,:) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,tl%np1)
          do q=1,qsize
             elem(ie)%state%Q(:,:,k,q)=elem(ie)%state%Qdp(:,:,k,q,np1_qdp)/dp_np1(:,:)
          enddo
       enddo
    enddo

    ! now we have:
    !   u(nm1)   dynamics at  t+dt_remap - 2*dt
    !   u(n0)    dynamics at  t+dt_remap - dt
    !   u(np1)   dynamics at  t+dt_remap
    !
    !   Q(1)   Q at t+dt_remap

    if (energy_fixer > 0) then
       call prim_energy_fixer(elem,hvcoord,hybrid,tl,nets,nete,nsubstep)
    endif

    ! =================================
    ! update dynamics time level pointers
    ! =================================
    call TimeLevel_update(tl,"leapfrog")

    ! now we have:
    !   u(nm1)   dynamics at  t+dt_remap - dt       (Robert-filtered)
    !   u(n0)    dynamics at  t+dt_remap
    !   u(np1)   undefined


  end subroutine prim_run_subcycle

  subroutine prim_step(elem, hybrid,nets,nete, dt, tl, hvcoord)
!
!   Take qsplit dynamics steps and one tracer step
!   for vertically lagrangian option, this subroutine does only the horizontal step
!
!   input:
!       tl%nm1   not used
!       tl%n0    data at time t
!       tl%np1   new values at t+dt_q
!
!   then we update timelevel pointers:
!       tl%nm1 = tl%n0
!       tl%n0  = tl%np1
!   so that:
!       tl%nm1   tracers:  t    dynamics:  t+(qsplit-1)*dt
!       tl%n0    time t + dt_q
!
!
    use hybvcoord_mod, only : hvcoord_t
    use time_mod, only : TimeLevel_t, timelevel_update, nsplit
    use control_mod, only: statefreq, ftype, qsplit, nu_p, rsplit
    use control_mod, only : tracer_transport_type
    use control_mod, only : tracer_grid_type, TRACER_GRIDTYPE_GLL
    use prim_advance_mod, only : prim_advance_exp
    use prim_advection_mod, only : prim_advec_tracers_remap, deriv
    use reduction_mod, only : parallelmax
    use time_mod,    only : time_at

    type (element_t) , intent(inout)        :: elem(:)

    type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)

    type (hvcoord_t), intent(in)      :: hvcoord         ! hybrid vertical coordinate struct

    integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete  ! ending thread element number   (private)
    real(kind=r8), intent(in)        :: dt  ! "timestep dependent" timestep
    type (TimeLevel_t), intent(inout)       :: tl
    real(kind=r8) :: st, st1, dp, dt_q
    integer :: ie, t, q,k,i,j,n, n_Q

    real (kind=r8)                          :: maxcflx, maxcfly

    real (kind=r8) :: dp_np1(np,np)

    dt_q = dt*qsplit

    ! ===============
    ! initialize mean flux accumulation variables and save some variables at n0
    ! for use by advection
    ! ===============
    do ie=nets,nete
      elem(ie)%derived%eta_dot_dpdn=0     ! mean vertical mass flux
      elem(ie)%derived%vn0=0              ! mean horizontal mass flux
      elem(ie)%derived%omega_p=0
      if (nu_p>0) then
         elem(ie)%derived%dpdiss_ave=0
         elem(ie)%derived%dpdiss_biharmonic=0
      endif

      if (rsplit==0) then
        ! save dp at time t for use in tracers
#if (defined COLUMN_OPENMP)
!$omp parallel do default(shared), private(k)
#endif
         do k=1,nlev
            elem(ie)%derived%dp(:,:,k)=&
                 ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                 ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,tl%n0)
         enddo
      else
         ! dp at time t:  use floating lagrangian levels:
         elem(ie)%derived%dp(:,:,:)=elem(ie)%state%dp3d(:,:,:,tl%n0)
      endif
    enddo

    ! ===============
    ! Dynamical Step
    ! ===============
    n_Q = tl%n0  ! n_Q = timelevel of FV tracers at time t.  need to save this
                 ! FV tracers still carry 3 timelevels
                 ! SE tracers only carry 2 timelevels
    call prim_advance_exp(elem, deriv(hybrid%ithr), hvcoord,   &
         hybrid, dt, tl, nets, nete)
    do n=2,qsplit
       call TimeLevel_update(tl,"leapfrog")
       call prim_advance_exp(elem, deriv(hybrid%ithr), hvcoord,   &
            hybrid, dt, tl, nets, nete)
       ! defer final timelevel update until after Q update.
    enddo
    ! current dynamics state variables:
    !    derived%dp              =  dp at start of timestep
    !    derived%vstar           =  velocity at start of tracer timestep
    !    derived%vn0             =  mean horiz. flux:   U*dp
    ! rsplit=0
    !        state%v(:,:,:,np1)      = velocity on reference levels
    !        state%ps_v(:,:,:,np1)   = ps
    ! rsplit>0
    !        state%v(:,:,:,np1)      = velocity on lagrangian levels 
    !        state%dp3d(:,:,:,np1)   = dp3d
    !


    ! ===============
    ! Tracer Advection.  
    ! in addition, this routine will apply the DSS to:
    !        derived%eta_dot_dpdn    =  mean vertical velocity (used for remap below)
    !        derived%omega           =
    ! Tracers are always vertically lagrangian.  
    ! For rsplit=0: 
    !   if tracer scheme needs v on lagrangian levels it has to vertically interpolate
    !   if tracer scheme needs dp3d, it needs to derive it from ps_v
    ! ===============
    if (tracer_grid_type == TRACER_GRIDTYPE_GLL) then
      call Prim_Advec_Tracers_remap(elem, deriv(hybrid%ithr),hvcoord,hybrid,&
           dt_q,tl,nets,nete)
    else
      stop "erro tracer"
    endif

  end subroutine prim_step


!=======================================================================================================!


  subroutine prim_finalize(hybrid)
    type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)

    ! ==========================
    ! end of the hybrid program
    ! ==========================
  end subroutine prim_finalize



!=======================================================================================================!
  subroutine prim_energy_fixer(elem,hvcoord,hybrid,tl,nets,nete,nsubstep)
!
! non-subcycle code:
!  Solution is given at times u(t-1),u(t),u(t+1)
!  E(n=1) = energy before dynamics
!  E(n=2) = energy after dynamics
!
!  fixer will add a constant to the temperature so E(n=2) = E(n=1)
!
    use parallel_mod, only: global_shared_buf, global_shared_sum
    use hybvcoord_mod, only : hvcoord_t
    use physical_constants, only : Cp
    use time_mod, only : timelevel_t
    use control_mod, only : use_cpstar, energy_fixer
    use hybvcoord_mod, only : hvcoord_t
    use global_norms_mod, only: wrap_repro_sum
    type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)
    integer :: t2,n,nets,nete
    type (element_t)     , intent(inout), target :: elem(:)
    type (hvcoord_t)                  :: hvcoord
    type (TimeLevel_t), intent(inout)       :: tl
    integer, intent(in)                    :: nsubstep

    integer :: ie,k,i,j,nmax
    real (kind=r8), dimension(np,np,nlev)  :: dp   ! delta pressure
    real (kind=r8), dimension(np,np,nlev)  :: sumlk
    real (kind=r8), pointer  :: PEner(:,:,:)
    real (kind=r8), dimension(np,np)  :: suml
    real (kind=r8) :: psum(nets:nete,4),psum_g(4),beta

    ! when forcing is applied during dynamics timstep, actual forcing is
    ! slightly different (about 0.1 W/m^2) then expected by the physics
    ! since u & T are changing while FU and FT are held constant.
    ! to correct for this, save compute de_from_forcing at step 1
    ! and then adjust by:  de_from_forcing_step1 - de_from_forcing_stepN
    real (kind=r8),save :: de_from_forcing_step1
    real (kind=r8)      :: de_from_forcing

    t2=tl%np1    ! timelevel for T
    if (use_cpstar /= 0 ) then
       call endrun('Energy fixer requires use_cpstar=0')
    endif


    psum = 0
    do ie=nets,nete

#if (defined COLUMN_OPENMP)
!$omp parallel do default(shared), private(k)
#endif
       do k=1,nlev
          dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,t2)
       enddo
       suml=0
#if (defined COLUMN_OPENMP)
!$omp parallel do default(shared), private(k, i, j)
#endif
       do k=1,nlev
          do i=1,np
          do j=1,np
                sumlk(i,j,k) = cp*dp(i,j,k)
          enddo
          enddo
       enddo
       suml=0
       do k=1,nlev
          do i=1,np
          do j=1,np
             suml(i,j) = suml(i,j) + sumlk(i,j,k)
          enddo
          enddo
       enddo
       PEner => elem(ie)%accum%PEner(:,:,:)

       ! psum(:,4) = energy before forcing
       ! psum(:,1) = energy after forcing, before dynamics
       ! psum(:,2) = energy after dynamics
       ! psum(:,3) = cp*dp (internal energy added is beta*psum(:,3))
       psum(ie,3) = psum(ie,3) + SUM(suml(:,:)*elem(ie)%spheremp(:,:))
       do n=1,2
          psum(ie,n) = psum(ie,n) + SUM(  elem(ie)%spheremp(:,:)*&
               (PEner(:,:,n) + &
               elem(ie)%accum%IEner(:,:,n) + &
               elem(ie)%accum%KEner(:,:,n) ) )
       enddo
    enddo

    nmax=3

    do ie=nets,nete
       do n=1,nmax
          global_shared_buf(ie,n) = psum(ie,n)
       enddo
    enddo
    call wrap_repro_sum(nvars=nmax, comm=hybrid%par%comm)
    do n=1,nmax
       psum_g(n) = global_shared_sum(n)
    enddo

    beta = ( psum_g(1)-psum_g(2) )/psum_g(3)

    ! apply fixer
    do ie=nets,nete
       elem(ie)%state%T(:,:,:,t2) =  elem(ie)%state%T(:,:,:,t2) + beta
    enddo
    end subroutine prim_energy_fixer
!=======================================================================================================!



    subroutine smooth_topo_datasets(phis,sghdyn,sgh30dyn,elem,hybrid,nets,nete)
    use control_mod, only : smooth_phis_numcycle,smooth_sgh_numcycle
    use hybrid_mod, only : hybrid_t
    use edgetype_mod, only : EdgeBuffer_t
    use edge_mod, only : edgevpack, edgevunpack
    use derivative_mod, only : derivative_t , laplace_sphere_wk
    use prim_advance_mod, only : smooth_phis
    use prim_advection_mod, only: deriv
    implicit none

    integer , intent(in) :: nets,nete
    real (kind=r8), intent(inout)   :: phis(np,np,nets:nete)
    real (kind=r8), intent(inout)   :: sghdyn(np,np,nets:nete)
    real (kind=r8), intent(inout)   :: sgh30dyn(np,np,nets:nete)
    type (hybrid_t)      , intent(in) :: hybrid
    type (element_t)     , intent(inout), target :: elem(:)
    ! local
    integer :: ie
    real (kind=r8) :: minf

    minf=-9e9
    if (hybrid%masterthread) &
       write(iulog,*) "Applying hyperviscosity smoother to PHIS"
    call smooth_phis(phis,elem,hybrid,deriv(hybrid%ithr),nets,nete,minf,smooth_phis_numcycle)

    minf=0
    if (hybrid%masterthread) &
       write(iulog,*) "Applying hyperviscosity smoother to SGH"
    call smooth_phis(sghdyn,elem,hybrid,deriv(hybrid%ithr),nets,nete,minf,smooth_sgh_numcycle)
    if (hybrid%masterthread) &
       write(iulog,*) "Applying hyperviscosity smoother to SGH30"
    call smooth_phis(sgh30dyn,elem,hybrid,deriv(hybrid%ithr),nets,nete,minf,smooth_sgh_numcycle)

    end subroutine smooth_topo_datasets

end module prim_driver_mod
