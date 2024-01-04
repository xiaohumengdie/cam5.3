module vertremap_mod

  use shr_kind_mod,   only: r8=>shr_kind_r8
  use vertremap_base, only: remap1
  use dimensions_mod, only: np,nlev,qsize
  use element_mod,    only: element_t
  use perf_mod,       only: t_startf, t_stopf  ! _EXTERNAL

  implicit none
  private
  public :: vertical_remap

contains

  subroutine vertical_remap(hybrid,elem,hvcoord,dt,np1,np1_qdp,nets,nete)

  ! This routine is called at the end of the vertically Lagrangian
  ! dynamics step to compute the vertical flux needed to get back
  ! to reference eta levels
  !
  ! input:
  !     derived%dp()  delta p on levels at beginning of timestep
  !     state%dp3d(np1)  delta p on levels at end of timestep
  ! output:
  !     state%ps_v(np1)          surface pressure at time np1
  !     derived%eta_dot_dpdn()   vertical flux from final Lagrangian
  !                              levels to reference eta levels
  !
  use hybvcoord_mod,  only: hvcoord_t
  use control_mod,    only: rsplit
  use hybrid_mod,     only: hybrid_t
  use cam_abortutils, only: endrun


  type (hybrid_t),  intent(in)    :: hybrid  ! distributed parallel structure (shared)
  type (element_t), intent(inout) :: elem(:)
  type (hvcoord_t)                :: hvcoord
  real (kind=r8)           :: dt

  integer :: ie,i,j,k,np1,nets,nete,np1_qdp
  real (kind=r8), dimension(np,np,nlev)  :: dp,dp_star
  real (kind=r8), dimension(np,np,nlev,2)  :: ttmp

  call t_startf('vertical_remap')

  ! reference levels:
  !   dp(k) = (hyai(k+1)-hyai(k))*ps0 + (hybi(k+1)-hybi(k))*ps_v(i,j)
  !   hybi(1)=0          pure pressure at top of atmosphere
  !   hyai(1)=ptop
  !   hyai(nlev+1) = 0   pure sigma at bottom
  !   hybi(nlev+1) = 1
  !
  ! sum over k=1,nlev
  !  sum(dp(k)) = (hyai(nlev+1)-hyai(1))*ps0 + (hybi(nlev+1)-hybi(1))*ps_v
  !             = -ps0 + ps_v
  !  ps_v =  ps0+sum(dp(k))
  !
  ! reference levels:
  !    dp(k) = (hyai(k+1)-hyai(k))*ps0 + (hybi(k+1)-hybi(k))*ps_v
  ! floating levels:
  !    dp_star(k) = dp(k) + dt_q*(eta_dot_dpdn(i,j,k+1) - eta_dot_dpdn(i,j,k) )
  ! hence:
  !    (dp_star(k)-dp(k))/dt_q = (eta_dot_dpdn(i,j,k+1) - eta_dot_dpdn(i,j,k) )
  !
  
  do ie=nets,nete
!     ! SET VERTICAL VELOCITY TO ZERO FOR DEBUGGING
!     elem(ie)%derived%eta_dot_dpdn(:,:,:)=0
     if (rsplit==0) then
        ! compute dp_star from eta_dot_dpdn():
        do k=1,nlev
           dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,np1)
           dp_star(:,:,k) = dp(:,:,k) + dt*(elem(ie)%derived%eta_dot_dpdn(:,:,k+1) -&
                elem(ie)%derived%eta_dot_dpdn(:,:,k))
        enddo
        if (minval(dp_star)<0) call endrun('negative layer thickness.  timestep or remap time too large')
     else
        !  REMAP u,v,T from levels in dp3d() to REF levels
        !
        ! update final ps_v
        elem(ie)%state%ps_v(:,:,np1) = hvcoord%hyai(1)*hvcoord%ps0 + &
             sum(elem(ie)%state%dp3d(:,:,:,np1),3)
        do k=1,nlev
           dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,np1)
           dp_star(:,:,k) = elem(ie)%state%dp3d(:,:,k,np1)
        enddo
        if (minval(dp_star)<0) call endrun('negative layer thickness.  timestep or remap time too large')

        ! remap the dynamics:
#undef REMAP_TE
#ifdef REMAP_TE
        ! remap u,v and cp*T + .5 u^2
        ttmp(:,:,:,1)=(elem(ie)%state%v(:,:,1,:,np1)**2 + &
             elem(ie)%state%v(:,:,2,:,np1)**2)/2 + &
             elem(ie)%state%t(:,:,:,np1)*cp
#else
        ttmp(:,:,:,1)=elem(ie)%state%t(:,:,:,np1)
#endif
        ttmp(:,:,:,1)=ttmp(:,:,:,1)*dp_star
        call remap1(ttmp,np,1,dp_star,dp)
        elem(ie)%state%t(:,:,:,np1)=ttmp(:,:,:,1)/dp

        ttmp(:,:,:,1)=elem(ie)%state%v(:,:,1,:,np1)*dp_star
        ttmp(:,:,:,2)=elem(ie)%state%v(:,:,2,:,np1)*dp_star
        call remap1(ttmp,np,2,dp_star,dp)
        elem(ie)%state%v(:,:,1,:,np1)=ttmp(:,:,:,1)/dp
        elem(ie)%state%v(:,:,2,:,np1)=ttmp(:,:,:,2)/dp
#ifdef REMAP_TE
        ! back out T from TE
        elem(ie)%state%t(:,:,:,np1) = &
             ( elem(ie)%state%t(:,:,:,np1) - ( (elem(ie)%state%v(:,:,1,:,np1)**2 + &
             elem(ie)%state%v(:,:,2,:,np1)**2)/2))/cp

#endif
     endif

     ! remap the tracers from lagrangian levels (dp_star)  to REF levels dp
     if (qsize>0) then
        call remap1(elem(ie)%state%Qdp(:,:,:,:,np1_qdp),np,qsize,dp_star,dp)
     endif
  enddo
  call t_stopf('vertical_remap')
  end subroutine vertical_remap

end module 
