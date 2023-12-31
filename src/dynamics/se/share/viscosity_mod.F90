module viscosity_mod
!
!  This module should be renamed "global_deriv_mod.F90"
! 
!  It is a collection of derivative operators that must be applied to the field 
!  over the sphere (as opposed to derivative operators that can be applied element 
!  by element)
!
!
  use shr_kind_mod,   only: r8=>shr_kind_r8
  use thread_mod
  use dimensions_mod, only: np, nlev,qsize,nelemd
  use hybrid_mod,     only: hybrid_t, config_thread_region
  use parallel_mod,   only: parallel_t
  use element_mod,    only: element_t
  use derivative_mod, only: derivative_t, laplace_sphere_wk, vlaplace_sphere_wk, vorticity_sphere, derivinit, divergence_sphere
  use edgetype_mod,   only: EdgeBuffer_t
  use edge_mod,       only: edgevpack, edgevunpack, edgevunpackmin, &
       edgevunpackmax, initEdgeBuffer, FreeEdgeBuffer, edgeSunpackmax, edgeSunpackmin, edgeSpack
  use bndry_mod,      only: bndry_exchange, bndry_exchange_start, bndry_exchange_finish
  use control_mod,    only: hypervis_scaling, nu, nu_div


  implicit none
  save

  public :: biharmonic_wk_scalar
#if 0
  public :: biharmonic_wk_omega
#endif
  public :: neighbor_minmax, neighbor_minmax_start,neighbor_minmax_finish
  public :: biharmonic_wk

  !
  ! compute vorticity/divergence and then project to make continious
  ! high-level routines uses only for I/O
  public :: compute_zeta_C0
  public :: compute_div_C0

  interface compute_zeta_C0
    module procedure compute_zeta_C0_hybrid       ! hybrid version
    module procedure compute_zeta_C0_par          ! single threaded
  end interface compute_zeta_C0
  interface compute_div_C0
    module procedure compute_div_C0_hybrid
    module procedure compute_div_C0_par
  end interface compute_div_C0

  public :: compute_zeta_C0_contra    ! for older versions of sweq which carry
  public :: compute_div_C0_contra     ! velocity around in contra-coordinates

  type (EdgeBuffer_t)          :: edge1

CONTAINS

subroutine biharmonic_wk_dp3d(elem,dptens,ptens,vtens,deriv,edge3,hybrid,nt,nets,nete)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! compute weak biharmonic operator
  !    input:  h,v (stored in elem()%, in lat-lon coordinates
  !    output: ptens,vtens  overwritten with weak biharmonic of h,v (output in lat-lon coordinates)
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  integer              , intent(in)  :: nt,nets,nete
  real (kind=r8), dimension(np,np,2,nlev,nets:nete)  :: vtens
  real (kind=r8), dimension(np,np,nlev,nets:nete) :: ptens,dptens
  type (EdgeBuffer_t)  , intent(inout) :: edge3
  type (derivative_t)  , intent(in) :: deriv
  
  ! local
  integer :: k,kptr,ie
  real (kind=r8), dimension(:,:), pointer :: rspheremv
  real (kind=r8), dimension(np,np) :: tmp
  real (kind=r8), dimension(np,np,2) :: v
  real (kind=r8) :: nu_ratio1, nu_ratio2
  logical var_coef1

   !if tensor hyperviscosity with tensor V is used, then biharmonic operator is (\grad\cdot V\grad) (\grad \cdot \grad) 
   !so tensor is only used on second call to laplace_sphere_wk
   var_coef1 = .true.
   if(hypervis_scaling > 0)    var_coef1 = .false.

   ! note: there is a scaling bug in the treatment of nu_div
   ! nu_ratio is applied twice, once in each laplace operator
   ! so in reality:   nu_div_actual = (nu_div/nu)**2 nu
   ! We should fix this, but it requires adjusting all CAM defaults
   nu_ratio1=1
   nu_ratio2=1
   if (nu_div/=nu) then
      if(hypervis_scaling /= 0) then
         ! we have a problem with the tensor in that we cant seperate
         ! div and curl components.  So we do, with tensor V:
         ! nu * (del V del ) * ( nu_ratio * grad(div) - curl(curl))
         nu_ratio1=(nu_div/nu)**2   ! preserve buggy scaling
         nu_ratio2=1
      else
         nu_ratio1=nu_div/nu
         nu_ratio2=nu_div/nu
      endif
   endif


   do ie=nets,nete

#if (defined COLUMN_OPENMP)
!$omp parallel do default(shared), private(k,tmp)
#endif
      do k=1,nlev
         tmp=elem(ie)%state%T(:,:,k,nt) 
         !ptens(:,:,k,ie)=laplace_sphere_wk(tmp,deriv,elem(ie),var_coef=var_coef1)
         call laplace_sphere_wk(tmp,deriv,elem(ie),ptens(:,:,k,ie),var_coef=var_coef1)
         tmp=elem(ie)%state%dp3d(:,:,k,nt) 
         !dptens(:,:,k,ie)=laplace_sphere_wk(tmp,deriv,elem(ie),var_coef=var_coef1)
         call laplace_sphere_wk(tmp,deriv,elem(ie),dptens(:,:,k,ie),var_coef=var_coef1)
         !vtens(:,:,:,k,ie)=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),&
         !     var_coef=var_coef1,nu_ratio=nu_ratio1)
         call vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),&
              vtens(:,:,:,k,ie), var_coef=var_coef1,nu_ratio=nu_ratio1)
      enddo
      kptr=0
      call edgeVpack(edge3, ptens(1,1,1,ie),nlev,kptr,ie)
      kptr=nlev
      call edgeVpack(edge3, vtens(1,1,1,1,ie),2*nlev,kptr,ie)
      kptr=3*nlev
      call edgeVpack(edge3, dptens(1,1,1,ie),nlev,kptr,ie)

   enddo
   
   call bndry_exchange(hybrid,edge3)
   
   do ie=nets,nete
      rspheremv     => elem(ie)%rspheremp(:,:)
      
      kptr=0
      call edgeVunpack(edge3, ptens(1,1,1,ie), nlev, kptr, ie)
      kptr=nlev
      call edgeVunpack(edge3, vtens(1,1,1,1,ie), 2*nlev, kptr, ie)
      kptr=3*nlev
      call edgeVunpack(edge3, dptens(1,1,1,ie), nlev, kptr, ie)
      
      ! apply inverse mass matrix, then apply laplace again
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,v,tmp)
#endif
      do k=1,nlev
         tmp(:,:)=rspheremv(:,:)*ptens(:,:,k,ie)
         !ptens(:,:,k,ie)=laplace_sphere_wk(tmp,deriv,elem(ie),var_coef=.true.)
         call laplace_sphere_wk(tmp,deriv,elem(ie),ptens(:,:,k,ie),var_coef=.true.)
         tmp(:,:)=rspheremv(:,:)*dptens(:,:,k,ie)
         !dptens(:,:,k,ie)=laplace_sphere_wk(tmp,deriv,elem(ie),var_coef=.true.)
         call laplace_sphere_wk(tmp,deriv,elem(ie),dptens(:,:,k,ie),var_coef=.true.)

         v(:,:,1)=rspheremv(:,:)*vtens(:,:,1,k,ie)
         v(:,:,2)=rspheremv(:,:)*vtens(:,:,2,k,ie)
         !vtens(:,:,:,k,ie)=vlaplace_sphere_wk(v(:,:,:),deriv,elem(ie),&
         !     var_coef=.true.,nu_ratio=nu_ratio2)
         call vlaplace_sphere_wk(v(:,:,:),deriv,elem(ie),&
              vtens(:,:,:,k,ie),var_coef=.true.,nu_ratio=nu_ratio2)

      enddo
   enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine biharmonic_wk_dp3d

subroutine biharmonic_wk(elem,pstens,ptens,vtens,deriv,edge3,hybrid,nt,nets,nete)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute weak biharmonic operator
!    input:  h,v (stored in elem()%, in lat-lon coordinates
!    output: ptens,vtens  overwritten with weak biharmonic of h,v (output in lat-lon coordinates)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(inout), target :: elem(:)
integer :: nt,nets,nete
real (kind=r8), dimension(np,np,2,nlev,nets:nete)  :: vtens
real (kind=r8), dimension(np,np,nlev,nets:nete) :: ptens
type (EdgeBuffer_t)  , intent(inout) :: edge3
type (derivative_t)  , intent(in) :: deriv
real (kind=r8), dimension(np,np,nets:nete) :: pstens

! local
integer :: k,kptr,i,j,ie,ic
real (kind=r8), dimension(:,:), pointer :: rspheremv
real (kind=r8), dimension(np,np) :: lap_ps
real (kind=r8), dimension(np,np,nlev) :: T
real (kind=r8), dimension(np,np,2) :: v
real (kind=r8) ::  nu_ratio1,nu_ratio2
logical var_coef1

   !if tensor hyperviscosity with tensor V is used, then biharmonic operator is (\grad\cdot V\grad) (\grad \cdot \grad) 
   !so tensor is only used on second call to laplace_sphere_wk
   var_coef1 = .true.
   if(hypervis_scaling > 0)  var_coef1= .false.

   ! note: there is a scaling bug in the treatment of nu_div
   ! nu_ratio is applied twice, once in each laplace operator
   ! so in reality:   nu_div_actual = (nu_div/nu)**2 nu
   ! We should fix this, but it requires adjusting all CAM defaults
   nu_ratio1=1
   nu_ratio2=1
   if (nu_div/=nu) then
      if(hypervis_scaling /= 0) then
         ! we have a problem with the tensor in that we cant seperate
         ! div and curl components.  So we do, with tensor V:
         ! nu * (del V del ) * ( nu_ratio * grad(div) - curl(curl))
         nu_ratio1=(nu_div/nu)**2   ! preserve buggy scaling
         nu_ratio2=1
      else
         nu_ratio1=nu_div/nu
         nu_ratio2=nu_div/nu
      endif
   endif


   do ie=nets,nete
      
      ! should filter lnps + PHI_s/RT?
      call laplace_sphere_wk(elem(ie)%state%ps_v(:,:,nt),deriv,elem(ie),pstens(:,:,ie),var_coef=var_coef1)
      
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k, j, i)
#endif
      do k=1,nlev
         do j=1,np
            do i=1,np
               T(i,j,k)=elem(ie)%state%T(i,j,k,nt) 
            enddo
         enddo
        
         !ptens(:,:,k,ie)=laplace_sphere_wk(T(:,:,k),deriv,elem(ie),var_coef=var_coef1)
         call laplace_sphere_wk(T(:,:,k),deriv,elem(ie),ptens(:,:,k,ie),var_coef=var_coef1)
         !vtens(:,:,:,k,ie)=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,&
         !     elem(ie),var_coef=var_coef1,nu_ratio=nu_ratio1)
         call vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,&
              elem(ie),vtens(:,:,:,k,ie),var_coef=var_coef1,nu_ratio=nu_ratio1)

      enddo
      kptr=0
      call edgeVpack(edge3, ptens(:,:,:,ie),nlev,kptr,ie)
      kptr=nlev
      call edgeVpack(edge3, vtens(:,:,1,:,ie),nlev,kptr,ie)
      kptr=2*nlev
      call edgeVpack(edge3, vtens(:,:,2,:,ie),nlev,kptr,ie)

      kptr=3*nlev
      call edgeVpack(edge3, pstens(:,:,ie),1,kptr,ie)
   enddo
   
   call bndry_exchange(hybrid,edge3)
   
   do ie=nets,nete
      rspheremv     => elem(ie)%rspheremp(:,:)
      
      kptr=0
      call edgeVunpack(edge3, ptens(:,:,:,ie), nlev, kptr, ie)
      kptr=nlev
      call edgeVunpack(edge3, vtens(:,:,1,:,ie), nlev, kptr, ie)
      kptr=2*nlev
      call edgeVunpack(edge3, vtens(:,:,1,:,ie), nlev, kptr, ie)
      
      ! apply inverse mass matrix, then apply laplace again
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k, j, i, v)
#endif
      do k=1,nlev
         do j=1,np
            do i=1,np
               T(i,j,k)=rspheremv(i,j)*ptens(i,j,k,ie)
               v(i,j,1)=rspheremv(i,j)*vtens(i,j,1,k,ie)
               v(i,j,2)=rspheremv(i,j)*vtens(i,j,2,k,ie)
            enddo
         enddo
         !ptens(:,:,k,ie)=laplace_sphere_wk(T(:,:,k),deriv,elem(ie),var_coef=.true.)
         call laplace_sphere_wk(T(:,:,k),deriv,elem(ie),ptens(:,:,k,ie),var_coef=.true.)
         !vtens(:,:,:,k,ie)=vlaplace_sphere_wk(v(:,:,:),deriv,elem(ie),var_coef=.true.,&
         !     nu_ratio=nu_ratio2)
         call vlaplace_sphere_wk(v(:,:,:),deriv,elem(ie),vtens(:,:,:,k,ie),var_coef=.true.,&
              nu_ratio=nu_ratio2)
      enddo
         
      kptr=3*nlev
      call edgeVunpack(edge3, pstens(:,:,ie), 1, kptr,ie) 
      ! apply inverse mass matrix, then apply laplace again
      lap_ps(:,:)=rspheremv(:,:)*pstens(:,:,ie)
      !pstens(:,:,ie)=laplace_sphere_wk(lap_ps,deriv,elem(ie),var_coef=.true.)
      call laplace_sphere_wk(lap_ps,deriv,elem(ie),pstens(:,:,ie),var_coef=.true.)

   enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine biharmonic_wk

#if 0
subroutine biharmonic_wk_omega(elem,ptens,deriv,edge3,hybrid,nets,nete,kbeg,kend)
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  integer              , intent(in)  :: nets,nete
  integer              , intent(in)  :: kbeg, kend
  real (kind=r8), dimension(np,np,nlev,nets:nete) :: ptens
  type (EdgeBuffer_t)  , intent(inout) :: edge3
  type (derivative_t)  , intent(in) :: deriv
  
  ! local
  integer :: i,j,k,kptr,ie,kblk
  real (kind=r8), dimension(:,:), pointer :: rspheremv
  real (kind=r8), dimension(np,np) :: tmp
  real (kind=r8), dimension(np,np) :: tmp2
  real (kind=r8), dimension(np,np,2) :: v
  real (kind=r8) :: nu_ratio1, nu_ratio2
  logical var_coef1
  
  kblk = kend - kbeg + 1
  
  !if tensor hyperviscosity with tensor V is used, then biharmonic operator is (\grad\cdot V\grad) (\grad \cdot \grad) 
  !so tensor is only used on second call to laplace_sphere_wk
  var_coef1 = .true.
  if(hypervis_scaling > 0)    var_coef1 = .false.
  
  nu_ratio1=1
  nu_ratio2=1
  
  do ie=nets,nete
    
    !$omp parallel do num_threads(vert_num_threads) private(k,tmp)
    do k=kbeg,kend
      tmp=elem(ie)%derived%omega(:,:,k) 
      call laplace_sphere_wk(tmp,deriv,elem(ie),ptens(:,:,k,ie),var_coef=var_coef1)
    enddo
    
    kptr = kbeg - 1
    call edgeVpack(edge3,ptens(:,:,kbeg:kend,ie),kblk,kptr,ie)
  enddo
  
  call bndry_exchange(hybrid,edge3,location='biharmonic_wk_omega')
  
  do ie=nets,nete
    rspheremv     => elem(ie)%rspheremp(:,:)
    
    kptr = kbeg - 1
    call edgeVunpack(edge3,ptens(:,:,kbeg:kend,ie),kblk,kptr,ie)
    
    ! apply inverse mass matrix, then apply laplace again
    !$omp parallel do num_threads(vert_num_threads) private(k,tmp)
    do k=kbeg,kend
      tmp(:,:)=rspheremv(:,:)*ptens(:,:,k,ie)
      call laplace_sphere_wk(tmp,deriv,elem(ie),ptens(:,:,k,ie),var_coef=.true.)
    enddo
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine biharmonic_wk_omega
#endif

subroutine biharmonic_wk_scalar(elem,qtens,deriv,edgeq,hybrid,nets,nete)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute weak biharmonic operator
!    input:  qtens = Q
!    output: qtens = weak biharmonic of Q
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(inout), target :: elem(:)
integer :: nets,nete
real (kind=r8), dimension(np,np,nlev,qsize,nets:nete) :: qtens
type (EdgeBuffer_t)  , intent(inout) :: edgeq
type (derivative_t)  , intent(in) :: deriv

! local
integer :: k,kptr,i,j,ie,ic,q
real (kind=r8), dimension(np,np) :: lap_p
logical var_coef1

   !if tensor hyperviscosity with tensor V is used, then biharmonic operator is (\grad\cdot V\grad) (\grad \cdot \grad) 
   !so tensor is only used on second call to laplace_sphere_wk
   var_coef1 = .true.
   if(hypervis_scaling > 0)    var_coef1 = .false.



   do ie=nets,nete
      do q=1,qsize      
         do k=1,nlev    !  Potential loop inversion (AAM)
           lap_p(:,:)=qtens(:,:,k,q,ie)
           call laplace_sphere_wk(lap_p,deriv,elem(ie),qtens(:,:,k,q,ie),var_coef=var_coef1)
         enddo
         kptr = nlev*(q-1) 
         call edgeVpack(edgeq, qtens(:,:,1:nlev,q,ie),nlev,kptr,ie)
      enddo
   enddo


   call bndry_exchange(hybrid,edgeq)
   
   do ie=nets,nete

      ! apply inverse mass matrix, then apply laplace again
      do q=1,qsize      
        kptr = nlev*(q-1) 
        call edgeVunpack(edgeq, qtens(:,:,1:nlev,q,ie),nlev,kptr,ie)
        do k=1,nlev    !  Potential loop inversion (AAM)
           lap_p(:,:)=elem(ie)%rspheremp(:,:)*qtens(:,:,k,q,ie)
           call laplace_sphere_wk(lap_p,deriv,elem(ie),qtens(:,:,k,q,ie),var_coef=.true.)
        enddo
      enddo
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine biharmonic_wk_scalar


subroutine make_C0(zeta,elem,hybrid,nets,nete)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! apply DSS (aka assembly procedure) to zeta.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
integer :: nets,nete
real (kind=r8), dimension(np,np,nlev,nets:nete) :: zeta

! local
integer :: k,i,j,ie,ic,kptr


    call initEdgeBuffer(hybrid%par,edge1,elem,nlev)

do ie=nets,nete
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
   do k=1,nlev
      zeta(:,:,k,ie)=zeta(:,:,k,ie)*elem(ie)%spheremp(:,:)
   enddo
   kptr=0
   call edgeVpack(edge1, zeta(1,1,1,ie),nlev,kptr,ie)
enddo
call bndry_exchange(hybrid,edge1)
do ie=nets,nete
   kptr=0
   call edgeVunpack(edge1, zeta(1,1,1,ie),nlev,kptr, ie)
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
   do k=1,nlev
      zeta(:,:,k,ie)=zeta(:,:,k,ie)*elem(ie)%rspheremp(:,:)
   enddo
enddo

call FreeEdgeBuffer(edge1) 

end subroutine


subroutine make_C0_vector(v,elem,hybrid,nets,nete)
#if 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! apply DSS to a velocity vector
! this is a low-performance routine used for I/O and analysis.
! no need to optimize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
integer :: nets,nete
real (kind=r8), dimension(np,np,2,nlev,nets:nete) :: v

! local
integer :: k,i,j,ie,ic,kptr
type (EdgeBuffer_t)          :: edge2


    call initEdgeBuffer(hybrid%par,edge2,elem,2*nlev)

do ie=nets,nete
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
   do k=1,nlev
      v(:,:,1,k,ie)=v(:,:,1,k,ie)*elem(ie)%spheremp(:,:)
      v(:,:,2,k,ie)=v(:,:,2,k,ie)*elem(ie)%spheremp(:,:)
   enddo
   kptr=0
   call edgeVpack(edge2, v(1,1,1,1,ie),2*nlev,kptr,ie)
enddo
call bndry_exchange(hybrid,edge2)
do ie=nets,nete
   kptr=0
   call edgeVunpack(edge2, v(1,1,1,1,ie),2*nlev,kptr,ie)
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
   do k=1,nlev
      v(:,:,1,k,ie)=v(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)
      v(:,:,2,k,ie)=v(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)
   enddo
enddo

call FreeEdgeBuffer(edge2) 
#endif
end subroutine






subroutine compute_zeta_C0_contra(zeta,elem,hybrid,nets,nete,nt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 vorticity.  That is, solve:  
!     < PHI, zeta > = <PHI, curl(elem%state%v >
!
!    input:  v (stored in elem()%, in contra-variant coordinates)
!    output: zeta(:,:,:,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
integer :: nt,nets,nete
real (kind=r8), dimension(np,np,nlev,nets:nete) :: zeta
real (kind=r8), dimension(np,np,2) :: ulatlon
real (kind=r8), dimension(np,np) :: v1,v2

! local
integer :: k,ie
type (derivative_t)          :: deriv

call derivinit(deriv)

do k=1,nlev
do ie=nets,nete
    v1 = elem(ie)%state%v(:,:,1,k,nt)
    v2 = elem(ie)%state%v(:,:,2,k,nt)
    ulatlon(:,:,1) = elem(ie)%D(:,:,1,1)*v1 + elem(ie)%D(:,:,1,2)*v2
    ulatlon(:,:,2) = elem(ie)%D(:,:,2,1)*v1 + elem(ie)%D(:,:,2,2)*v2
   call vorticity_sphere(ulatlon,deriv,elem(ie),zeta(:,:,k,ie))
enddo
enddo

call make_C0(zeta,elem,hybrid,nets,nete)

end subroutine



subroutine compute_div_C0_contra(zeta,elem,hybrid,nets,nete,nt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 divergence. That is, solve:  
!     < PHI, zeta > = <PHI, div(elem%state%v >
!
!    input:  v (stored in elem()%, in contra-variant coordinates)
!    output: zeta(:,:,:,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
integer :: nt,nets,nete
real (kind=r8), dimension(np,np,nlev,nets:nete) :: zeta
real (kind=r8), dimension(np,np,2) :: ulatlon
real (kind=r8), dimension(np,np) :: v1,v2

! local
integer :: k,ie
type (derivative_t)          :: deriv

call derivinit(deriv)

do k=1,nlev
do ie=nets,nete
    v1 = elem(ie)%state%v(:,:,1,k,nt)
    v2 = elem(ie)%state%v(:,:,2,k,nt)
    ulatlon(:,:,1) = elem(ie)%D(:,:,1,1)*v1 + elem(ie)%D(:,:,1,2)*v2
    ulatlon(:,:,2) = elem(ie)%D(:,:,2,1)*v1 + elem(ie)%D(:,:,2,2)*v2
   call divergence_sphere(ulatlon,deriv,elem(ie),zeta(:,:,k,ie))
enddo
enddo

call make_C0(zeta,elem,hybrid,nets,nete)

end subroutine

subroutine compute_zeta_C0_par(zeta,elem,par,nt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 vorticity.  That is, solve:  
!     < PHI, zeta > = <PHI, curl(elem%state%v >
!
!    input:  v (stored in elem()%, in lat-lon coordinates)
!    output: zeta(:,:,:,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type (parallel_t) :: par
type (element_t)     , intent(in), target :: elem(:)
real (kind=r8), dimension(np,np,nlev,nelemd) :: zeta
integer :: nt

! local
type (hybrid_t)              :: hybrid
integer :: k,i,j,ie,ic
type (derivative_t)          :: deriv

! single thread
hybrid = config_thread_region(par,'serial')

call compute_zeta_C0_hybrid(zeta,elem,hybrid,1,nelemd,nt)

end subroutine


subroutine compute_div_C0_par(zeta,elem,par,nt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 divergence. That is, solve:  
!     < PHI, zeta > = <PHI, div(elem%state%v >
!
!    input:  v (stored in elem()%, in lat-lon coordinates)
!    output: zeta(:,:,:,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (parallel_t) :: par
type (element_t)     , intent(in), target :: elem(:)
real (kind=r8), dimension(np,np,nlev,nelemd) :: zeta
integer :: nt

! local
type (hybrid_t)              :: hybrid
integer :: k,i,j,ie,ic
type (derivative_t)          :: deriv

! single thread
hybrid = config_thread_region(par,'serial')

call compute_div_C0_hybrid(zeta,elem,hybrid,1,nelemd,nt)

end subroutine



subroutine compute_zeta_C0_hybrid(zeta,elem,hybrid,nets,nete,nt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 vorticity.  That is, solve:  
!     < PHI, zeta > = <PHI, curl(elem%state%v >
!
!    input:  v (stored in elem()%, in lat-lon coordinates)
!    output: zeta(:,:,:,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
integer :: nt,nets,nete
real (kind=r8), dimension(np,np,nlev,nets:nete) :: zeta

! local
integer :: k,i,j,ie,ic
type (derivative_t)          :: deriv

call derivinit(deriv)

do ie=nets,nete
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
do k=1,nlev
   call vorticity_sphere(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),zeta(:,:,k,ie))
enddo
enddo

call make_C0(zeta,elem,hybrid,nets,nete)

end subroutine


subroutine compute_div_C0_hybrid(zeta,elem,hybrid,nets,nete,nt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 divergence. That is, solve:  
!     < PHI, zeta > = <PHI, div(elem%state%v >
!
!    input:  v (stored in elem()%, in lat-lon coordinates)
!    output: zeta(:,:,:,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
integer :: nt,nets,nete
real (kind=r8), dimension(np,np,nlev,nets:nete) :: zeta

! local
integer :: k,i,j,ie,ic
type (derivative_t)          :: deriv

call derivinit(deriv)

do ie=nets,nete
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
do k=1,nlev
   call divergence_sphere(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),zeta(:,:,k,ie))
enddo
enddo

call make_C0(zeta,elem,hybrid,nets,nete)

end subroutine








subroutine neighbor_minmax(hybrid,edgeMinMax,nets,nete,min_neigh,max_neigh)
 
   type (hybrid_t)      , intent(in) :: hybrid
   type (EdgeBuffer_t)  , intent(inout) :: edgeMinMax
   integer :: nets,nete
   real (kind=r8) :: min_neigh(nlev,qsize,nets:nete)
   real (kind=r8) :: max_neigh(nlev,qsize,nets:nete)

   ! local 
   integer:: ie, q, k
   integer:: kblk,qblk,kptr

   kblk = nlev  ! calculate size of the block of vertical levels
   qblk = qsize   ! calculate size of the block of tracers

   do ie=nets,nete
      do q = 1, qsize
         kptr = nlev*(q - 1)
         call  edgeSpack(edgeMinMax,min_neigh(:,q,ie),kblk,kptr,ie)
         kptr = qsize*nlev + nlev*(q - 1)
         call  edgeSpack(edgeMinMax,max_neigh(:,q,ie),kblk,kptr,ie)
      enddo
   enddo
   
   call bndry_exchange(hybrid,edgeMinMax,location='neighbor_minmax')

   do ie=nets,nete
      do q=1,qsize
         kptr = nlev*(q - 1)
         call  edgeSunpackMIN(edgeMinMax,min_neigh(:,q,ie),kblk,kptr,ie)
         kptr = qsize*nlev + nlev*(q - 1) 
         call  edgeSunpackMAX(edgeMinMax,max_neigh(:,q,ie),kblk,kptr,ie)
         do k=1,nlev
            min_neigh(k,q,ie) = max(min_neigh(k,q,ie),0d0)
         enddo
      enddo
   enddo

end subroutine neighbor_minmax
  

subroutine neighbor_minmax_start(hybrid,edgeMinMax,nets,nete,min_neigh,max_neigh)

   type (hybrid_t)      , intent(in) :: hybrid
   type (EdgeBuffer_t)  , intent(inout) :: edgeMinMax
   integer :: nets,nete
   real (kind=r8) :: min_neigh(nlev,qsize,nets:nete)
   real (kind=r8) :: max_neigh(nlev,qsize,nets:nete)
   integer :: kblk, qblk
   integer :: kbeg, kend, qbeg, qend

   ! local 
   integer :: ie,q, k,kptr

   !call get_loop_ranges(hybrid,kbeg=kbeg,kend=kend,qbeg=qbeg,qend=qend)
   kbeg = 1
   kend = nlev
   qbeg = 1
   qend = qsize

   kblk = kend - kbeg + 1   ! calculate size of the block of vertical levels
   qblk = qend - qbeg + 1   ! calculate size of the block of tracers

   do ie=nets,nete
      do q=qbeg, qend
         kptr = nlev*(q - 1) + kbeg - 1
         call  edgeSpack(edgeMinMax,min_neigh(kbeg:kend,q,ie),kblk,kptr,ie)
         kptr = qsize*nlev + nlev*(q - 1) + kbeg - 1
         call  edgeSpack(edgeMinMax,max_neigh(kbeg:kend,q,ie),kblk,kptr,ie)
      enddo
   enddo

   call bndry_exchange_start(hybrid,edgeMinMax,location='viscosity_mod:882')

end subroutine neighbor_minmax_start

subroutine neighbor_minmax_finish(hybrid,edgeMinMax,nets,nete,min_neigh,max_neigh)

   type (hybrid_t)      , intent(in) :: hybrid
   type (EdgeBuffer_t)  , intent(inout) :: edgeMinMax
   integer :: nets,nete
   real (kind=r8) :: min_neigh(nlev,qsize,nets:nete)
   real (kind=r8) :: max_neigh(nlev,qsize,nets:nete)
   integer :: kblk, qblk
   integer :: ie,q, k,kptr
   integer :: kbeg, kend, qbeg, qend

   !call get_loop_ranges(hybrid,kbeg=kbeg,kend=kend,qbeg=qbeg,qend=qend)
   qbeg=1
   qend=qsize
   kbeg=1
   kend=nlev

   kblk = kend - kbeg + 1   ! calculate size of the block of vertical levels
   qblk = qend - qbeg + 1   ! calculate size of the block of tracers

   call bndry_exchange_finish(hybrid,edgeMinMax,location='viscosity_mod:902')

   do ie=nets,nete
      do q=qbeg, qend
         kptr = nlev*(q - 1) + kbeg - 1
         call  edgeSunpackMIN(edgeMinMax,min_neigh(kbeg:kend,q,ie),kblk,kptr,ie)
         kptr = qsize*nlev + nlev*(q - 1) + kbeg - 1
         call  edgeSunpackMAX(edgeMinMax,max_neigh(kbeg:kend,q,ie),kblk,kptr,ie)
         do k=kbeg,kend
            min_neigh(k,q,ie) = max(min_neigh(k,q,ie),0.0_r8)
         enddo
      enddo
   enddo

end subroutine neighbor_minmax_finish

end module viscosity_mod
