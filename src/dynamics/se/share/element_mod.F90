module element_mod

  use shr_kind_mod,           only: r8=>shr_kind_r8, i8=>shr_kind_i8
  use coordinate_systems_mod, only: spherical_polar_t, cartesian2D_t, cartesian3D_t, distance
  use dimensions_mod,         only: np, npsq, nlev, nlevp, qsize_d, max_neigh_edges
  use edgetype_mod,           only: edgedescriptor_t, rotation_t
  use gridgraph_mod,          only: gridvertex_t

  implicit none
  private
  integer, public, parameter :: timelevels = 3


! =========== PRIMITIVE-EQUATION DATA-STRUCTURES =====================

  type, public :: elem_state_t

    ! prognostic variables for preqx solver

    ! prognostics must match those in prim_restart_mod.F90
    ! vertically-lagrangian code advects dp3d instead of ps_v
    ! tracers Q, Qdp always use 2 level time scheme

    real (kind=r8) :: v     (np,np,2,nlev,timelevels)            ! velocity                           1
    real (kind=r8) :: T     (np,np,nlev,timelevels)              ! temperature                        2
    real (kind=r8) :: dp3d  (np,np,nlev,timelevels)              ! delta p on levels                  8
    real (kind=r8) :: lnps  (np,np,timelevels)                   ! log surface pressure               3
    real (kind=r8) :: ps_v  (np,np,timelevels)                   ! surface pressure                   4
    real (kind=r8) :: phis  (np,np)                              ! surface geopotential (prescribed)  5
    real (kind=r8) :: Q     (np,np,nlev,qsize_d)                 ! Tracer concentration               6
    real (kind=r8) :: Qdp   (np,np,nlev,qsize_d,2)               ! Tracer mass                        7
  end type elem_state_t

  integer, public, parameter :: StateComponents=8! num prognistics variables (for prim_restart_mod.F90)

  !___________________________________________________________________
  type, public :: derived_state_t
     !
     ! storage for subcycling tracers/dynamics
     !
    real (kind=r8) :: vn0  (np,np,2,nlev)                      ! velocity for SE tracer advection
    real (kind=r8) :: vstar(np,np,2,nlev)                      ! velocity on Lagrangian surfaces
    real (kind=r8) :: dpdiss_biharmonic(np,np,nlev)            ! mean dp dissipation tendency, if nu_p>0
    real (kind=r8) :: dpdiss_ave(np,np,nlev)                   ! mean dp used to compute psdiss_tens

    ! diagnostics for explicit timestep
    real (kind=r8) :: phi(np,np,nlev)                          ! geopotential
    real (kind=r8) :: omega_p(np,np,nlev)                      ! vertical tendency (derived)
    real (kind=r8) :: eta_dot_dpdn(np,np,nlevp)                ! mean vertical flux from dynamics

    ! semi-implicit diagnostics: computed in explict-component, reused in Helmholtz-component.
    real (kind=r8) :: grad_lnps(np,np,2)                       ! gradient of log surface pressure
    real (kind=r8) :: zeta(np,np,nlev)                         ! relative vorticity
    real (kind=r8) :: div(np,np,nlev,timelevels)               ! divergence

    ! tracer advection fields used for consistency and limiters
    real (kind=r8) :: dp(np,np,nlev)                           ! for dp_tracers at physics timestep
    real (kind=r8) :: divdp(np,np,nlev)                        ! divergence of dp
    real (kind=r8) :: divdp_proj(np,np,nlev)                   ! DSSed divdp

    ! forcing terms for CAM
    real (kind=r8) :: FQ(np,np,nlev,qsize_d, 1)                ! tracer forcing
    real (kind=r8) :: FM(np,np,2,nlev, 1)                      ! momentum forcing
    real (kind=r8) :: FT(np,np,nlev, 1)                        ! temperature forcing
    real (kind=r8) :: omega_prescribed(np,np,nlev)             ! prescribed vertical tendency

    real (kind=r8) :: pecnd(np,np,nlev)                        ! pressure perturbation from condensate
    real (kind=r8) :: FQps(np,np,timelevels)                   ! forcing of FQ on ps_v

  end type derived_state_t

  !___________________________________________________________________
  type, public :: elem_accum_t


    ! the "4" timelevels represents data computed at:
    !  1  t-.5
    !  2  t+.5   after dynamics
    !  3  t+.5   after forcing
    !  4  t+.5   after Robert
    ! after calling TimeLevelUpdate, all times above decrease by 1.0

    real (kind=r8) :: KEner(np,np,4)
    real (kind=r8) :: PEner(np,np,4)
    real (kind=r8) :: IEner(np,np,4)
    real (kind=r8) :: IEner_wet(np,np,4)
    real (kind=r8) :: Qvar(np,np,qsize_d,4)                    ! Q variance at half time levels
    real (kind=r8) :: Qmass(np,np,qsize_d,4)                   ! Q mass at half time levels
    real (kind=r8) :: Q1mass(np,np,qsize_d)                    ! Q mass at full time levels

  end type elem_accum_t


! ============= DATA-STRUCTURES COMMON TO ALL SOLVERS ================

  type, public :: index_t
     integer :: ia(npsq),ja(npsq)
     integer :: is,ie
     integer :: NumUniquePts
     integer :: UniquePtOffset
  end type index_t

  !___________________________________________________________________
  type, public :: element_t
     integer :: LocalId
     integer :: GlobalId

     ! Coordinate values of element points
     type (spherical_polar_t) :: spherep(np,np)                       ! Spherical coords of GLL points

     ! Equ-angular gnomonic projection coordinates
     type (cartesian2D_t)     :: cartp(np,np)                         ! gnomonic coords of GLL points
     type (cartesian2D_t)     :: corners(4)                           ! gnomonic coords of element corners
     real (kind=r8)    :: u2qmap(4,2)                          ! bilinear map from ref element to quad in cubedsphere coordinates
                                                                      ! SHOULD BE REMOVED
     ! 3D cartesian coordinates
     type (cartesian3D_t)     :: corners3D(4)

     ! Element diagnostics
     real (kind=r8)    :: area                                 ! Area of element
     real (kind=r8)    :: normDinv                             ! some type of norm of Dinv used for CFL
     real (kind=r8)    :: dx_short                             ! short length scale in km
     real (kind=r8)    :: dx_long                              ! long length scale in km

     real (kind=r8)    :: variable_hyperviscosity(np,np)       ! hyperviscosity based on above
     real (kind=r8)    :: hv_courant                           ! hyperviscosity courant number
     real (kind=r8)    :: tensorVisc(np,np,2,2)                !og, matrix V for tensor viscosity

     ! Edge connectivity information
!     integer :: node_numbers(4)
!     integer :: node_multiplicity(4)                 ! number of elements sharing corner node

     type (GridVertex_t)      :: vertex                               ! element grid vertex information
     type (EdgeDescriptor_t)  :: desc

     type (elem_state_t)      :: state

     type (derived_state_t)   :: derived
     type (elem_accum_t)       :: accum
     ! Metric terms
     real (kind=r8)    :: met(np,np,2,2)                       ! metric tensor on velocity and pressure grid
     real (kind=r8)    :: metinv(np,np,2,2)                    ! metric tensor on velocity and pressure grid
     real (kind=r8)    :: metdet(np,np)                        ! g = SQRT(det(g_ij)) on velocity and pressure grid
     real (kind=r8)    :: rmetdet(np,np)                       ! 1/metdet on velocity pressure grid
     real (kind=r8)    :: D(np,np,2,2)                         ! Map covariant field on cube to vector field on the sphere
     real (kind=r8)    :: Dinv(np,np,2,2)                      ! Map vector field on the sphere to covariant v on cube

     ! Convert vector fields from spherical to rectangular components
     ! The transpose of this operation is its pseudoinverse.
     real (kind=r8)    :: vec_sphere2cart(np,np,3,2)

     ! Mass matrix terms for an element on a cube face
     real (kind=r8)    :: mp(np,np)                            ! mass matrix on v and p grid
     real (kind=r8)    :: rmp(np,np)                           ! inverse mass matrix on v and p grid

     ! Mass matrix terms for an element on the sphere
     ! This mass matrix is used when solving the equations in weak form
     ! with the natural (surface area of the sphere) inner product
     real (kind=r8)    :: spheremp(np,np)                      ! mass matrix on v and p grid
     real (kind=r8)    :: rspheremp(np,np)                     ! inverse mass matrix on v and p grid

     integer(i8) :: gdofP(np,np)                     ! global degree of freedom (P-grid)

     real (kind=r8)    :: fcor(np,np)                          ! Coreolis term

     type (index_t) :: idxP
     type (index_t),pointer :: idxV
     integer :: FaceNum

     ! force element_t to be a multiple of 8 bytes.
     ! on BGP, code will crash (signal 7, or signal 15) if 8 byte alignment is off
     ! check core file for:
     ! core.63:Generated by interrupt..(Alignment Exception DEAR=0xa1ef671c ESR=0x01800000 CCR0=0x4800a002)
     integer :: dummy
  end type element_t

  !___________________________________________________________________
  public :: element_coordinates
  public :: element_var_coordinates
  public :: element_var_coordinates3D
  public :: GetColumnIdP,GetColumnIdV
  public :: allocate_element_desc

contains

! ===================== ELEMENT_MOD METHODS ==========================

  function GetColumnIdP(elem,i,j) result(col_id)

    ! Get unique identifier for a Physics column on the P-grid

    type(element_t), intent(in) :: elem
    integer, intent(in) :: i,j
    integer :: col_id
    col_id = elem%gdofP(i,j)
  end function GetColumnIdP

  !___________________________________________________________________
  function GetColumnIdV(elem,i,j) result(col_id)

    !  Get unique identifier for a Physics column on the V-grid

    type(element_t), intent(in) :: elem
    integer, intent(in) :: i,j
    integer :: col_id
    col_id = elem%gdofP(i,j)
  end function GetColumnIdV

  !___________________________________________________________________
  function element_coordinates(start,end,points) result(cart)

    ! Initialize 2D rectilinear element colocation points

    type (cartesian2D_t), intent(in) :: start
    type (cartesian2D_t), intent(in) :: end
    real(r8),             intent(in) :: points(:)
    type (cartesian2D_t)             :: cart(SIZE(points),SIZE(points))

    type (cartesian2D_t) :: length, centroid
    real(r8)             :: y
    integer              :: i,j

    length%x   = 0.50D0*(end%x-start%x)
    length%y   = 0.50D0*(end%y-start%y)
    centroid%x = 0.50D0*(end%x+start%x)
    centroid%y = 0.50D0*(end%y+start%y)
    do j=1,SIZE(points)
       y = centroid%y + length%y*points(j)
       do i=1,SIZE(points)
          cart(i,j)%x = centroid%x + length%x*points(i)
          cart(i,j)%y = y
       end do
    end do
  end function element_coordinates

  !___________________________________________________________________
  function element_var_coordinates(c,points) result(cart)

    type (cartesian2D_t), intent(in) :: c(4)
    real(r8),             intent(in) :: points(:)
    type (cartesian2D_t)             :: cart(SIZE(points),SIZE(points))

    real(r8) :: p(size(points))
    real(r8) :: q(size(points))
    integer  :: i,j

    p(:) = (1.0D0-points(:))/2.0D0
    q(:) = (1.0D0+points(:))/2.0D0

    do j=1,SIZE(points)
       do i=1,SIZE(points)
          cart(i,j)%x = p(i)*p(j)*c(1)%x &
                      + q(i)*p(j)*c(2)%x &
                      + q(i)*q(j)*c(3)%x &
                      + p(i)*q(j)*c(4)%x
          cart(i,j)%y = p(i)*p(j)*c(1)%y &
                      + q(i)*p(j)*c(2)%y &
                      + q(i)*q(j)*c(3)%y &
                      + p(i)*q(j)*c(4)%y
       end do
    end do
  end function element_var_coordinates

  !___________________________________________________________________
  function element_var_coordinates3d(c,points) result(cart)

    type(cartesian3D_t), intent(in) :: c(4)
    real(r8), intent(in) :: points(:)

    type(cartesian3D_t)  :: cart(SIZE(points),SIZE(points))

    real(r8) :: p(size(points))
    real(r8) :: q(size(points)), r
    integer  :: i,j

    p(:) = (1.0D0-points(:))/2.0D0
    q(:) = (1.0D0+points(:))/2.0D0

    do j=1,SIZE(points)
       do i=1,SIZE(points)
          cart(i,j)%x = p(i)*p(j)*c(1)%x &
                      + q(i)*p(j)*c(2)%x &
                      + q(i)*q(j)*c(3)%x &
                      + p(i)*q(j)*c(4)%x
          cart(i,j)%y = p(i)*p(j)*c(1)%y &
                      + q(i)*p(j)*c(2)%y &
                      + q(i)*q(j)*c(3)%y &
                      + p(i)*q(j)*c(4)%y
          cart(i,j)%z = p(i)*p(j)*c(1)%z &
                      + q(i)*p(j)*c(2)%z &
                      + q(i)*q(j)*c(3)%z &
                      + p(i)*q(j)*c(4)%z

          ! project back to sphere:
          r = distance(cart(i,j))
          cart(i,j)%x = cart(i,j)%x/r
          cart(i,j)%y = cart(i,j)%y/r
          cart(i,j)%z = cart(i,j)%z/r
       end do
    end do
  end function element_var_coordinates3d

  !___________________________________________________________________
  subroutine allocate_element_desc(elem)

    type (element_t), intent(inout)   :: elem(:)
    integer                           :: num, j,i

    num = SIZE(elem)

    do j=1,num
       allocate(elem(j)%desc%putmapP(max_neigh_edges))
       allocate(elem(j)%desc%getmapP(max_neigh_edges))
       allocate(elem(j)%desc%putmapP_ghost(max_neigh_edges))
       allocate(elem(j)%desc%getmapP_ghost(max_neigh_edges))
       allocate(elem(j)%desc%putmapS(max_neigh_edges))
       allocate(elem(j)%desc%getmapS(max_neigh_edges))
       allocate(elem(j)%desc%reverse(max_neigh_edges))
       allocate(elem(j)%desc%globalID(max_neigh_edges))
       allocate(elem(j)%desc%loc2buf(max_neigh_edges))
       do i=1,max_neigh_edges
          elem(j)%desc%loc2buf(i)=i
          elem(j)%desc%globalID(i)=-1
       enddo

    end do
  end subroutine allocate_element_desc


end module element_mod
