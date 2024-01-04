! This module contains constants and namelist variables used through out the model
! to avoid circular dependancies please do not 'use' any further modules here.
!
module control_mod
  use shr_kind_mod,     only: r8=>shr_kind_r8

  integer, public, parameter :: MAX_STRING_LEN=240
  integer, public, parameter :: MAX_FILE_LEN=240
!  character(len=MAX_STRING_LEN)    , public :: integration    ! time integration (only one currently supported is "explicit")

! Tracer transport type
! We potentially have five types of tracer advection. However, not all of them
! may be chosen at runtime due to compile-type restrictions on arrays
  integer, public, parameter :: TRACERTRANSPORT_SE_GLL          = 1
  integer, public, parameter :: TRACERTRANSPORT_SEMILAGRANG_GLL = 2
  integer, public            :: tracer_transport_type = TRACERTRANSPORT_SE_GLL
  integer, public, parameter :: TRACER_GRIDTYPE_GLL = 11
  integer, public            :: tracer_grid_type = TRACER_GRIDTYPE_GLL


!shallow water advection tests:
!kmass points to a level with density.  other levels contain test tracers

  integer, public  :: tstep_type= 0                           ! 0 = leapfrog
                                                              ! 1 = RK (foward-in-time)
  integer, public  :: rk_stage_user  = 0                      ! number of RK stages to use
  integer, public  :: ftype = 0                                ! Forcing Type
                                                               ! ftype = 0  HOMME ApplyColumn() type forcing process split
                                                               ! ftype = -1   ignore forcing  (used for testing energy balance)
  integer, public  :: use_cpstar=0                             ! use cp or cp* in T equation                               
  integer, public  :: energy_fixer = 0    !  -1: No fixer, use non-staggered formula
                                          !   0: No Fixer, use staggered in time formula
                                          !       (only for leapfrog)
                                          !   1 or 4:  Enable fixer, non-staggered formula
                                              
  integer, public :: qsplit = 1           ! ratio of dynamics tsteps to tracer tsteps
  integer, public :: rsplit = 0           ! for vertically lagrangian dynamics, apply remap
                                          ! every rsplit tracer timesteps
  integer, public :: physics = 0          ! Defines if the program is to use its own physics (HOMME standalone), valid values 1,2,3
                                          ! physics = 0, no physics
                                          ! physics = 1, Use physics


! vert_remap_q_alg:    0  default value, Zerroukat monotonic splines
!                      1  PPM vertical remap with mirroring at the boundaries
!                         (solid wall bc's, high-order throughout)
!                      2  PPM vertical remap without mirroring at the boundaries
!                         (no bc's enforced, first-order at two cells bordering top and bottom boundaries)
  integer, public :: vert_remap_q_alg = 0


 integer, public :: cubed_sphere_map = -1  ! -1 = chosen at run time
                                           !  0 = equi-angle Gnomonic (default)
                                           !  1 = equi-spaced Gnomonic (not yet coded)
                                           !  2 = element-local projection  (for var-res)
                                           !  3 = parametric (not yet coded)

!tolerance to define smth small, was introduced for lim 8 in 2d and 3d
  real (kind=r8), public, parameter :: tol_limiter=1.0e-13_r8

  integer              , public :: limiter_option = 0

  integer              , public :: partmethod     ! partition methods
  character(len=MAX_STRING_LEN)    , public :: topology       ! options: "cube" is supported
  integer              , public :: tasknum
  integer              , public :: remapfreq      ! remap frequency of synopsis of system state (steps)
  character(len=MAX_STRING_LEN) :: remap_type     ! selected remapping option
  integer              , public :: statefreq      ! output frequency of synopsis of system state (steps)
  integer              , public :: runtype
  integer              , public :: timerdetail
  integer              , public :: numnodes
  integer              , public :: multilevel

  character(len=MAX_STRING_LEN)    , public :: columnpackage
  character(len=MAX_STRING_LEN)    , public :: moisture


  integer              , public :: fine_ne = -1              ! set for refined exodus meshes (variable viscosity)
  real (kind=r8), public :: max_hypervis_courant = 1d99 ! upper bound for Courant number
                                                               ! (only used for variable viscosity, recommend 1.9 in namelist)
  real (kind=r8), public :: nu      = 7.0D5           ! viscosity (momentum equ)
  real (kind=r8), public :: nu_div  = -1              ! viscsoity (momentum equ, div component)
  real (kind=r8), public :: nu_s    = -1              ! default = nu   T equ. viscosity
  real (kind=r8), public :: nu_q    = -1              ! default = nu   tracer viscosity
  real (kind=r8), public :: nu_p    = 0.0D5           ! default = 0    ps equ. viscosity
  real (kind=r8), public :: nu_top  = 0.0D5           ! top-of-the-model viscosity
  integer, public :: hypervis_subcycle=1    ! number of subcycles for hyper viscsosity timestep
  integer, public :: hypervis_subcycle_q=1  ! number of subcycles for hyper viscsosity timestep on TRACERS
  integer, public :: hypervis_order=0     ! laplace**hypervis_order.  0=not used  1=regular viscosity, 2=grad**4
  integer, public :: psurf_vis = 0        ! 0 = use laplace on eta surfaces
                                          ! 1 = use (approx.) laplace on p surfaces

  real (kind=r8), public :: hypervis_power=0     ! if not 0, use variable hyperviscosity based on element area
  real (kind=r8), public :: hypervis_scaling=0      ! use tensor hyperviscosity

!
!three types of hyper viscosity are supported right now:
! (1) const hv:    nu * del^2 del^2
! (2) scalar hv:   nu(lat,lon) * del^2 del^2
! (3) tensor hv,   nu * ( \div * tensor * \grad ) * del^2
!
! (1) default:  hypervis_power=0, hypervis_scaling=0
! (2) Original version for var-res grids. (M. Levy)
!            scalar coefficient within each element
!            hypervisc_scaling=0
!            set hypervis_power>0 and set fine_ne, max_hypervis_courant
! (3) tensor HV var-res grids
!            tensor within each element:
!            set hypervis_scaling > 0 (typical values would be 3.2 or 4.0)
!            hypervis_power=0
!            (\div * tensor * \grad) operator uses cartesian laplace
!


  ! hyperviscosity parameters used for smoothing topography
  integer, public :: smooth_phis_numcycle = 0   ! 0 = disable
  integer, public :: smooth_sgh_numcycle = 0   ! 0 = disabled
  real (kind=r8), public :: smooth_phis_nudt = 0

  real (kind=r8), public :: initial_total_mass = 0    ! initial perturbation in JW test case
  integer, public, parameter :: west  = 1
  integer, public, parameter :: east  = 2
  integer, public, parameter :: south = 3
  integer, public, parameter :: north = 4

  integer, public, parameter :: swest = 5
  integer, public, parameter :: seast = 6
  integer, public, parameter :: nwest = 7
  integer, public, parameter :: neast = 8

end module control_mod
