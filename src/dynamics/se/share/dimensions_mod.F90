module dimensions_mod
  use shr_kind_mod, only: r8=>shr_kind_r8
  use constituents, only: qsize_d=>pcnst ! _EXTERNAL

  implicit none
  private

! set MAX number of tracers.  actual number of tracers is a run time argument  
  integer, parameter, public :: np = NP
  integer, parameter, public :: nc  = NC

  integer         :: qsize = 0

  integer, parameter, public :: nhe=1        !Max. Courant number
  integer, parameter, public :: nhr=2        !halo width needed for reconstruction - phl
!  integer, parameter, public :: ns=3         !quadratic halo interpolation - recommended setting for nc=3
!  integer, parameter, public :: ns=4         !cubic halo interpolation     - recommended setting for nc=4
  integer, parameter, public :: ns=NC

  !nhc determines width of halo exchanged with neighboring elements
  integer, parameter, public :: nhc = nhr+(nhe-1)+(ns-MOD(ns,2))/2

  integer, parameter, public :: npsq = np*np
  integer, parameter, public :: nlev=PLEV
  integer, parameter, public :: nlevp=nlev+1

  !default for non-refined mesh (note that these are *not* parameters now)
  integer, public  :: max_elements_attached_to_node = 4
  integer, public  :: s_nv = 6
  integer, public  :: max_corner_elem               = 1 !max_elements_attached_to_node-3
  integer, public  :: max_neigh_edges               = 8 !4 + 4*max_corner_elem


  public :: qsize,qsize_d

  integer, public :: ne
  integer, public :: nelem       ! total number of elements
  integer, public :: nelemd      ! number of elements per MPI task
  integer, public :: nelemdmax   ! max number of elements on any MPI task
  integer, public :: nPhysProc                          ! This is the number of physics processors/ per dynamics processor
  integer, public :: nnodes,npart,nmpi_per_node
  integer, public :: GlobalUniqueCols



  public :: set_mesh_dimensions

contains

  subroutine set_mesh_dimensions()

    ! new "params"
    max_elements_attached_to_node = 7  ! variable resolution
    s_nv = 2*max_elements_attached_to_node 

    !recalculate these
    max_corner_elem               = max_elements_attached_to_node-3
    max_neigh_edges               = 4 + 4*max_corner_elem


  end subroutine set_mesh_dimensions


end module dimensions_mod

