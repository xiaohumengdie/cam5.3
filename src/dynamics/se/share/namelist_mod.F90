module namelist_mod
  use cam_logfile,    only: iulog
  use params_mod,     only: recursive, sfcurve
  use shr_kind_mod,   only: r8=>shr_kind_r8
  use spmd_utils,     only: mpi_integer, mpi_real8, mpi_character, mpi_logical
  use control_mod,    only: MAX_STRING_LEN, MAX_FILE_LEN
  use control_mod,    only: partmethod, multilevel, numnodes, tasknum, remapfreq
  use control_mod,    only: nu_top, runtype, cubed_sphere_map, moisture, initial_total_mass
  use control_mod,    only: columnpackage, hypervis_power, hypervis_scaling
  use thread_mod,     only: omp_get_max_threads, max_num_threads, horz_num_threads, vert_num_threads, tracer_num_threads
  use dimensions_mod, only: ne, np, nnodes, nmpi_per_node, npart, set_mesh_dimensions
  use cam_abortutils, only: endrun
  use parallel_mod,   only: parallel_t, partitionfornodes, useframes
  use interpolate_mod,only: vector_uvars, vector_vvars, max_vecvars, interpolate_analysis, replace_vec_by_vordiv
  use interpolate_mod,only: set_interp_parameter, get_interp_parameter

!=============================================================================!
  implicit none
  private
!
! This module should contain no global data and should only be 'use'd to
!    call one of the public interfaces below
!
  public :: homme_set_defaults
  public :: homme_postprocess_namelist
  public :: read_analysis_nl

 contains

  ! ============================================
  ! homme_set_defaults:
  !
  !  Set default values for namelist variables
  !
  ! ============================================
  subroutine homme_set_defaults()
    npart               = 1
    useframes           = 0
    multilevel          = 1
    numnodes            = -1
    runtype             = 0
    remapfreq           = 240
    tasknum             =-1
    moisture            = "dry"
    columnpackage       = "none"
    nu_top              = 0
    initial_total_mass  = 0
    ne                  = 0

  end subroutine homme_set_defaults

  subroutine read_analysis_nl(par, NLFileName)
    use units, only : getunit, freeunit
    use mesh_mod, only : MeshOpen
    character(len=*), intent(in) :: NLFilename  ! Namelist filename
    type (parallel_t), intent(in) ::  par
    integer :: interp_nlat, interp_nlon, interp_gridtype, interp_type
    integer :: i, ii, j
    integer  :: ierr
    character(len=80) :: errstr, arg
    integer :: unitn
    character(len=*), parameter ::  subname = "homme:read_analysis_nl"


    namelist /analysis_nl/    &
        interp_nlat,          &
        interp_nlon,          &
        interp_gridtype,      &
        interp_type,          &
        interpolate_analysis


    if (par%masterproc) then

!      Default interpolation grid  (0 = auto compute based on ne,nv)  interpolation is off by default
#ifdef PIO_INTERP
       interpolate_analysis=.true.
#else
       interpolate_analysis=.false.
#endif
       interp_nlat =  0
       interp_nlon = 0
       interp_gridtype = 2
       interp_type = 0
       replace_vec_by_vordiv(:)=.false.
       vector_uvars(:)=''
       vector_vvars(:)=''
       vector_uvars(1:10) = (/'U       ','UBOT    ','U200    ','U250    ',&
            'U850    ','FU      ','CONVU   ','DIFFU   ','UTGWORO ','UFLX    '/)
       vector_vvars(1:10) = (/'V       ','VBOT    ','V200    ','V250    ',&
            'V850    ','FV      ','CONVV   ','DIFFV   ','VTGWORO ','VFLX    '/)

       write(iulog,*)"reading analysis namelist..."
       unitn=getunit()
       open( unitn, file=trim(nlfilename), status='old' )
       ierr = 1
       do while ( ierr > 0 )
          read (unitn,analysis_nl,iostat=ierr)
          if (ierr < 0) then
             print *,'analysis_nl namelist read returns an'// &
                  ' end of file or end of record condition, ignoring.'
          end if
       end do
       close( unitn )
       call freeunit( unitn )
    end if


    call MPI_barrier(par%comm,ierr)

    call MPI_bcast(interpolate_analysis, 7,MPI_logical,par%root,par%comm,ierr)
    call MPI_bcast(interp_nlat , 1,mpi_integer,par%root,par%comm,ierr)
    call MPI_bcast(interp_nlon , 1,mpi_integer,par%root,par%comm,ierr)
    call MPI_bcast(interp_gridtype , 1,mpi_integer,par%root,par%comm,ierr)
    call MPI_bcast(interp_type , 1,mpi_integer,par%root,par%comm,ierr)

    call MPI_bcast(replace_vec_by_vordiv ,MAX_VECVARS ,MPI_logical,par%root,par%comm,ierr)
    call MPI_bcast(vector_uvars ,10*MAX_VECVARS ,mpi_character,par%root,par%comm,ierr)
    call MPI_bcast(vector_vvars ,10*MAX_VECVARS ,mpi_character,par%root,par%comm,ierr)

    call set_interp_parameter('gridtype',interp_gridtype)
    call set_interp_parameter("itype",interp_type)
    if(any(interpolate_analysis)) then
       if (interp_nlon==0 .or. interp_nlat==0) then
          ! compute interpolation grid based on number of points around equator
          call set_interp_parameter('auto',4*ne*(np-1))
          interp_nlat = get_interp_parameter('nlat')
          interp_nlon = get_interp_parameter('nlon')
       else
          call set_interp_parameter('nlat',interp_nlat)
          call set_interp_parameter('nlon',interp_nlon)
       endif
    endif

    if (par%masterproc) then
       write(iulog,*)"done reading namelist..."

       write(iulog,*)" analysis interpolation = ", interpolate_analysis
       if(any(interpolate_analysis)) then
          write(iulog,*)" analysis interp nlat = ",interp_nlat
          write(iulog,*)" analysis interp nlon = ",interp_nlon
          write(iulog,*)" analysis interp gridtype = ",interp_gridtype
          write(iulog,*)" analysis interpolation type = ",interp_type
       end if
    endif

  end subroutine read_analysis_nl

  subroutine homme_postprocess_namelist(mesh_file, par)
    use mesh_mod,        only: MeshOpen
!    use dimensions_mod,  only: ntrac
    ! Dummy arguments
    character(len=*),  intent(in) :: mesh_file
    type (parallel_t), intent(in) :: par

    ! Local variable
    character(len=*), parameter :: subname = 'HOMME_POSTPROCESS_NAMELIST: '

!    if(par%masterproc) then
!      write(iulog, *) subname, 'omp_get_max_threads() = ', max_num_threads
!    end if
!
!    if((vert_num_threads > 1) .and. (limiter_option .ne. 8)) then
!       if(par%masterproc) then
!         write(iulog, *) subname, 'WARNING: vertical threading on supported for limiter_option != 8 '
!       end if
!       vert_num_threads = 1
!    endif

    if (ne /= 0) then
!      if (mesh_file /= "none" .and. mesh_file /= "/dev/null") then
      if (mesh_file /= "none") then
        if (par%masterproc) then
          write(iulog, *) subname, "mesh_file:", trim(mesh_file),        &
               " and ne:",ne," are both sepcified in the input file."
          write(iulog,*) "            Specify one or the other, but not both."
        end if
        call endrun(subname//"Do not specify ne if using a mesh file input.")
      end if
    end if
    if (par%masterproc) then
      write(iulog,*) subname, "Mesh File:", trim(mesh_file)
    end if
    if (ne == 0) then
      if (par%masterproc) then
        write (iulog,*) subname, "Opening Mesh File:", trim(mesh_file)
      end if
      call set_mesh_dimensions()
      call MeshOpen(mesh_file, par)
    end if

    ! set map
    if (cubed_sphere_map < 0) then
      if (ne == 0) then
        cubed_sphere_map = 2  ! element_local for var-res grids
      else
        cubed_sphere_map = 0  ! default is equi-angle gnomonic
      end if
    end if

    if (par%masterproc) then
      write (iulog,*) subname, "Reference element projection: cubed_sphere_map=",cubed_sphere_map
    end if

    !logic around different hyperviscosity options
    if (hypervis_power /= 0) then
      if (hypervis_scaling /= 0) then
        if (par%masterproc) then
          write(iulog, *) subname, 'Both hypervis_power and hypervis_scaling are nonzero.'
          write(iulog, *) '        (1) Set hypervis_power=1, hypervis_scaling=0 for HV based on an element area.'
          write(iulog, *) '        (2) Set hypervis_power=0 and hypervis_scaling=1 for HV based on a tensor.'
          write(iulog, *) '        (3) Set hypervis_power=0 and hypervis_scaling=0 for constant HV.'
        end if
        call endrun(subname//"ERROR: hypervis_power>0 and hypervis_scaling>0")
      end if
    end if

    if (multilevel <= 0) then
      nmpi_per_node = 1
    end if

    nnodes = npart / nmpi_per_node

    if((numnodes > 0) .and. (multilevel == 1)) then
      nnodes = numnodes
      nmpi_per_node = npart/nnodes
    end if

    ! ====================================================================
    !  Do not perform node level partitioning if you are only on one node
    ! ====================================================================
    if((nnodes .eq. 1) .and. PartitionForNodes) then
      PartitionForNodes = .FALSE.
    end if

  end subroutine homme_postprocess_namelist
end module namelist_mod
