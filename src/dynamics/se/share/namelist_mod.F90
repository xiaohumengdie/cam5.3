module namelist_mod
  !-----------------
  use cam_logfile,    only: iulog
  !-----------------
  use params_mod,     only: recursive, sfcurve
  !-----------------
  use cube_mod, only : rotate_grid
  use shr_kind_mod,   only: r8=>shr_kind_r8
  use spmd_utils,     only: mpi_integer, mpi_real8, mpi_character, mpi_logical
  !-----------------
  use control_mod,    only:   &
       MAX_STRING_LEN,&
       MAX_FILE_LEN,  &
       partmethod,    &       ! Mesh partitioning method (METIS)
       topology,      &       ! Mesh topology
       multilevel,    &
       numnodes,      &
       tasknum,       &       ! used dg model in AIX machine
       remapfreq,     &       ! number of steps per remapping call
       remap_type,    &       ! selected remapping option
       statefreq,     &       ! number of steps per printstate call
       runtype,       &
       tstep_type, &
       cubed_sphere_map, &
       qsplit, &
       rsplit, &
       physics, &
       rk_stage_user, &
       ftype,        &
       energy_fixer,        &
       limiter_option, &
       fine_ne,       &
       max_hypervis_courant, &
       nu,            &
       nu_s,          &
       nu_q,          &
       nu_div,          &
       nu_p,          &
       nu_top,        &
       hypervis_scaling,   & ! use tensor HV instead of scalar coefficient
       psurf_vis,    &
       hypervis_order,    &
       hypervis_power,    &
       hypervis_subcycle, &
       hypervis_subcycle_q, &
       smooth_phis_numcycle, &
       smooth_sgh_numcycle, &
       smooth_phis_nudt, &
       initial_total_mass, &  ! set > 0 to set the initial_total_mass
       columnpackage, &
       moisture,      &
       vert_remap_q_alg
      

  !-----------------
  use thread_mod, only : omp_get_max_threads, max_num_threads, horz_num_threads, vert_num_threads, tracer_num_threads
  !-----------------
  use dimensions_mod, only : ne, np, nnodes, nmpi_per_node, npart, qsize, qsize_d, set_mesh_dimensions
  !-----------------
  use time_mod, only : nsplit, smooth, phys_tscale
  !-----------------
  use cam_abortutils, only: endrun
  use parallel_mod,   only: parallel_t, partitionfornodes, useframes
  !-----------------

  use interpolate_mod, only : vector_uvars, vector_vvars, max_vecvars, interpolate_analysis, replace_vec_by_vordiv
  use interpolate_mod, only : set_interp_parameter, get_interp_parameter

!=============================================================================!
  implicit none
  private
!
! This module should contain no global data and should only be 'use'd to
!    call one of the public interfaces below
!
  public :: readnl

 contains

  ! ============================================
  ! readnl:
  !
  !  Read in the namelists...
  !
  ! ============================================
  subroutine readnl(par, NLFileName)
    use units, only : getunit, freeunit
    use mesh_mod, only : MeshOpen
    character(len=*), intent(in) :: NLFilename  ! Namelist filename
    type (parallel_t), intent(in) ::  par
    character(len=MAX_FILE_LEN) :: mesh_file
    integer :: se_ftype, se_limiter_option
    integer :: se_phys_tscale, se_nsplit
    integer :: interp_nlat, interp_nlon, interp_gridtype, interp_type
    integer :: i, ii, j
    integer  :: ierr
    character(len=80) :: errstr, arg
    real(kind=r8) :: dt_max
    character(len=MAX_STRING_LEN) :: se_topology
    integer :: se_partmethod
    integer :: se_ne
    integer :: unitn
    character(len=*), parameter ::  subname = "homme:namelist_mod"
! These items are only here to keep readnl from crashing. Remove when possible
    integer :: se_fv_nphys
    character(len=80)  :: se_write_phys_grid
    character(len=256) :: se_phys_grid_file
    ! ============================================
    ! Namelists
    ! ============================================

    namelist /ctl_nl/ PARTMETHOD,    &       ! Mesh partitioning method (METIS)
                      TOPOLOGY,      &       ! Mesh topology
                     se_partmethod,    &
                     se_topology,      &
                     se_ne,            &
                     se_limiter_option, &
                     npart,         &
                     multilevel,    &
                     useframes,     &
                     numnodes,      &
                     ne,            &       ! element resolution factor
                     tasknum,       &
                     remapfreq,     &       ! number of steps per remapping call
                     remap_type,    &       ! selected remapping option
                     statefreq,     &       ! number of steps per printstate call
                     tstep_type, &
                     cubed_sphere_map, &
                     qsplit, &
                     rsplit, &
                     physics, &             ! The type of physics, 0=none, 1=multicloud or 2= emanuel.
                     rk_stage_user, &
                     se_ftype,        &       ! forcing type
                     energy_fixer,        &       ! forcing type
                     fine_ne,       &
                     max_hypervis_courant, &
                     nu,            &
                     nu_s,          &
                     nu_q,          &
                     nu_div,          &
                     nu_p,          &
                     nu_top,        &
                     psurf_vis,    &
                     hypervis_order,    &
                     hypervis_power,    &
                     hypervis_subcycle, &
                     hypervis_subcycle_q, &
                     hypervis_scaling, &
                     smooth_phis_numcycle, &
                     smooth_sgh_numcycle, &
                     smooth_phis_nudt, &
                     initial_total_mass, &
                     rotate_grid,   &
                     mesh_file,     &               ! Name of mesh file
                     vert_remap_q_alg

    namelist  /ctl_nl/ SE_NSPLIT,  &       ! number of dynamics steps per physics timestep
                       se_phys_tscale, &
! These items are only here to keep readnl from crashing. Remove when possible
                       se_fv_nphys,    &      ! Linear size of FV physics grid
                       se_write_phys_grid, &  ! Write physics grid file if .true.
                       se_phys_grid_file      ! Physics grid filename

    namelist /analysis_nl/    &
        interp_nlat,          &
        interp_nlon,          &
        interp_gridtype,      &
        interp_type,          &
        interpolate_analysis

!=======================================================================================================!
    ! ==========================
    ! Set the default partmethod
    ! ==========================

    PARTMETHOD    = RECURSIVE
    npart         = 1
    useframes     = 0
    multilevel    = 1
    ! set all CAM defaults
    ! CAM requires forward-in-time, subcycled dynamics
    ! RK2 3 stage tracers, sign-preserving conservative
    tstep_type              = 1      ! forward-in-time RK methods
    qsplit=4; rk_stage_user=3
    se_limiter_option=4
    se_ftype = 2
    energy_fixer = -1      ! no fixer, non-staggered-in-time formulas
    se_partmethod = -1
    se_ne       = -1
    se_topology = 'none'
    se_phys_tscale=0
    se_nsplit = 1
    qsize = qsize_d
    numnodes      = -1
    runtype       = 0
    statefreq     = 1
    remapfreq     = 240
    remap_type    = "parabolic"
    tasknum       =-1
    moisture      = "dry"
    columnpackage = "none"
    nu_top=0
    initial_total_mass=0
    mesh_file='none'
    ne              = 0


    ! =======================
    ! Read namelist variables
    ! =======================


    if (par%masterproc) then

       write(iulog,*)"reading ctl namelist..."
       unitn=getunit()
       open( unitn, file=trim(nlfilename), status='old' )
       ierr = 1
       do while ( ierr /= 0 )
          read (unitn,ctl_nl,iostat=ierr)
          if (ierr < 0) then
            write(6,*) 'ierr =',ierr
             call endrun( subname//':: namelist read returns an'// &
                  ' end of file or end of record condition' )
          end if
       end do
       close( unitn )
       call freeunit( unitn )
#ifndef _USEMETIS
      !=================================
      ! override the selected partition
      ! method and set it to SFCURVE
      !=================================
      PARTMETHOD = SFCURVE
#endif

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

    if(se_partmethod /= -1) partmethod = se_partmethod
    if(se_ne /= -1) ne = se_ne
    if(se_topology .ne. 'none') topology = se_topology

    call MPI_barrier(par%comm,ierr)

    npart  = par%nprocs

    ! =====================================
    !  Spread the namelist stuff around
    ! =====================================

    call MPI_bcast(PARTMETHOD ,1,mpi_integer,par%root,par%comm,ierr)
    call MPI_bcast(TOPOLOGY     ,MAX_STRING_LEN,mpi_character  ,par%root,par%comm,ierr)
    call MPI_bcast(tasknum ,1,mpi_integer,par%root,par%comm,ierr)

    call MPI_bcast( ne        ,1,mpi_integer,par%root,par%comm,ierr)
    call MPI_bcast(qsize     ,1,mpi_integer,par%root,par%comm,ierr)


    call MPI_bcast(remapfreq ,1,mpi_integer,par%root,par%comm,ierr)
    call MPI_bcast(remap_type, MAX_STRING_LEN, mpi_character, par%root, par%comm, ierr)
    call MPI_bcast(statefreq ,1,mpi_integer,par%root,par%comm,ierr)
    call MPI_bcast(multilevel ,1,mpi_integer,par%root,par%comm,ierr)
    call MPI_bcast(useframes ,1,mpi_integer,par%root,par%comm,ierr)
    call MPI_bcast(runtype   ,1,mpi_integer,par%root,par%comm,ierr)
    phys_tscale = se_phys_tscale
    limiter_option  = se_limiter_option
    nsplit = se_nsplit
    call MPI_bcast(smooth    ,1,MPI_real8   ,par%root,par%comm,ierr)
    call MPI_bcast(phys_tscale,1,MPI_real8   ,par%root,par%comm,ierr)
    call MPI_bcast(NSPLIT,1,mpi_integer,par%root,par%comm,ierr)
    call MPI_bcast(limiter_option  ,1,mpi_integer   ,par%root,par%comm,ierr)
    call MPI_bcast(se_ftype     ,1,mpi_integer   ,par%root,par%comm,ierr)
    call MPI_bcast(energy_fixer,1,mpi_integer   ,par%root,par%comm,ierr)
    call MPI_bcast(vert_remap_q_alg,1,mpi_integer   ,par%root,par%comm,ierr)

    call MPI_bcast(fine_ne    ,1,mpi_integer,par%root,par%comm,ierr)
    call MPI_bcast(max_hypervis_courant,1,MPI_real8   ,par%root,par%comm,ierr)
    call MPI_bcast(nu         ,1,MPI_real8   ,par%root,par%comm,ierr)
    call MPI_bcast(nu_s         ,1,MPI_real8   ,par%root,par%comm,ierr)
    call MPI_bcast(nu_q         ,1,MPI_real8   ,par%root,par%comm,ierr)
    call MPI_bcast(nu_div       ,1,MPI_real8   ,par%root,par%comm,ierr)
    call MPI_bcast(nu_p         ,1,MPI_real8   ,par%root,par%comm,ierr)
    call MPI_bcast(nu_top   ,1,MPI_real8   ,par%root,par%comm,ierr)
    call MPI_bcast(psurf_vis,1,mpi_integer   ,par%root,par%comm,ierr)
    call MPI_bcast(hypervis_order,1,mpi_integer   ,par%root,par%comm,ierr)
    call MPI_bcast(hypervis_power,1,MPI_real8   ,par%root,par%comm,ierr)
    call MPI_bcast(hypervis_scaling,1,MPI_real8   ,par%root,par%comm,ierr)
    call MPI_bcast(hypervis_subcycle,1,mpi_integer   ,par%root,par%comm,ierr)
    call MPI_bcast(hypervis_subcycle_q,1,mpi_integer   ,par%root,par%comm,ierr)
    call MPI_bcast(smooth_phis_numcycle,1,mpi_integer   ,par%root,par%comm,ierr)
    call MPI_bcast(smooth_sgh_numcycle,1,mpi_integer   ,par%root,par%comm,ierr)
    call MPI_bcast(smooth_phis_nudt,1,MPI_real8   ,par%root,par%comm,ierr)
    call MPI_bcast(initial_total_mass ,1,MPI_real8   ,par%root,par%comm,ierr)
    call MPI_bcast(rotate_grid   ,1,MPI_real8   ,par%root,par%comm,ierr)
    call MPI_bcast(mesh_file,MAX_FILE_LEN,mpi_character ,par%root,par%comm,ierr)
    call MPI_bcast(tstep_type,1,mpi_integer ,par%root,par%comm,ierr)
    call MPI_bcast(cubed_sphere_map,1,mpi_integer ,par%root,par%comm,ierr)
    call MPI_bcast(qsplit,1,mpi_integer ,par%root,par%comm,ierr)
    call MPI_bcast(rsplit,1,mpi_integer ,par%root,par%comm,ierr)
    call MPI_bcast(physics,1,mpi_integer ,par%root,par%comm,ierr)
    call MPI_bcast(rk_stage_user,1,mpi_integer ,par%root,par%comm,ierr)
    call MPI_bcast(moisture,MAX_STRING_LEN,mpi_character ,par%root,par%comm,ierr)
    call MPI_bcast(columnpackage,MAX_STRING_LEN,mpi_character,par%root,par%comm,ierr)


    if (mesh_file /= "none" .AND. ne /=0) then
      write (*,*) "namelist_mod: mesh_file:",mesh_file, &
                  " and ne:",ne," are both sepcified in the input file."
      write (*,*) "Specify one or the other, but not both."
      call endrun("Do not specify ne if using a mesh file input.")
    end if
    if (par%masterproc) write (iulog,*) "Mesh File:", trim(mesh_file)
    if (ne.eq.0) then
       if (par%masterproc) write (iulog,*) "Opening Mesh File:", trim(mesh_file)
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

    ftype = se_ftype

    rk_stage_user=3  ! 3d PRIM code only supports 3 stage RK tracer advection


    ! CHECK phys timescale, requires se_ftype=0 (pure tendencies for forcing)
    if (phys_tscale/=0) then
       if (ftype>0) call endrun('user specified se_phys_tscale requires se_ftype<=0')
    endif
    if (limiter_option==8 .or. limiter_option==84) then
       if (hypervis_subcycle_q/=1) then
          call endrun('limiter 8,84 requires hypervis_subcycle_q=1')
       endif
    endif

!=======================================================================================================!
    nmpi_per_node=1
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

    ! some default diffusion coefficiets
    if(nu_s<0) nu_s=nu
    if(nu_q<0) nu_q=nu
    if(nu_div<0) nu_div=nu


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

    if (par%masterproc) then
       write(iulog,*)"done reading namelist..."


       write(iulog,*)"readnl: topology      = ",TRIM( TOPOLOGY )

       write(iulog,*)"readnl: ne,np         = ",NE,np
       write(iulog,*)"readnl: partmethod    = ",PARTMETHOD
       write(iulog,*)'readnl: nmpi_per_node = ',nmpi_per_node
       write(iulog,*)'readnl: multilevel    = ',multilevel
       write(iulog,*)'readnl: useframes     = ',useframes
       write(iulog,*)'readnl: nnodes        = ',nnodes
       write(iulog,*)'readnl: npart         = ',npart

       print *
       write(iulog,*)"readnl: tstep_type    = ",tstep_type
       write(iulog,*)"readnl: vert_remap_q_alg  = ",vert_remap_q_alg
       write(iulog,*)"readnl: se_nsplit         = ", NSPLIT
       write(iulog,*)"readnl: se_ftype          = ",ftype
       write(iulog,*)"readnl: se_limiter_option = ",limiter_option
       write(iulog,*)"readnl: qsplit        = ",qsplit
       write(iulog,*)"readnl: vertical remap frequency rsplit (0=disabled): ",rsplit
       write(iulog,*)"readnl: physics       = ",physics

       write(iulog,*)"readnl: energy_fixer  = ",energy_fixer
       write(iulog,*)"readnl: runtype       = ",runtype


       if (hypervis_power /= 0)then
          write(iulog,*)"Variable scalar hyperviscosity: hypervis_power=",hypervis_power
          write(iulog,*)"max_hypervis_courant = ", max_hypervis_courant
          write(iulog,*)"Equivalent ne in fine region = ", fine_ne
       elseif(hypervis_scaling /=0)then
          write(iulog,*)"Tensor hyperviscosity:  hypervis_scaling=",hypervis_scaling
       else
          write(iulog,*)"Constant (hyper)viscosity used."
       endif

       write(iulog,*)"hypervis_subcycle, hypervis_subcycle_q = ",&
            hypervis_subcycle,hypervis_subcycle_q
       !write(iulog,*)"psurf_vis: ",psurf_vis
       write(iulog,'(a,2e9.2)')"viscosity:  nu (vor/div) = ",nu,nu_div
       write(iulog,'(a,2e9.2)')"viscosity:  nu_s      = ",nu_s
       write(iulog,'(a,2e9.2)')"viscosity:  nu_q      = ",nu_q
       write(iulog,'(a,2e9.2)')"viscosity:  nu_p      = ",nu_p
       write(iulog,'(a,2e9.2)')"viscosity:  nu_top      = ",nu_top
       write(iulog,*)"PHIS smoothing:  ",smooth_phis_numcycle,smooth_phis_nudt
       write(iulog,*)"SGH  smoothing:  ",smooth_sgh_numcycle

       if(initial_total_mass>0) then
          write(iulog,*) "initial_total_mass = ",initial_total_mass
       end if

       write(iulog,*)" analysis interpolation = ", interpolate_analysis
       if(any(interpolate_analysis)) then
          write(iulog,*)" analysis interp nlat = ",interp_nlat
          write(iulog,*)" analysis interp nlon = ",interp_nlon
          write(iulog,*)" analysis interp gridtype = ",interp_gridtype
          write(iulog,*)" analysis interpolation type = ",interp_type
       end if
!=======================================================================================================!
    endif

  end subroutine readnl

end module namelist_mod
