Module dyn_comp

  use shr_kind_mod,           only: r8=>shr_kind_r8, shr_kind_cl
  use element_mod, only : element_t, elem_state_t
  use time_mod, only : TimeLevel_t, nsplit
  use hybvcoord_mod, only : hvcoord_t
  use hybrid_mod, only : hybrid_t
  use perf_mod, only: t_startf, t_stopf
  use cam_logfile, only : iulog
  use time_manager, only: is_first_step
  use spmd_utils,  only : iam, npes_cam => npes
  use pio,         only: file_desc_t
  use time_mod, only : phys_tscale

  implicit none
  private
  save

  ! PUBLIC MEMBER FUNCTIONS:
  public dyn_init1, dyn_init2, dyn_run, dyn_final

  ! PUBLIC DATA MEMBERS:
  public dyn_import_t, dyn_export_t


  type (TimeLevel_t)   , public :: TimeLevel     ! main time level struct (used by tracers)

  type dyn_import_t
     type (element_t), pointer :: elem(:) => null()
  end type dyn_import_t

  type dyn_export_t
     type (element_t), pointer :: elem(:) => null()
  end type dyn_export_t
  type (hvcoord_t), public  :: hvcoord

  integer, parameter  ::  DYN_RUN_SUCCESS           = 0
  integer, parameter  ::  DYN_RUN_FAILURE           = -1

  ! !DESCRIPTION: This module implements the SE Dynamical Core as
  !               an ESMF gridded component.  It is specific to SE
  !               and does not use ESMF.
  !
  ! \paragraph{Overview}
  !
  !   This module contains an ESMF wrapper for the SE
  !   Dynamical Core used in the Community Atmospheric Model. 
  !
  ! !REVISION HISTORY:
  !
  !  JPE  06.05.31:  created
  !
  !----------------------------------------------------------------------

  ! Enumeration of DYNAMICS_IN_COUPLINGS


  logical, parameter         :: DEBUG = .true.

  real(r8), parameter        :: ONE    = 1.0_r8

  character(*), parameter, public :: MODULE_NAME = "dyn_comp"
  character(*), parameter, public :: VERSION     = "$Id$" 

  ! Frontogenesis indices
  integer, public :: frontgf_idx = -1
  integer, public :: frontga_idx = -1

!===============================================================================
contains
!===============================================================================

subroutine dyn_readnl(NLFileName)

   use namelist_utils, only: find_group_name
   use namelist_mod,   only: homme_set_defaults, homme_postprocess_namelist, read_analysis_nl
   use units,          only: getunit, freeunit
   use spmd_utils,     only: masterproc, masterprocid, mpicom, npes
   use spmd_utils,     only: mpi_real8, mpi_integer, mpi_character, mpi_logical

   use control_mod,    only: TRACERTRANSPORT_SE_GLL, tracer_transport_type
   use control_mod,    only: TRACER_GRIDTYPE_GLL, tracer_grid_type
   use control_mod,    only: energy_fixer, hypervis_order, hypervis_subcycle
   use control_mod,    only: hypervis_subcycle_q, statefreq
   use control_mod,    only: nu, nu_div, nu_p, nu_q, nu_top, qsplit, rsplit, nu_s
   use control_mod,    only: vert_remap_q_alg, tstep_type, rk_stage_user
   use control_mod,    only: ftype, limiter_option, partmethod
   use control_mod,    only: topology, numnodes, tasknum, moisture
   use control_mod,    only: columnpackage, remapfreq, remap_type
   use control_mod,    only: initial_total_mass
   use control_mod,    only: fine_ne, hypervis_power, hypervis_scaling
   use control_mod,    only: max_hypervis_courant
    use cam_abortutils, only: endrun
    use dimensions_mod, only: qsize, qsize_d, npsq, ne, npart
    use constituents,   only: pcnst
    use params_mod,     only: SFCURVE
    use parallel_mod,   only: par, initmpi
!!XXgoldyXX: v For future CSLAM / physgrid commit
!    use dp_grids,       only: fv_nphys, fv_nphys2, nphys_pts, write_phys_grid, phys_grid_file
!!XXgoldyXX: ^ For future CSLAM / physgrid commit

    ! Dummy argument
    character(len=*), intent(in) :: NLFileName

    ! Local variables
    integer                      :: unitn, ierr
   ! SE Namelist variables
!    integer                      :: se_fine_ne
    integer                      :: se_ftype
    integer                      :: se_hypervis_order
    integer                      :: se_hypervis_subcycle
    integer                      :: se_hypervis_subcycle_q
    integer                      :: se_limiter_option
!    real(r8)                     :: se_max_hypervis_courant
    character(len=SHR_KIND_CL)   :: se_mesh_file
    integer                      :: se_ne
    integer                      :: se_npes
    integer                      :: se_nsplit
    real(r8)                     :: se_nu
    real(r8)                     :: se_nu_div
    real(r8)                     :: se_nu_p
    real(r8)                     :: se_nu_q
    real(r8)                     :: se_nu_top
    integer                      :: se_qsplit
!    logical                      :: se_refined_mesh
    integer                      :: se_rsplit
    integer                      :: se_statefreq
    integer                      :: se_tstep_type
    integer                      :: se_phys_tscale
    integer                      :: se_vert_remap_q_alg
!!XXgoldyXX: v For future CSLAM / physgrid commit
!    character(len=METHOD_LEN)     :: se_tracer_transport_method
!    character(len=METHOD_LEN)     :: se_cslam_ideal_test
!    character(len=METHOD_LEN)     :: se_cslam_test_type
!    character(len=METHOD_LEN)     :: se_write_phys_grid
!    character(len=shr_kind_cl)    :: se_phys_grid_file
!    integer                       :: se_fv_nphys = 0
!!XXgoldyXX: ^ For future CSLAM / physgrid commit

    namelist /dyn_se_inparm/      &
         se_ftype,                & ! forcing type
         se_hypervis_order,       &
         se_hypervis_subcycle,    &
         se_hypervis_subcycle_q,  &
         se_limiter_option,       &
         se_ne,                   &
         se_npes,                 &
         se_nsplit,               & ! # of dynamics steps per physics timestep
         se_nu,                   &
         se_nu_div,               &
         se_nu_p,                 &
         se_nu_q,                 &
         se_nu_top,               &
         se_qsplit,               &
         se_rsplit,               &
         se_statefreq,            & ! number of steps per printstate call
         se_tstep_type,           &
         se_phys_tscale,          &
         se_vert_remap_q_alg

    !--------------------------------------------------------------------------

   ! namelist default values should be in namelist (but you know users . . .)
    ! NB: Of course, these should keep up with what is in namelist_defaults ...
    se_mesh_file            = 'none'
    se_ftype                = 2
    se_hypervis_order       = 2
    se_hypervis_subcycle    = 3
    se_hypervis_subcycle_q  = 1
    se_limiter_option       = 4
    se_ne                   = -1
    se_npes                 = npes
    se_nsplit               = 1
    se_nu                   = 1.0e15_r8
    se_nu_div               = 2.5e15_r8
    se_nu_p                 = 1.0e15_r8
    se_nu_q                 = -1.0_r8
    se_nu_top               = 2.5e5_r8
    se_qsplit               = 4
    se_rsplit               = 3
    se_statefreq            = 1
    se_tstep_type           = 1
    se_phys_tscale          = 0
    se_vert_remap_q_alg     = 1

    ! Read the namelist (dyn_se_inparm)
    call MPI_barrier(mpicom, ierr)
    if (masterproc) then
      write(iulog, *) "dyn_readnl: reading dyn_se_inparm namelist..."
      unitn = getunit()
      open( unitn, file=trim(NLFileName), status='old' )
      call find_group_name(unitn, 'dyn_se_inparm', status=ierr)
      if (ierr == 0) then
        read(unitn, dyn_se_inparm, iostat=ierr)
        if (ierr /= 0) then
          call endrun('dyn_readnl: ERROR reading dyn_se_inparm namelist')
        end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

    ! Broadcast namelist values to all PEs
    call MPI_bcast(se_ftype, 1, mpi_integer, masterprocid, mpicom, ierr)
    call MPI_bcast(se_hypervis_order, 1, mpi_integer, masterprocid, mpicom, ierr)
    call MPI_bcast(se_hypervis_subcycle, 1, mpi_integer, masterprocid, mpicom, ierr)
    call MPI_bcast(se_hypervis_subcycle_q, 1, mpi_integer, masterprocid, mpicom, ierr)
    call MPI_bcast(se_limiter_option, 1, mpi_integer, masterprocid, mpicom, ierr)
    call MPI_bcast(se_ne, 1, mpi_integer, masterprocid, mpicom, ierr)
    call MPI_bcast(se_npes, 1, mpi_integer, masterprocid, mpicom, ierr)
    call MPI_bcast(se_nsplit, 1, mpi_integer, masterprocid, mpicom, ierr)
    call MPI_bcast(se_nu, 1, mpi_real8, masterprocid, mpicom, ierr)
    call MPI_bcast(se_nu_div, 1, mpi_real8, masterprocid, mpicom, ierr)
    call MPI_bcast(se_nu_p, 1, mpi_real8, masterprocid, mpicom, ierr)
    call MPI_bcast(se_nu_q, 1, mpi_real8, masterprocid, mpicom, ierr)
    call MPI_bcast(se_nu_top, 1, mpi_real8, masterprocid, mpicom, ierr)
    call MPI_bcast(se_qsplit, 1, mpi_integer, masterprocid, mpicom, ierr)
    call MPI_bcast(se_rsplit, 1, mpi_integer, masterprocid, mpicom, ierr)
    call MPI_bcast(se_statefreq, 1, mpi_integer, masterprocid, mpicom, ierr)
    call MPI_bcast(se_tstep_type, 1, mpi_integer, masterprocid, mpicom, ierr)
    call MPI_bcast(se_phys_tscale, 1 ,mpi_real8, masterprocid, mpicom, ierr)
    call MPI_bcast(se_vert_remap_q_alg, 1, mpi_integer, masterprocid, mpicom, ierr)

    phys_tscale = se_phys_tscale

    ! Initialize the SE structure that holds the MPI decomposition information
    if (se_npes <= 0) then
      se_npes = npes
    end if
    par = initmpi(se_npes)

    ! Fix up unresolved default values
    ! default diffusion coefficiets
    if (se_nu_q < 0) then
      se_nu_q = se_nu
    end if
    if (se_nu_div < 0) then
      se_nu_div = se_nu
    end if

    ! Set HOMME defaults
    call homme_set_defaults()

    call read_analysis_nl(par, NLFileName)
    ! Set HOMME variables not in CAM's namelist but with different CAM defaults
    partmethod           = SFCURVE
    npart                = se_npes
    energy_fixer         = -1      ! no fixer, non-staggered-in-time formulas
    ! CAM requires forward-in-time, subcycled dynamics
    ! RK2 3 stage tracers, sign-preserving conservative
    rk_stage_user        = 3
!    integration          = 'explicit'
    qsize                = qsize_d
    topology             = "cube"
    ! Finally, set the HOMME variables which have different names
!    fine_ne              = se_fine_ne
    ftype                = se_ftype
    hypervis_order       = se_hypervis_order
    hypervis_subcycle    = se_hypervis_subcycle
    hypervis_subcycle_q  = se_hypervis_subcycle_q
    limiter_option       = se_limiter_option
!    max_hypervis_courant = se_max_hypervis_courant
    if(se_ne /= -1)   ne = se_ne
!    ne                   = se_ne
    nsplit               = se_nsplit
    nu                   = se_nu
    nu_div               = se_nu_div
    nu_p                 = se_nu_p
    nu_q                 = se_nu_q
    nu_top               = se_nu_top
    qsplit               = se_qsplit
    rsplit               = se_rsplit
    statefreq            = se_statefreq
    tstep_type           = se_tstep_type
    vert_remap_q_alg     = se_vert_remap_q_alg

    if(nu_s<0) nu_s=nu
    if(nu_q<0) nu_q=nu
    if(nu_div<0) nu_div=nu
    call homme_postprocess_namelist(se_mesh_file, par)

    if (phys_tscale/=0) then
       if (ftype>0) call endrun('user specified se_phys_tscale requires se_ftype<=0')
    endif

    if (masterproc) then
      write(iulog, '(a,i0)') 'dyn_readnl: se_ftype = ',se_ftype
      write(iulog, '(a,i0)') 'dyn_readnl: se_hypervis_order = ',se_hypervis_order
      write(iulog, '(a,i0)') 'dyn_readnl: se_hypervis_subcycle = ',se_hypervis_subcycle
      write(iulog, '(a,i0)') 'dyn_readnl: se_hypervis_subcycle_q = ',se_hypervis_subcycle_q
      write(iulog, '(a,i0)') 'dyn_readnl: se_limiter_option = ',se_limiter_option
!      if (.not. se_refined_mesh) then
!        write(iulog, '(a,i0)') 'dyn_readnl: se_ne = ',se_ne
!      end if
      write(iulog, '(a,i0)') 'dyn_readnl: se_npes = ',se_npes
      write(iulog, '(a,i0)') 'dyn_readnl: se_nsplit = ',se_nsplit
      write(iulog, '(a,e9.2)') 'dyn_readnl: se_nu = ',se_nu
      write(iulog, '(a,e9.2)') 'dyn_readnl: se_nu_div = ',se_nu_div
      write(iulog, '(a,e9.2)') 'dyn_readnl: se_nu_p = ',se_nu_p
      write(iulog, '(a,e9.2)') 'dyn_readnl: se_nu_q = ',se_nu_q
      write(iulog, '(a,e9.2)') 'dyn_readnl: se_nu_top = ',se_nu_top
      write(iulog, '(a,i0)') 'dyn_readnl: se_qsplit = ',se_qsplit
      write(iulog, '(a,i0)') 'dyn_readnl: se_rsplit = ',se_rsplit
      write(iulog, '(a,i0)') 'dyn_readnl: se_statefreq = ',se_statefreq
      write(iulog, '(a,i0)') 'dyn_readnl: se_tstep_type = ',se_tstep_type
      write(iulog, '(a,i0)') 'dyn_readnl: se_vert_remap_q_alg = ',se_vert_remap_q_alg
!      if (se_refined_mesh) then
!        write(iulog, *) 'dyn_readnl: Refined mesh simulation'
!        write(iulog, *) 'dyn_readnl: se_mesh_file = ',trim(se_mesh_file)
!        write(iulog, '(a,e11.4)') 'dyn_readnl: se_max_hypervis_courant = ',se_max_hypervis_courant
!      end if
      if (limiter_option==8 .or. limiter_option==84) then
         if (hypervis_subcycle_q/=1) then
            call endrun('limiter 8,84 requires hypervis_subcycle_q=1')
         endif
      endif

    end if

  end subroutine dyn_readnl

  subroutine dyn_init1(fh, NLFileName, dyn_in, dyn_out)

  ! Initialize the dynamical core

    use pio,             only: file_desc_t
    use hycoef,          only: hycoef_init
    use ref_pres,        only: ref_pres_init

    use dyn_grid,        only: dyn_grid_init, elem, get_dyn_grid_parm
    use dyn_grid,        only: set_horiz_grid_cnt_d

    use spmd_utils,      only: mpi_integer, mpicom, mpi_logical
    use native_mapping,  only: create_native_mapping_files, native_mapping_readnl
    use time_manager,    only: get_nstep, dtime

    use dimensions_mod,  only: globaluniquecols, nelem, nelemd, nelemdmax
    use parallel_mod,    only: par
    use prim_init,       only: prim_init1
    use thread_mod,      only: horz_num_threads
    use control_mod,     only: runtype, qsplit, rsplit
    use time_mod,        only: tstep

    use phys_control,    only: use_gw_front, use_gw_front_igw
    use physics_buffer,  only: pbuf_add_field, dtype_r8
    use ppgrid,          only: pcols, pver
    use pmgrid,              only: dyndecomp_set
    use rgrid,               only: fullgrid
    use interpolate_mod,     only: interpolate_analysis

    ! PARAMETERS:
    type(file_desc_t),   intent(in)  :: fh       ! PIO file handle for initial or restart file
    character(len=*),    intent(in)  :: NLFileName
    type (dyn_import_t), intent(OUT) :: dyn_in
    type (dyn_export_t), intent(OUT) :: dyn_out

#ifdef _OPENMP    
    integer omp_get_num_threads
#endif
    integer :: neltmp(3)
    logical :: nellogtmp(7)
    integer :: npes_se
    integer :: nthreads
    !----------------------------------------------------------------------

    if (use_gw_front .or. use_gw_front_igw) then
       call pbuf_add_field("FRONTGF", "global", dtype_r8, (/pcols,pver/), &
            frontgf_idx)
       call pbuf_add_field("FRONTGA", "global", dtype_r8, (/pcols,pver/), &
            frontga_idx)
    end if

    ! Initialize dynamics grid
    call dyn_grid_init()

    ! Read the SE specific part of the namelist
    call dyn_readnl(NLFileName)

    ! override the setting in the SE namelist, it's redundant anyway
    if (.not. is_first_step()) runtype = 1

    ! Initialize hybrid coordinate arrays.
    call hycoef_init(fh)

    ! Initialize physics grid reference pressures (needed by initialize_radbuffer)
    call ref_pres_init()

    ! legacy reduced grid code -- should be removed
    fullgrid=.true.

#ifdef _OPENMP    
    if(par%masterproc) then
       nthreads = omp_get_num_threads()
       write(iulog,*) " "
       write(iulog,*) "dyn_init1: number of OpenMP threads = ", nthreads
       write(iulog,*) " "
    endif
#if defined (COLUMN_OPENMP)
    if(par%masterproc) then
       write(iulog,*) " "
       write(iulog,*) "dyn_init1: using OpenMP within element instead of across elements"
       write(iulog,*) " "
    endif
#endif
#else
    if(par%masterproc) then
       write(iulog,*) " "
       write(iulog,*) "dyn_init1: openmp not activated"
       write(iulog,*) " "
    endif
#endif
    if(iam < par%nprocs) then
       call prim_init1(elem,par,TimeLevel)

       dyn_in%elem => elem
       dyn_out%elem => elem
    
       call set_horiz_grid_cnt_d(GlobalUniqueCols)

       neltmp(1) = nelemdmax
       neltmp(2) = nelem
       neltmp(3) = get_dyn_grid_parm('plon')
       nellogtmp(1:7) = interpolate_analysis(1:7)
    else
       nelemd = 0
       neltmp(1) = 0
       neltmp(2) = 0
       neltmp(3) = 0
       nellogtmp(1:7) = .true.
    endif

    dyndecomp_set = .true.
    if (par%nprocs .lt. npes_cam) then
! Broadcast quantities to auxiliary processes
#ifdef SPMD
       call mpibcast(neltmp, 3, mpi_integer, 0, mpicom)
       call mpibcast(nellogtmp, 7, mpi_logical, 0, mpicom)
#endif
      if (iam .ge. par%nprocs) then
        nelemdmax = neltmp(1)
        nelem     = neltmp(2)
        call set_horiz_grid_cnt_d(neltmp(3))
      interpolate_analysis(1:7) = nellogtmp(1:7)
     endif
    endif

    !
    ! This subroutine creates mapping files using SE basis functions if requested
    !
    call native_mapping_readnl(NLFileName)
    call create_native_mapping_files( par, elem,'native')
    call create_native_mapping_files( par, elem,'bilin')

    ! Dynamics timestep
    !
    !  Note: dtime = progress made in one timestep.  value in namelist
    !        dtime = the frequency at which physics is called
    !        tstep = the dynamics timestep:  
    !

    if (rsplit==0) then
       ! non-lagrangian code
       tstep = dtime/real(nsplit*qsplit,r8)
       TimeLevel%nstep = get_nstep()*nsplit*qsplit
   else
      ! lagrangian code
       tstep = dtime/real(nsplit*qsplit*rsplit,r8)
       TimeLevel%nstep = get_nstep()*nsplit*qsplit*rsplit
    endif

    ! initial SE (subcycled) nstep
    TimeLevel%nstep0 = 0

  end subroutine dyn_init1


  subroutine dyn_init2(dyn_in)
    use dimensions_mod,   only: nlev, nelemd
    use prim_driver_mod,  only: prim_init2
    use prim_si_ref_mod,  only: prim_set_mass
    use hybrid_mod,       only: config_thread_region, get_loop_ranges, init_loop_ranges
    use hycoef,           only: hyam, hybm, hyai, hybi, ps0
    use parallel_mod,     only: par
    use time_mod,         only: time_at
    use control_mod,      only: moisture, runtype
    use hybrid_mod,       only: config_thread_region, get_loop_ranges, init_loop_ranges
    use cam_control_mod,  only: aqua_planet, ideal_phys, adiabatic
    use comsrf,           only: landm, sgh, sgh30
    use nctopo_util_mod,  only: nctopo_util_driver
    use cam_instance,     only: inst_index
    use thread_mod,       only: horz_num_threads

    type (dyn_import_t), intent(inout) :: dyn_in

    type(element_t),    pointer :: elem(:)

    integer :: ithr, nets, nete, ie, k
    real(r8), parameter :: Tinit=300.0_r8
    real(r8) :: dyn_ps0
    type(hybrid_t) :: hybrid

    elem  => dyn_in%elem

    dyn_ps0=ps0
    hvcoord%hyam=hyam
    hvcoord%hyai=hyai
    hvcoord%hybm=hybm
    hvcoord%hybi=hybi
    hvcoord%ps0=dyn_ps0  
    do k=1,nlev
       hvcoord%hybd(k) = hvcoord%hybi(k+1) - hvcoord%hybi(k)
    end do
    call init_loop_ranges(nelemd)
    if(iam < par%nprocs) then


       if ( horz_num_threads == 1) then
         hybrid = config_thread_region(par,'serial')
       else
         hybrid = config_thread_region(par,'horizontal')
       endif
       call get_loop_ranges(hybrid,ibeg=nets,iend=nete)
!       write(*,200) nets, nete
!200 format(2x, 2i4, ' nets, nete - used init')
       moisture='moist'

       if(adiabatic) then
          moisture='dry'
          if(runtype == 0) then
             do ie=nets,nete
                elem(ie)%state%q(:,:,:,:)=0.0_r8
                elem(ie)%derived%fq(:,:,:,:,:)=0.0_r8
             end do
          end if
       else if(ideal_phys) then
          moisture='dry'
          if(runtype == 0) then
             do ie=nets,nete
                elem(ie)%state%lnps(:,:,:) =LOG(dyn_ps0)

                elem(ie)%state%ps_v(:,:,:) =dyn_ps0

                elem(ie)%state%phis(:,:)=0.0_r8

                elem(ie)%state%T(:,:,:,:) =Tinit

                elem(ie)%state%v(:,:,:,:,:) =0.0_r8

                elem(ie)%state%q(:,:,:,:)=0.0_r8

             end do
          end if
       else if(aqua_planet .and. runtype==0)  then
          do ie=nets,nete
             !          elem(ie)%state%lnps(:,:,:) =LOG(dyn_ps0)
             !          elem(ie)%state%ps_v(:,:,:) =dyn_ps0
             elem(ie)%state%phis(:,:)=0.0_r8
          end do
          if(allocated(landm)) landm=0.0_r8
          if(allocated(sgh)) sgh=0.0_r8
          if(allocated(sgh30)) sgh30=0.0_r8
       end if

       do ie=nets,nete
          elem(ie)%derived%FM=0.0_r8
          elem(ie)%derived%FT=0.0_r8
          elem(ie)%derived%FQ=0.0_r8
       end do

       ! scale PS to achieve prescribed dry mass
       if (runtype == 0) then
          ! new run, scale mass to value given in namelist, if needed
          call prim_set_mass(elem, TimeLevel,hybrid,hvcoord,nets,nete)
       endif
       call prim_init2(elem,hybrid,nets,nete, TimeLevel, hvcoord)
       !
       ! This subroutine is used to create nc_topo files, if requested
       ! 
       call nctopo_util_driver(elem,hybrid,nets,nete)
    end if

    if (inst_index == 1) then
       call write_grid_mapping(par, elem)
    end if

  end subroutine dyn_init2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !-----------------------------------------------------------------------
  !BOP
  ! !ROUTINE:  RUN --- Driver for the 
  !
  ! !INTERFACE:
  subroutine dyn_run( dyn_state, rc )

    ! !USES:
    use parallel_mod,     only: par
    use prim_driver_mod,  only: prim_run_subcycle
    use dimensions_mod,   only: nlev, nelemd
    use thread_mod,       only: horz_num_threads
    use time_mod,         only: tstep
    use hybrid_mod,       only: init_loop_ranges, get_loop_ranges, config_thread_region
!    use perf_mod, only: t_startf, t_stopf
    implicit none


    type (dyn_export_t), intent(inout)       :: dyn_state   !  container
    type(hybrid_t) :: hybrid

    integer, intent(out)               :: rc      ! Return code
    integer ::  n
    integer :: nets, nete, ithr
    integer :: ie

    ! !DESCRIPTION:
    !
    call init_loop_ranges(nelemd)
    if(iam < par%nprocs) then
       if (horz_num_threads == 1) then
         hybrid = config_thread_region(par,'serial')
       else
         hybrid = config_thread_region(par,'horizontal')
       endif
       call get_loop_ranges(hybrid,ibeg=nets,iend=nete)
!       write(*,200) nets, nete
!200 format(2x, 2i4, ' nets, nete - dyn_run')

       do n=1, nsplit
          ! forward-in-time RK, with subcycling
          call prim_run_subcycle(dyn_state%elem,hybrid,nets,nete,&
               tstep, TimeLevel, hvcoord, n)
       end do
    end if
    rc = DYN_RUN_SUCCESS

    !EOC
  end subroutine dyn_run
  !-----------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dyn_final(DYN_STATE, RESTART_FILE)

    type (elem_state_t), target     :: DYN_STATE
    character(LEN=*)   , intent(IN) :: RESTART_FILE



  end subroutine dyn_final



  subroutine write_grid_mapping(par, elem)
    use parallel_mod,     only: parallel_t
    use element_mod, only : element_t
    use cam_pio_utils, only : cam_pio_createfile, pio_subsystem
    use pio, only : file_desc_t, pio_def_dim, var_desc_t, pio_int, pio_def_var, &
         pio_enddef, pio_closefile, pio_initdecomp, io_desc_t, pio_write_darray, &
         pio_freedecomp, pio_setdebuglevel
    use dimensions_mod, only : np, nelem, nelemd
    use dof_mod, only : createmetadata

    type(parallel_t) :: par
    type(element_t) :: elem(:)
    type(file_desc_t) :: nc
    type(var_desc_t) :: vid
    type(io_desc_t) :: iodesc
    integer :: dim1, dim2, ierr, i, j, ie, cc, base, ii, jj
    integer, parameter :: npm12 = (np-1)*(np-1)
    integer :: subelement_corners(npm12*nelemd,4)
    integer :: dof(npm12*nelemd*4)


    ! Create a CS grid mapping file for postprocessing tools

       ! write meta data for physics on GLL nodes
       call cam_pio_createfile(nc, 'SEMapping.nc', 0)
   
       ierr = pio_def_dim(nc, 'ncenters', npm12*nelem, dim1)
       ierr = pio_def_dim(nc, 'ncorners', 4, dim2)
       ierr = pio_def_var(nc, 'element_corners', PIO_INT, (/dim1,dim2/),vid)
    
       ierr = pio_enddef(nc)
       call createmetadata(par, elem, subelement_corners)

       jj=0
       do cc=0,3
          do ie=1,nelemd
             base = ((elem(ie)%globalid-1)+cc*nelem)*npm12
             ii=0
             do j=1,np-1
                do i=1,np-1
                   ii=ii+1
                   jj=jj+1
                   dof(jj) = base+ii
                end do
             end do
          end do
       end do

       call pio_initdecomp(pio_subsystem, pio_int, (/nelem*npm12,4/), dof, iodesc)

       call pio_write_darray(nc, vid, iodesc, reshape(subelement_corners,(/nelemd*npm12*4/)), ierr)
       
       call pio_closefile(nc)

       call pio_freedecomp(pio_subsystem, iodesc)

  end subroutine write_grid_mapping

end module dyn_comp
