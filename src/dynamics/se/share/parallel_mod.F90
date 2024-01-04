module parallel_mod
  ! ---------------------------
  use shr_kind_mod,   only: r8=>shr_kind_r8
  ! ---------------------------
  use dimensions_mod, only : nmpi_per_node, nlev, qsize_d
  ! ---------------------------
  use spmd_utils,     only: MPI_STATUS_SIZE, MPI_MAX_ERROR_STRING, MPI_TAG_UB
  use cam_logfile,  only: iulog

  implicit none
  private

  integer,  public, parameter   :: ORDERED         = 1
  integer,  public, parameter   :: FAST            = 2
  integer,  public, parameter   :: BNDRY_TAG_BASE  = 0
  integer,  public, parameter   :: THREAD_TAG_BITS = 9
  integer,  public, parameter   :: MAX_ACTIVE_MSG = (MPI_TAG_UB/2**THREAD_TAG_BITS) - 1
  integer,  public, parameter   :: HME_status_size = MPI_STATUS_SIZE

  integer,  public, parameter   :: HME_BNDRY_P2P   = 1
  integer,  public, parameter   :: HME_BNDRY_MASHM = 2
  integer,  public, parameter   :: HME_BNDRY_A2A   = 3
  integer,  public, parameter   :: HME_BNDRY_A2AO  = 4

  integer,  public, parameter   :: nrepro_vars = MAX(10, nlev*qsize_d)
  integer,  public, parameter   :: HME_MPATTERN_P  = 101
  integer,  public, parameter   :: HME_MPATTERN_S  = 102
  integer,  public, parameter   :: HME_MPATTERN_G  = 103

  integer,  public              :: MaxNumberFrames
  integer,  public              :: numframes
  integer,  public              :: useframes
  logical,  public              :: PartitionForNodes
  logical,  public              :: PartitionForFrames
  integer,  public              :: MPIreal_t,MPIinteger_t,MPIChar_t,MPILogical_t
  integer,  public              :: MPIaddr_kind
  integer,  public              :: iam

  ! Namelist-selectable type of boundary comms (AUTO,P2P,A2A,MASHM)
  integer,  public              :: boundaryCommMethod

  integer,  public, allocatable :: status(:,:)
  integer,  public, allocatable :: Rrequest(:)
  integer,  public, allocatable :: Srequest(:)

  real(r8), public, allocatable :: FrameWeight(:)
  integer,  public, allocatable :: FrameIndex(:)
  integer,  public, allocatable :: FrameCount(:)
  integer,  public              :: nComPoints
  integer,  public              :: nPackPoints

  real(r8), public, allocatable :: global_shared_buf(:,:)
  real(r8), public              :: global_shared_sum(nrepro_vars)

  ! ==================================================
  ! Define type parallel_t for distributed memory info
  ! ==================================================

  integer, parameter :: ncomponents=1

  type, public :: parallel_t
    integer :: rank                       ! local rank
    integer :: root                       ! local root
    integer :: nprocs                     ! number of processes in group
    integer :: comm                       ! communicator
    integer :: intercomm                  ! inter communicator list
    integer :: intracomm                  ! intra node communicator
    integer :: intracommsize              ! number of MPI ranks in intra communicator
    integer :: intracommrank              ! rank in intra communicator
    integer :: commGraphFull
    integer :: commGraphInter
    integer :: commGraphIntra
    integer :: groupGraphFull
    logical :: masterproc
  end type

  type (parallel_t), public :: par ! info for distributed memory programming

  ! ===================================================
  ! Module Interfaces
  ! ===================================================

  public :: initmpi
  public :: syncmp
  public :: copy_par

  interface assignment ( = )
    module procedure copy_par
  end interface

CONTAINS

! ================================================
!   copy_par: copy constructor for parallel_t type
!
!
!   Overload assignment operator for parallel_t
! ================================================

  subroutine copy_par(par2,par1)
    type(parallel_t), intent(out) :: par2
    type(parallel_t), intent(in)  :: par1

    par2%rank       = par1%rank
    par2%root       = par1%root
    par2%nprocs     = par1%nprocs
    par2%comm       = par1%comm
    par2%intercomm  = par1%intercomm
    par2%commGraphFull   = par1%commGraphFull
    par2%commGraphInter  = par1%commGraphInter
    par2%commGraphIntra  = par1%commGraphIntra
    par2%groupGraphFull  = par1%groupGraphFull
    par2%masterproc = par1%masterproc

  end subroutine copy_par

! ================================================
!  initmpi:
!  Initializes the parallel (message passing)
!  environment, returns a parallel_t structure..
! ================================================

  function initmpi(npes_in) result(par)
    use spmd_utils, only : mpicom
    integer, intent(in), optional ::  npes_in
    type (parallel_t) par

#ifdef _MPI

#include <mpif.h>
#ifdef _AIX
    integer                              :: ii         
    character(len=2)                                    :: cfn
#endif

    integer              :: ierr,tmp
    integer              :: FrameNumber
    logical :: running   ! state of MPI at beginning of initmpi call
    character(len=MPI_MAX_PROCESSOR_NAME)               :: my_name
    character(len=MPI_MAX_PROCESSOR_NAME), allocatable  :: the_names(:)

    integer, allocatable :: tarray(:)
    integer              :: namelen, i
    integer              :: color
    integer :: iam_cam, npes_cam
    integer :: npes_homme
    !================================================
    !     Basic MPI initialization
    ! ================================================

    call MPI_initialized(running, ierr)

    if (.not.running) then
       call MPI_init(ierr)
    end if

    par%root          = 0
    par%intercomm = 0
    
#ifdef CAM
    if(present(npes_in)) then
       npes_homme=npes_in
    else
       npes_homme=npes_cam
    end if
    call MPI_comm_size(mpicom,npes_cam,ierr)
    call MPI_comm_rank(mpicom,iam_cam,ierr)
    color = iam_cam/npes_homme
    call mpi_comm_split(mpicom, color, iam_cam, par%comm, ierr)
#else
    par%comm     = MPI_COMM_WORLD
#endif
    call MPI_comm_rank(par%comm,par%rank,ierr)
    call MPI_comm_size(par%comm,par%nprocs,ierr)

    par%masterproc = .FALSE.
    if(par%rank .eq. par%root) par%masterproc = .TRUE.
    if (par%masterproc) write(iulog,*)'number of MPI processes: ',par%nprocs
           
    if (MPI_DOUBLE_PRECISION==20 .and. MPI_REAL8==18) then
       ! LAM MPI defined MPI_REAL8 differently from MPI_DOUBLE_PRECISION
       ! and LAM MPI's allreduce does not accept on MPI_REAL8
       MPIreal_t    = MPI_DOUBLE_PRECISION
    else
       MPIreal_t    = MPI_REAL8
    endif
    MPIinteger_t = MPI_INTEGER
    MPIchar_t    = MPI_CHARACTER 
    MPILogical_t = MPI_LOGICAL
    MPIaddr_kind = MPI_ADDRESS_KIND

    ! ================================================ 
    !  Determine where this MPI process is running 
    !   then use this information to determined the 
    !   number of MPI processes per node    
    ! ================================================ 

    my_name(:) = ''
    call MPI_Get_Processor_Name(my_name,namelen,ierr)

    allocate(the_names(par%nprocs))
    do i=1,par%nprocs
       the_names(i)(:) =  ''
    enddo
    ! ================================================ 
    !   Collect all the machine names 
    ! ================================================ 
    call MPI_Allgather(my_name,MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER, &
           the_names,MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,par%comm,ierr)

    ! ======================================================================
    !   Calculate how many other MPI processes are on my node 
    ! ======================================================================
    nmpi_per_node = 0
    do i=1,par%nprocs
      if( TRIM(ADJUSTL(my_name)) .eq. TRIM(ADJUSTL(the_names(i)))   ) then 
        nmpi_per_node = nmpi_per_node + 1
      endif
    enddo

    ! =======================================================================
    !  Verify that everybody agrees on this number otherwise do not do 
    !  the multi-level partitioning
    ! =======================================================================
    call MPI_Allreduce(nmpi_per_node,tmp,1,MPIinteger_t,MPI_BAND,par%comm,ierr)
    if(tmp .ne. nmpi_per_node) then 
      if (par%masterproc) write(iulog,*)'initmpi:  disagrement accross nodes for nmpi_per_node'
      nmpi_per_node = 1
      PartitionForNodes=.FALSE.
    else
      PartitionForNodes=.TRUE.
    endif


#ifdef _AIX
    PartitionForFrames=.FALSE.
    if((my_name(1:2) .eq. 'bv') .or.  &   ! Bluvista
       (my_name(1:2) .eq. 'bl') .or.  &    ! Blueice
       (my_name(1:2) .eq. 'bh') .or.  &    ! Blue Horizon
       (my_name(1:2) .eq. 'bs') .or.  &    ! BlueSky
       (my_name(1:2) .eq. 's0')       &    ! Seaborg
      ) then

       ! ================================================================================
       ! Note: the frame based optimization is only supported on blackforest or babyblue
       ! ================================================================================
       cfn = my_name(3:4)
       read(cfn,'(i2)') FrameNumber

       ! ======================================================
       ! Make sure that the system does not have too may frames 
       ! ======================================================
       call MPI_Allreduce(FrameNumber,MaxNumberFrames,1,MPIinteger_t,MPI_MAX,par%comm,ierr)
       MaxNumberFrames=MaxNumberFrames+1
       allocate(FrameCount(MaxNumberFrames))
       allocate(tarray(MaxNumberFrames))

       call MPI_Allreduce(useframes,tmp,1,MPIinteger_t,MPI_BAND,par%comm,ierr)
       if(tmp .ne. useframes) then 
          if (par%masterpoc) write(iulog,*) "initmpi:  disagreement accross nodes for useframes"
          PartitionForFrames=.FALSE.
       endif

       if(PartitionForFrames) then 
         tarray(:) = 0
         tarray(FrameNumber+1) = 1
         call MPI_Allreduce(tarray,FrameCount,MaxNumberFrames,MPIinteger_t,MPI_SUM,par%comm,ierr)
         if(par%masterproc)  write(iulog,*)'initmpi: FrameCount : ',FrameCount
           numFrames = COUNT(FrameCount .ne. 0)
           allocate(FrameWeight(numFrames))
           allocate(FrameIndex(numFrames))

         ii=1
         do i=1,MaxNumberFrames
           if(FrameCount(i) .ne. 0) then 
             FrameWeight(ii) = real(FrameCount(i),kind=4)/real(par%nprocs,kind=4)
             FrameIndex(ii)  = i
             ii=ii+1
           endif
         enddo
       endif

      FrameCount(:)=FrameCount(:)/nmpi_per_node
      ! ==========================================
      ! We are not running on more than one frame 
      ! ==========================================
      if(numFrames .eq. 1)  PartitionForFrames=.FALSE.

    endif

    write(iulog,*) 'initmpi: mpi task ',par%rank,': ',nmpi_per_node,' task(s) on node ',my_name(1:namelen),  &
                   'on frame # ',FrameNumber 
#endif
    if(PartitionForFrames) then 
      if(par%masterproc) write(iulog,*)'initmpi: FrameWeight: ',FrameWeight
      if(par%masterproc) write(iulog,*)'initmpi: FrameIndex: ',FrameIndex
    endif

    deallocate(the_names)
 
#else
    par%root          =  0 
    par%rank          =  0
    par%nprocs        =  1
    par%comm          = -1
    par%intercomm     = -1
    par%masterproc    = .TRUE.
    nmpi_per_node     =  2
    PartitionForNodes = .TRUE.
#endif
    !===================================================
    !  Kind of lame but set this variable to be 1 based 
    !===================================================
    iam = par%rank+1

  end function initmpi

  ! =====================================
  ! syncmp:
  !
  ! sychronize message passing domains
  !
  ! =====================================
  subroutine syncmp(par)
    use cam_abortutils, only: endrun
    use spmd_utils,     only: MPI_MAX_ERROR_STRING, MPI_ERROR

    type (parallel_t), intent(in)       :: par

    integer                             :: errorcode, errorlen, ierr
    character(len=MPI_MAX_ERROR_STRING) :: errorstring

    call MPI_barrier(par%comm, ierr)

    if(ierr == MPI_ERROR) then
      errorcode = ierr
      call MPI_Error_String(errorcode, errorstring, errorlen, ierr)
      call endrun(errorstring)
    end if
  end subroutine syncmp

end module parallel_mod
