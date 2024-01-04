module bndry_mod
  use shr_kind_mod,   only: r8=>shr_kind_r8, i8=>shr_kind_i8
  use parallel_mod,   only: HME_BNDRY_A2A, HME_BNDRY_A2AO
  use thread_mod,     only: omp_in_parallel, omp_get_thread_num
  use gbarrier_mod,   only: gbarrier
  use cam_abortutils, only: endrun
  use cam_logfile,     only: iulog


  implicit none
  private

  interface bndry_exchange
     module procedure bndry_exchange_threaded
     module procedure bndry_exchange_nonthreaded
     module procedure long_bndry_exchange_nonth
  end interface
  public :: bndry_exchange


  interface bndry_exchange_start
     module procedure bndry_exchange_threaded_start
     module procedure bndry_exchange_nonthreaded_start
  end interface
  public :: bndry_exchange_start

  interface bndry_exchange_finish
     module procedure bndry_exchange_threaded_finish
     module procedure bndry_exchange_nonthreaded_finish
  end interface
  public :: bndry_exchange_finish



contains

  subroutine bndry_exchange_a2a(par,nthreads,ithr,buffer,location)
    use edgetype_mod,  only: Edgebuffer_t
    use schedtype_mod, only: schedule_t, cycle_t, schedule
    use thread_mod,    only: omp_in_parallel, omp_get_thread_num
    use perf_mod,      only: t_startf, t_stopf
    use spmd_utils,    only: mpi_real8, mpi_success, mpi_status_size
    use parallel_mod,  only: parallel_t

    type (parallel_t)            :: par
    integer, intent(in)          :: nthreads
    integer                      :: ithr ! The OpenMP thread ID
    type (EdgeBuffer_t)          :: buffer
    character(len=*),  optional  :: location

    type (Schedule_t), pointer   :: pSchedule
    type (Cycle_t),    pointer   :: pCycle
    integer                      :: icycle,ierr
    integer                      :: length
    integer                      :: iptr,source,nlyr
    integer                      :: nSendCycles,nRecvCycles
    integer                      :: errorcode,errorlen
    character*(80)               :: errorstring
    character(len=*),  parameter :: subname = 'bndry_exchange_a2a'
    character(len=80)            :: locstring


    integer                      :: i,j
    integer                      :: request
    integer                      :: lstatus(MPI_status_size)

! Neighborhood collectives are only in MPI3 and up
#ifdef SPMD
#ifdef _MPI3

   if(ithr == 0) then

      call MPI_Ineighbor_Alltoallv(buffer%buf,buffer%scountsFull,buffer%sdisplsFull,Mpi_real8, &
                     buffer%receive,buffer%rcountsFull,buffer%rdisplsFull,Mpi_real8,par%commGraphFull,request,ierr)
      if(ierr .ne. MPI_SUCCESS) then
        errorcode=ierr
        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
        write(iulog,*) subname,': Error after call to MPI_Ineighbor_alltoallv: ',errorstring
      endif

      if(present(location)) then
        locstring = TRIM(subname) // ': ' // TRIM(location)
      else
        locstring = TRIM(subname)
      endif
      ! location 1 for copyBuffer
      call copyBuffer(nthreads,ithr,buffer,locstring)

      call MPI_wait(request,lstatus,ierr)
   else

      if(present(location)) then
        locstring = TRIM(subname) // ': ' // TRIM(location)
      else
        locstring = TRIM(subname)
      endif
      call copyBuffer(nthreads,ithr,buffer,locstring)

   endif
#else
    call endrun('bndry_exchange_a2a requires MPI-3 feature support')
#endif
#endif

  end subroutine bndry_exchange_a2a

  subroutine copyBuffer(nthreads,ithr,buffer,location)
    use edgetype_mod, only : Edgebuffer_t
    integer :: nthreads
    integer :: ithr
    type (EdgeBuffer_t)          :: buffer
    character(len=80)            :: location
    logical ::  ompThreadMissmatch
    integer lenMovePtr, iptr,length,i,j

    ompThreadMissmatch = .false.
    lenMovePtr = size(buffer%moveptr)
    if ( lenMOveptr .ne. nthreads) then
      ompthreadMissmatch = .true.
      write(*,30) TRIM(location), lenMoveptr, nthreads
    endif

    if (.not. ompthreadMissmatch) then
      iptr   = buffer%moveptr(ithr+1)
      length = buffer%moveLength(ithr+1)
      if(length>0) then
        do i=0,length-1
           buffer%receive(iptr+i) = buffer%buf(iptr+i)
        enddo
      endif
    else if(ompthreadMissmatch .and. ithr == 0) then
       do j=1,lenMovePtr
          iptr   = buffer%moveptr(j)
          length = buffer%moveLength(j)
          if(length>0) then
             do i=0,length-1
                buffer%receive(iptr+i) = buffer%buf(iptr+i)
             enddo
          endif
       enddo
    endif
30  format(a,'Potential perf issue: ',a,'LenMoveptr,nthreads: ',2(i3))
  end subroutine copyBuffer

  subroutine bndry_exchange_a2ao(par,nthreads,ithr,buffer,location)
    use edgetype_mod, only : Edgebuffer_t
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use thread_mod, only : omp_in_parallel, omp_get_thread_num
    use perf_mod, only : t_startf, t_stopf
    use spmd_utils,   only: mpi_real8, mpi_success, mpi_status_size
    use parallel_mod, only: parallel_t
    use perf_mod, only : t_startf, t_stopf

    type (parallel_t)                 :: par
    integer, intent(in)               :: nthreads
    integer                           :: ithr  ! The OpenMP thread ID
    type (EdgeBuffer_t)               :: buffer
    character(len=*), optional        :: location

    integer                           :: ierr
    integer                           :: errorcode,errorlen
    character(len=80)                 :: errorstring
    character(len=*),       parameter :: subname = 'bndry_exchange_a2ao'
    character(len=80)                 :: locstring

    integer :: requestIntra,requestInter
    integer :: lstatus(MPI_status_size)

! Neighborhood collectives are only in MPI3 and up
#ifdef SPMD
#ifdef _MPI3

   if(ithr == 0) then

      call t_startf('bndry_a2ao')
      ! Start Inter-node communication
      call MPI_Ineighbor_Alltoallv(buffer%buf,buffer%scountsInter,buffer%sdisplsInter,MPI_real8, &
           buffer%receive,buffer%rcountsInter,buffer%rdisplsInter,MPI_real8,par%commGraphInter,requestInter,ierr)
      if(ierr .ne. MPI_SUCCESS) then
        errorcode=ierr
        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
        write(iulog,*) subname,': Error after call to MPI_Ineighbor_alltoallv: ',errorstring
      endif
      ! Start Intra-node communication
      call MPI_Ineighbor_Alltoallv(buffer%buf,buffer%scountsIntra,buffer%sdisplsIntra,MPI_real8, &
           buffer%receive,buffer%rcountsIntra,buffer%rdisplsIntra,MPI_real8,par%commGraphIntra,requestIntra,ierr)
      if(ierr .ne. MPI_SUCCESS) then
        errorcode=ierr
        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
        write(iulog,*) subname,': Error after call to MPI_Ineighbor_alltoallv: ',errorstring
      endif

      if(present(location)) then
        locstring = TRIM(subname) // ': ' // TRIM(location)
      else
        locstring = TRIM(subname)
      endif
      ! Finish the Intra-node communication
      call MPI_wait(requestIntra,lstatus,ierr)

      ! location 3 for copyBuffer
      call copyBuffer(nthreads,ithr,buffer,locstring)

      ! Finish the Inter-node communication
      call MPI_wait(requestInter,lstatus,ierr)
      call t_stopf('bndry_a2ao')

   else

      if(present(location)) then
        locstring = TRIM(subname) // ': ' // TRIM(location)
      else
        locstring = TRIM(subname)
      endif
      !Copy buffer for ithr!=0
      call copyBuffer(nthreads,ithr,buffer,locstring)

   endif
#else
    call endrun('bndry_exchange_a2ao requires MPI-3 feature support')
#endif
#endif

  end subroutine bndry_exchange_a2ao

  subroutine bndry_exchange_p2p(par,nthreads,ithr,buffer,location)
    use edgetype_mod,  only: Edgebuffer_t
    use schedtype_mod, only: schedule_t, cycle_t, schedule
    use thread_mod,    only: omp_in_parallel, omp_get_thread_num
    use spmd_utils,    only: mpi_real8, mpi_success
    use parallel_mod,  only: parallel_t
    use perf_mod,      only: t_startf, t_stopf

    type (parallel_t)                 :: par
    integer, intent(in)               :: nthreads
    integer                           :: ithr
    type (EdgeBuffer_t)               :: buffer
    character(len=*), optional        :: location

    type (Schedule_t),pointer         :: pSchedule
    type (Cycle_t),pointer            :: pCycle
    integer                           :: dest,length,tag
    integer                           :: icycle,ierr
    integer                           :: iptr,source,nlyr
    integer                           :: nSendCycles,nRecvCycles
    integer                           :: errorcode,errorlen
    character*(80)                    :: errorstring
    character(len=*), parameter       :: subname = 'bndry_exchange_p2p'
    character(len=80)                 :: locstring
    logical, parameter :: Debug=.FALSE.

    integer                           :: i,j
    logical :: ompthreadMissmatch
    integer :: lenMovePtr

    pSchedule => Schedule(1)
    nlyr = buffer%nlyr
    ompthreadMissmatch = .FALSE.

    lenMovePtr = size(buffer%moveptr)

  if(ithr == 0) then
    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles


    !==================================================
    !  Fire off the sends
    !==================================================

    do icycle=1,nSendCycles
       pCycle      => pSchedule%SendCycle(icycle)
       dest            = pCycle%dest - 1
       length          = buffer%scountsFull(icycle)
       tag             = buffer%tag
       iptr            = buffer%sdisplsFull(icycle) + 1
       if(Debug) write(iulog,*) subname,': MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
       call MPI_Isend(buffer%buf(iptr),length,Mpi_real8,dest,tag,par%comm,buffer%Srequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          write(iulog,*) subname,': Error after call to MPI_Isend: ',errorstring
       endif
    end do    ! icycle

    !==================================================
    !  Post the Receives
    !==================================================
    do icycle=1,nRecvCycles
       pCycle         => pSchedule%RecvCycle(icycle)
       source          = pCycle%source - 1
       length          = buffer%rcountsFull(icycle)
       tag             = buffer%tag
       iptr            = buffer%rdisplsFull(icycle) + 1
       if(Debug) write(iulog,*) subname,': MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
       call MPI_Irecv(buffer%receive(iptr),length,Mpi_real8, &
            source,tag,par%comm,buffer%Rrequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          write(iulog,*) subname,': Error after call to MPI_Irecv: ',errorstring
       endif
    end do    ! icycle
    if(present(location)) then
      locstring = TRIM(subname) // ': ' // TRIM(location)
    else
      locstring = TRIM(subname)
    endif
    call copyBuffer(nthreads,ithr,buffer,locstring)
    if (nSendCycles>0) call MPI_Waitall(nSendCycles,buffer%Srequest,buffer%status,ierr)
    if (nRecvCycles>0) call MPI_Waitall(nRecvCycles,buffer%Rrequest,buffer%status,ierr)
  else
    if(present(location)) then
      locstring = TRIM(subname) // ': ' // TRIM(location)
    else
      locstring = TRIM(subname)
    endif
    call copyBuffer(nthreads,ithr,buffer,locstring)
  endif

  end subroutine bndry_exchange_p2p

  subroutine bndry_exchange_p2p_start(par,nthreads,ithr,buffer,location)

    use edgetype_mod,  only: Edgebuffer_t
    use schedtype_mod, only: schedule_t, cycle_t, schedule
    use thread_mod,    only: omp_in_parallel, omp_get_thread_num
    use spmd_utils,    only: mpi_real8, mpi_success
    use parallel_mod,  only: parallel_t

    type (parallel_t)                 :: par
    integer, intent(in)               :: nthreads
    integer                           :: ithr
    type (EdgeBuffer_t)               :: buffer
    character (len=*), optional :: location

    type (Schedule_t),pointer         :: pSchedule
    type (Cycle_t),pointer            :: pCycle
    integer                           :: dest,length,tag
    integer                           :: icycle,ierr
    integer                           :: iptr,source,nlyr
    integer                           :: nSendCycles,nRecvCycles
    integer                           :: errorcode,errorlen
    character*(80)                    :: errorstring
    character(len=*),       parameter :: subname = 'bndry_exchange_p2p_start'
    logical,                parameter :: Debug=.FALSE.

    integer                           :: i,j, lenMovePtr
    logical :: ompthreadMissmatch

    pSchedule => Schedule(1)
    nlyr = buffer%nlyr
    ompthreadMissmatch = .FALSE.

    lenMovePtr = size(buffer%moveptr)

  if(ithr == 0) then
    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles

    !==================================================
    !  Fire off the sends
    !==================================================

    do icycle=1,nSendCycles
       pCycle      => pSchedule%SendCycle(icycle)
       dest            = pCycle%dest - 1
       length          = buffer%scountsFull(icycle)
       tag             = buffer%tag
       iptr            = buffer%sdisplsFull(icycle) + 1
       if(Debug) write(iulog,*) subname,': MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
       call MPI_Isend(buffer%buf(iptr),length,Mpi_real8,dest,tag,par%comm,buffer%Srequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          write(iulog,*) subname,': Error after call to MPI_Isend: ',errorstring
       endif
    end do    ! icycle

    !==================================================
    !  Post the Receives
    !==================================================
    do icycle=1,nRecvCycles
       pCycle         => pSchedule%RecvCycle(icycle)
       source          = pCycle%source - 1
       length          = buffer%rcountsFull(icycle)
       tag             = buffer%tag
       iptr            = buffer%rdisplsFull(icycle) + 1
       if(Debug) write(iulog,*) subname,': MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
       call MPI_Irecv(buffer%receive(iptr),length,Mpi_real8, &
            source,tag,par%comm,buffer%Rrequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          write(iulog,*) subname,': Error after call to MPI_Irecv: ',errorstring
       endif
    end do    ! icycle
  endif

  end subroutine bndry_exchange_p2p_start

  subroutine bndry_exchange_p2p_finish(par,nthreads,ithr,buffer,location)
    use edgetype_mod,  only: Edgebuffer_t
    use schedtype_mod, only: schedule_t, cycle_t, schedule
    use thread_mod,    only: omp_in_parallel, omp_get_thread_num
    use parallel_mod,  only: parallel_t
    use perf_mod,      only: t_startf, t_stopf


    type (parallel_t)            :: par
    integer, intent(in)          :: nthreads
    integer                      :: ithr
    type (EdgeBuffer_t)          :: buffer
    character(len=*),  optional  :: location

    type (Schedule_t), pointer   :: pSchedule
    type (Cycle_t),    pointer   :: pCycle
    integer                      :: dest,length,tag
    integer                      :: icycle,ierr
    integer                      :: iptr,source,nlyr
    integer                      :: nSendCycles,nRecvCycles
    integer                      :: errorcode,errorlen
    character*(80)               :: errorstring
    character(len=*),  parameter :: subname = 'bndry_exchange_p2p_finish'
    character(len=80)            :: locstring
    logical,           parameter :: Debug=.FALSE.

    integer                      :: i,j
    logical                      :: ompthreadMissmatch
    integer                      :: lenMovePtr


  pSchedule => Schedule(1)
  nlyr = buffer%nlyr
  if(present(location)) then
    locstring = TRIM(subname) // ': ' // TRIM(location)
  else
    locstring = TRIM(subname)
  endif
  call copyBuffer(nthreads,ithr,buffer,locstring)

  if(ithr == 0) then

    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles

    if (nSendCycles>0) call MPI_Waitall(nSendCycles,buffer%Srequest,buffer%status,ierr)
    if (nRecvCycles>0) call MPI_Waitall(nRecvCycles,buffer%Rrequest,buffer%status,ierr)

  endif

  end subroutine bndry_exchange_p2p_finish

  subroutine long_bndry_exchange_nonth(par,buffer)
    use edgetype_mod,  only: LongEdgebuffer_t
    use schedtype_mod, only: schedule_t, cycle_t, schedule
    use thread_mod,    only: omp_in_parallel
    use parallel_mod,  only: parallel_t, status, srequest, rrequest
    use spmd_utils,    only: mpi_integer, mpi_success

    type (parallel_t)            :: par
    type (LongEdgeBuffer_t)      :: buffer

    type (Schedule_t), pointer   :: pSchedule
    type (Cycle_t),    pointer   :: pCycle
    integer                      :: dest,length,tag
    integer                      :: icycle,ierr
    integer                      :: iptr,source,nlyr
    integer                      :: nSendCycles,nRecvCycles
    integer                      :: errorcode,errorlen
    character*(80)               :: errorstring
    character(len=*),  parameter :: subname = 'long_bndry_exchange_nonth'
    logical,           parameter :: Debug=.FALSE.

    integer                      :: i

#ifdef SPMD
    if(omp_in_parallel()) then
       print *,subname,': Warning you are calling a non-thread safe'
       print *,'         routine inside a threaded region....     '
       print *,'                Results are not predictable!!            '
    endif


    ! Setup the pointer to proper Schedule
    pSchedule => Schedule(1)
    nlyr = buffer%nlyr

    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles


    !==================================================
    !  Fire off the sends
    !==================================================

    do icycle=1,nSendCycles
       pCycle      => pSchedule%SendCycle(icycle)
       dest            = pCycle%dest - 1
       length      = nlyr * pCycle%lengthP
       tag             = pCycle%tag
       iptr            = pCycle%ptrP

       call MPI_Isend(buffer%buf(1,iptr),length,Mpi_integer,dest,tag,par%comm,Srequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          write(iulog,*) subname,': Error after call to MPI_Isend: ',errorstring
       endif
    end do    ! icycle

    !==================================================
    !  Post the Receives
    !==================================================
    do icycle=1,nRecvCycles
       pCycle         => pSchedule%RecvCycle(icycle)
       source          = pCycle%source - 1
       length      = nlyr * pCycle%lengthP
       tag             = pCycle%tag
       iptr            = pCycle%ptrP

       call MPI_Irecv(buffer%receive(1,iptr),length,Mpi_integer, &
            source,tag,par%comm,Rrequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          write(iulog,*) subname,': Error after call to MPI_Irecv: ',errorstring
       endif
    end do    ! icycle


    !==================================================
    !  Wait for all the receives to complete
    !==================================================

    if (nSendCycles>0) call MPI_Waitall(nSendCycles,Srequest,status,ierr)
    if (nRecvCycles>0) call MPI_Waitall(nRecvCycles,Rrequest,status,ierr)
    do icycle=1,nRecvCycles
       pCycle         => pSchedule%RecvCycle(icycle)
       length             = pCycle%lengthP
       iptr            = pCycle%ptrP
       do i=0,length-1
          buffer%buf(1:nlyr,iptr+i) = buffer%receive(1:nlyr,iptr+i)
       enddo
    end do   ! icycle

#endif

  end subroutine long_bndry_exchange_nonth
  !********************************************************************************
  !
  !********************************************************************************


 subroutine bndry_exchange_threaded(hybrid,buffer,location)
    use hybrid_mod, only : hybrid_t
    use edgetype_mod, only : Edgebuffer_t
    use perf_mod, only: t_startf, t_stopf, t_adj_detailf
    implicit none

    type (hybrid_t)             :: hybrid
    type (EdgeBuffer_t)         :: buffer
    character(len=*), optional  :: location

    character(len=*), parameter :: subname = 'bndry_exchange_threaded'

    call gbarrier(buffer%gbarrier, hybrid%ithr)
    if(buffer%bndry_type == HME_BNDRY_A2A) then
       call bndry_exchange_a2a(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_A2AO) then
       call bndry_exchange_a2ao(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    else
       call bndry_exchange_p2p(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    endif
    call gbarrier(buffer%gbarrier, hybrid%ithr)

 end subroutine bndry_exchange_threaded

 subroutine bndry_exchange_threaded_start(hybrid,buffer,location)
    use hybrid_mod, only : hybrid_t
    use edgetype_mod, only : Edgebuffer_t
    use perf_mod, only: t_startf, t_stopf, t_adj_detailf
    implicit none

    type (hybrid_t)             :: hybrid
    type (EdgeBuffer_t)         :: buffer
    character(len=*), optional  :: location

    character(len=*), parameter :: subname = 'bndry_exchange_threaded_start'

    call gbarrier(buffer%gbarrier, hybrid%ithr)
    call bndry_exchange_p2p_start(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)

 end subroutine bndry_exchange_threaded_start

 subroutine bndry_exchange_threaded_finish(hybrid,buffer,location)
    use hybrid_mod, only : hybrid_t
    use edgetype_mod, only : Edgebuffer_t
    use perf_mod, only: t_startf, t_stopf, t_adj_detailf
    implicit none

    type (hybrid_t)             :: hybrid
    type (EdgeBuffer_t)         :: buffer
    character(len=*), optional  :: location

    character(len=*), parameter :: subname = 'bndry_exchange_threaded_finish'

    call bndry_exchange_p2p_finish(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    call gbarrier(buffer%gbarrier, hybrid%ithr)

 end subroutine bndry_exchange_threaded_finish

 subroutine bndry_exchange_nonthreaded(par,buffer,location)
    use parallel_mod, only : parallel_t
    use edgetype_mod, only : Edgebuffer_t
    implicit none

    type (parallel_t)           :: par
    type (EdgeBuffer_t)         :: buffer
    character(len=*), optional  :: location

    integer                     :: ithr
    integer                     :: nthreads
!    integer                     :: localsense
    character(len=*), parameter :: subname = 'bndry_exchange_nonthreaded'

    ithr=0
    nthreads = 1
    if(buffer%bndry_type == HME_BNDRY_A2A) then
       call bndry_exchange_a2a(par,nthreads,ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_A2AO) then
       call bndry_exchange_a2ao(par,nthreads,ithr,buffer,location)
    else
       call bndry_exchange_p2p(par,nthreads,ithr,buffer,location)
    endif

  end subroutine bndry_exchange_nonthreaded

 subroutine bndry_exchange_nonthreaded_start(par,buffer,location)
    use parallel_mod, only : parallel_t
    use edgetype_mod, only : Edgebuffer_t
    implicit none

    type (parallel_t)           :: par
    type (EdgeBuffer_t)         :: buffer
    character (len=*), optional :: location

    integer                     :: ithr
    integer                     :: nthreads
    character(len=*), parameter :: subname = 'bndry_exchange_nonthreaded_start'

    ithr=0
    nthreads=1
    call bndry_exchange_p2p_start(par,nthreads,ithr,buffer,location)

 end subroutine bndry_exchange_nonthreaded_start

 subroutine bndry_exchange_nonthreaded_finish(par,buffer,location)
    use parallel_mod, only : parallel_t
    use edgetype_mod, only : Edgebuffer_t
    implicit none

    type (parallel_t)                 :: par
    integer                           :: ithr
    type (EdgeBuffer_t)               :: buffer
    character (len=*), optional :: location
    integer :: nthreads

    character(len=*), parameter :: subname = 'bndry_exchange_nonthreaded_finish'

    ithr=0
    nthreads=1
    call bndry_exchange_p2p_finish(par,nthreads,ithr,buffer,location)

  end subroutine bndry_exchange_nonthreaded_finish

end module bndry_mod
