#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module bndry_mod
  use parallel_mod, only : abortmp,iam, HME_BNDRY_P2P, HME_BNDRY_MASHM, HME_BNDRY_A2A, HME_BNDRY_A2AO, &
                           HME_BNDRY_GET1, HME_BNDRY_GET2, HME_BNDRY_PUT1, HME_BNDRY_PUT2
  use thread_mod, only : omp_in_parallel, omp_get_thread_num
  use gbarrier_mod, only: gbarrier
  use perf_mod, only : t_startf, t_stopf ! EXTERNAL
  implicit none
  private
#ifdef _MPI
#include <mpif.h>
#endif

  public :: bndry_exchangeV
  public :: bndry_exchangeS
  public :: bndry_exchangeS_start
  public :: bndry_exchangeS_finish



  interface bndry_exchangeV
     module procedure bndry_exchangeV_threaded
     module procedure bndry_exchangeV_nonthreaded
     module procedure long_bndry_exchangeV_nonth
  end interface

  interface bndry_exchangeS
     module procedure bndry_exchangeS_threaded 
     module procedure bndry_exchangeS_nonthreaded
  end interface

  interface bndry_exchangeS_finish
     module procedure bndry_exchangeS_threaded_finish
     module procedure bndry_exchangeS_nonthreaded_finish
  end interface

  interface bndry_exchangeS_start
     module procedure bndry_exchangeS_threaded_start
     module procedure bndry_exchangeS_nonthreaded_start
  end interface

contains 

  subroutine bndry_exchange_a2a(par,nthreads,ithr,buffer,location)
    use kinds, only : log_kind
    use edgetype_mod, only : Edgebuffer_t
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use thread_mod, only : omp_in_parallel, omp_get_thread_num
#ifdef _MPI
    use parallel_mod, only : parallel_t, abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success, mpi_version,hme_status_size
#else
    use parallel_mod, only : parallel_t, abortmp
#endif
    type (parallel_t)                 :: par
    integer, intent(in)               :: nthreads
    integer                           :: ithr  ! The OpenMP thread ID
    type (EdgeBuffer_t)               :: buffer
    character(len=*), optional        :: location
!
!    type (Schedule_t),pointer         :: pSchedule
!    type (Cycle_t),pointer            :: pCycle
!    integer                           :: icycle,ierr
!    integer                           :: length
!    integer                           :: iptr,source,nlyr
!    integer                           :: nSendCycles,nRecvCycles
    integer                           :: ierr
    integer                           :: errorcode,errorlen
    character(len=80)                 :: errorstring
    character(len=*),       parameter :: subname = 'bndry_exchange_a2a'
    character(len=80)                 :: locstring
!
!    logical(kind=log_kind), parameter :: Debug=.FALSE.
!    logical :: ompthreadMissmatch
!
!    integer                           :: i,j
!    integer :: lenMovePtr
    integer :: request
    integer :: lstatus(HME_status_size)

! Neighborhood collectives are only in MPI3 and up
#ifdef _MPI3
   if(ithr == 0) then 

      call MPI_Ineighbor_Alltoallv(buffer%buf,buffer%scountsFull,buffer%sdisplsFull,MPIreal_t, &
                     buffer%receive,buffer%rcountsFull,buffer%rdisplsFull,MPIreal_t,par%commGraphFull,request,ierr)
      if(ierr .ne. MPI_SUCCESS) then
        errorcode=ierr
        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
        print *,subname,': Error after call to MPI_neighbor_alltoallv: ',errorstring
      endif

      locstring = TRIM(subname) // ': ' // TRIM(location)
      ! location 1 for copyBuffer
      call copyBuffer(nthreads,ithr,buffer,locstring)

      call MPI_wait(request,lstatus,ierr)
   else

      locstring = TRIM(subname) // ': ' // TRIM(location)
      call copyBuffer(nthreads,ithr,buffer,locstring)

   endif
#else
    call abortmp('bndry_exchange_a2a requires MPI-3 feature support')
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
30  format(a,'Potential performance issue: ',a,'LenMoveptr,nthreads: ',2(i3))
  end subroutine copyBuffer

  subroutine bndry_exchange_put1(par,nthreads,ithr,buffer,location)

    use kinds, only : log_kind,long_kind
    use edgetype_mod, only : Edgebuffer_t
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use thread_mod, only : omp_in_parallel, omp_get_thread_num
#ifdef _MPI
    use parallel_mod, only : parallel_t, abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success, mpi_version, mpi_lock_shared, hme_status_size
#else
    use parallel_mod, only : parallel_t, abortmp
#endif
    type (parallel_t)                 :: par
    integer, intent(in)               :: nthreads
    integer                           :: ithr  ! The OpenMP thread ID
    type (EdgeBuffer_t)               :: buffer
    character(len=*), optional        :: location

    type (Schedule_t),pointer         :: pSchedule
    type (Cycle_t),pointer            :: pCycle
    integer                           :: icycle,ierr
    integer                           :: errorcode,errorlen
    character(len=80)                 :: errorstring
    character(len=*),       parameter :: subname = 'bndry_exchange_put1'
    character(len=80)                 :: locstring
!
    logical(kind=log_kind), parameter :: Debug=.TRUE.

    integer :: requestIntra,requestInter
    integer :: lstatus(HME_status_size)
    integer :: iptr,tag,length,source,dest,nRecvCycles, nSendCycles
    integer(kind=long_kind) :: disp

   pSchedule => schedule(1)
   if(ithr == 0) then 

    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles

    call MPI_Win_fence(0,buffer%win,ierr)
    do icycle=1,nSendCycles
       pCycle         => pSchedule%SendCycle(icycle)
       dest            = pCycle%dest - 1
       length          = buffer%scountsFull(icycle)
       tag             = buffer%tag
       iptr            = buffer%sdisplsFull(icycle)
       ! disp            = 0  !FIXME: This needs to point to the correct spot in the remove ranks memory
       disp            = INT(buffer%putDisplsFull(icycle),kind=long_kind)  !FIXME: This needs to point to the correct spot in the remove ranks memory
!       if(Debug .and. (source == 103 .or. iam == 104)) print *,'IAM: ',iam, subname,': MPI_Irecv: SRC:',source,'LENGTH:',length,'PTR: ',iptr, 'TAG: ',tag
!       call MPI_Win_lock(MPI_LOCK_SHARED,source,0,buffer%win,ierr)
!       if(Debug) print *,'IAM: ', iam, 'After MPI_win_lock'
       call MPI_Put(buffer%buf(iptr+1),length,MPIreal_t,dest,disp,length,MPIreal_t,buffer%win,ierr)
!       if(Debug) print *,'IAM: ', iam, 'After MPI_Get'
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,subname,': Error after call to MPI_Get: ',errorstring
       endif
    end do    ! icycle

    ! location 3 for copyBuffer
    call copyBuffer(nthreads,ithr,buffer,locstring)

    call MPI_Win_fence(0,buffer%win,ierr)

   else

    locstring = TRIM(subname) // ': ' // TRIM(location) 
    !Copy buffer for ithr!=0
    call copyBuffer(nthreads,ithr,buffer,locstring)

   endif

!   stop 'At end of bndry_exchange_put1'

  end subroutine bndry_exchange_put1

  subroutine bndry_exchange_put2(par,nthreads,ithr,buffer,location)
    use kinds, only : log_kind,long_kind
    use edgetype_mod, only : Edgebuffer_t
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use thread_mod, only : omp_in_parallel, omp_get_thread_num
#ifdef _MPI
    use parallel_mod, only : parallel_t, abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success, mpi_version, mpi_lock_shared, hme_status_size
#else
    use parallel_mod, only : parallel_t, abortmp
#endif
    type (parallel_t)                 :: par
    integer, intent(in)               :: nthreads
    integer                           :: ithr  ! The OpenMP thread ID
    type (EdgeBuffer_t)               :: buffer
    character(len=*), optional        :: location

    type (Schedule_t),pointer         :: pSchedule
    type (Cycle_t),pointer            :: pCycle
    integer                           :: icycle,ierr
    integer                           :: errorcode,errorlen
    character(len=80)                 :: errorstring
    character(len=*),       parameter :: subname = 'bndry_exchange_put2'
    character(len=80)                 :: locstring
!
    logical(kind=log_kind), parameter :: Debug=.TRUE.

    integer :: requestIntra,requestInter
    integer :: lstatus(HME_status_size)
    integer :: iptr,tag,length,source,dest, nRecvCycles, nSendCycles
    integer(kind=long_kind) :: disp

   pSchedule => schedule(1)
   if(ithr == 0) then 

    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles

    call MPI_Win_post(par%groupGraphFull,0,buffer%win,ierr)
    !==================================================
    ! Get the data from the remote process
    !==================================================
    call MPI_win_start(par%groupGraphFull,0,buffer%win,ierr)
    do icycle=1,nSendCycles
       pCycle         => pSchedule%SendCycle(icycle)
       dest            = pCycle%dest - 1
       length          = buffer%scountsFull(icycle)
       tag             = buffer%tag
       iptr            = buffer%sdisplsFull(icycle)
       disp            = INT(buffer%putDisplsFull(icycle),kind=long_kind)
       call MPI_Put(buffer%buf(iptr+1),length,MPIreal_t,dest,disp,length,MPIreal_t,buffer%win,ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,subname,': Error after call to MPI_Get: ',errorstring
       endif
    end do    ! icycle

    ! location 3 for copyBuffer
    call copyBuffer(nthreads,ithr,buffer,locstring)

    call MPI_Win_complete(buffer%win,ierr)
    call MPI_Win_wait(buffer%win,ierr)

   else

    locstring = TRIM(subname) // ': ' // TRIM(location) 
    !Copy buffer for ithr!=0
    call copyBuffer(nthreads,ithr,buffer,locstring)

   endif

  end subroutine bndry_exchange_put2

  subroutine bndry_exchange_get1(par,nthreads,ithr,buffer,location)
    use kinds, only : log_kind,long_kind
    use edgetype_mod, only : Edgebuffer_t
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use thread_mod, only : omp_in_parallel, omp_get_thread_num
#ifdef _MPI
    use parallel_mod, only : parallel_t, abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success, mpi_version, mpi_lock_shared, hme_status_size
#else
    use parallel_mod, only : parallel_t, abortmp
#endif
    type (parallel_t)                 :: par
    integer, intent(in)               :: nthreads
    integer                           :: ithr  ! The OpenMP thread ID
    type (EdgeBuffer_t)               :: buffer
    character(len=*), optional        :: location

    type (Schedule_t),pointer         :: pSchedule
    type (Cycle_t),pointer            :: pCycle
    integer                           :: icycle,ierr
    integer                           :: errorcode,errorlen
    character(len=80)                 :: errorstring
    character(len=*),       parameter :: subname = 'bndry_exchange_get1'
    character(len=80)                 :: locstring
!
    logical(kind=log_kind), parameter :: Debug=.TRUE.

    integer :: requestIntra,requestInter
    integer :: lstatus(HME_status_size)
    integer :: iptr,tag,length,source,nRecvCycles, nSendCycles
    integer(kind=long_kind) :: disp

   pSchedule => schedule(1)
   if(ithr == 0) then 

    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles

    call MPI_Win_fence(0,buffer%win,ierr)
!    call MPI_Win_post(par%commGraphfull,0,buffer%win,ierr)
    !==================================================
    ! Get the data from the remote process
    !==================================================
!    call MPI_win_start(par%commGraphfull,0,buffer%win,ierr)
    do icycle=1,nRecvCycles
       pCycle         => pSchedule%RecvCycle(icycle)
       source          = pCycle%source - 1
       length          = buffer%rcountsFull(icycle)
       tag             = buffer%tag
       iptr            = buffer%rdisplsFull(icycle)
       disp            = INT(buffer%getDisplsFull(icycle),kind=long_kind)
       call MPI_Get(buffer%receive(iptr+1),length,MPIreal_t,source,disp,length,MPIreal_t,buffer%win,ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,subname,': Error after call to MPI_Get: ',errorstring
       endif
    end do    ! icycle

    ! location 3 for copyBuffer
    call copyBuffer(nthreads,ithr,buffer,locstring)

    call MPI_Win_fence(0,buffer%win,ierr)

   else

    locstring = TRIM(subname) // ': ' // TRIM(location) 
    !Copy buffer for ithr!=0
    call copyBuffer(nthreads,ithr,buffer,locstring)

   endif

!   stop 'At end of bndry_exchange_get1'

  end subroutine bndry_exchange_get1

  subroutine bndry_exchange_get2(par,nthreads,ithr,buffer,location)
    use kinds, only : log_kind,long_kind
    use edgetype_mod, only : Edgebuffer_t
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use thread_mod, only : omp_in_parallel, omp_get_thread_num
#ifdef _MPI
    use parallel_mod, only : parallel_t, abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success, mpi_version, mpi_lock_shared, hme_status_size
#else
    use parallel_mod, only : parallel_t, abortmp
#endif
    type (parallel_t)                 :: par
    integer, intent(in)               :: nthreads
    integer                           :: ithr  ! The OpenMP thread ID
    type (EdgeBuffer_t)               :: buffer
    character(len=*), optional        :: location

    type (Schedule_t),pointer         :: pSchedule
    type (Cycle_t),pointer            :: pCycle
    integer                           :: icycle,ierr
    integer                           :: errorcode,errorlen
    character(len=80)                 :: errorstring
    character(len=*),       parameter :: subname = 'bndry_exchange_get2'
    character(len=80)                 :: locstring
!
    logical(kind=log_kind), parameter :: Debug=.TRUE.

    integer :: requestIntra,requestInter
    integer :: lstatus(HME_status_size)
    integer :: iptr,tag,length,source,nRecvCycles, nSendCycles
    integer(kind=long_kind) :: disp

   pSchedule => schedule(1)
   if(ithr == 0) then 

    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles

!    call MPI_Win_fence(0,buffer%win,ierr)
    call MPI_Win_post(par%groupGraphFull,0,buffer%win,ierr)
    !==================================================
    ! Get the data from the remote process
    !==================================================
    call MPI_win_start(par%groupGraphFull,0,buffer%win,ierr)
    do icycle=1,nRecvCycles
       pCycle         => pSchedule%RecvCycle(icycle)
       source          = pCycle%source - 1
       length          = buffer%rcountsFull(icycle)
       tag             = buffer%tag
       iptr            = buffer%rdisplsFull(icycle)
       ! disp            = 0  !FIXME: This needs to point to the correct spot in the remove ranks memory
       disp            = INT(buffer%getDisplsFull(icycle),kind=long_kind)  !FIXME: This needs to point to the correct spot in the remove ranks memory
!       if(Debug .and. (source == 103 .or. iam == 104)) print *,'IAM: ',iam, subname,': MPI_Irecv: SRC:',source,'LENGTH:',length,'PTR: ',iptr, 'TAG: ',tag
!       call MPI_Win_lock(MPI_LOCK_SHARED,source,0,buffer%win,ierr)
!       if(Debug) print *,'IAM: ', iam, 'After MPI_win_lock'
       call MPI_Get(buffer%receive(iptr+1),length,MPIreal_t,source,disp,length,MPIreal_t,buffer%win,ierr)
!       if(Debug) print *,'IAM: ', iam, 'After MPI_Get'
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,subname,': Error after call to MPI_Get: ',errorstring
       endif
!       call MPI_win_unlock(source,buffer%win,ierr)
!       if(Debug) print *,'IAM: ', iam, 'After MPI_win_unlock'
    end do    ! icycle

!    if(Debug) print *,'IAM: ', iam, 'Before copyBuffer'
    ! location 3 for copyBuffer
    call copyBuffer(nthreads,ithr,buffer,locstring)

    call MPI_Win_complete(buffer%win,ierr)
    call MPI_Win_wait(buffer%win,ierr)

   else

    locstring = TRIM(subname) // ': ' // TRIM(location) 
    !Copy buffer for ithr!=0
    call copyBuffer(nthreads,ithr,buffer,locstring)

   endif

!   stop 'At end of bndry_exchange_get2'

  end subroutine bndry_exchange_get2

  subroutine bndry_exchange_a2ao(par,nthreads,ithr,buffer,location)
    use kinds, only : log_kind
    use edgetype_mod, only : Edgebuffer_t
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use thread_mod, only : omp_in_parallel, omp_get_thread_num
#ifdef _MPI
    use parallel_mod, only : parallel_t, abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success, mpi_version,hme_status_size
#else
    use parallel_mod, only : parallel_t, abortmp
#endif
    type (parallel_t)                 :: par
    integer, intent(in)               :: nthreads
    integer                           :: ithr  ! The OpenMP thread ID
    type (EdgeBuffer_t)               :: buffer
    character(len=*), optional        :: location

!    type (Schedule_t),pointer         :: pSchedule
!    type (Cycle_t),pointer            :: pCycle
!    integer                           :: icycle,ierr
!    integer                           :: length
!    integer                           :: iptr,source,nlyr
!    integer                           :: nSendCycles,nRecvCycles
    integer                           :: ierr
    integer                           :: errorcode,errorlen
    character(len=80)                 :: errorstring
    character(len=*),       parameter :: subname = 'bndry_exchange_a2ao'
    character(len=80)                 :: locstring
!
!    logical(kind=log_kind), parameter :: Debug=.FALSE.
!    logical :: ompthreadMissmatch
!
!    integer                           :: i,j
!    integer :: lenMovePtr
    integer :: requestIntra,requestInter
    integer :: lstatus(HME_status_size)

! Neighborhood collectives are only in MPI3 and up
#ifdef _MPI3

   if(ithr == 0) then 

      ! Start Inter-node communication
      call MPI_Ineighbor_Alltoallv(buffer%buf,buffer%scountsInter,buffer%sdisplsInter,MPIreal_t, &
                     buffer%receive,buffer%rcountsInter,buffer%rdisplsInter,MPIreal_t,par%commGraphInter,requestInter,ierr)
      if(ierr .ne. MPI_SUCCESS) then
        errorcode=ierr
        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
        print *,subname,': Error after call to MPI_neighbor_alltoallv: ',errorstring
      endif
      locstring = TRIM(subname) // ': ' // TRIM(location) 

      ! Start Intra-node communication
      call MPI_Ineighbor_Alltoallv(buffer%buf,buffer%scountsIntra,buffer%sdisplsIntra,MPIreal_t, &
                     buffer%receive,buffer%rcountsIntra,buffer%rdisplsIntra,MPIreal_t,par%commGraphIntra,requestIntra,ierr)
      if(ierr .ne. MPI_SUCCESS) then
        errorcode=ierr
        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
        print *,subname,': Error after call to MPI_neighbor_alltoallv: ',errorstring
      endif

      ! Finish the Intra-node communication 
      call MPI_wait(requestIntra,lstatus,ierr)

      ! location 3 for copyBuffer
      call copyBuffer(nthreads,ithr,buffer,locstring)

      ! Finish the Inter-node communication 
      call MPI_wait(requestInter,lstatus,ierr)

   else

      locstring = TRIM(subname) // ': ' // TRIM(location) 
      !Copy buffer for ithr!=0
      call copyBuffer(nthreads,ithr,buffer,locstring)

   endif
#else
    call abortmp('bndry_exchange_a2ao requires MPI-3 feature support')
#endif

  end subroutine bndry_exchange_a2ao

  subroutine bndry_exchangeV_p2p(par,nthreads,ithr,buffer,location)
    use kinds, only : log_kind
    use edgetype_mod, only : Edgebuffer_t
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use thread_mod, only : omp_in_parallel, omp_get_thread_num
#ifdef _MPI
    use parallel_mod, only : parallel_t, abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success
#else
    use parallel_mod, only : parallel_t, abortmp
#endif
    type (parallel_t)                 :: par
    integer, intent(in)               :: nthreads
    integer                           :: ithr  ! The OpenMP thread ID
    type (EdgeBuffer_t)               :: buffer
    character(len=*), optional        :: location

    type (Schedule_t),pointer         :: pSchedule
    type (Cycle_t),pointer            :: pCycle
    integer                           :: dest,length,tag
    integer                           :: icycle,ierr
    integer                           :: iptr,source,nlyr
!    integer                           :: imptr
    integer                           :: nSendCycles,nRecvCycles
    integer                           :: errorcode,errorlen
    character*(80)                    :: errorstring
    character(len=*),       parameter :: subname = 'bndry_exchangeV_p2p'
    character(len=80)                 :: locstring
!
    logical(kind=log_kind), parameter :: Debug=.FALSE.
!    logical :: ompthreadMissmatch
!
!    integer                           :: i,j
!    integer :: lenMovePtr

    pSchedule => Schedule(1)
!    nlyr = buffer%nlyr
!    ompthreadMissmatch = .false. 

!    lenMovePtr = size(buffer%moveptr)

  if(ithr == 0) then 
    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles

    !==================================================
    !  Post the Receives 
    !==================================================
    do icycle=1,nRecvCycles
       pCycle         => pSchedule%RecvCycle(icycle)
       source          = pCycle%source - 1
       length          = buffer%rcountsFull(icycle)
       tag             = buffer%tag
       iptr            = buffer%rdisplsFull(icycle)
       if(Debug) print *,'IAM: ',iam, subname,': MPI_Irecv: SRC:',source,'LENGTH:',length,'PTR: ',iptr, 'TAG: ',tag
       call MPI_Irecv(buffer%receive(iptr+1),length,MPIreal_t, &
            source,tag,par%comm,buffer%Rrequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,subname,': Error after call to MPI_Irecv: ',errorstring
       endif
    end do    ! icycle

    !==================================================
    !  Fire off the sends
    !==================================================
    do icycle=1,nSendCycles
       pCycle         => pSchedule%SendCycle(icycle)
       dest            = pCycle%dest - 1
       length          = buffer%scountsFull(icycle)
       tag             = buffer%tag
       iptr            = buffer%sdisplsFull(icycle)
       if(Debug) print *,'IAM: ',iam, subname,': MPI_Isend: DEST:',dest,'LENGTH:',length,'PTR: ',iptr, 'TAG: ',tag
       call MPI_Isend(buffer%buf(iptr+1),length,MPIreal_t,dest,tag,par%comm,buffer%Srequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,subname,': Error after call to MPI_Isend: ',errorstring
       endif
    end do    ! icycle
    
      locstring = TRIM(subname) // ': ' // TRIM(location) 
      call copyBuffer(nthreads,ithr,buffer,locstring)

      call MPI_Waitall(nSendCycles,buffer%Srequest,buffer%status,ierr)
      call MPI_Waitall(nRecvCycles,buffer%Rrequest,buffer%status,ierr)
  else  ! else for non-master threads

      locstring = TRIM(subname) // ': ' // TRIM(location) 
      !Copy buffer for ithr!=0
      call copyBuffer(nthreads,ithr,buffer,locstring)

  endif

  end subroutine bndry_exchangeV_p2p

  subroutine bndry_exchangeS_p2p(par,nthreads,ithr,buffer,location)
    use kinds, only : log_kind
    use edgetype_mod, only : Edgebuffer_t
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use thread_mod, only : omp_in_parallel, omp_get_thread_num
#ifdef _MPI
    use parallel_mod, only : parallel_t, abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success
#else
    use parallel_mod, only : parallel_t, abortmp
#endif

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
    character(len=*), parameter       :: subname = 'bndry_exchangeS_p2p'
    character(len=80)                 :: locstring
    logical(kind=log_kind), parameter :: Debug=.FALSE.

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
       if(Debug) print *,subname,': MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
       call MPI_Isend(buffer%buf(iptr),length,MPIreal_t,dest,tag,par%comm,buffer%Srequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,subname,': Error after call to MPI_Isend: ',errorstring
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
       if(Debug) print *,subname,': MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
       call MPI_Irecv(buffer%receive(iptr),length,MPIreal_t, &
            source,tag,par%comm,buffer%Rrequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,subname,': Error after call to MPI_Irecv: ',errorstring
       endif
    end do    ! icycle

    locstring = TRIM(subname) // ': ' // TRIM(location)
    call copyBuffer(nthreads,ithr,buffer,locstring)
    call MPI_Waitall(nSendCycles,buffer%Srequest,buffer%status,ierr)
    call MPI_Waitall(nRecvCycles,buffer%Rrequest,buffer%status,ierr)
  else
    locstring = TRIM(subname) // ': ' // TRIM(location)
    call copyBuffer(nthreads,ithr,buffer,locstring)
  endif

  end subroutine bndry_exchangeS_p2p

  subroutine bndry_exchangeS_p2p_start(par,nthreads,ithr,buffer,location)

    use kinds, only : log_kind
    use edgetype_mod, only : Edgebuffer_t
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use thread_mod, only : omp_in_parallel, omp_get_thread_num
#ifdef _MPI
    use parallel_mod, only : parallel_t, abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success
#else
    use parallel_mod, only : parallel_t, abortmp
#endif
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
    character(len=*),       parameter :: subname = 'bndry_exchangeS_p2p_start'
    logical(kind=log_kind), parameter :: Debug=.FALSE.

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
       if(Debug) print *,subname,': MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
       call MPI_Isend(buffer%buf(iptr),length,MPIreal_t,dest,tag,par%comm,buffer%Srequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,subname,': Error after call to MPI_Isend: ',errorstring
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
       if(Debug) print *,subname,': MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
       call MPI_Irecv(buffer%receive(iptr),length,MPIreal_t, &
            source,tag,par%comm,buffer%Rrequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,subname,': Error after call to MPI_Irecv: ',errorstring
       endif
    end do    ! icycle
  endif
    
  end subroutine bndry_exchangeS_p2p_start

  subroutine bndry_exchangeS_p2p_finish(par,nthreads,ithr,buffer,location)
    use kinds, only : log_kind
    use edgetype_mod, only : Edgebuffer_t
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use thread_mod, only : omp_in_parallel, omp_get_thread_num
#ifdef _MPI
    use parallel_mod, only : parallel_t, abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success
#else
    use parallel_mod, only : parallel_t, abortmp
#endif
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
    character(len=*),       parameter :: subname = 'bndry_exchangeS_p2p_finish'
    character(len=80)                 :: locstring
    logical(kind=log_kind), parameter :: Debug=.FALSE.

    integer                           :: i,j
    logical :: ompthreadMissmatch
    integer        :: lenMovePtr

    pSchedule => Schedule(1)
    nlyr = buffer%nlyr

  locstring = TRIM(subname) // ': ' // TRIM(location)
  call copyBuffer(nthreads,ithr,buffer,locstring)

  if(ithr == 0) then 

    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles

    call MPI_Waitall(nSendCycles,buffer%Srequest,buffer%status,ierr)
    call MPI_Waitall(nRecvCycles,buffer%Rrequest,buffer%status,ierr)
  endif

  end subroutine bndry_exchangeS_p2p_finish

  subroutine long_bndry_exchangeV_nonth(par,buffer)
    use kinds, only : log_kind
    use edgetype_mod, only : LongEdgebuffer_t
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use thread_mod, only : omp_in_parallel
#ifdef _MPI
    use parallel_mod, only : parallel_t, abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success
#else
    use parallel_mod, only : parallel_t, abortmp
#endif
    type (parallel_t)                 :: par
    type (LongEdgeBuffer_t)           :: buffer

    type (Schedule_t),pointer         :: pSchedule
    type (Cycle_t),pointer            :: pCycle
    integer                           :: dest,length,tag
    integer                           :: icycle,ierr
    integer                           :: iptr,source,nlyr
    integer                           :: nSendCycles,nRecvCycles
    integer                           :: errorcode,errorlen
    character*(80)                    :: errorstring
    character(len=*),       parameter :: subname = 'long_bndry_exchangeV_nonth'
    logical(kind=log_kind), parameter :: Debug=.FALSE.

    integer                           :: i

#ifdef _MPI
    if(omp_in_parallel()) then
       print *,subname,': Warning you are calling a non-thread safe'
       print *,'		 routine inside a threaded region....     '
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
       !DBG if(Debug) print *,subname,': MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
       call MPI_Isend(buffer%buf(1,iptr),length,MPIinteger_t,dest,tag,par%comm,Srequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,subname,': Error after call to MPI_Isend: ',errorstring
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
       !DBG if(Debug) print *,subname,': MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
       call MPI_Irecv(buffer%receive(1,iptr),length,MPIinteger_t, &
            source,tag,par%comm,Rrequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,subname,': Error after call to MPI_Irecv: ',errorstring
       endif
    end do    ! icycle


    !==================================================
    !  Wait for all the receives to complete
    !==================================================

    call MPI_Waitall(nSendCycles,Srequest,status,ierr)
    call MPI_Waitall(nRecvCycles,Rrequest,status,ierr)
    do icycle=1,nRecvCycles
       pCycle         => pSchedule%RecvCycle(icycle)
       length             = pCycle%lengthP
       iptr            = pCycle%ptrP
       do i=0,length-1
          buffer%buf(1:nlyr,iptr+i) = buffer%receive(1:nlyr,iptr+i)
       enddo
    end do   ! icycle

#endif

  end subroutine long_bndry_exchangeV_nonth
  !********************************************************************************
  !
  !********************************************************************************
 subroutine bndry_exchangeV_threaded(hybrid,buffer,location)
    use hybrid_mod, only : hybrid_t, PrintHybrid
    use edgetype_mod, only : Edgebuffer_t
    implicit none

    type (hybrid_t)               :: hybrid
    type (EdgeBuffer_t)           :: buffer
    character(len=*),   parameter :: subname = 'bndry_exchangeV_threaded'
    character(len=*), optional        :: location
    integer :: localsense

!    call t_adj_detailf(+2)
    call gbarrier(buffer%gbarrier, hybrid%ithr)
    if(buffer%bndry_type == HME_BNDRY_A2A) then 
       call bndry_exchange_a2a(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    else if(buffer%bndry_type == HME_BNDRY_A2AO) then 
       call bndry_exchange_a2ao(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    else if(buffer%bndry_type == HME_BNDRY_GET1) then 
       call bndry_exchange_get1(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    else if(buffer%bndry_type == HME_BNDRY_GET2) then 
       call bndry_exchange_get2(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    else if(buffer%bndry_type == HME_BNDRY_PUT1) then 
       call bndry_exchange_put1(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    else if(buffer%bndry_type == HME_BNDRY_PUT2) then 
       call bndry_exchange_put2(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    else
       call bndry_exchangeV_p2p(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    endif
    call gbarrier(buffer%gbarrier, hybrid%ithr)
!    call t_adj_detailf(-2)

  end subroutine bndry_exchangeV_threaded

  subroutine bndry_exchangeV_nonthreaded(par,buffer,location)
    use parallel_mod, only : parallel_t
    use edgetype_mod, only : Edgebuffer_t
    implicit none

    type (parallel_t)           :: par
    type (EdgeBuffer_t)         :: buffer
    character(len=*), optional  :: location

    ! local
    character(len=*), parameter :: subname = 'bndry_exchangeV_nonthreaded'
    integer                     :: ithr
    integer                     :: nthreads
    integer                     :: localsense

!    call t_adj_detailf(+2)
    ithr=0
    nthreads = 1
!    localsense = 0
    if(buffer%bndry_type == HME_BNDRY_A2A) then 
       call bndry_exchange_a2a(par,nthreads,ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_A2AO) then 
       call bndry_exchange_a2ao(par,nthreads,ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_GET1) then 
       call bndry_exchange_get1(par,nthreads,ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_GET2) then 
       call bndry_exchange_get2(par,nthreads,ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_PUT1) then 
       call bndry_exchange_put1(par,nthreads,ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_PUT2) then 
       call bndry_exchange_put2(par,nthreads,ithr,buffer,location)
    else
       call bndry_exchangeV_p2p(par,nthreads,ithr,buffer,location)
    endif
!    call t_adj_detailf(-2)

  end subroutine bndry_exchangeV_nonthreaded

 subroutine bndry_exchangeS_threaded(hybrid,buffer,location)
    use hybrid_mod, only : hybrid_t
    use edgetype_mod, only : Edgebuffer_t
    implicit none

    type (hybrid_t)             :: hybrid
    type (EdgeBuffer_t)         :: buffer
    character(len=*), optional  :: location

    character(len=*), parameter :: subname = 'bndry_exchangeS_threaded'

!    call t_adj_detailf(+2)
    call gbarrier(buffer%gbarrier, hybrid%ithr)
    if(buffer%bndry_type == HME_BNDRY_A2A) then
       call bndry_exchange_a2a(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_A2AO) then 
       call bndry_exchange_a2ao(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_GET1) then 
       call bndry_exchange_get1(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_GET2) then 
       call bndry_exchange_get2(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_PUT1) then 
       call bndry_exchange_put1(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_PUT2) then 
       call bndry_exchange_put2(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    else
       call bndry_exchangeS_p2p(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    endif
    call gbarrier(buffer%gbarrier, hybrid%ithr)
!    call t_adj_detailf(-2)

 end subroutine bndry_exchangeS_threaded

 subroutine bndry_exchangeS_threaded_start(hybrid,buffer,location)
    use hybrid_mod, only : hybrid_t
    use edgetype_mod, only : Edgebuffer_t
    implicit none

    type (hybrid_t)             :: hybrid
    type (EdgeBuffer_t)         :: buffer
    character(len=*), optional  :: location

    character(len=*), parameter :: subname = 'bndry_exchangeS_threaded_start'

!    call t_adj_detailf(+2)
    call gbarrier(buffer%gbarrier, hybrid%ithr)
    call bndry_exchangeS_p2p_start(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
!    call t_adj_detailf(-2)

 end subroutine bndry_exchangeS_threaded_start

 subroutine bndry_exchangeS_threaded_finish(hybrid,buffer,location)
    use hybrid_mod, only : hybrid_t
    use edgetype_mod, only : Edgebuffer_t
    implicit none

    type (hybrid_t)             :: hybrid
    type (EdgeBuffer_t)         :: buffer
    character(len=*), optional  :: location

    character(len=*), parameter :: subname = 'bndry_exchangeS_threaded_finish'

!    call t_adj_detailf(+2)
    call bndry_exchangeS_p2p_finish(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    call gbarrier(buffer%gbarrier, hybrid%ithr)
!    call t_adj_detailf(-2)

 end subroutine bndry_exchangeS_threaded_finish

 subroutine bndry_exchangeS_nonthreaded(par,buffer,location)
    use parallel_mod, only : parallel_t
    use edgetype_mod, only : Edgebuffer_t
    implicit none

    type (parallel_t)           :: par
    type (EdgeBuffer_t)         :: buffer
    character(len=*), optional  :: location

    integer                     :: ithr
    integer                     :: nthreads
!    integer                     :: localsense
    character(len=*), parameter :: subname = 'bndry_exchangeS_nonthreaded'

!    call t_adj_detailf(+2)
    ithr=0
    nthreads = 1
    if(buffer%bndry_type == HME_BNDRY_A2A) then 
       call bndry_exchange_a2a(par,nthreads,ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_A2AO) then 
       call bndry_exchange_a2ao(par,nthreads,ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_GET1) then 
       call bndry_exchange_get1(par,nthreads,ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_GET2) then 
       call bndry_exchange_get2(par,nthreads,ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_PUT1) then 
       call bndry_exchange_put1(par,nthreads,ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_PUT2) then 
       call bndry_exchange_put2(par,nthreads,ithr,buffer,location)
    else
       call bndry_exchangeS_p2p(par,nthreads,ithr,buffer,location)
    endif
!    call t_adj_detailf(-2)

  end subroutine bndry_exchangeS_nonthreaded

 subroutine bndry_exchangeS_nonthreaded_start(par,buffer,location)
    use parallel_mod, only : parallel_t
    use edgetype_mod, only : Edgebuffer_t
    implicit none

    type (parallel_t)           :: par
    type (EdgeBuffer_t)         :: buffer
    character (len=*), optional :: location

    integer                     :: ithr
    integer                     :: nthreads
    character(len=*), parameter :: subname = 'bndry_exchangeS_nonthreaded_start'

!    call t_adj_detailf(+2)
    ithr=0
    nthreads=1
    call bndry_exchangeS_p2p_start(par,nthreads,ithr,buffer,location)
!    call t_adj_detailf(-2)

 end subroutine bndry_exchangeS_nonthreaded_start

 subroutine bndry_exchangeS_nonthreaded_finish(par,buffer,location)
    use parallel_mod, only : parallel_t
    use edgetype_mod, only : Edgebuffer_t
    implicit none

    type (parallel_t)                 :: par
    integer                           :: ithr
    type (EdgeBuffer_t)               :: buffer
    character (len=*), optional :: location
    integer :: nthreads

    character(len=*), parameter :: subname = 'bndry_exchangeS_nonthreaded_finish'

!    call t_adj_detailf(+2)
    ithr=0
    nthreads=1
    call bndry_exchangeS_p2p_finish(par,nthreads,ithr,buffer,location)
!    call t_adj_detailf(-2)

  end subroutine bndry_exchangeS_nonthreaded_finish

end module bndry_mod
