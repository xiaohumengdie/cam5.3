module reduction_mod
  use shr_kind_mod,   only: r8=>shr_kind_r8
  use spmd_utils,     only: mpi_sum, mpi_min, mpi_max, mpi_real8, mpi_integer
  use spmd_utils,     only: mpi_success
  use cam_abortutils, only: endrun

  implicit none
  private

  type, public :: ReductionBuffer_int_1d_t
     integer, dimension(:), pointer :: buf
     integer :: len=0
     integer :: ctr
  end type ReductionBuffer_int_1d_t

  type, public :: ReductionBuffer_r_1d_t
     real (kind=r8), dimension(:), pointer :: buf
     integer :: len=0
     integer :: ctr
  end type ReductionBuffer_r_1d_t

  type, public :: ReductionBuffer_ordered_1d_t
     real (kind=r8), dimension(:,:),pointer :: buf
     integer :: len=0
     integer :: ctr
  end type ReductionBuffer_ordered_1d_t

  public :: ParallelMin
  public :: ParallelMax

  type (ReductionBuffer_int_1d_t), public :: red_max_int
  type (ReductionBuffer_int_1d_t), public :: red_sum_int
  type (ReductionBuffer_r_1d_t),   public :: red_flops,red_timer
  type (ReductionBuffer_r_1d_t),   public :: red_max
  type (ReductionBuffer_r_1d_t),   public :: red_min
  type (ReductionBuffer_r_1d_t),   public :: red_sum
#ifndef Darwin
  SAVE red_max_int, red_sum_int, red_flops, red_max, red_min, red_sum, red_timer
#endif

  interface ParallelMin
     module procedure ParallelMin1d
     module procedure ParallelMin0d
  end interface
  interface ParallelMax
     module procedure ParallelMax1d_int
     module procedure ParallelMax2d_int
     module procedure ParallelMax1d
     module procedure ParallelMax0d
     module procedure ParallelMax0d_int
  end interface

  interface pmax_mt
     module procedure pmax_mt_int_1d
     module procedure pmax_mt_r_1d
  end interface

  interface pmin_mt
     module procedure pmin_mt_r_1d
  end interface

  interface InitReductionBuffer
     module procedure InitReductionBuffer_int_1d
     module procedure InitReductionBuffer_r_1d
     module procedure InitReductionBuffer_ordered_1d
  end interface

  public :: InitReductionBuffer
  public :: pmax_mt, pmin_mt
  public :: ElementSum_1d

contains

  function ParallelMin1d(data,hybrid) result(pmin)
    use hybrid_mod, only : hybrid_t

    real(kind=r8),   intent(in) :: data(:)
    type (hybrid_t), intent(in) :: hybrid
    real(kind=r8)               :: pmin

    real(kind=r8)               :: tmp(1)


    tmp(1) = MINVAL(data)
    call pmin_mt(red_min,tmp,1,hybrid)
    pmin = red_min%buf(1)

  end function ParallelMin1d

  function ParallelMin0d(data,hybrid) result(pmin)
    use hybrid_mod, only : hybrid_t
    implicit none
    real(kind=r8), intent(in)    :: data
    type (hybrid_t),      intent(in)    :: hybrid
    real(kind=r8)                :: pmin
    real(kind=r8)                :: tmp(1)
    tmp(1) = data
    call pmin_mt(red_min,tmp,1,hybrid)
    pmin = red_min%buf(1)

  end function ParallelMin0d
  !==================================================
  function ParallelMax2d_int(data, n, m, hybrid) result(pmax)
    use hybrid_mod, only : hybrid_t
    implicit none
    integer, intent(in)                 :: n,m
    integer, intent(in), dimension(n,m) :: data
    type (hybrid_t), intent(in)         :: hybrid
    integer, dimension(n,m)             :: pmax
    integer, dimension(n*m)             :: tmp
    integer :: i,j
    do i=1,n 
      do j=1,m
        tmp(i+(j-1)*n) = data(i,j)
      enddo 
    enddo 
    call pmax_mt(red_max_int,tmp,n*m,hybrid)
    do i=1,n 
      do j=1,m
        pmax(i,j) = red_max_int%buf(i+(j-1)*n) 
      enddo 
    enddo 
  end function ParallelMax2d_int

  function ParallelMax1d_int(data, len, hybrid) result(pmax)
    use hybrid_mod, only : hybrid_t
    implicit none
    integer, intent(in)                 :: len
    integer, intent(in), dimension(len) :: data
    type (hybrid_t), intent(in)         :: hybrid
    integer, dimension(len)             :: pmax, tmp

    tmp = data(:)
    call pmax_mt(red_max_int,tmp,len,hybrid)
    pmax(:) = red_max_int%buf(1:len)

  end function ParallelMax1d_int
  function ParallelMax1d(data,hybrid) result(pmax)
    use hybrid_mod, only : hybrid_t
    implicit none
    real(kind=r8), intent(in)    :: data(:)
    type (hybrid_t),      intent(in)    :: hybrid
    real(kind=r8)                :: pmax

    real(kind=r8)                :: tmp(1)


    tmp(1) = MAXVAL(data)
    call pmax_mt(red_max,tmp,1,hybrid)
    pmax = red_max%buf(1)

  end function ParallelMax1d
  function ParallelMax0d(data,hybrid) result(pmax)
    use hybrid_mod, only : hybrid_t
    implicit none
    real(kind=r8), intent(in)    :: data
    type (hybrid_t),      intent(in)    :: hybrid
    real(kind=r8)                :: pmax
    real(kind=r8)                :: tmp(1)

    tmp(1)=data

    call pmax_mt(red_max,tmp,1,hybrid)
    pmax = red_max%buf(1)

  end function ParallelMax0d
  function ParallelMax0d_int(data,hybrid) result(pmax)
    use hybrid_mod, only : hybrid_t
    implicit none
    integer             , intent(in)    :: data
    type (hybrid_t),      intent(in)    :: hybrid
    integer                             :: pmax
    integer                             :: tmp(1)

    tmp(1)=data

    call pmax_mt(red_max_int,tmp,1,hybrid)
    pmax = red_max_int%buf(1)

  end function ParallelMax0d_int
  !==================================================
  subroutine InitReductionBuffer_int_1d(red,len)
    use thread_mod, only: omp_get_num_threads
    integer, intent(in)           :: len
    type (ReductionBuffer_int_1d_t),intent(out) :: red

    if (omp_get_num_threads()>1) then
       call endrun("Error: attempt to allocate reduction buffer in threaded region")
    endif

    ! if buffer is already allocated and large enough, do nothing
    if (len > red%len) then
       !buffer is too small, or has not yet been allocated
       if (red%len>0) deallocate(red%buf)
       red%len  = len
       allocate(red%buf(len))
       red%buf  = 0
       red%ctr  = 0
    endif

  end subroutine InitReductionBuffer_int_1d
  !****************************************************************
  subroutine InitReductionBuffer_r_1d(red,len)
    use thread_mod, only: omp_get_num_threads
    integer, intent(in)           :: len
    type (ReductionBuffer_r_1d_t),intent(out) :: red

    if (omp_get_num_threads()>1) then
       call endrun("Error: attempt to allocate reduction buffer in threaded region")
    endif

    if (len > red%len) then
       if (red%len>0) deallocate(red%buf)
       red%len  = len
       allocate(red%buf(len))
       red%buf  = 0.0_R8
       red%ctr  = 0
    endif
  end subroutine InitReductionBuffer_r_1d
  !****************************************************************
  subroutine InitReductionBuffer_ordered_1d(red,len,nthread)
    use thread_mod, only: omp_get_num_threads
    integer, intent(in)           :: len
    integer, intent(in)           :: nthread
    type (ReductionBuffer_ordered_1d_t),intent(out) :: red

    if (omp_get_num_threads()>1) then
       call endrun("Error: attempt to allocate reduction buffer in threaded region")
    endif

    if (len > red%len) then
       if (red%len>0) deallocate(red%buf)
       red%len  = len
       allocate(red%buf(len,nthread+1))
       red%buf  = 0.0_R8
       red%ctr  = 0
    endif
  end subroutine InitReductionBuffer_ordered_1d

  ! =======================================
  ! pmax_mt:
  !
  ! thread safe, parallel reduce maximum
  ! of a one dimensional reduction vector
  ! =======================================

  subroutine pmax_mt_int_1d(red,redp,len,hybrid)
    use hybrid_mod, only : hybrid_t

    type (ReductionBuffer_int_1d_t)   :: red       ! shared memory reduction buffer struct
    integer,               intent(in) :: len       ! buffer length
    integer, intent(inout)            :: redp(len) ! thread private vector of partial sum
    type (hybrid_t),       intent(in) :: hybrid    ! parallel handle

    ! Local variables
#ifdef _MPI
    integer ierr
#endif

    integer  :: k
    if (len>red%len) then
      call endrun('ERROR: threadsafe reduction buffer too small')
    end if

    if (red%ctr == 0) red%buf(1:len)= -9999
    if (red%ctr < hybrid%NThreads) then
       do k=1,len
          red%buf(k)=MAX(red%buf(k),redp(k))
       enddo
       red%ctr=red%ctr+1
    end if
    if (red%ctr == hybrid%NThreads) red%ctr=0
#ifdef _MPI
    if (hybrid%ithr==0) then

       call MPI_Allreduce(red%buf(1),redp,len,Mpi_integer, &
            MPI_MAX,hybrid%par%comm,ierr)

       red%buf(1:len)=redp(1:len)
    end if
#endif

  end subroutine pmax_mt_int_1d
  
  subroutine pmax_mt_r_1d(red,redp,len,hybrid)
    use hybrid_mod, only : hybrid_t

    type (ReductionBuffer_r_1d_t)     :: red     ! shared memory reduction buffer struct
    real (kind=r8), intent(inout) :: redp(:) ! thread private vector of partial sum
    integer,               intent(in) :: len     ! buffer length
    type (hybrid_t),       intent(in) :: hybrid  ! parallel handle

    ! Local variables
#ifdef _MPI
    integer ierr
#endif

    integer  :: k
    if (len>red%len) then
      call endrun('ERROR: threadsafe reduction buffer too small')
    end if

    if (red%ctr == 0) red%buf(1:len)= -9.11e30_r8
    if (red%ctr < hybrid%NThreads) then
       do k=1,len
          red%buf(k)=MAX(red%buf(k),redp(k))
       enddo
       red%ctr=red%ctr+1
    end if
    if (red%ctr == hybrid%NThreads) red%ctr=0
#ifdef _MPI
    if (hybrid%ithr==0) then

       call MPI_Allreduce(red%buf(1),redp,len,Mpi_real8, &
            MPI_MAX,hybrid%par%comm,ierr)

       red%buf(1:len)=redp(1:len)
    end if
#endif

  end subroutine pmax_mt_r_1d

  ! =======================================
  ! pmin_mt:
  !
  ! thread safe, parallel reduce maximum
  ! of a one dimensional reduction vector
  ! =======================================

  subroutine pmin_mt_r_1d(red,redp,len,hybrid)
    use hybrid_mod, only : hybrid_t

    type (ReductionBuffer_r_1d_t)     :: red     ! shared memory reduction buffer struct
    real (kind=r8), intent(inout) :: redp(:) ! thread private vector of partial sum
    integer,               intent(in) :: len     ! buffer length
    type (hybrid_t),       intent(in) :: hybrid  ! parallel handle

    ! Local variables

#ifdef _MPI
    integer :: ierr
#endif
    integer :: k

    if (len>red%len) then
      call endrun('ERROR: threadsafe reduction buffer too small')
    end if

    if (red%ctr == 0) red%buf(1:len)= 9.11e30_r8
    if (red%ctr < hybrid%NThreads) then
       do k=1,len
          red%buf(k)=MIN(red%buf(k),redp(k))
       enddo
       red%ctr=red%ctr+1
    end if
    if (red%ctr == hybrid%NThreads) red%ctr=0
#ifdef _MPI
    if (hybrid%ithr==0) then

       call MPI_Allreduce(red%buf(1),redp,len,Mpi_real8, &
            MPI_MIN,hybrid%par%comm,ierr)

       red%buf(1:len)=redp(1:len)
    end if
#endif

  end subroutine pmin_mt_r_1d

  subroutine ElementSum_1d(res,variable,type,hybrid)
    use hybrid_mod,     only: hybrid_t
    use dimensions_mod, only: nelem
    use parallel_mod,   only: ORDERED

    ! ==========================
    !     Arguments
    ! ==========================
    real(kind=r8),   intent(out) :: res
    real(kind=r8),   intent(in)  :: variable(:)
    integer,         intent(in)  :: type
    type (hybrid_t), intent(in)  :: hybrid 

    ! ==========================
    !       Local Variables
    ! ==========================

    !
    ! Note this is a real kludge here since it may be used for
    !  arrays of size other then nelem
    !

#ifdef _MPI
    integer           :: errorcode,errorlen
    character(len=80) :: errorstring

    real(kind=r8)     :: local_sum
    integer           :: ierr
#else
    integer           :: i
#endif

#ifdef _MPI
    if(hybrid%ithr == 0) then 
          local_sum=SUM(variable)
          call MPI_Barrier(hybrid%par%comm,ierr)

          call MPI_Allreduce(local_sum,res,1,Mpi_real8, &
               MPI_SUM,hybrid%par%comm,ierr)
          if(ierr .ne. MPI_SUCCESS) then 
             errorcode=ierr
             call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
             print *,'ElementSum_1d: Error after call to MPI_Allreduce: ',errorstring
          endif
    endif
#else
    if(hybrid%ithr == 0) then 
       if(type == ORDERED) then
          ! ===========================
          !  Perform the ordererd sum
          ! ===========================
          res = 0.0_r8
          do i=1,nelem
             res = res + variable(i)
          enddo
       else
          res=SUM(variable)
       endif
    endif
#endif

  end subroutine ElementSum_1d

end module reduction_mod
