#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module filter_mod
  use kinds, only : real_kind
  use dimensions_mod, only : np
  use perf_mod, only: t_startf, t_stopf, t_barrierf ! _EXTERNAL
  implicit none
  private

  real (kind=real_kind), parameter :: one=1.0D0, zero=0.0D0
  type, public :: filter_t
     real (kind=real_kind) :: FmatP(np,np)    ! Filter matrix
     character(len=8)      :: type    ! filter name: bv,fm, etc...
  end type filter_t

  public :: taylor_filter_create
  public :: fm_filter_create

  public :: filter_P
  public :: fm_transfer
  public :: bv_transfer

  interface filter_matrix_create
     module procedure filter_matrix_create_inv
     module procedure filter_matrix_create_noinv
  end interface

  private :: bvsigma
  private :: filter_matrix_create
contains
  ! ==================================
  ! taylor_filter_create:
  !
  ! Create a Taylor filter (requires
  ! communication) and use transfer
  ! function T...
  ! ==================================
  function taylor_filter_create(Tp, mu, gquadp) result(flt)
    use kinds, only : real_kind
    use quadrature_mod, only : quadrature_t, legendre, quad_norm
    real (kind=real_kind), dimension(:) :: Tp          ! transfer fn       
    real (kind=real_kind), intent(in)   :: mu             ! filter viscosity
    type (quadrature_t)  , intent(in)   :: gquadp  ! quadrature points and wts

    type (filter_t)                   :: flt

    ! Local variables 

    real (kind=real_kind), dimension(:), allocatable   :: gammaP ! normalization factor
    real (kind=real_kind), dimension(:,:), allocatable :: legP    ! legendre polynomials

    real (kind=real_kind), dimension(:,:), allocatable :: AP     ! spectral-physical 
    real (kind=real_kind), dimension(:,:), allocatable :: AinvP   ! physical-spectral


    integer                 :: npts
    integer                 :: n,j
    logical, parameter      :: Debug = .FALSE.

    call t_startf('taylor_filter_create')
    allocate(gammaP(np))

    allocate(legP(np,np))
    allocate(AP(np,np))

    allocate(AinvP(np,np))

    gammaP=quad_norm(gquadP,np)

    do j=1,np
       if(Debug) print *,'taylor_filter_create: gammaP(j) := ',gammaP(j)
       legP(:,j)=legendre(gquadP%points(j),np-1)
    end do

    ! -----------------------------------------
    ! compute T(k), the spectral coefficent 
    ! transfer function for the B-V filter
    ! -----------------------------------------

    AP = legP

    do n=1,np
       do j=1,np
          AinvP(j,n)=legP(n,j)*gquadP%weights(j)/gammaP(n)
       end do
    end do

    call filter_matrix_create(AP,Tp,AinvP,flt%FmatP,np,mu)

    flt%type = "bv"

    deallocate(gammaP)
    deallocate(legP)
    deallocate(AP)
    deallocate(AinvP)
    call t_stopf('taylor_filter_create')

  end function taylor_filter_create

  ! ==================================
  ! fm_filter_create:
  !
  ! Fischer-Mullen Filter Constructor
  ! with transfer function T...
  ! ==================================

  function fm_filter_create(Tp, mu, gquadP) result(flt)
    use quadrature_mod, only : quadrature_t, legendre

    real (kind=real_kind),   dimension(:) :: Tp    ! transfer fn
    real (kind=real_kind), intent(in) :: mu       ! filter viscosity
    type (quadrature_t)  , intent(in) :: gquadP    ! quadrature points/wts

    type (filter_t)                   :: flt

    ! Local variables 


    real (kind=real_kind), dimension(:,:), allocatable :: legP,legV    ! legendre polynomials

    real (kind=real_kind),   dimension(:,:), allocatable :: AP,AV      ! spectral-physical 

    integer                 :: k,n,j

    call t_startf('fm_filter_create')

    allocate(legP(np,np))

    allocate(AP(np,np))

    !JMD     allocate(flt%Fmat(npts,npts))
    !JMD     allocate(flt%tempt(npts,npts))

    do j=1,np
       legP(:,j)=legendre(gquadP%points(j),np-1)
       AP(1,j)=legP(1,j)
       AP(2,j)=legP(2,j)
       do n=3,np
          AP(n,j)=legP(n,j)-legP(n-2,j)
       end do
    end do


    call filter_matrix_create(AP,Tp,flt%FmatP,np,mu)

    flt%type = "fm"

    deallocate(legP)
    deallocate(AP)
    call t_stopf('fm_filter_create')

  end function fm_filter_create

  ! ==================================================
  ! This routing builds a Fischer-Mullen transfer 
  ! function that looks like:
  !
  !
  !     ^
  !  T  |
  !     |                 |
  !  1  |__________      _v_
  !     |          -_     
  !     |            \  wght
  !     |             \  ___
  !     |             |   ^
  !  0  |-------------|---|>
  !               ^   ^ 
  !     0         kc  npts  k-->
  !
  !  Where kc := npts-kcut is the point below which T = 1.
  ! ==================================================

  function fm_transfer(kcut,wght,npts) result(T)

    integer               :: kcut
    real (kind=real_kind) :: wght
    integer               :: npts
    real (kind=real_kind) :: T(npts)

    ! Local variables

    integer :: k
    integer :: kc
    real (kind=real_kind) :: cut
    real (kind=real_kind) :: amp

    call t_startf('fm_transfer')
    kc  = npts-kcut
    cut = kcut

    do k=1,kc
       T(k)=one
    end do

    do k=kc+1,npts
       amp = wght*(k-kc)*(k-kc)/(cut*cut)    ! quadratic growth
       T(k) = one - amp
    end do
    call t_stopf('fm_transfer')

  end function fm_transfer

  ! ===========================================
  ! bv_transfer:
  ! 
  ! compute Boyd-Vandeven transfer 
  ! function for a spectral element 
  ! with npts degrees of freedom...
  ! ===========================================

  function bv_transfer(p,s,npts) result(T)

    real (kind=real_kind), intent(in) :: p     ! order of sigmoid
    real (kind=real_kind), intent(in) :: s     ! scale of filter
    integer,               intent(in) :: npts  ! number of points
    real (kind=real_kind) :: T(npts)

    integer k
    real (kind=real_kind) :: arg
    real (kind=real_kind) :: rat

    call t_startf('bv_transfer')
    do k=1,npts
       rat = REAL(k,kind=real_kind)/REAL(npts,kind=real_kind)
       if (rat<s) then
          T(k)=one
       else
          arg = (rat-s)/(one-s)
          !
          !JPE: Need to assure that arg <= 1.0 
          !
          T(k)= bvsigma(p,min(arg,one))
       end if
    end do
    call t_stopf('bv_transfer')

  end function bv_transfer

  ! ===========================================
  ! bvsigma:
  !
  ! Boyd - Vandeven sigma function
  ! see p.98 of Taylor, Tribbia and Iskandarani
  ! ===========================================

  function bvsigma(p,x) result(sigma)

    ! CAM code is required to use this function, to handle differences in
    ! Fortran 2008 support and compiler extensions.
    use shr_spfn_mod, only: erfc => shr_spfn_erfc

    real (kind=real_kind), intent(in) :: p
    real (kind=real_kind), intent(in) :: x
    real (kind=real_kind)             :: sigma

    ! Local variables

    real (kind=real_kind) :: om
    real (kind=real_kind) :: xfac
    real (kind=real_kind) :: arg

    call t_startf('bvsigma')

    om=ABS(x)-0.5D0

    if (x==zero) then
       sigma=one
    else if (x==one) then
       sigma=zero
    else
       xfac=one
       if (om /= zero) xfac = SQRT(-LOG((one-4.0D0*om**2))/(4.0D0*om**2))
       arg = 2.0D0*SQRT(p)*om*xfac
       sigma = 0.50D0*erfc(arg)
    end if
    call t_stopf('bvsigma')

  end function bvsigma

  ! ===========================================
  ! filter:
  !
  ! Apply a Filter Matrix
  ! in both directions (x and y) of a 2-d domain.
  !
  ! ===========================================

  subroutine filter_P(p,flt) 
    use kinds, only : int_kind
    type (filter_t)                    :: flt
    real(kind=real_kind),intent(inout) :: p(np,np)

    ! Local

    integer  :: i,j,l
    integer(kind=int_kind), parameter :: unroll = 2
    real(kind=real_kind) :: temptp(np,np) 
    real(kind=real_kind) :: sumx00,sumx01
    real(kind=real_kind) :: sumx10,sumx11


       do j=1,np
          do l=1,np
             sumx00=0.0d0
!DIR$ UNROLL(NP)
             do i=1,np
		sumx00 = sumx00 + flt%FmatP(i,l  )*p(i,j  )
             enddo
             temptP(j  ,l  ) = sumx00	
          enddo
       enddo
       do j=1,np
	  do i=1,np
	     sumx00=0.0d0
!DIR$ UNROLL(NP)
	     do l=1,np
		sumx00 = sumx00 +  flt%FmatP(l,j  )*temptP(l,i  )
	     enddo
	     p(i  ,j  ) = sumx00
          enddo
       enddo

  end subroutine filter_P
  

  ! =================================================================================
  !
  ! filter_matrix_create_inv:
  !
  ! This routing builds a 1D filter matrix, F, from a diagonal transfer function T 
  ! and the transform matrix A that maps a field u into spectral coefficients u_c,
  ! and its inverse A^-1, which is provided, via
  !    
  !   F = (1-mu)*I + mu*A.T.A^-1
  !
  !  where mu is the viscosity term, or weight of the filtering. 
  !
  ! ==================================================================================

  subroutine filter_matrix_create_inv(A,T,Ainv,F,npts,mu)

    integer                           :: npts      
    real (kind=real_kind), intent(in) :: A(npts,npts)
    real (kind=real_kind), intent(in) :: T(npts)
    real (kind=real_kind), intent(in) :: Ainv(npts,npts)
    real (kind=real_kind), intent(in) :: mu

    real (kind=real_kind) :: F(npts,npts)

    ! Local variables

    integer j,k,n
    real (kind=real_kind) :: fsum

    do k=1,npts  
       do j=1,npts
          fsum=zero
          do n=1,npts
             fsum = fsum + A(n,j)*T(n)*Ainv(k,n)
          end do
          if (j == k) then
             F(k,j) = (one-mu) + mu*fsum
          else
             F(k,j) = mu*fsum
          end if
       end do
    end do

  end subroutine filter_matrix_create_inv

  ! =================================================================================
  !
  ! filter_matrix_create_noinv:
  !
  ! This routing builds a 1D filter matrix, F, from a diagonal transfer function T 
  ! and the transform matrix A that maps a field u into spectral coefficients u_c,
  ! and its inverse A^-1, which is computed from A, via
  !    
  !   F = (1-mu)*I + mu*A.T.A^-1
  !
  !  where mu is the viscosity term, or weight of the filtering. 
  !
  ! =================================================================================

  subroutine filter_matrix_create_noinv(A,T,F,npts,mu)

    integer                           :: npts      
    real (kind=real_kind), intent(in) :: A(npts,npts)
    real (kind=real_kind), intent(in) :: T(npts)
    real (kind=real_kind), intent(in) :: mu

    real (kind=real_kind) :: F(npts,npts)

    ! Local variables

    integer               :: n,j,k
    integer               :: ipiv(npts)
    real (kind=real_kind) :: Ainv(npts,npts)
    real (kind=real_kind) :: fsum
    integer               :: ierr

    Ainv = A
    ierr=gaujordf(Ainv,npts,npts,ipiv)

    do k=1,npts  
       do j=1,npts
          fsum=zero
          do n=1,npts
             fsum = fsum + A(n,j)*T(n)*Ainv(k,n)
          end do
          if (j == k) then
             F(k,j) = (one-mu) + mu*fsum
          else
             F(k,j) = mu*fsum
          end if
       end do
    end do


  end subroutine filter_matrix_create_noinv

  ! ========================================================
  !
  ! Gauss-Jordan matrix inversion with full pivoting
  !
  ! Num. Rec. p. 30, 2nd Ed., Fortran
  !
  ! appropriated from code supplied by Paul Fischer and
  ! converted to F90 by RDL
  !
  ! a     is an m x n real matrix to be inverted in place
  ! ipiv  is the pivot vector
  !
  ! ========================================================

  function gaujordf(a,m,n,ipiv) result(ierr)

    integer, intent(in)   :: m,n
    real (kind=real_kind) :: a(m,n)
    integer               :: ipiv(n)
    integer               :: ierr

    integer :: indr(m),indc(n)
    integer :: i,j,k
    integer :: ir,jc

    real (kind=real_kind) :: amx,piv,tmp,work
    real (kind=real_kind) :: rmult(m)
    real (kind=real_kind) :: eps

    ierr = 0
    eps = 1.0D-14
    !
    do k=1,m
       ipiv(k)=0
    end do
    !
    do k=1,m
       amx=zero
       !
       !    Pivot search
       !
       do i=1,m
          if (ipiv(i) /= 1) then
             do j=1,m
                if (ipiv(j) == 0) then
                   if (abs(a(i,j)) >= amx) then
                      amx = ABS(a(i,j))
                      ir  = i
                      jc  = j
                   end if
                else if (ipiv(j) > 1) then
                   ierr = -ipiv(j)
                   return
                end if
             end do
          end if
       end do

       ipiv(jc) = ipiv(jc) + 1
       !
       !    Swap rows
       ! 
       if (ir /= jc) then
          do j=1,n
             tmp  = a(ir,j)
             a(ir,j) = a(jc,j)
             a(jc,j) = tmp
          end do
       end if

       indr(k)=ir
       indc(k)=jc

       !     write(6 ,*) k,' Piv:',jc,a(jc,jc)
       !       write(28,*) k,' Piv:',jc,a(jc,jc)

       if (abs(a(jc,jc)) < eps) then
          write(6 ,*) 'Gauss Jordan Pivot too small:',jc,a(jc,jc)
          !       write(28,*) 'small Gauss Jordan Piv:',jc,a(jc,jc)
          ierr = jc
          return
       end if

       piv = one/a(jc,jc)
       a(jc,jc)=one

       do j=1,n
          a(jc,j) = a(jc,j)*piv
       end do

       do j=1,n
          work    = a(jc,j)
          a(jc,j) = a(1 ,j)
          a(1 ,j) = work
       end do

       do i=2,m
          rmult(i) = a(i,jc)
          a(i,jc)  = zero
       end do

       do j=1,n
          do i=2,m
             a(i,j) = a(i,j) - rmult(i)*a(1,j)
          end do
       end do

       do j=1,n
          work    = a(jc,j)
          a(jc,j) = a(1 ,j)
          a(1 ,j) = work
       end do

    end do

    !
    !     Unscramble matrix
    !

    do j=m,1,-1
       if (indr(j) /= indc(j)) then
          do i=1,m
             tmp=a(i,indr(j))
             a(i,indr(j))=a(i,indc(j))
             a(i,indc(j))=tmp
          end do
       end if
    end do
  end function gaujordf

end module filter_mod
