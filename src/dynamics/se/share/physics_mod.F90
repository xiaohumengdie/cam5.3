module physics_mod
  
  use shr_kind_mod,   only: r8=>shr_kind_r8
  use physconst,      only: cpair, Cpwater_vapor => cpwv, Rwater_vapor => rh2o, Rgas => rair
  use dimensions_mod, only : np, nlev
  implicit none
  
  private
  
  public :: Virtual_Temperature
  public :: Virtual_Specific_Heat

 interface Virtual_Temperature
    module procedure Virtual_Temperature1d
    module procedure Virtual_Temperature3d
 end interface


contains
  
  !===========================
  !
  ! For help or information:
  ! 
  ! Amik St-Cyr
  ! 
  ! e-mail: amik@ucar.edu
  !
  !===========================
 
  !================================
  ! For reference see Emanuel 1994 
  !================================
  
  function Virtual_Temperature1d(Tin,rin) result(Tv)
    
    real (kind=r8),intent(in) :: Tin
    real (kind=r8),intent(in) :: rin
    real (kind=r8)            :: Tv

!    Tv = Tin*(1_r8 + rin/Rd_on_Rv)/(1_r8 + rin)

    Tv = Tin*(1_r8 + (Rwater_vapor/Rgas - 1.0_r8)*rin)


  end function Virtual_Temperature1d

  function Virtual_Temperature3d(T,Q) result(T_v)
    real (kind=r8),intent(in) :: T(np,np,nlev)
    real (kind=r8),intent(in) :: Q(np,np,nlev)
    real (kind=r8) :: T_v(np,np,nlev)
    integer :: i, j, k

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j)
#endif
    do k=1,nlev
       do j=1,np
          do i=1,np
             T_v(i,j,k) = Virtual_Temperature1d(T(i,j,k), Q(i,j,k))
          end do
       end do
    end do
  end function Virtual_Temperature3d

  function Virtual_Specific_Heat(rin) result(Cp_star)
    
    real (kind=r8),intent(in) :: rin
    real (kind=r8)            :: Cp_star
 
    Cp_star = cpair*(1.0_r8 + (Cpwater_vapor/cpair - 1.0_r8)*rin)
   
  end function Virtual_Specific_Heat
     
  end module physics_mod
