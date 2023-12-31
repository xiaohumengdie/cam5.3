!-----------------------------------------------------------------------
! $Id: array_index.F90 6559 2013-10-15 20:59:28Z janhft@uwm.edu $
!-----------------------------------------------------------------------
module array_index

! Description:
!   Contains indices to variables in larger arrays.
!   Note that the 'ii' is necessary because 'i' is used in
!   statistics to track locations in the zt/zm/sfc derived types.

! References:
!   None
!-----------------------------------------------------------------------
  implicit none

  ! Variables
  ! Microphysics mixing ratios
  integer, public :: &
    iirrainm,    & ! Rain water mixing ratio  [kg/kg]
    iirsnowm,    & ! Snow mixing ratio        [kg/kg]
    iiricem,     & ! Ice mixing ratio         [kg/kg]
    iirgraupelm    ! Graupel mixing ratio     [kg/kg]
!$omp threadprivate(iirrainm, iirsnowm, iiricem, iirgraupelm)

  ! Microphysics concentrations
  integer, public :: &
    iiNrm,       & ! Rain drop concentration                       [num/kg]
    iiNsnowm,    & ! Snow concentration                            [num/kg]
    iiNim,       & ! Ice concentration                             [num/kg]
    iiNgraupelm, & ! Graupel concentration                         [num/kg]
    iiNcnm,      & ! Cloud nuclei concentration                    [num/kg]
    iiNcm          ! Cloud droplet concentration (not part of PDF) [num/kg]
!$omp threadprivate(iiNrm, iiNsnowm, iiNim, iiNgraupelm, iiNcnm,iiNcm)

  ! Scalar quantities
  integer, public :: & 
    iisclr_rt, iisclr_thl, iisclr_CO2, & ! [kg/kg]/[K]/[1e6 mol/mol]
    iiedsclr_rt, iiedsclr_thl, iiedsclr_CO2 ! "    "
!$omp threadprivate(iisclr_rt, iisclr_thl, iisclr_CO2, &
!$omp   iiedsclr_rt, iiedsclr_thl, iiedsclr_CO2)

  private ! Default Scope

end module array_index
!-----------------------------------------------------------------------
