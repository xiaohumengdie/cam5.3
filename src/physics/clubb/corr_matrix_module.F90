!$Id: corr_matrix_module.F90 6631 2013-12-10 17:53:20Z raut@uwm.edu $
!-------------------------------------------------------------------------------
module corr_matrix_module

  use clubb_precision, only: &
      core_rknd

  implicit none

  ! Latin hypercube indices / Correlation array indices
  integer, public :: &
    iiPDF_s_mellor = -1, &
    iiPDF_t_mellor = -1, &
    iiPDF_w        = -1
!$omp threadprivate(iiPDF_s_mellor, iiPDF_t_mellor, iiPDF_w)

  integer, public :: &
   iiPDF_rrain    = -1, &
   iiPDF_rsnow    = -1, &
   iiPDF_rice     = -1, &
   iiPDF_rgraupel = -1
!$omp threadprivate(iiPDF_rrain, iiPDF_rsnow, iiPDF_rice, iiPDF_rgraupel)

  integer, public :: &
   iiPDF_Nr       = -1, &
   iiPDF_Nsnow    = -1, &
   iiPDF_Ni       = -1, &
   iiPDF_Ngraupel = -1, &
   iiPDF_Ncn      = -1
!$omp threadprivate(iiPDF_Nr, iiPDF_Nsnow, iiPDF_Ni, iiPDF_Ngraupel, iiPDF_Ncn)

  integer, parameter, public :: &
    d_var_total = 12 ! Size of the default correlation arrays

  integer, public :: &
    d_variables
!$omp threadprivate(d_variables)

  real( kind = core_rknd ), public, dimension(:), allocatable :: &
    xp2_on_xm2_array_cloud, &
    xp2_on_xm2_array_below

  real( kind = core_rknd ), public, dimension(:,:), allocatable :: &
    corr_array_cloud, &
    corr_array_below
!$omp threadprivate(xp2_on_xm2_array_cloud, xp2_on_xm2_array_below, &
!$omp   corr_array_cloud, corr_array_below)

  real( kind = core_rknd ), public, dimension(:,:), allocatable :: &
      corr_array_cloud_def, &
      corr_array_below_def
!$omp threadprivate( corr_array_cloud_def, corr_array_below_def )


  private

  public :: read_correlation_matrix, setup_pdf_indices, setup_corr_varnce_array, &
            cleanup_corr_matrix_arrays, hm_idx, init_clubb_arrays

  private :: get_corr_var_index, return_pdf_index, def_corr_idx


  contains

  !-----------------------------------------------------------------------------
  subroutine init_default_corr_arrays(  ) 

    ! Description:
    !   Initializes the default correlation arrays with correlations from the 
    !   arm_97 case.
    !-----------------------------------------------------------------------------
  
    implicit none

    ! This "renaming" is used to shorten the matrix declarations below.
    integer, parameter :: c = core_rknd

    ! ---- Begin Code ----
 
    allocate( corr_array_cloud_def(d_var_total,d_var_total) )
    allocate( corr_array_below_def(d_var_total,d_var_total) )

    corr_array_cloud_def = reshape( &

(/ 1._c, .3_c, .09_c , .09_c , .242_c , .285_c , -.08_c , .28_c , .06_c , .04_c , 0._c, 0._c, &! s
   0._c, 1._c, .027_c, .027_c, .0726_c, .0855_c, -.024_c, .084_c, .018_c, .012_c, 0._c, 0._c, &! t
   0._c, 0._c, 1._c  , .34_c , 0._c   , 0._c   ,  .44_c , .55_c , .65_c , .73_c , 0._c, 0._c, &! w
   0._c, 0._c, 0._c  , 1._c  , 0._c   , 0._c   ,  .39_c , .29_c , .14_c , .21_c , 0._c, 0._c, &! Ncn
   0._c, 0._c, 0._c  , 0._c  , 1._c   , .768_c ,  0._c  , 0._c  , 0._c  , 0._c  , 0._c, 0._c, &! rr
   0._c, 0._c, 0._c  , 0._c  , 0._c   , 1._c   ,  0._c  , 0._c  , 0._c  , 0._c  , 0._c, 0._c, &! Nr
   0._c, 0._c, 0._c  , 0._c  , 0._c   , 0._c   ,  1._c  , .77_c , .29_c , .49_c , 0._c, 0._c, &! ri
   0._c, 0._c, 0._c  , 0._c  , 0._c   , 0._c   ,  0._c  , 1._c  , .43_c , .60_c , 0._c, 0._c, &! Ni
   0._c, 0._c, 0._c  , 0._c  , 0._c   , 0._c   ,  0._c  , 0._c  , 1._c  , .95_c , 0._c, 0._c, &! rs
   0._c, 0._c, 0._c  , 0._c  , 0._c   , 0._c   ,  0._c  , 0._c  , 0._c  , 1._c  , 0._c, 0._c, &! Ns
   0._c, 0._c, 0._c  , 0._c  , 0._c   , 0._c   ,  0._c  , 0._c  , 0._c  , 0._c  , 1._c, 0._c, &! rg
   0._c, 0._c, 0._c  , 0._c  , 0._c   , 0._c   ,  0._c  , 0._c  , 0._c  , 0._c  , 0._c, 1._c/),&!Ng

    shape(corr_array_cloud_def))
!  s     t     w       Ncn     rr       Nr        ri      Ni      rs      Ns      rg    Ng 

    corr_array_cloud_def = transpose( corr_array_cloud_def )


    corr_array_below_def = reshape( &

(/ 1._c, .3_c, .09_c , .09_c , .056_c , .015_c , -.08_c , .28_c , .06_c , .04_c , 0._c, 0._c, &! s
   0._c, 1._c, .027_c, .027_c, .168_c , .0045_c, -.024_c, .084_c, .018_c, .012_c, 0._c, 0._c, &! t
   0._c, 0._c, 1._c  , .34_c , 0._c   , 0._c   , .44_c  , .55_c , .65_c , .73_c , 0._c, 0._c, &! w
   0._c, 0._c, 0._c  , 1._c  , 0._c   , 0._c   , .39_c  , .29_c , .14_c , .21_c , 0._c, 0._c, &! Ncn
   0._c, 0._c, 0._c  , 0._c  , 1._c   , .886_c , 0._c   , 0._c  , 0._c  , 0._c  , 0._c, 0._c, &! rr
   0._c, 0._c, 0._c  , 0._c  , 0._c   , 1._c   , 0._c   , 0._c  , 0._c  , 0._c  , 0._c, 0._c, &! Nr
   0._c, 0._c, 0._c  , 0._c  , 0._c   , 0._c   , 1._c   , .77_c , .29_c , .49_c , 0._c, 0._c, &! ri
   0._c, 0._c, 0._c  , 0._c  , 0._c   , 0._c   , 0._c   , 1._c  , .43_c , .60_c , 0._c, 0._c, &! Ni
   0._c, 0._c, 0._c  , 0._c  , 0._c   , 0._c   , 0._c   , 0._c  , 1._c  , .95_c , 0._c, 0._c, &! rs
   0._c, 0._c, 0._c  , 0._c  , 0._c   , 0._c   , 0._c   , 0._c  , 0._c  , 1._c  , 0._c, 0._c, &! Ns
   0._c, 0._c, 0._c  , 0._c  , 0._c   , 0._c   , 0._c   , 0._c  , 0._c  , 0._c  , 1._c, 0._c, &! rg
   0._c, 0._c, 0._c  , 0._c  , 0._c   , 0._c   , 0._c   , 0._c  , 0._c  , 0._c  , 0._c, 1._c/),&!Ng

    shape(corr_array_below_def))
!  s     t     w       Ncn     rr       Nr       ri       Ni      rs      Ns      rg    Ng 


    corr_array_below_def = transpose( corr_array_below_def )

  end subroutine init_default_corr_arrays

  !-----------------------------------------------------------------------------
  pure function def_corr_idx( iiPDF_x ) result(ii_def_corr)

    ! Description:
    !   Map from a iiPDF index to the corresponding index in the default 
    !   correlation arrays.
    !-----------------------------------------------------------------------------

    implicit none

    ! Constant Parameters

    ! Indices that represent the order in the default corr arrays
    ! (s, t, w, Ncn, rr, Nr, ri, Ni, rs, Ns, rg, Ng)
    integer, parameter :: &
    ii_s = 1, &
    ii_t = 2, &
    ii_w = 3, &
    ii_Ncn = 4, &
    ii_rr = 5, &
    ii_Nr = 6, &
    ii_ri = 7, &
    ii_Ni = 8, &
    ii_rs = 9, &
    ii_Ns = 10, &
    ii_rg = 11, &
    ii_Ng = 12

    ! Input Variables

    integer, intent(in) :: iiPDF_x

    ! Return Variable

    integer :: ii_def_corr

    ! ---- Begin Code ----

    ii_def_corr = -1

      if (iiPDF_x == iiPDF_s_mellor) then
         ii_def_corr = ii_s

      elseif (iiPDF_x == iiPDF_t_mellor) then
        ii_def_corr = ii_t

      elseif (iiPDF_x == iiPDF_w) then
        ii_def_corr = ii_w

      elseif (iiPDF_x == iiPDF_Ncn) then
        ii_def_corr = ii_Ncn

      elseif (iiPDF_x == iiPDF_rrain) then
        ii_def_corr = ii_rr

      elseif (iiPDF_x == iiPDF_Nr) then
        ii_def_corr = ii_Nr

      elseif (iiPDF_x == iiPDF_rice) then
        ii_def_corr = ii_ri

      elseif (iiPDF_x == iiPDF_Ni) then
        ii_def_corr = ii_Ni

      elseif (iiPDF_x == iiPDF_rsnow) then
        ii_def_corr = ii_rs

      elseif (iiPDF_x == iiPDF_Nsnow) then
        ii_def_corr = ii_Ns

      elseif (iiPDF_x == iiPDF_rgraupel) then
        ii_def_corr = ii_rg

      elseif (iiPDF_x == iiPDF_Ngraupel) then
        ii_def_corr = ii_Ng

      endif
  end function def_corr_idx

  !-----------------------------------------------------------------------------
  subroutine set_corr_arrays_to_default(  ) 

    ! Description:
    !   If there are no corr_array.in files for the current case, default 
    !   correlations from arm_97 are used (e.g. in WRF_CLUBB). 
    !-----------------------------------------------------------------------------
  
    use constants_clubb, only: &
        zero, &
        one

    implicit none

    ! Local Variables
    integer :: i, j ! Loop iterators


    ! ---- Begin Code ----

    corr_array_cloud = zero
    corr_array_below = zero

    do i = 1, d_variables
       corr_array_cloud(i,i) = one
       corr_array_below(i,i) = one
    enddo

    do i = 1, d_variables-1
       do j = i+1, d_variables
          if ( def_corr_idx(i) > def_corr_idx(j) ) then
             corr_array_cloud(j, i) = corr_array_cloud_def(def_corr_idx(j), def_corr_idx(i))
             corr_array_below(j, i) = corr_array_below_def(def_corr_idx(j), def_corr_idx(i))
          else
             corr_array_cloud(j, i) = corr_array_cloud_def(def_corr_idx(i), def_corr_idx(j))
             corr_array_below(j, i) = corr_array_below_def(def_corr_idx(i), def_corr_idx(j))
          endif
       enddo
    enddo

  end subroutine set_corr_arrays_to_default


  !-----------------------------------------------------------------------------
  subroutine read_correlation_matrix( iunit, input_file, d_variables, &
                                      corr_array )

    ! Description:
    !   Reads a correlation variance array from a file and stores it in an array.
    !-----------------------------------------------------------------------------

    use input_reader, only: &
      one_dim_read_var, & ! Variable(s)
      read_one_dim_file, deallocate_one_dim_vars, count_columns ! Procedure(s)

    use matrix_operations, only: set_lower_triangular_matrix ! Procedure(s)

    use constants_clubb, only: fstderr ! Variable(s)

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: &
      iunit, &    ! File I/O unit
      d_variables ! number of variables in the array

    character(len=*), intent(in) :: input_file ! Path to the file

    ! Input/Output Variable(s)
    real( kind = core_rknd ), dimension(d_variables,d_variables), intent(inout) :: &
      corr_array ! Correlation variance array

    ! Local Variable(s)

    type(one_dim_read_var), allocatable, dimension(:) :: &
      retVars ! stores the variables read in from the corr_varnce.in file

    integer ::   &
      var_index1,    & ! variable index
      var_index2,    & ! variable index
      nCols,         & ! the number of columns in the file
      i, j         ! Loop index


    !--------------------------- BEGIN CODE -------------------------

    nCols = count_columns( iunit, input_file )

    ! Allocate all arrays based on d_variables
    allocate( retVars(1:nCols) )

    ! Initializing to zero means that correlations we don't have
    ! (e.g. Nc and any variable other than s_mellor ) are assumed to be 0.
    corr_array(:,:) = 0.0_core_rknd

    ! Set main diagonal to 1
    do i=1, d_variables
      corr_array(i,i) = 1.0_core_rknd
    end do

    ! Read the values from the specified file
    call read_one_dim_file( iunit, nCols, input_file, retVars )

    if( size( retVars(1)%values ) /= nCols ) then
      write(fstderr, *) "Correlation matrix must have an equal number of rows and cols in file ", &
            input_file
      stop "Bad data in correlation file."
    end if

    ! Start at 2 because the first index is always just 1.0 in the first row
    ! and the rest of the rows are ignored
    do i=2, nCols
      var_index1 = get_corr_var_index( retVars(i)%name )
      if( var_index1 > -1 ) then
        do j=1, (i-1)
          var_index2 = get_corr_var_index( retVars(j)%name )
          if( var_index2 > -1 ) then
            call set_lower_triangular_matrix &
                 ( d_variables, var_index1, var_index2, retVars(i)%values(j), &
                   corr_array )
          end if
        end do
      end if
    end do

    call deallocate_one_dim_vars( nCols, retVars )

    return
  end subroutine read_correlation_matrix

  !--------------------------------------------------------------------------
  function get_corr_var_index( var_name ) result( i )

    ! Definition:
    !   Returns the index for a variable based on its name.
    !--------------------------------------------------------------------------

    implicit none

    character(len=*), intent(in) :: var_name ! The name of the variable

    ! Output variable
    integer :: i

    !------------------ BEGIN CODE -----------------------------
    i = -1

    select case( trim(var_name) )

    case( "s" )
      i = iiPDF_s_mellor

    case( "t" )
      i = iiPDF_t_mellor

    case( "w" )
      i = iiPDF_w

    case( "Ncn" )
      i = iiPDF_Ncn

    case( "rrain" )
      i = iiPDF_rrain

    case( "Nr" )
      i = iiPDF_Nr

    case( "rice" )
      i = iiPDF_rice

    case( "Ni" )
      i = iiPDF_Ni

    case( "rsnow" )
      i = iiPDF_rsnow

    case( "Nsnow" )
      i = iiPDF_Nsnow

    end select

    return

  end function get_corr_var_index

  !-----------------------------------------------------------------------
  subroutine setup_pdf_indices( hydromet_dim, iirrainm, iiNrm, &
                                iiricem, iiNim, iirsnowm, iiNsnowm, &
                                l_ice_micro )

  ! Description:
  !
  !   Setup for the iiPDF indices. These indices are used to address s, t, w
  !   and the hydrometeors in the mean/stdev/corr arrays
  !
  ! References:
  !
  !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      hydromet_dim    ! Total number of hydrometeor species.

    integer, intent(in) :: &
      iirrainm, & ! Index of rain water mixing ratio
      iiNrm,    & ! Index of rain drop concentration
      iiricem,  & ! Index of ice mixing ratio
      iiNim,    & ! Index of ice crystal concentration
      iirsnowm, & ! Index of snow mixing ratio
      iiNsnowm    ! Index of snow concentration

    logical, intent(in) :: &
      l_ice_micro  ! Whether the microphysics scheme will do ice

    ! Local Variables
    integer :: &
      pdf_count, & ! Count number of PDF variables
      i            ! Hydrometeor loop index

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    iiPDF_s_mellor = 1 ! Extended liquid water mixing ratio
    iiPDF_t_mellor = 2 ! 't' orthogonal to 's'
    iiPDF_w        = 3 ! vertical velocity
    iiPDF_Ncn      = 4 ! Cloud nuclei concentration or extended Nc.

    pdf_count = iiPDF_Ncn

    ! Loop over hydrometeors.
    ! Hydrometeor indices in the PDF arrays should be in the same order as
    ! found in the hydrometeor arrays.
    if ( hydromet_dim > 0 ) then

       do i = 1, hydromet_dim, 1

          if ( i == iirrainm ) then
             pdf_count = pdf_count + 1
             iiPDF_rrain = pdf_count
          endif

          if ( i == iiNrm ) then
             pdf_count = pdf_count + 1
             iiPDF_Nr = pdf_count
          endif

          if ( l_ice_micro ) then

             if ( i == iiricem ) then
                pdf_count = pdf_count + 1
                iiPDF_rice = pdf_count
             endif

             if ( i == iiNim ) then
                pdf_count = pdf_count + 1
                iiPDF_Ni = pdf_count
             endif

             if ( i == iirsnowm ) then
                pdf_count = pdf_count + 1
                iiPDF_rsnow = pdf_count
             endif

             if ( i == iiNsnowm ) then
                pdf_count = pdf_count + 1
                iiPDF_Nsnow = pdf_count
             endif

          else

             iiPDF_rice = -1
             iiPDF_Ni = -1
             iiPDF_rsnow = -1
             iiPDF_Nsnow = -1

          endif ! l_ice_micro

       enddo ! i = 1, hydromet_dim, 1

    endif ! hydromet_dim > 0

    ! Disabled until we have values for the correlations of graupel and
    ! other variates in the latin hypercube sampling.
    iiPDF_rgraupel = -1
    iiPDF_Ngraupel = -1

    d_variables = pdf_count


    return

  end subroutine setup_pdf_indices
  !-----------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  subroutine return_pdf_index( hydromet_index, pdf_count, pdf_index )

  ! Description:
  !   Set the Latin hypercube variable index if the hydrometeor exists
  ! References:
  !   None
  !-------------------------------------------------------------------------

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      hydromet_index

    ! Input/Output Variables
    integer, intent(inout) :: &
      pdf_count

    ! Output Variables
    integer, intent(out) :: &
      pdf_index

    ! ---- Begin Code ----

    if ( hydromet_index > 0 ) then
      pdf_count = pdf_count + 1
      pdf_index = pdf_count
    else
      pdf_index = -1
    end if

    return
  end subroutine return_pdf_index

!===============================================================================
  subroutine setup_corr_varnce_array( input_file_cloud, input_file_below, &
                                      iunit )

! Description:
!   Setup an array with the x'^2/xm^2 variables on the diagonal and the other
!   elements to be correlations between various variables.

! References:
!   None.
!-------------------------------------------------------------------------------

    use parameters_microphys, only: &
      rrp2_on_rrm2_cloud,   & ! Variables
      Nrp2_on_Nrm2_cloud,   &
      Ncp2_on_Ncm2_cloud => Ncnp2_on_Ncnm2_cloud, &
      rrp2_on_rrm2_below,   &
      Nrp2_on_Nrm2_below,   &
      Ncp2_on_Ncm2_below => Ncnp2_on_Ncnm2_below

    use parameters_microphys, only: &
      rsnowp2_on_rsnowm2_cloud, & ! Variables
      Nsnowp2_on_Nsnowm2_cloud, &
      ricep2_on_ricem2_cloud, &
      Nicep2_on_Nicem2_cloud, &
      rsnowp2_on_rsnowm2_below, &
      Nsnowp2_on_Nsnowm2_below, &
      ricep2_on_ricem2_below, &
      Nicep2_on_Nicem2_below

    use parameters_microphys, only: &
      l_fix_s_t_correlations ! Variable(s)

!   use matrix_operations, only: print_lower_triangular_matrix ! Procedure(s)

    use constants_clubb, only: &
      fstdout, & ! Constant(s)
      fstderr, &
      zero

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: max, epsilon, trim

    ! Input Variables
    integer, intent(in) :: &
      iunit ! The file unit

    character(len=*), intent(in) :: &
      input_file_cloud, & ! Path to the in cloud correlation file
      input_file_below    ! Path to the out of cloud correlation file

    ! Local variables
    integer :: iiPDF_Nc

    character(len=1) :: response
    logical :: l_warning, corr_file_exist
    integer :: i

    ! ---- Begin Code ----

    iiPDF_Nc = iiPDF_Ncn

    allocate( corr_array_cloud(d_variables,d_variables) )
    allocate( corr_array_below(d_variables,d_variables) )

    allocate( xp2_on_xm2_array_cloud(d_variables) )
    allocate( xp2_on_xm2_array_below(d_variables) )

    xp2_on_xm2_array_cloud(:) = zero
    xp2_on_xm2_array_below(:) = zero

    ! corr_file_exist is true if the *_corr_array_cloud.in file exists
    ! Note: It is assumed that if the *_corr_array_cloud.in file exists
    !       then *_corr_array_below.in also exists
    inquire( file = input_file_cloud, exist = corr_file_exist )

    if ( corr_file_exist ) then

       call read_correlation_matrix( iunit, trim( input_file_cloud ), d_variables, & ! In
                                     corr_array_cloud ) ! Out

       call read_correlation_matrix( iunit, trim( input_file_below ), d_variables, & ! In
                                     corr_array_below ) ! Out

    else ! Read in default correlation matrices (from arm_97_corr_array_cloud)

       write(fstderr,*) "Warning: "//trim( input_file_cloud )//" was not found! " // &
                        "The default correlation arrays (hardwired from arm_97) will be used."

       call init_default_corr_arrays( )

       call set_corr_arrays_to_default( )

    endif

    ! Sanity check to avoid confusing non-convergence results.
    if ( .not. l_fix_s_t_correlations .and. iiPDF_Nc > 0 ) then
      l_warning = .false.
      do i = 1, d_variables
        if ( ( corr_array_cloud(i,iiPDF_Nc) /= zero .or.  &
               corr_array_below(i,iiPDF_Nc) /= zero ) .and. &
             i /= iiPDF_Nc ) then
          l_warning = .true.
        end if
      end do ! 1..d_variables
      if ( l_warning ) then
        write(fstderr,*) "Warning: the specified correlations for s and Nc are non-zero."
        write(fstderr,*) "The latin hypercube code will not converge to the analytic solution "// &
          "using these settings."
        write(fstderr,'(A)',advance='no') "Continue? "
        read(*,*) response
        if ( response(1:1) /= 'y' .and. response(1:1) /= 'Y' ) then
           stop "Exiting..."
        end if
      end if
    end if ! l_fix_s_t_correlations

    if ( iiPDF_Nc > 0 ) then
      xp2_on_xm2_array_cloud(iiPDF_Nc) = Ncp2_on_Ncm2_cloud
    end if

    if ( iiPDF_rrain > 0 ) then
      xp2_on_xm2_array_cloud(iiPDF_rrain) = rrp2_on_rrm2_cloud
      if ( iiPDF_Nr > 0 ) then
        xp2_on_xm2_array_cloud(iiPDF_Nr) = Nrp2_on_Nrm2_cloud
      end if ! iiPDF_Nr > 0
    end if ! iiPDF_rrain > 0

    if ( iiPDF_rsnow > 0 ) then
      xp2_on_xm2_array_cloud(iiPDF_rsnow) = rsnowp2_on_rsnowm2_cloud


      if ( iiPDF_Nsnow > 0 ) then
        xp2_on_xm2_array_cloud(iiPDF_Nsnow) = Nsnowp2_on_Nsnowm2_cloud


      end if ! iiPDF_Nsnow > 0
    end if ! iiPDF_rsnow > 0

    if ( iiPDF_rice > 0 ) then
      xp2_on_xm2_array_cloud(iiPDF_rice) = ricep2_on_ricem2_cloud


      if ( iiPDF_Ni > 0 ) then
        xp2_on_xm2_array_cloud(iiPDF_Ni) = Nicep2_on_Nicem2_cloud

      end if ! iiPDF_Ni > 0
    end if ! iiPDF_rice > 0

    ! Sampling for graupel (disabled)
    if ( iiPDF_rgraupel > 0 ) then
      xp2_on_xm2_array_cloud(iiPDF_rgraupel) = -999._core_rknd


      if ( iiPDF_Ngraupel > 0 ) then
        xp2_on_xm2_array_cloud(iiPDF_Ngraupel) = -999._core_rknd


      end if ! iiPDF_Ngraupel > 0
    end if ! iiPDF_rgraupel > 0

    if ( iiPDF_Nc > 0 ) then
      ! The epsilon is a kluge to prevent a singular matrix in generate_lh_sample
      xp2_on_xm2_array_below(iiPDF_Nc) = &
        max( Ncp2_on_Ncm2_below, epsilon( Ncp2_on_Ncm2_below ) )

    end if

    if ( iiPDF_rrain > 0 ) then
      xp2_on_xm2_array_below(iiPDF_rrain) = rrp2_on_rrm2_below



      if ( iiPDF_Nr > 0 ) then
        xp2_on_xm2_array_below(iiPDF_Nr) = Nrp2_on_Nrm2_below


      end if ! iiPDF_Nr > 0
    end if ! iiPDF_rrain > 0

    if ( iiPDF_rsnow > 0 ) then
      xp2_on_xm2_array_below(iiPDF_rsnow) = rsnowp2_on_rsnowm2_below


      if ( iiPDF_Nsnow > 0 ) then
        xp2_on_xm2_array_below(iiPDF_Nsnow) = Nsnowp2_on_Nsnowm2_below

      end if ! iiPDF_Nsnow > 0
    end if ! iiPDF_rsnow > 0

    if ( iiPDF_rice > 0 ) then
      xp2_on_xm2_array_below(iiPDF_rice) = ricep2_on_ricem2_below


      if ( iiPDF_Ni > 0 ) then
        xp2_on_xm2_array_below(iiPDF_Ni) =  Nicep2_on_Nicem2_below
      end if ! iiPDF_Ni > 0

    end if ! iiPDF_rice > 0

    if ( iiPDF_rgraupel > 0 ) then
      xp2_on_xm2_array_below(iiPDF_rgraupel) = -999._core_rknd


      if ( iiPDF_Ngraupel > 0 ) then
        xp2_on_xm2_array_below(iiPDF_Ngraupel) = -999._core_rknd


      end if ! iiPDF_Ngraupel > 0
    end if ! iiPDF_rgraupel > 0

    return
  end subroutine setup_corr_varnce_array

  !-----------------------------------------------------------------------------
  subroutine cleanup_corr_matrix_arrays( )

    ! Description:
    !   De-allocate latin hypercube arrays
    ! References:
    !   None
    !---------------------------------------------------------------------------
    implicit none

    ! External
    intrinsic :: allocated

    ! ---- Begin Code ----

    if ( allocated( corr_array_cloud ) ) then
      deallocate( corr_array_cloud )
    end if

    if ( allocated( corr_array_below ) ) then
      deallocate( corr_array_below )
    end if

    if ( allocated( xp2_on_xm2_array_cloud ) ) then
      deallocate( xp2_on_xm2_array_cloud )
    end if

    if ( allocated( xp2_on_xm2_array_below ) ) then
      deallocate( xp2_on_xm2_array_below )
    end if

    if ( allocated( corr_array_cloud_def ) ) then
      deallocate( corr_array_cloud_def )
    end if

    if ( allocated( corr_array_below_def ) ) then
      deallocate( corr_array_below_def )
    end if

    return
  end subroutine cleanup_corr_matrix_arrays

 !-----------------------------------------------------------------------------
  subroutine init_clubb_arrays( hydromet_dim, iirrainm, iiNrm, iirsnowm, & ! Variables
                                iiricem, iiNcm, iiNsnowm, iiNim, &
                                l_ice_micro, iunit )

    ! Description: This subroutine sets up arrays that are necessary for WRF.
    !   
    !-----------------------------------------------------------------------------

    implicit none

    ! Constant Parameter(s)
    character(len=*), parameter :: &
      LH_file_path_below = "./clubb_corr_array_below.in", &
      LH_file_path_cloud = "./clubb_corr_array_cloud.in"    

    ! Input Variables
    integer, intent(in) :: &
      hydromet_dim, &
      iirrainm, &
      iiNrm, &
      iirsnowm, & 
      iiricem, & 
      iiNcm, & 
      iiNsnowm, & 
      iiNim, &
      iunit

    logical, intent(in) :: l_ice_micro

    ! ---- Begin Code ----

    call setup_pdf_indices( hydromet_dim, iirrainm, iiNrm, &
                            iiricem, iiNim, iirsnowm, iiNsnowm, &
                            l_ice_micro )

    ! Setup the arrays and indices containing the correlations, etc.
    call setup_corr_varnce_array( LH_file_path_cloud, LH_file_path_below, iunit )

    return    

  end subroutine init_clubb_arrays

  !-----------------------------------------------------------------------
  function hm_idx(iiPDF_idx) result(ii_idx)
  ! Description:
  ! Returns the position of a certain hydrometeor within the hydromet arrays
  ! according to its iiPDF index.

  ! References:
  !
  !-----------------------------------------------------------------------

    use array_index, only: &
        iiNcnm, &
        iirrainm, &
        iiNrm, &
        iirsnowm, &
        iiNsnowm, &
        iiricem, &
        iiNim, &
        iirgraupelm, &
        iiNgraupelm

      implicit none

    ! Input Variables
    integer, intent(in) :: iiPDF_idx

    ! Return Variable
    integer :: ii_idx

    ! Local Variables

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    if ( iiPDF_idx == iiPDF_Ncn ) then
       ii_idx = iiNcnm
    endif

    if ( iiPDF_idx == iiPDF_rrain ) then
       ii_idx = iirrainm
    endif

    if ( iiPDF_idx == iiPDF_Nr ) then
       ii_idx = iiNrm
    endif

    if ( iiPDF_idx == iiPDF_rsnow ) then
       ii_idx = iirsnowm
    endif

    if ( iiPDF_idx == iiPDF_Nsnow ) then
       ii_idx = iiNsnowm
    endif

    if ( iiPDF_idx == iiPDF_rice ) then
       ii_idx = iiricem
    endif

    if ( iiPDF_idx == iiPDF_Ni ) then
       ii_idx = iiNim
    endif

    if ( iiPDF_idx == iiPDF_rgraupel ) then
       ii_idx = iirgraupelm
    endif

    if ( iiPDF_idx == iiPDF_Ngraupel ) then
       ii_idx = iiNgraupelm
    endif

    return

  end function hm_idx
  !-----------------------------------------------------------------------

end module corr_matrix_module
