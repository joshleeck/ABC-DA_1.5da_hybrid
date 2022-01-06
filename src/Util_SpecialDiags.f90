PROGRAM Util_SpecialDiags

!*****************************************************
!*   Code to make the following disgnostics          *
!*                                                   *
!* 1. Correlations between control parameters (real sp)
!* 2. Correlations between control parameters (spec sp)
!* 3. Correlations between control variables         *
!* 4. Ensemble average KE spectrum                   *
!* 5. Scale dependent balance diagnostics            *
!* 6. Variances of total, balanced, and unbalanced scaled dens
!*                                                   *
!*   Ross Bannister, NCEO, April/May 2020            *
!*                                                   *
!*****************************************************


USE DefConsTypes, ONLY :         &
    ZREAL8,                      &
    pi,                          &
    nlongs,                      &
    nlevs,                       &
    dx,                          &
    CVT_type,                    &
    ABC_type,                    &
    CV_type,                     &
    datadirCVT,                  &
    CVT_file,                    &
    NEnsMems,                    &
    Nlats,                       &
    Nens,                        &
    datadirABCperts,             &
    dims_type,                   &
    datadirABCfcs,               &
    fft_wsave_x,                 &
    fft_work_x,                  &
    LevMeanBalr


IMPLICIT NONE

! Declare variables
!==========================
INTEGER, PARAMETER        :: nscales = 20
LOGICAL                   :: DiagSwitch(6)
CHARACTER(LEN=320)        :: CVT_filename
INTEGER                   :: Neffmems, ens, item, mem, lat, orig_CVT_vert_opt_sym
REAL(ZREAL8)              :: Neffmems_r, Nens_r, Nlongs_r
CHARACTER(LEN=320)        :: meanfile
CHARACTER(LEN=320)        :: ABCfile
TYPE(dims_type)           :: dims
TYPE(CVT_type)            :: CVT
TYPE(ABC_type)            :: ABC_pert, ABC_mean, ABC_state
TYPE(CV_type)             :: ControlVar, ControlVarS
TYPE(CV_type)             :: Ctrl_Cor_with_1_inner, Ctrl_Cor_with_2_inner, Ctrl_Cor_with_3_inner
TYPE(CV_type)             :: Ctrl_Cor_with_4_inner, Ctrl_Cor_with_5_inner
TYPE(CV_type)             :: Ctrl_Cor_with_1_outer, Ctrl_Cor_with_2_outer, Ctrl_Cor_with_3_outer
TYPE(CV_type)             :: Ctrl_Cor_with_4_outer, Ctrl_Cor_with_5_outer
TYPE(CV_type)             :: StdDevs
REAL(ZREAL8), ALLOCATABLE :: data_fl(:), EnsAvKE(:,:)
INTEGER                   :: x, z, k, real_index, imag_index, last_index, scale
REAL(ZREAL8)              :: spatialscale
REAL(ZREAL8)              :: INT_HF
REAL(ZREAL8), ALLOCATABLE :: GeoCorr(:,:)
REAL(ZREAL8), ALLOCATABLE :: LBCorr(:,:)
REAL(ZREAL8), ALLOCATABLE :: HydCorr(:,:)
REAL(ZREAL8), ALLOCATABLE :: GeoCorr_outer(:,:)
REAL(ZREAL8), ALLOCATABLE :: LBCorr_outer(:,:)
REAL(ZREAL8), ALLOCATABLE :: HydCorr_outer(:,:)
REAL(ZREAL8), ALLOCATABLE :: GeoMean(:,:,:)
REAL(ZREAL8), ALLOCATABLE :: LBMean(:,:,:)
REAL(ZREAL8), ALLOCATABLE :: HydMean(:,:,:)
REAL(ZREAL8), ALLOCATABLE :: Geo(:,:,:)
REAL(ZREAL8), ALLOCATABLE :: LB(:,:,:)
REAL(ZREAL8), ALLOCATABLE :: Hyd(:,:,:)
REAL(ZREAL8), ALLOCATABLE :: GeoVar(:,:,:)
REAL(ZREAL8), ALLOCATABLE :: LBVar(:,:,:)
REAL(ZREAL8), ALLOCATABLE :: HydVar(:,:,:)
REAL(ZREAL8), ALLOCATABLE :: GeoCorr_inner(:,:)
REAL(ZREAL8), ALLOCATABLE :: LBCorr_inner(:,:)
REAL(ZREAL8), ALLOCATABLE :: HydCorr_inner(:,:)
INTEGER,      ALLOCATABLE :: Scales(:)


PRINT*, '*************************************************************************'
PRINT*, 'Running Util_SpecialDiags'
PRINT*, '*************************************************************************'

! Read namelist
CALL SetOptions

! Decide which diagnostics to compute
! -----------------------------------
! Real space correlations between control parameters
DiagSwitch(1) = .TRUE.
! Spectral space correlations between control parameters
DiagSwitch(2) = .TRUE.
! Actual control variables
DiagSwitch(3) = .TRUE.
! Compute average KE spectrum
DiagSwitch(4) = .TRUE.
! Find strength of balance correlation with scale
DiagSwitch(5) = .TRUE.
! Variances of scaled density (total, balanced, and unbalanced)
DiagSwitch(6) = .TRUE.



IF (DiagSwitch(1) .OR. DiagSwitch(2) .OR. DiagSwitch(3) .OR. DiagSwitch(6)) THEN
  PRINT *, 'Reading options from CVT file'
  CALL Initialise_CVT (CVT)
  CVT_filename = TRIM(datadirCVT) // '/' // TRIM(CVT_file)
  CALL Read_Covs (CVT_filename, CVT,           &
                  .TRUE.,                      &
                  .TRUE., .TRUE.,              &
                  .TRUE., .TRUE.)
  PRINT *, '-- done'
END IF



CALL Initialise_dims (dims)
Neffmems   = NEnsMems * Nlats
Neffmems_r = REAL(Neffmems)
Nens_r     = REAL(Nens)
Nlongs_r   = REAL(Nlongs)



IF (Diagswitch(1)) THEN
  PRINT *, '========================================================================='
  PRINT *, 'Performing real-space correlations between control parameters'
  PRINT *, '========================================================================='

  CALL Initialise_model_vars (ABC_mean, .FALSE.)
  CALL Initialise_model_vars (ABC_pert, .FALSE.)
  CALL Initialise_CVs(ControlVar, .FALSE.)
  CALL Initialise_CVs(StdDevs, .FALSE.)

  ! Initialise the correlations to zero
  CALL Initialise_CVs(Ctrl_Cor_with_1_outer, .FALSE.)
  CALL Initialise_CVs(Ctrl_Cor_with_2_outer, .FALSE.)
  CALL Initialise_CVs(Ctrl_Cor_with_3_outer, .FALSE.)
  CALL Initialise_CVs(Ctrl_Cor_with_4_outer, .FALSE.)
  CALL Initialise_CVs(Ctrl_Cor_with_5_outer, .FALSE.)


  PRINT *, 'There are ', Neffmems, ' effective ensemble members.'
  DO ens = 1, Nens

    ! Initialise the correlations to zero
    CALL Initialise_CVs(Ctrl_Cor_with_1_inner, .FALSE.)
    CALL Initialise_CVs(Ctrl_Cor_with_2_inner, .FALSE.)
    CALL Initialise_CVs(Ctrl_Cor_with_3_inner, .FALSE.)
    CALL Initialise_CVs(Ctrl_Cor_with_4_inner, .FALSE.)
    CALL Initialise_CVs(Ctrl_Cor_with_5_inner, .FALSE.)

    ! Read-in the ensemble mean for this ensemble
    WRITE (meanfile, '(A,A,I0.3,A)') TRIM(datadirABCperts), '/MeanABC', ens, '.nc'
    PRINT *, 'Reading mean state ', TRIM (meanfile)
    CALL Read_state_2d (meanfile, ABC_mean, dims, 1, .TRUE.)
    PRINT *, '-- done'

    ! Read-in the ABC perts and convert to parameter perts
    item = 0
    DO mem = 1, NEnsMems
      DO lat = 1, Nlats
        item = item + 1

        ! Read-in a pert found from stage 2
        WRITE (ABCfile, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirABCperts), '/PertABC_Ens', ens, '_Item', item, '.nc'
        PRINT *, 'Reading file ', TRIM(ABCfile)
        CALL Read_state_2d (ABCfile, ABC_pert, dims, 1, .TRUE.)
        PRINT *, '-- done'
        IF ((ens == 1) .AND. (item == 1)) THEN
          CALL Set_ht_dep_cons (dims)
        END IF

        ! Pass through the inverse parameter transform
        PRINT *, 'Performing inverse parameter transform'
        CALL U_p_inv (ABC_mean, ControlVar, ABC_pert, CVT % CVT_order,                        &
                      CVT % CVT_param_opt_gb, CVT % CVT_param_opt_hb, CVT % CVT_param_opt_ab,   &
                      CVT % CVT_param_opt_reg, LevMeanBalr,                                     &
                      CVT % Regression(1:nlevs,1:nlevs), dims,                                  &
                      ((ens == 1) .AND. (item ==1)),  &  ! Output diagnostics only if this is met
                      '.')                               ! For diagnostic location
        PRINT *, 'Done'


        ! ========================================================================
        ! Contribute to the covariance between control parameters in real-space
        PRINT *, 'Computing covariance contribution'
        ! psi with psi
        Ctrl_Cor_with_1_inner % v1(1:nlongs, 1:nlevs) = Ctrl_Cor_with_1_inner % v1(1:nlongs, 1:nlevs) +  &
          ControlVar % v1(1:nlongs, 1:nlevs) * ControlVar % v1(1:nlongs, 1:nlevs)
        ! psi with chi
        Ctrl_Cor_with_1_inner % v2(1:nlongs, 1:nlevs) = Ctrl_Cor_with_1_inner % v2(1:nlongs, 1:nlevs) +  &
          ControlVar % v1(1:nlongs, 1:nlevs) * ControlVar % v2(1:nlongs, 1:nlevs)
        ! psi with (unbalanced) r
        Ctrl_Cor_with_1_inner % v3(1:nlongs, 1:nlevs) = Ctrl_Cor_with_1_inner % v3(1:nlongs, 1:nlevs) +  &
          ControlVar % v1(1:nlongs, 1:nlevs) * ControlVar % v3(1:nlongs, 1:nlevs)
        ! psi with (unbalanced) b
        Ctrl_Cor_with_1_inner % v4(1:nlongs, 1:nlevs) = Ctrl_Cor_with_1_inner % v4(1:nlongs, 1:nlevs) +  &
          ControlVar % v1(1:nlongs, 1:nlevs) * ControlVar % v4(1:nlongs, 1:nlevs)
        ! psi with (unbalanced) w
        Ctrl_Cor_with_1_inner % v5(1:nlongs, 1:nlevs) = Ctrl_Cor_with_1_inner % v5(1:nlongs, 1:nlevs) +  &
          ControlVar % v1(1:nlongs, 1:nlevs) * ControlVar % v5(1:nlongs, 1:nlevs)

        ! chi with psi
        Ctrl_Cor_with_2_inner % v1(1:nlongs, 1:nlevs) = Ctrl_Cor_with_2_inner % v1(1:nlongs, 1:nlevs) +  &
          ControlVar % v2(1:nlongs, 1:nlevs) * ControlVar % v1(1:nlongs, 1:nlevs)
        ! chi with chi
        Ctrl_Cor_with_2_inner % v2(1:nlongs, 1:nlevs) = Ctrl_Cor_with_2_inner % v2(1:nlongs, 1:nlevs) +  &
          ControlVar % v2(1:nlongs, 1:nlevs) * ControlVar % v2(1:nlongs, 1:nlevs)
        ! chi with (unbalanced) r
        Ctrl_Cor_with_2_inner % v3(1:nlongs, 1:nlevs) = Ctrl_Cor_with_2_inner % v3(1:nlongs, 1:nlevs) +  &
          ControlVar % v2(1:nlongs, 1:nlevs) * ControlVar % v3(1:nlongs, 1:nlevs)
        ! chi with (unbalanced) b
        Ctrl_Cor_with_2_inner % v4(1:nlongs, 1:nlevs) = Ctrl_Cor_with_2_inner % v4(1:nlongs, 1:nlevs) +  &
          ControlVar % v2(1:nlongs, 1:nlevs) * ControlVar % v4(1:nlongs, 1:nlevs)
        ! chi with (unbalanced) w
        Ctrl_Cor_with_2_inner % v5(1:nlongs, 1:nlevs) = Ctrl_Cor_with_2_inner % v5(1:nlongs, 1:nlevs) +  &
          ControlVar % v2(1:nlongs, 1:nlevs) * ControlVar % v5(1:nlongs, 1:nlevs)

        ! (unbalanced) r with psi
        Ctrl_Cor_with_3_inner % v1(1:nlongs, 1:nlevs) = Ctrl_Cor_with_3_inner % v1(1:nlongs, 1:nlevs) +  &
          ControlVar % v3(1:nlongs, 1:nlevs) * ControlVar % v1(1:nlongs, 1:nlevs)
        ! (unbalanced) r with chi
        Ctrl_Cor_with_3_inner % v2(1:nlongs, 1:nlevs) = Ctrl_Cor_with_3_inner % v2(1:nlongs, 1:nlevs) +  &
          ControlVar % v3(1:nlongs, 1:nlevs) * ControlVar % v2(1:nlongs, 1:nlevs)
        ! (unbalanced) r with (unbalanced) r
        Ctrl_Cor_with_3_inner % v3(1:nlongs, 1:nlevs) = Ctrl_Cor_with_3_inner % v3(1:nlongs, 1:nlevs) +  &
          ControlVar % v3(1:nlongs, 1:nlevs) * ControlVar % v3(1:nlongs, 1:nlevs)
        ! (unbalanced) r with (unbalanced) b
        Ctrl_Cor_with_3_inner % v4(1:nlongs, 1:nlevs) = Ctrl_Cor_with_3_inner % v4(1:nlongs, 1:nlevs) +  &
          ControlVar % v3(1:nlongs, 1:nlevs) * ControlVar % v4(1:nlongs, 1:nlevs)
        ! (unbalanced) r with (unbalanced) w
        Ctrl_Cor_with_3_inner % v5(1:nlongs, 1:nlevs) = Ctrl_Cor_with_3_inner % v5(1:nlongs, 1:nlevs) +  &
          ControlVar % v3(1:nlongs, 1:nlevs) * ControlVar % v5(1:nlongs, 1:nlevs)

        ! (unbalanced) b with psi
        Ctrl_Cor_with_4_inner % v1(1:nlongs, 1:nlevs) = Ctrl_Cor_with_4_inner % v1(1:nlongs, 1:nlevs) +  &
          ControlVar % v4(1:nlongs, 1:nlevs) * ControlVar % v1(1:nlongs, 1:nlevs)
        ! (unbalanced) b with chi
        Ctrl_Cor_with_4_inner % v2(1:nlongs, 1:nlevs) = Ctrl_Cor_with_4_inner % v2(1:nlongs, 1:nlevs) +  &
          ControlVar % v4(1:nlongs, 1:nlevs) * ControlVar % v2(1:nlongs, 1:nlevs)
        ! (unbalanced) b with (unbalanced) r
        Ctrl_Cor_with_4_inner % v3(1:nlongs, 1:nlevs) = Ctrl_Cor_with_4_inner % v3(1:nlongs, 1:nlevs) +  &
          ControlVar % v4(1:nlongs, 1:nlevs) * ControlVar % v3(1:nlongs, 1:nlevs)
        ! (unbalanced) b with (unbalanced) b
        Ctrl_Cor_with_4_inner % v4(1:nlongs, 1:nlevs) = Ctrl_Cor_with_4_inner % v4(1:nlongs, 1:nlevs) +  &
          ControlVar % v4(1:nlongs, 1:nlevs) * ControlVar % v4(1:nlongs, 1:nlevs)
        ! (unbalanced) b with (unbalanced) w
        Ctrl_Cor_with_4_inner % v5(1:nlongs, 1:nlevs) = Ctrl_Cor_with_4_inner % v5(1:nlongs, 1:nlevs) +  &
          ControlVar % v4(1:nlongs, 1:nlevs) * ControlVar % v5(1:nlongs, 1:nlevs)

        ! (unbalanced) w with psi
        Ctrl_Cor_with_5_inner % v1(1:nlongs, 1:nlevs) = Ctrl_Cor_with_5_inner % v1(1:nlongs, 1:nlevs) +  &
          ControlVar % v5(1:nlongs, 1:nlevs) * ControlVar % v1(1:nlongs, 1:nlevs)
        ! (unbalanced) w with chi
        Ctrl_Cor_with_5_inner % v2(1:nlongs, 1:nlevs) = Ctrl_Cor_with_5_inner % v2(1:nlongs, 1:nlevs) +  &
          ControlVar % v5(1:nlongs, 1:nlevs) * ControlVar % v2(1:nlongs, 1:nlevs)
        ! (unbalanced) w with (unbalanced) r
        Ctrl_Cor_with_5_inner % v3(1:nlongs, 1:nlevs) = Ctrl_Cor_with_5_inner % v3(1:nlongs, 1:nlevs) +  &
          ControlVar % v5(1:nlongs, 1:nlevs) * ControlVar % v3(1:nlongs, 1:nlevs)
        ! (unbalanced) w with (unbalanced) b
        Ctrl_Cor_with_5_inner % v4(1:nlongs, 1:nlevs) = Ctrl_Cor_with_5_inner % v4(1:nlongs, 1:nlevs) +  &
          ControlVar % v5(1:nlongs, 1:nlevs) * ControlVar % v4(1:nlongs, 1:nlevs)
        ! (unbalanced) w with (unbalanced) w
        Ctrl_Cor_with_5_inner % v5(1:nlongs, 1:nlevs) = Ctrl_Cor_with_5_inner % v5(1:nlongs, 1:nlevs) +  &
          ControlVar % v5(1:nlongs, 1:nlevs) * ControlVar % v5(1:nlongs, 1:nlevs)
        PRINT *, 'Done'
        ! ========================================================================

      END DO
    END DO

    ! Normalise to give covariances
    PRINT *, 'Normalising'
    CALL Div_CV_cons(Ctrl_Cor_with_1_inner, Neffmems_r)
    CALL Div_CV_cons(Ctrl_Cor_with_2_inner, Neffmems_r)
    CALL Div_CV_cons(Ctrl_Cor_with_3_inner, Neffmems_r)
    CALL Div_CV_cons(Ctrl_Cor_with_4_inner, Neffmems_r)
    CALL Div_CV_cons(Ctrl_Cor_with_5_inner, Neffmems_r)
    PRINT *, 'Done'

    ! Compute the standard deviations
    PRINT *, 'Computing standard deviations'
    StdDevs % v1(1:nlongs, 1:nlevs) = SQRT(Ctrl_Cor_with_1_inner % v1(1:nlongs, 1:nlevs))
    StdDevs % v2(1:nlongs, 1:nlevs) = SQRT(Ctrl_Cor_with_2_inner % v2(1:nlongs, 1:nlevs))
    StdDevs % v3(1:nlongs, 1:nlevs) = SQRT(Ctrl_Cor_with_3_inner % v3(1:nlongs, 1:nlevs))
    StdDevs % v4(1:nlongs, 1:nlevs) = SQRT(Ctrl_Cor_with_4_inner % v4(1:nlongs, 1:nlevs))
    StdDevs % v5(1:nlongs, 1:nlevs) = SQRT(Ctrl_Cor_with_5_inner % v5(1:nlongs, 1:nlevs))
    PRINT *, 'Done'

    ! Compute correlations
    PRINT *, 'Computing correlations'
    Ctrl_Cor_with_1_inner % v1(1:nlongs, 1:nlevs) = Ctrl_Cor_with_1_inner % v1(1:nlongs, 1:nlevs) / &
      ( StdDevs % v1(1:nlongs, 1:nlevs) * StdDevs % v1(1:nlongs, 1:nlevs) )
    Ctrl_Cor_with_1_inner % v2(1:nlongs, 1:nlevs) = Ctrl_Cor_with_1_inner % v2(1:nlongs, 1:nlevs) / &
      ( StdDevs % v1(1:nlongs, 1:nlevs) * StdDevs % v2(1:nlongs, 1:nlevs) )
    Ctrl_Cor_with_1_inner % v3(1:nlongs, 1:nlevs) = Ctrl_Cor_with_1_inner % v3(1:nlongs, 1:nlevs) / &
      ( StdDevs % v1(1:nlongs, 1:nlevs) * StdDevs % v3(1:nlongs, 1:nlevs) )
    Ctrl_Cor_with_1_inner % v4(1:nlongs, 1:nlevs) = Ctrl_Cor_with_1_inner % v4(1:nlongs, 1:nlevs) / &
      ( StdDevs % v1(1:nlongs, 1:nlevs) * StdDevs % v4(1:nlongs, 1:nlevs) )
    Ctrl_Cor_with_1_inner % v5(1:nlongs, 1:nlevs) = Ctrl_Cor_with_1_inner % v5(1:nlongs, 1:nlevs) / &
      ( StdDevs % v1(1:nlongs, 1:nlevs) * StdDevs % v5(1:nlongs, 1:nlevs) )

    Ctrl_Cor_with_2_inner % v1(1:nlongs, 1:nlevs) = Ctrl_Cor_with_2_inner % v1(1:nlongs, 1:nlevs) / &
      ( StdDevs % v2(1:nlongs, 1:nlevs) * StdDevs % v1(1:nlongs, 1:nlevs) )
    Ctrl_Cor_with_2_inner % v2(1:nlongs, 1:nlevs) = Ctrl_Cor_with_2_inner % v2(1:nlongs, 1:nlevs) / &
      ( StdDevs % v2(1:nlongs, 1:nlevs) * StdDevs % v2(1:nlongs, 1:nlevs) )
    Ctrl_Cor_with_2_inner % v3(1:nlongs, 1:nlevs) = Ctrl_Cor_with_2_inner % v3(1:nlongs, 1:nlevs) / &
      ( StdDevs % v2(1:nlongs, 1:nlevs) * StdDevs % v3(1:nlongs, 1:nlevs) )
    Ctrl_Cor_with_2_inner % v4(1:nlongs, 1:nlevs) = Ctrl_Cor_with_2_inner % v4(1:nlongs, 1:nlevs) / &
      ( StdDevs % v2(1:nlongs, 1:nlevs) * StdDevs % v4(1:nlongs, 1:nlevs) )
    Ctrl_Cor_with_2_inner % v5(1:nlongs, 1:nlevs) = Ctrl_Cor_with_2_inner % v5(1:nlongs, 1:nlevs) / &
      ( StdDevs % v2(1:nlongs, 1:nlevs) * StdDevs % v5(1:nlongs, 1:nlevs) )

    Ctrl_Cor_with_3_inner % v1(1:nlongs, 1:nlevs) = Ctrl_Cor_with_3_inner % v1(1:nlongs, 1:nlevs) / &
      ( StdDevs % v3(1:nlongs, 1:nlevs) * StdDevs % v1(1:nlongs, 1:nlevs) )
    Ctrl_Cor_with_3_inner % v2(1:nlongs, 1:nlevs) = Ctrl_Cor_with_3_inner % v2(1:nlongs, 1:nlevs) / &
      ( StdDevs % v3(1:nlongs, 1:nlevs) * StdDevs % v2(1:nlongs, 1:nlevs) )
    Ctrl_Cor_with_3_inner % v3(1:nlongs, 1:nlevs) = Ctrl_Cor_with_3_inner % v3(1:nlongs, 1:nlevs) / &
      ( StdDevs % v3(1:nlongs, 1:nlevs) * StdDevs % v3(1:nlongs, 1:nlevs) )
    Ctrl_Cor_with_3_inner % v4(1:nlongs, 1:nlevs) = Ctrl_Cor_with_3_inner % v4(1:nlongs, 1:nlevs) / &
      ( StdDevs % v3(1:nlongs, 1:nlevs) * StdDevs % v4(1:nlongs, 1:nlevs) )
    Ctrl_Cor_with_3_inner % v5(1:nlongs, 1:nlevs) = Ctrl_Cor_with_3_inner % v5(1:nlongs, 1:nlevs) / &
      ( StdDevs % v3(1:nlongs, 1:nlevs) * StdDevs % v5(1:nlongs, 1:nlevs) )

    Ctrl_Cor_with_4_inner % v1(1:nlongs, 1:nlevs) = Ctrl_Cor_with_4_inner % v1(1:nlongs, 1:nlevs) / &
      ( StdDevs % v4(1:nlongs, 1:nlevs) * StdDevs % v1(1:nlongs, 1:nlevs) )
    Ctrl_Cor_with_4_inner % v2(1:nlongs, 1:nlevs) = Ctrl_Cor_with_4_inner % v2(1:nlongs, 1:nlevs) / &
      ( StdDevs % v4(1:nlongs, 1:nlevs) * StdDevs % v2(1:nlongs, 1:nlevs) )
    Ctrl_Cor_with_4_inner % v3(1:nlongs, 1:nlevs) = Ctrl_Cor_with_4_inner % v3(1:nlongs, 1:nlevs) / &
      ( StdDevs % v4(1:nlongs, 1:nlevs) * StdDevs % v3(1:nlongs, 1:nlevs) )
    Ctrl_Cor_with_4_inner % v4(1:nlongs, 1:nlevs) = Ctrl_Cor_with_4_inner % v4(1:nlongs, 1:nlevs) / &
      ( StdDevs % v4(1:nlongs, 1:nlevs) * StdDevs % v4(1:nlongs, 1:nlevs) )
    Ctrl_Cor_with_4_inner % v5(1:nlongs, 1:nlevs) = Ctrl_Cor_with_4_inner % v5(1:nlongs, 1:nlevs) / &
      ( StdDevs % v4(1:nlongs, 1:nlevs) * StdDevs % v5(1:nlongs, 1:nlevs) )

    Ctrl_Cor_with_5_inner % v1(1:nlongs, 1:nlevs) = Ctrl_Cor_with_5_inner % v1(1:nlongs, 1:nlevs) / &
      ( StdDevs % v5(1:nlongs, 1:nlevs) * StdDevs % v1(1:nlongs, 1:nlevs) )
    Ctrl_Cor_with_5_inner % v2(1:nlongs, 1:nlevs) = Ctrl_Cor_with_5_inner % v2(1:nlongs, 1:nlevs) / &
      ( StdDevs % v5(1:nlongs, 1:nlevs) * StdDevs % v2(1:nlongs, 1:nlevs) )
    Ctrl_Cor_with_5_inner % v3(1:nlongs, 1:nlevs) = Ctrl_Cor_with_5_inner % v3(1:nlongs, 1:nlevs) / &
      ( StdDevs % v5(1:nlongs, 1:nlevs) * StdDevs % v3(1:nlongs, 1:nlevs) )
    Ctrl_Cor_with_5_inner % v4(1:nlongs, 1:nlevs) = Ctrl_Cor_with_5_inner % v4(1:nlongs, 1:nlevs) / &
      ( StdDevs % v5(1:nlongs, 1:nlevs) * StdDevs % v4(1:nlongs, 1:nlevs) )
    Ctrl_Cor_with_5_inner % v5(1:nlongs, 1:nlevs) = Ctrl_Cor_with_5_inner % v5(1:nlongs, 1:nlevs) / &
      ( StdDevs % v5(1:nlongs, 1:nlevs) * StdDevs % v5(1:nlongs, 1:nlevs) )
    PRINT *, 'Done'


    ! Add on to the outer for averages
    PRINT *, 'Averaging over ensembles'
    CALL Add_CVs (Ctrl_Cor_with_1_outer, Ctrl_Cor_with_1_inner)
    CALL Add_CVs (Ctrl_Cor_with_2_outer, Ctrl_Cor_with_2_inner)
    CALL Add_CVs (Ctrl_Cor_with_3_outer, Ctrl_Cor_with_3_inner)
    CALL Add_CVs (Ctrl_Cor_with_4_outer, Ctrl_Cor_with_4_inner)
    CALL Add_CVs (Ctrl_Cor_with_5_outer, Ctrl_Cor_with_5_inner)
    PRINT *, 'Done'


  END DO


  ! Normalise correlations
  PRINT *, 'Normalising'
  CALL Div_CV_cons(Ctrl_Cor_with_1_outer, Nens_r)
  CALL Div_CV_cons(Ctrl_Cor_with_2_outer, Nens_r)
  CALL Div_CV_cons(Ctrl_Cor_with_3_outer, Nens_r)
  CALL Div_CV_cons(Ctrl_Cor_with_4_outer, Nens_r)
  CALL Div_CV_cons(Ctrl_Cor_with_5_outer, Nens_r)
  PRINT *, 'Done'


  ! Write out the correlations
  PRINT *, 'Outputting the correlations and standard deviations'
  CALL Write_CV ('RealSpWith_psi.nc', Ctrl_Cor_with_1_outer, 1, CVT,    &
                 dims % longs_v(1:nlongs), dims % half_levs(1:nlevs))
  CALL Write_CV ('RealSpWith_chi.nc', Ctrl_Cor_with_2_outer, 1, CVT,    &
                 dims % longs_v(1:nlongs), dims % half_levs(1:nlevs))
  CALL Write_CV ('RealSpWith_ur.nc',  Ctrl_Cor_with_3_outer, 1, CVT,    &
                 dims % longs_v(1:nlongs), dims % half_levs(1:nlevs))
  CALL Write_CV ('RealSpWith_ub.nc',  Ctrl_Cor_with_4_outer, 1, CVT,    &
                 dims % longs_v(1:nlongs), dims % half_levs(1:nlevs))
  CALL Write_CV ('RealSpWith_uw.nc',  Ctrl_Cor_with_5_outer, 1, CVT,    &
                 dims % longs_v(1:nlongs), dims % half_levs(1:nlevs))

  CALL Write_CV ('StdDevs.nc', StdDevs, 1, CVT,                         &
                 dims % longs_v(1:nlongs), dims % half_levs(1:nlevs))
  PRINT *, 'Done'

  ! Tidy up
  CALL Deallocate_model_vars (ABC_mean)
  CALL Deallocate_model_vars (ABC_pert)
  CALL Deallocate_CVs(ControlVar)
  CALL Deallocate_CVs(StdDevs)
  CALL Deallocate_CVs(Ctrl_Cor_with_1_outer)
  CALL Deallocate_CVs(Ctrl_Cor_with_2_outer)
  CALL Deallocate_CVs(Ctrl_Cor_with_3_outer)
  CALL Deallocate_CVs(Ctrl_Cor_with_4_outer)
  CALL Deallocate_CVs(Ctrl_Cor_with_5_outer)
  CALL Deallocate_CVs(Ctrl_Cor_with_1_inner)
  CALL Deallocate_CVs(Ctrl_Cor_with_2_inner)
  CALL Deallocate_CVs(Ctrl_Cor_with_3_inner)
  CALL Deallocate_CVs(Ctrl_Cor_with_4_inner)
  CALL Deallocate_CVs(Ctrl_Cor_with_5_inner)

END IF






IF (DiagSwitch(2)) THEN
  PRINT *, '========================================================================='
  PRINT *, 'Performing spectral-space correlations between control parameters'
  PRINT *, '========================================================================='

  CALL Initialise_model_vars (ABC_mean, .FALSE.)
  CALL Initialise_model_vars (ABC_pert, .FALSE.)
  CALL Initialise_CVs(ControlVar, .FALSE.)
  CALL Initialise_CVs(ControlVarS, .FALSE.)
  CALL Initialise_CVs(StdDevs, .FALSE.)


  ! Initialise the correlations to zero
  CALL Initialise_CVs(Ctrl_Cor_with_1_outer, .FALSE.)
  CALL Initialise_CVs(Ctrl_Cor_with_2_outer, .FALSE.)
  CALL Initialise_CVs(Ctrl_Cor_with_3_outer, .FALSE.)
  CALL Initialise_CVs(Ctrl_Cor_with_4_outer, .FALSE.)
  CALL Initialise_CVs(Ctrl_Cor_with_5_outer, .FALSE.)


  PRINT *, 'There are ', Neffmems, ' effective ensemble members.'
  DO ens = 1, Nens

    ! Initialise the correlations to zero
    CALL Initialise_CVs(Ctrl_Cor_with_1_inner, .FALSE.)
    CALL Initialise_CVs(Ctrl_Cor_with_2_inner, .FALSE.)
    CALL Initialise_CVs(Ctrl_Cor_with_3_inner, .FALSE.)
    CALL Initialise_CVs(Ctrl_Cor_with_4_inner, .FALSE.)
    CALL Initialise_CVs(Ctrl_Cor_with_5_inner, .FALSE.)

    ! Read-in the ensemble mean for this ensemble
    WRITE (meanfile, '(A,A,I0.3,A)') TRIM(datadirABCperts), '/MeanABC', ens, '.nc'
    PRINT *, 'Reading mean state ', TRIM (meanfile)
    CALL Read_state_2d (meanfile, ABC_mean, dims, 1, .TRUE.)
    PRINT *, '-- done'

    ! Read-in the ABC perts and convert to parameter perts
    item = 0
    DO mem = 1, NEnsMems
      DO lat = 1, Nlats
        item = item + 1

        ! Read-in a pert found from stage 2
        WRITE (ABCfile, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirABCperts), '/PertABC_Ens', ens, '_Item', item, '.nc'
        PRINT *, 'Reading file ', TRIM(ABCfile)
        CALL Read_state_2d (ABCfile, ABC_pert, dims, 1, .TRUE.)
        PRINT *, '-- done'
        IF ((ens == 1) .AND. (item == 1)) THEN
          CALL Set_ht_dep_cons (dims)
        END IF

        ! Pass through the inverse parameter transform
        PRINT *, 'Performing inverse parameter transform'
        CALL U_p_inv (ABC_mean, ControlVar, ABC_pert, CVT % CVT_order,                        &
                      CVT % CVT_param_opt_gb, CVT % CVT_param_opt_hb, CVT % CVT_param_opt_ab,   &
                      CVT % CVT_param_opt_reg, LevMeanBalr,                                     &
                      CVT % Regression(1:nlevs,1:nlevs), dims,                                  &
                      ((ens == 1) .AND. (item ==1)),  &  ! Output diagnostics only if this is met
                      '.')                               ! For diagnostic location
        PRINT *, 'Done'


        ! ========================================================================
        ! Fourier transform the control parameters
        PRINT *, 'Performing FT'
        CALL fft_real2spec ( ControlVar % v1(1:nlongs, 1:nlevs), ControlVarS % v1(1:nlongs, 1:nlevs) )
        CALL fft_real2spec ( ControlVar % v2(1:nlongs, 1:nlevs), ControlVarS % v2(1:nlongs, 1:nlevs) )
        CALL fft_real2spec ( ControlVar % v3(1:nlongs, 1:nlevs), ControlVarS % v3(1:nlongs, 1:nlevs) )
        CALL fft_real2spec ( ControlVar % v4(1:nlongs, 1:nlevs), ControlVarS % v4(1:nlongs, 1:nlevs) )
        CALL fft_real2spec ( ControlVar % v5(1:nlongs, 1:nlevs), ControlVarS % v5(1:nlongs, 1:nlevs) )
        PRINT *, 'Done'


        ! ========================================================================
        ! Contribute to the covariance between control parameters in spectral-space
        PRINT *, 'Computing covariance contribution'

        ! Deal with the largest scale first
        ! ---------------------------------
        ! psi with psi
        Ctrl_Cor_with_1_inner % v1(1, 1:nlevs) = Ctrl_Cor_with_1_inner % v1(1, 1:nlevs) +  &
          ControlVarS % v1(1, 1:nlevs) * ControlVarS % v1(1, 1:nlevs)
        ! psi with chi
        Ctrl_Cor_with_1_inner % v2(1, 1:nlevs) = Ctrl_Cor_with_1_inner % v2(1, 1:nlevs) +  &
          ControlVarS % v1(1, 1:nlevs) * ControlVarS % v2(1, 1:nlevs)
        ! psi with (unbalanced) r
        Ctrl_Cor_with_1_inner % v3(1, 1:nlevs) = Ctrl_Cor_with_1_inner % v3(1, 1:nlevs) +  &
          ControlVarS % v1(1, 1:nlevs) * ControlVarS % v3(1, 1:nlevs)
        ! psi with (unbalanced) b
        Ctrl_Cor_with_1_inner % v4(1, 1:nlevs) = Ctrl_Cor_with_1_inner % v4(1, 1:nlevs) +  &
          ControlVarS % v1(1, 1:nlevs) * ControlVarS % v4(1, 1:nlevs)
        ! psi with (unbalanced) w
        Ctrl_Cor_with_1_inner % v5(1, 1:nlevs) = Ctrl_Cor_with_1_inner % v5(1, 1:nlevs) +  &
          ControlVarS % v1(1, 1:nlevs) * ControlVarS % v5(1, 1:nlevs)

        ! chi with psi
        Ctrl_Cor_with_2_inner % v1(1, 1:nlevs) = Ctrl_Cor_with_2_inner % v1(1, 1:nlevs) +  &
          ControlVarS % v2(1, 1:nlevs) * ControlVarS % v1(1, 1:nlevs)
        ! chi with chi
        Ctrl_Cor_with_2_inner % v2(1, 1:nlevs) = Ctrl_Cor_with_2_inner % v2(1, 1:nlevs) +  &
          ControlVarS % v2(1, 1:nlevs) * ControlVarS % v2(1, 1:nlevs)
        ! chi with (unbalanced) r
        Ctrl_Cor_with_2_inner % v3(1, 1:nlevs) = Ctrl_Cor_with_2_inner % v3(1, 1:nlevs) +  &
          ControlVarS % v2(1, 1:nlevs) * ControlVarS % v3(1, 1:nlevs)
        ! chi with (unbalanced) b
        Ctrl_Cor_with_2_inner % v4(1, 1:nlevs) = Ctrl_Cor_with_2_inner % v4(1, 1:nlevs) +  &
          ControlVarS % v2(1, 1:nlevs) * ControlVarS % v4(1, 1:nlevs)
        ! chi with (unbalanced) w
        Ctrl_Cor_with_2_inner % v5(1, 1:nlevs) = Ctrl_Cor_with_2_inner % v5(1, 1:nlevs) +  &
          ControlVarS % v2(1, 1:nlevs) * ControlVarS % v5(1, 1:nlevs)

        ! (unbalanced) r with psi
        Ctrl_Cor_with_3_inner % v1(1, 1:nlevs) = Ctrl_Cor_with_3_inner % v1(1, 1:nlevs) +  &
          ControlVarS % v3(1, 1:nlevs) * ControlVarS % v1(1, 1:nlevs)
        ! (unbalanced) r with chi
        Ctrl_Cor_with_3_inner % v2(1, 1:nlevs) = Ctrl_Cor_with_3_inner % v2(1, 1:nlevs) +  &
          ControlVarS % v3(1, 1:nlevs) * ControlVarS % v2(1, 1:nlevs)
        ! (unbalanced) r with (unbalanced) r
        Ctrl_Cor_with_3_inner % v3(1, 1:nlevs) = Ctrl_Cor_with_3_inner % v3(1, 1:nlevs) +  &
          ControlVarS % v3(1, 1:nlevs) * ControlVarS % v3(1, 1:nlevs)
        ! (unbalanced) r with (unbalanced) b
        Ctrl_Cor_with_3_inner % v4(1, 1:nlevs) = Ctrl_Cor_with_3_inner % v4(1, 1:nlevs) +  &
          ControlVarS % v3(1, 1:nlevs) * ControlVarS % v4(1, 1:nlevs)
        ! (unbalanced) r with (unbalanced) w
        Ctrl_Cor_with_3_inner % v5(1, 1:nlevs) = Ctrl_Cor_with_3_inner % v5(1, 1:nlevs) +  &
          ControlVarS % v3(1, 1:nlevs) * ControlVarS % v5(1, 1:nlevs)

        ! (unbalanced) b with psi
        Ctrl_Cor_with_4_inner % v1(1, 1:nlevs) = Ctrl_Cor_with_4_inner % v1(1, 1:nlevs) +  &
          ControlVarS % v4(1, 1:nlevs) * ControlVarS % v1(1, 1:nlevs)
        ! (unbalanced) b with chi
        Ctrl_Cor_with_4_inner % v2(1, 1:nlevs) = Ctrl_Cor_with_4_inner % v2(1, 1:nlevs) +  &
          ControlVarS % v4(1, 1:nlevs) * ControlVarS % v2(1, 1:nlevs)
        ! (unbalanced) b with (unbalanced) r
        Ctrl_Cor_with_4_inner % v3(1, 1:nlevs) = Ctrl_Cor_with_4_inner % v3(1, 1:nlevs) +  &
          ControlVarS % v4(1, 1:nlevs) * ControlVarS % v3(1, 1:nlevs)
        ! (unbalanced) b with (unbalanced) b
        Ctrl_Cor_with_4_inner % v4(1, 1:nlevs) = Ctrl_Cor_with_4_inner % v4(1, 1:nlevs) +  &
          ControlVarS % v4(1, 1:nlevs) * ControlVarS % v4(1, 1:nlevs)
        ! (unbalanced) b with (unbalanced) w
        Ctrl_Cor_with_4_inner % v5(1, 1:nlevs) = Ctrl_Cor_with_4_inner % v5(1, 1:nlevs) +  &
          ControlVarS % v4(1, 1:nlevs) * ControlVarS % v5(1, 1:nlevs)

        ! (unbalanced) w with psi
        Ctrl_Cor_with_5_inner % v1(1, 1:nlevs) = Ctrl_Cor_with_5_inner % v1(1, 1:nlevs) +  &
          ControlVarS % v5(1, 1:nlevs) * ControlVarS % v1(1, 1:nlevs)
        ! (unbalanced) w with chi
        Ctrl_Cor_with_5_inner % v2(1, 1:nlevs) = Ctrl_Cor_with_5_inner % v2(1, 1:nlevs) +  &
          ControlVarS % v5(1, 1:nlevs) * ControlVarS % v2(1, 1:nlevs)
        ! (unbalanced) w with (unbalanced) r
        Ctrl_Cor_with_5_inner % v3(1, 1:nlevs) = Ctrl_Cor_with_5_inner % v3(1, 1:nlevs) +  &
          ControlVarS % v5(1, 1:nlevs) * ControlVarS % v3(1, 1:nlevs)
        ! (unbalanced) w with (unbalanced) b
        Ctrl_Cor_with_5_inner % v4(1, 1:nlevs) = Ctrl_Cor_with_5_inner % v4(1, 1:nlevs) +  &
          ControlVarS % v5(1, 1:nlevs) * ControlVarS % v4(1, 1:nlevs)
        ! (unbalanced) w with (unbalanced) w
        Ctrl_Cor_with_5_inner % v5(1, 1:nlevs) = Ctrl_Cor_with_5_inner % v5(1, 1:nlevs) +  &
          ControlVarS % v5(1, 1:nlevs) * ControlVarS % v5(1, 1:nlevs)

        ! Deal with the bulk of the scales
        ! --------------------------------
        DO k = 2, nlongs/2
          real_index = 2*k-2
          imag_index = 2*k-1

          ! psi with psi
          Ctrl_Cor_with_1_inner % v1(k, 1:nlevs) = Ctrl_Cor_with_1_inner % v1(k, 1:nlevs) +  &
            ControlVarS % v1(real_index, 1:nlevs) * ControlVarS % v1(real_index, 1:nlevs) +  &
            ControlVarS % v1(imag_index, 1:nlevs) * ControlVarS % v1(imag_index, 1:nlevs)
          ! psi with chi
          Ctrl_Cor_with_1_inner % v2(k, 1:nlevs) = Ctrl_Cor_with_1_inner % v2(k, 1:nlevs) +  &
            ControlVarS % v1(real_index, 1:nlevs) * ControlVarS % v2(real_index, 1:nlevs) +  &
            ControlVarS % v1(imag_index, 1:nlevs) * ControlVarS % v2(imag_index, 1:nlevs)
          ! psi with (unbalanced) r
          Ctrl_Cor_with_1_inner % v3(k, 1:nlevs) = Ctrl_Cor_with_1_inner % v3(k, 1:nlevs) +  &
            ControlVarS % v1(real_index, 1:nlevs) * ControlVarS % v3(real_index, 1:nlevs) +  &
            ControlVarS % v1(imag_index, 1:nlevs) * ControlVarS % v3(imag_index, 1:nlevs)
          ! psi with (unbalanced) b
          Ctrl_Cor_with_1_inner % v4(k, 1:nlevs) = Ctrl_Cor_with_1_inner % v4(k, 1:nlevs) +  &
            ControlVarS % v1(real_index, 1:nlevs) * ControlVarS % v4(real_index, 1:nlevs) +  &
            ControlVarS % v1(imag_index, 1:nlevs) * ControlVarS % v4(imag_index, 1:nlevs)
          ! psi with (unbalanced) w
          Ctrl_Cor_with_1_inner % v5(k, 1:nlevs) = Ctrl_Cor_with_1_inner % v5(k, 1:nlevs) +  &
            ControlVarS % v1(real_index, 1:nlevs) * ControlVarS % v5(real_index, 1:nlevs) +  &
            ControlVarS % v1(imag_index, 1:nlevs) * ControlVarS % v5(imag_index, 1:nlevs)

          ! chi with psi
          Ctrl_Cor_with_2_inner % v1(k, 1:nlevs) = Ctrl_Cor_with_2_inner % v1(k, 1:nlevs) +  &
            ControlVarS % v2(real_index, 1:nlevs) * ControlVarS % v1(real_index, 1:nlevs) +  &
            ControlVarS % v2(imag_index, 1:nlevs) * ControlVarS % v1(imag_index, 1:nlevs)
          ! chi with chi
          Ctrl_Cor_with_2_inner % v2(k, 1:nlevs) = Ctrl_Cor_with_2_inner % v2(k, 1:nlevs) +  &
            ControlVarS % v2(real_index, 1:nlevs) * ControlVarS % v2(real_index, 1:nlevs) +  &
            ControlVarS % v2(imag_index, 1:nlevs) * ControlVarS % v2(imag_index, 1:nlevs)
          ! chi with (unbalanced) r
          Ctrl_Cor_with_2_inner % v3(k, 1:nlevs) = Ctrl_Cor_with_2_inner % v3(k, 1:nlevs) +  &
            ControlVarS % v2(real_index, 1:nlevs) * ControlVarS % v3(real_index, 1:nlevs) +  &
            ControlVarS % v2(imag_index, 1:nlevs) * ControlVarS % v3(imag_index, 1:nlevs)
          ! chi with (unbalanced) b
          Ctrl_Cor_with_2_inner % v4(k, 1:nlevs) = Ctrl_Cor_with_2_inner % v4(k, 1:nlevs) +  &
            ControlVarS % v2(real_index, 1:nlevs) * ControlVarS % v4(real_index, 1:nlevs) +  &
            ControlVarS % v2(imag_index, 1:nlevs) * ControlVarS % v4(imag_index, 1:nlevs)
          ! chi with (unbalanced) w
          Ctrl_Cor_with_2_inner % v5(k, 1:nlevs) = Ctrl_Cor_with_2_inner % v5(k, 1:nlevs) +  &
            ControlVarS % v2(real_index, 1:nlevs) * ControlVarS % v5(real_index, 1:nlevs) +  &
            ControlVarS % v2(imag_index, 1:nlevs) * ControlVarS % v5(imag_index, 1:nlevs)

          ! (unbalanced) r with psi
          Ctrl_Cor_with_3_inner % v1(k, 1:nlevs) = Ctrl_Cor_with_3_inner % v1(k, 1:nlevs) +  &
            ControlVarS % v3(real_index, 1:nlevs) * ControlVarS % v1(real_index, 1:nlevs) +  &
            ControlVarS % v3(imag_index, 1:nlevs) * ControlVarS % v1(imag_index, 1:nlevs)
          ! (unbalanced) r with chi
          Ctrl_Cor_with_3_inner % v2(k, 1:nlevs) = Ctrl_Cor_with_3_inner % v2(k, 1:nlevs) +  &
            ControlVarS % v3(real_index, 1:nlevs) * ControlVarS % v2(real_index, 1:nlevs) +  &
            ControlVarS % v3(imag_index, 1:nlevs) * ControlVarS % v2(imag_index, 1:nlevs)
          ! (unbalanced) r with (unbalanced) r
          Ctrl_Cor_with_3_inner % v3(k, 1:nlevs) = Ctrl_Cor_with_3_inner % v3(k, 1:nlevs) +  &
            ControlVarS % v3(real_index, 1:nlevs) * ControlVarS % v3(real_index, 1:nlevs) +  &
            ControlVarS % v3(imag_index, 1:nlevs) * ControlVarS % v3(imag_index, 1:nlevs)
          ! (unbalanced) r with (unbalanced) b
          Ctrl_Cor_with_3_inner % v4(k, 1:nlevs) = Ctrl_Cor_with_3_inner % v4(k, 1:nlevs) +  &
            ControlVarS % v3(real_index, 1:nlevs) * ControlVarS % v4(real_index, 1:nlevs) +  &
            ControlVarS % v3(imag_index, 1:nlevs) * ControlVarS % v4(imag_index, 1:nlevs)
          ! (unbalanced) r with (unbalanced) w
          Ctrl_Cor_with_3_inner % v5(k, 1:nlevs) = Ctrl_Cor_with_3_inner % v5(k, 1:nlevs) +  &
            ControlVarS % v3(real_index, 1:nlevs) * ControlVarS % v5(real_index, 1:nlevs) +  &
            ControlVarS % v3(imag_index, 1:nlevs) * ControlVarS % v5(imag_index, 1:nlevs)

          ! (unbalanced) b with psi
          Ctrl_Cor_with_4_inner % v1(k, 1:nlevs) = Ctrl_Cor_with_4_inner % v1(k, 1:nlevs) +  &
            ControlVarS % v4(real_index, 1:nlevs) * ControlVarS % v1(real_index, 1:nlevs) +  &
            ControlVarS % v4(imag_index, 1:nlevs) * ControlVarS % v1(imag_index, 1:nlevs)
          ! (unbalanced) b with chi
          Ctrl_Cor_with_4_inner % v2(k, 1:nlevs) = Ctrl_Cor_with_4_inner % v2(k, 1:nlevs) +  &
            ControlVarS % v4(real_index, 1:nlevs) * ControlVarS % v2(real_index, 1:nlevs) +  &
            ControlVarS % v4(imag_index, 1:nlevs) * ControlVarS % v2(imag_index, 1:nlevs)
          ! (unbalanced) b with (unbalanced) r
          Ctrl_Cor_with_4_inner % v3(k, 1:nlevs) = Ctrl_Cor_with_4_inner % v3(k, 1:nlevs) +  &
            ControlVarS % v4(real_index, 1:nlevs) * ControlVarS % v3(real_index, 1:nlevs) +  &
            ControlVarS % v4(imag_index, 1:nlevs) * ControlVarS % v3(imag_index, 1:nlevs)
          ! (unbalanced) b with (unbalanced) b
          Ctrl_Cor_with_4_inner % v4(k, 1:nlevs) = Ctrl_Cor_with_4_inner % v4(k, 1:nlevs) +  &
            ControlVarS % v4(real_index, 1:nlevs) * ControlVarS % v4(real_index, 1:nlevs) +  &
            ControlVarS % v4(imag_index, 1:nlevs) * ControlVarS % v4(imag_index, 1:nlevs)
          ! (unbalanced) b with (unbalanced) w
          Ctrl_Cor_with_4_inner % v5(k, 1:nlevs) = Ctrl_Cor_with_4_inner % v5(k, 1:nlevs) +  &
            ControlVarS % v4(real_index, 1:nlevs) * ControlVarS % v5(real_index, 1:nlevs) +  &
            ControlVarS % v4(imag_index, 1:nlevs) * ControlVarS % v5(imag_index, 1:nlevs)

          ! (unbalanced) w with psi
          Ctrl_Cor_with_5_inner % v1(k, 1:nlevs) = Ctrl_Cor_with_5_inner % v1(k, 1:nlevs) +  &
            ControlVarS % v5(real_index, 1:nlevs) * ControlVarS % v1(real_index, 1:nlevs) +  &
            ControlVarS % v5(imag_index, 1:nlevs) * ControlVarS % v1(imag_index, 1:nlevs)
          ! (unbalanced) w with chi
          Ctrl_Cor_with_5_inner % v2(k, 1:nlevs) = Ctrl_Cor_with_5_inner % v2(k, 1:nlevs) +  &
            ControlVarS % v5(real_index, 1:nlevs) * ControlVarS % v2(real_index, 1:nlevs) +  &
            ControlVarS % v5(imag_index, 1:nlevs) * ControlVarS % v2(imag_index, 1:nlevs)
          ! (unbalanced) w with (unbalanced) r
          Ctrl_Cor_with_5_inner % v3(k, 1:nlevs) = Ctrl_Cor_with_5_inner % v3(k, 1:nlevs) +  &
            ControlVarS % v5(real_index, 1:nlevs) * ControlVarS % v3(real_index, 1:nlevs) +  &
            ControlVarS % v5(imag_index, 1:nlevs) * ControlVarS % v3(imag_index, 1:nlevs)
          ! (unbalanced) w with (unbalanced) b
          Ctrl_Cor_with_5_inner % v4(k, 1:nlevs) = Ctrl_Cor_with_5_inner % v4(k, 1:nlevs) +  &
            ControlVarS % v5(real_index, 1:nlevs) * ControlVarS % v4(real_index, 1:nlevs) +  &
            ControlVarS % v5(imag_index, 1:nlevs) * ControlVarS % v4(imag_index, 1:nlevs)
          ! (unbalanced) w with (unbalanced) w
          Ctrl_Cor_with_5_inner % v5(k, 1:nlevs) = Ctrl_Cor_with_5_inner % v5(k, 1:nlevs) +  &
            ControlVarS % v5(real_index, 1:nlevs) * ControlVarS % v5(real_index, 1:nlevs) +  &
            ControlVarS % v5(imag_index, 1:nlevs) * ControlVarS % v5(imag_index, 1:nlevs)
        END DO

        ! Deal with the smallest scale
        ! ----------------------------
        last_index = nlongs/2+1

        ! psi with psi
        Ctrl_Cor_with_1_inner % v1(last_index, 1:nlevs) = Ctrl_Cor_with_1_inner % v1(last_index, 1:nlevs) +  &
          ControlVarS % v1(nlongs, 1:nlevs) * ControlVarS % v1(nlongs, 1:nlevs)
        ! psi with chi
        Ctrl_Cor_with_1_inner % v2(last_index, 1:nlevs) = Ctrl_Cor_with_1_inner % v2(last_index, 1:nlevs) +  &
          ControlVarS % v1(nlongs, 1:nlevs) * ControlVarS % v2(nlongs, 1:nlevs)
        ! psi with (unbalanced) r
        Ctrl_Cor_with_1_inner % v3(last_index, 1:nlevs) = Ctrl_Cor_with_1_inner % v3(last_index, 1:nlevs) +  &
          ControlVarS % v1(nlongs, 1:nlevs) * ControlVarS % v3(nlongs, 1:nlevs)
        ! psi with (unbalanced) b
        Ctrl_Cor_with_1_inner % v4(last_index, 1:nlevs) = Ctrl_Cor_with_1_inner % v4(last_index, 1:nlevs) +  &
          ControlVarS % v1(nlongs, 1:nlevs) * ControlVarS % v4(nlongs, 1:nlevs)
        ! psi with (unbalanced) w
        Ctrl_Cor_with_1_inner % v5(last_index, 1:nlevs) = Ctrl_Cor_with_1_inner % v5(last_index, 1:nlevs) +  &
          ControlVarS % v1(nlongs, 1:nlevs) * ControlVarS % v5(nlongs, 1:nlevs)

        ! chi with psi
        Ctrl_Cor_with_2_inner % v1(last_index, 1:nlevs) = Ctrl_Cor_with_2_inner % v1(last_index, 1:nlevs) +  &
          ControlVarS % v2(nlongs, 1:nlevs) * ControlVarS % v1(nlongs, 1:nlevs)
        ! chi with chi
        Ctrl_Cor_with_2_inner % v2(last_index, 1:nlevs) = Ctrl_Cor_with_2_inner % v2(last_index, 1:nlevs) +  &
          ControlVarS % v2(nlongs, 1:nlevs) * ControlVarS % v2(nlongs, 1:nlevs)
        ! chi with (unbalanced) r
        Ctrl_Cor_with_2_inner % v3(last_index, 1:nlevs) = Ctrl_Cor_with_2_inner % v3(last_index, 1:nlevs) +  &
          ControlVarS % v2(nlongs, 1:nlevs) * ControlVarS % v3(nlongs, 1:nlevs)
        ! chi with (unbalanced) b
        Ctrl_Cor_with_2_inner % v4(last_index, 1:nlevs) = Ctrl_Cor_with_2_inner % v4(last_index, 1:nlevs) +  &
          ControlVarS % v2(nlongs, 1:nlevs) * ControlVarS % v4(nlongs, 1:nlevs)
        ! chi with (unbalanced) w
        Ctrl_Cor_with_2_inner % v5(last_index, 1:nlevs) = Ctrl_Cor_with_2_inner % v5(last_index, 1:nlevs) +  &
          ControlVarS % v2(nlongs, 1:nlevs) * ControlVarS % v5(nlongs, 1:nlevs)

        ! (unbalanced) r with psi
        Ctrl_Cor_with_3_inner % v1(last_index, 1:nlevs) = Ctrl_Cor_with_3_inner % v1(last_index, 1:nlevs) +  &
          ControlVarS % v3(nlongs, 1:nlevs) * ControlVarS % v1(nlongs, 1:nlevs)
        ! (unbalanced) r with chi
        Ctrl_Cor_with_3_inner % v2(last_index, 1:nlevs) = Ctrl_Cor_with_3_inner % v2(last_index, 1:nlevs) +  &
          ControlVarS % v3(nlongs, 1:nlevs) * ControlVarS % v2(nlongs, 1:nlevs)
        ! (unbalanced) r with (unbalanced) r
        Ctrl_Cor_with_3_inner % v3(last_index, 1:nlevs) = Ctrl_Cor_with_3_inner % v3(last_index, 1:nlevs) +  &
          ControlVarS % v3(nlongs, 1:nlevs) * ControlVarS % v3(nlongs, 1:nlevs)
        ! (unbalanced) r with (unbalanced) b
        Ctrl_Cor_with_3_inner % v4(last_index, 1:nlevs) = Ctrl_Cor_with_3_inner % v4(last_index, 1:nlevs) +  &
          ControlVarS % v3(nlongs, 1:nlevs) * ControlVarS % v4(nlongs, 1:nlevs)
        ! (unbalanced) r with (unbalanced) w
        Ctrl_Cor_with_3_inner % v5(last_index, 1:nlevs) = Ctrl_Cor_with_3_inner % v5(last_index, 1:nlevs) +  &
          ControlVarS % v3(nlongs, 1:nlevs) * ControlVarS % v5(nlongs, 1:nlevs)

        ! (unbalanced) b with psi
        Ctrl_Cor_with_4_inner % v1(last_index, 1:nlevs) = Ctrl_Cor_with_4_inner % v1(last_index, 1:nlevs) +  &
          ControlVarS % v4(nlongs, 1:nlevs) * ControlVarS % v1(nlongs, 1:nlevs)
        ! (unbalanced) b with chi
        Ctrl_Cor_with_4_inner % v2(last_index, 1:nlevs) = Ctrl_Cor_with_4_inner % v2(last_index, 1:nlevs) +  &
          ControlVarS % v4(nlongs, 1:nlevs) * ControlVarS % v2(nlongs, 1:nlevs)
        ! (unbalanced) b with (unbalanced) r
        Ctrl_Cor_with_4_inner % v3(last_index, 1:nlevs) = Ctrl_Cor_with_4_inner % v3(last_index, 1:nlevs) +  &
          ControlVarS % v4(nlongs, 1:nlevs) * ControlVarS % v3(nlongs, 1:nlevs)
        ! (unbalanced) b with (unbalanced) b
        Ctrl_Cor_with_4_inner % v4(last_index, 1:nlevs) = Ctrl_Cor_with_4_inner % v4(last_index, 1:nlevs) +  &
          ControlVarS % v4(nlongs, 1:nlevs) * ControlVarS % v4(nlongs, 1:nlevs)
        ! (unbalanced) b with (unbalanced) w
        Ctrl_Cor_with_4_inner % v5(last_index, 1:nlevs) = Ctrl_Cor_with_4_inner % v5(last_index, 1:nlevs) +  &
          ControlVarS % v4(nlongs, 1:nlevs) * ControlVarS % v5(nlongs, 1:nlevs)

        ! (unbalanced) w with psi
        Ctrl_Cor_with_5_inner % v1(last_index, 1:nlevs) = Ctrl_Cor_with_5_inner % v1(last_index, 1:nlevs) +  &
          ControlVarS % v5(nlongs, 1:nlevs) * ControlVarS % v1(nlongs, 1:nlevs)
        ! (unbalanced) w with chi
        Ctrl_Cor_with_5_inner % v2(last_index, 1:nlevs) = Ctrl_Cor_with_5_inner % v2(last_index, 1:nlevs) +  &
          ControlVarS % v5(nlongs, 1:nlevs) * ControlVarS % v2(nlongs, 1:nlevs)
        ! (unbalanced) w with (unbalanced) r
        Ctrl_Cor_with_5_inner % v3(last_index, 1:nlevs) = Ctrl_Cor_with_5_inner % v3(last_index, 1:nlevs) +  &
          ControlVarS % v5(nlongs, 1:nlevs) * ControlVarS % v3(nlongs, 1:nlevs)
        ! (unbalanced) w with (unbalanced) b
        Ctrl_Cor_with_5_inner % v4(last_index, 1:nlevs) = Ctrl_Cor_with_5_inner % v4(last_index, 1:nlevs) +  &
          ControlVarS % v5(nlongs, 1:nlevs) * ControlVarS % v4(nlongs, 1:nlevs)
        ! (unbalanced) w with (unbalanced) w
        Ctrl_Cor_with_5_inner % v5(last_index, 1:nlevs) = Ctrl_Cor_with_5_inner % v5(last_index, 1:nlevs) +  &
          ControlVarS % v5(nlongs, 1:nlevs) * ControlVarS % v5(nlongs, 1:nlevs)
        PRINT *, 'Done'
        ! ========================================================================

      END DO
    END DO

    ! Normalise to give covariances
    PRINT *, 'Normalising'
    CALL Div_CV_cons(Ctrl_Cor_with_1_inner, Neffmems_r)
    CALL Div_CV_cons(Ctrl_Cor_with_2_inner, Neffmems_r)
    CALL Div_CV_cons(Ctrl_Cor_with_3_inner, Neffmems_r)
    CALL Div_CV_cons(Ctrl_Cor_with_4_inner, Neffmems_r)
    CALL Div_CV_cons(Ctrl_Cor_with_5_inner, Neffmems_r)
    PRINT *, 'Done'

    ! Compute the standard deviations
    PRINT *, 'Computing standard deviations'
    StdDevs % v1(1:nlongs/2+1, 1:nlevs) = SQRT(Ctrl_Cor_with_1_inner % v1(1:nlongs/2+1, 1:nlevs))
    StdDevs % v2(1:nlongs/2+1, 1:nlevs) = SQRT(Ctrl_Cor_with_2_inner % v2(1:nlongs/2+1, 1:nlevs))
    StdDevs % v3(1:nlongs/2+1, 1:nlevs) = SQRT(Ctrl_Cor_with_3_inner % v3(1:nlongs/2+1, 1:nlevs))
    StdDevs % v4(1:nlongs/2+1, 1:nlevs) = SQRT(Ctrl_Cor_with_4_inner % v4(1:nlongs/2+1, 1:nlevs))
    StdDevs % v5(1:nlongs/2+1, 1:nlevs) = SQRT(Ctrl_Cor_with_5_inner % v5(1:nlongs/2+1, 1:nlevs))
    PRINT *, 'Done'

    ! Compute correlations
    PRINT *, 'Computing correlations'
    Ctrl_Cor_with_1_inner % v1(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_1_inner % v1(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v1(1:nlongs/2+1, 1:nlevs) * StdDevs % v1(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_1_inner % v2(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_1_inner % v2(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v1(1:nlongs/2+1, 1:nlevs) * StdDevs % v2(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_1_inner % v3(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_1_inner % v3(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v1(1:nlongs/2+1, 1:nlevs) * StdDevs % v3(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_1_inner % v4(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_1_inner % v4(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v1(1:nlongs/2+1, 1:nlevs) * StdDevs % v4(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_1_inner % v5(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_1_inner % v5(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v1(1:nlongs/2+1, 1:nlevs) * StdDevs % v5(1:nlongs/2+1, 1:nlevs) )

    Ctrl_Cor_with_2_inner % v1(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_2_inner % v1(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v2(1:nlongs/2+1, 1:nlevs) * StdDevs % v1(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_2_inner % v2(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_2_inner % v2(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v2(1:nlongs/2+1, 1:nlevs) * StdDevs % v2(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_2_inner % v3(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_2_inner % v3(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v2(1:nlongs/2+1, 1:nlevs) * StdDevs % v3(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_2_inner % v4(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_2_inner % v4(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v2(1:nlongs/2+1, 1:nlevs) * StdDevs % v4(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_2_inner % v5(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_2_inner % v5(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v2(1:nlongs/2+1, 1:nlevs) * StdDevs % v5(1:nlongs/2+1, 1:nlevs) )

    Ctrl_Cor_with_3_inner % v1(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_3_inner % v1(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v3(1:nlongs/2+1, 1:nlevs) * StdDevs % v1(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_3_inner % v2(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_3_inner % v2(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v3(1:nlongs/2+1, 1:nlevs) * StdDevs % v2(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_3_inner % v3(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_3_inner % v3(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v3(1:nlongs/2+1, 1:nlevs) * StdDevs % v3(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_3_inner % v4(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_3_inner % v4(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v3(1:nlongs/2+1, 1:nlevs) * StdDevs % v4(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_3_inner % v5(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_3_inner % v5(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v3(1:nlongs/2+1, 1:nlevs) * StdDevs % v5(1:nlongs/2+1, 1:nlevs) )

    Ctrl_Cor_with_4_inner % v1(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_4_inner % v1(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v4(1:nlongs/2+1, 1:nlevs) * StdDevs % v1(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_4_inner % v2(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_4_inner % v2(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v4(1:nlongs/2+1, 1:nlevs) * StdDevs % v2(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_4_inner % v3(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_4_inner % v3(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v4(1:nlongs/2+1, 1:nlevs) * StdDevs % v3(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_4_inner % v4(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_4_inner % v4(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v4(1:nlongs/2+1, 1:nlevs) * StdDevs % v4(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_4_inner % v5(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_4_inner % v5(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v4(1:nlongs/2+1, 1:nlevs) * StdDevs % v5(1:nlongs/2+1, 1:nlevs) )

    Ctrl_Cor_with_5_inner % v1(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_5_inner % v1(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v5(1:nlongs/2+1, 1:nlevs) * StdDevs % v1(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_5_inner % v2(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_5_inner % v2(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v5(1:nlongs/2+1, 1:nlevs) * StdDevs % v2(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_5_inner % v3(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_5_inner % v3(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v5(1:nlongs/2+1, 1:nlevs) * StdDevs % v3(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_5_inner % v4(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_5_inner % v4(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v5(1:nlongs/2+1, 1:nlevs) * StdDevs % v4(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_5_inner % v5(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_5_inner % v5(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v5(1:nlongs/2+1, 1:nlevs) * StdDevs % v5(1:nlongs/2+1, 1:nlevs) )
    PRINT *, 'Done'


    ! Add on to the outer for averages
    PRINT *, 'Averaging over ensembles'
    CALL Add_CVs (Ctrl_Cor_with_1_outer, Ctrl_Cor_with_1_inner)
    CALL Add_CVs (Ctrl_Cor_with_2_outer, Ctrl_Cor_with_2_inner)
    CALL Add_CVs (Ctrl_Cor_with_3_outer, Ctrl_Cor_with_3_inner)
    CALL Add_CVs (Ctrl_Cor_with_4_outer, Ctrl_Cor_with_4_inner)
    CALL Add_CVs (Ctrl_Cor_with_5_outer, Ctrl_Cor_with_5_inner)
    PRINT *, 'Done'


  END DO


  ! Normalise correlations
  PRINT *, 'Normalising'
  CALL Div_CV_cons(Ctrl_Cor_with_1_outer, Nens_r)
  CALL Div_CV_cons(Ctrl_Cor_with_2_outer, Nens_r)
  CALL Div_CV_cons(Ctrl_Cor_with_3_outer, Nens_r)
  CALL Div_CV_cons(Ctrl_Cor_with_4_outer, Nens_r)
  CALL Div_CV_cons(Ctrl_Cor_with_5_outer, Nens_r)
  PRINT *, 'Done'


  ! Write out the correlations

  PRINT *, 'Outputting the correlations'
  CALL Write_one_field ('SpecSpCor_psi_psi.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_1_outer % v1(1:nlongs/2+1,1:nlevs), &
                        'SpecSpCor_psi_psi')
  CALL Write_one_field ('SpecSpCor_psi_chi.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_1_outer % v2(1:nlongs/2+1,1:nlevs), &
                        'SpecSpCor_psi_chi')
  CALL Write_one_field ('SpecSpCor_psi_ru.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_1_outer % v3(1:nlongs/2+1,1:nlevs), &
                        'SpecSpCor_psi_ru')
  CALL Write_one_field ('SpecSpCor_psi_bu.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_1_outer % v4(1:nlongs/2+1,1:nlevs), &
                        'SpecSpCor_psi_bu')
  CALL Write_one_field ('SpecSpCor_psi_w.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_1_outer % v5(1:nlongs/2+1,1:nlevs), &
                        'SpecSpCor_psi_w')

  CALL Write_one_field ('SpecSpCor_chi_psi.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_2_outer % v1(1:nlongs/2+1,1:nlevs), &
                        'SpecSpCor_chi_psi')
  CALL Write_one_field ('SpecSpCor_chi_chi.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_2_outer % v2(1:nlongs/2+1,1:nlevs), &
                        'SpecSpCor_chi_chi')
  CALL Write_one_field ('SpecSpCor_chi_ru.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_2_outer % v3(1:nlongs/2+1,1:nlevs), &
                        'SpecSpCor_chi_ru')
  CALL Write_one_field ('SpecSpCor_chi_bu.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_2_outer % v4(1:nlongs/2+1,1:nlevs), &
                        'SpecSpCor_chi_bu')
  CALL Write_one_field ('SpecSpCor_chi_w.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_2_outer % v5(1:nlongs/2+1,1:nlevs), &
                        'SpecSpCor_chi_w')

  CALL Write_one_field ('SpecSpCor_ru_psi.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_3_outer % v1(1:nlongs/2+1,1:nlevs), &
                        'SpecSpCor_ru_psi')
  CALL Write_one_field ('SpecSpCor_ru_chi.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_3_outer % v2(1:nlongs/2+1,1:nlevs), &
                        'SpecSpCor_ru_chi')
  CALL Write_one_field ('SpecSpCor_ru_ru.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_3_outer % v3(1:nlongs/2+1,1:nlevs), &
                        'SpecSpCor_ru_ru')
  CALL Write_one_field ('SpecSpCor_ru_bu.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_3_outer % v4(1:nlongs/2+1,1:nlevs), &
                        'SpecSpCor_ru_bu')
  CALL Write_one_field ('SpecSpCor_ru_w.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_3_outer % v5(1:nlongs/2+1,1:nlevs), &
                        'SpecSpCor_ru_w')

  CALL Write_one_field ('SpecSpCor_bu_psi.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_4_outer % v1(1:nlongs/2+1,1:nlevs), &
                        'SpecSpCor_bu_psi')
  CALL Write_one_field ('SpecSpCor_bu_chi.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_4_outer % v2(1:nlongs/2+1,1:nlevs), &
                        'SpecSpCor_bu_chi')
  CALL Write_one_field ('SpecSpCor_bu_ru.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_4_outer % v3(1:nlongs/2+1,1:nlevs), &
                        'SpecSpCor_bu_ru')
  CALL Write_one_field ('SpecSpCor_bu_bu.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_4_outer % v4(1:nlongs/2+1,1:nlevs), &
                        'SpecSpCor_bu_bu')
  CALL Write_one_field ('SpecSpCor_bu_w.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_4_outer % v5(1:nlongs/2+1,1:nlevs), &
                        'SpecSpCor_bu_w')

  CALL Write_one_field ('SpecSpCor_w_psi.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_5_outer % v1(1:nlongs/2+1,1:nlevs), &
                        'SpecSpCor_w_psi')
  CALL Write_one_field ('SpecSpCor_w_chi.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_5_outer % v2(1:nlongs/2+1,1:nlevs), &
                        'SpecSpCor_w_chi')
  CALL Write_one_field ('SpecSpCor_w_ru.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_5_outer % v3(1:nlongs/2+1,1:nlevs), &
                        'SpecSpCor_w_ru')
  CALL Write_one_field ('SpecSpCor_w_bu.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_5_outer % v4(1:nlongs/2+1,1:nlevs), &
                        'SpecSpCor_w_bu')
  CALL Write_one_field ('SpecSpCor_w_w.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_5_outer % v5(1:nlongs/2+1,1:nlevs), &
                        'SpecSpCor_w_w')


  ! Write out the standard deviations
  CALL Write_one_field ('SpecSpstddev_psi.nc', nlongs/2+1, nlevs, StdDevs % v1(1:nlongs/2+1,1:nlevs), &
                        'SpecSpstddev_psi')
  CALL Write_one_field ('SpecSpstddev_chi.nc', nlongs/2+1, nlevs, StdDevs % v2(1:nlongs/2+1,1:nlevs), &
                        'SpecSpstddev_chi')
  CALL Write_one_field ('SpecSpstddev_ru.nc', nlongs/2+1, nlevs, StdDevs % v3(1:nlongs/2+1,1:nlevs), &
                        'SpecSpstddev_ru')
  CALL Write_one_field ('SpecSpstddev_bu.nc', nlongs/2+1, nlevs, StdDevs % v4(1:nlongs/2+1,1:nlevs), &
                        'SpecSpstddev_bu')
  CALL Write_one_field ('SpecSpstddev_w.nc', nlongs/2+1, nlevs, StdDevs % v5(1:nlongs/2+1,1:nlevs), &
                        'SpecSpstddev_w')

  PRINT *, 'Done'

  ! Tidy up
  CALL Deallocate_model_vars (ABC_mean)
  CALL Deallocate_model_vars (ABC_pert)
  CALL Deallocate_CVs(ControlVar)
  CALL Deallocate_CVs(ControlVarS)
  CALL Deallocate_CVs(StdDevs)
  CALL Deallocate_CVs(Ctrl_Cor_with_1_outer)
  CALL Deallocate_CVs(Ctrl_Cor_with_2_outer)
  CALL Deallocate_CVs(Ctrl_Cor_with_3_outer)
  CALL Deallocate_CVs(Ctrl_Cor_with_4_outer)
  CALL Deallocate_CVs(Ctrl_Cor_with_5_outer)
  CALL Deallocate_CVs(Ctrl_Cor_with_1_inner)
  CALL Deallocate_CVs(Ctrl_Cor_with_2_inner)
  CALL Deallocate_CVs(Ctrl_Cor_with_3_inner)
  CALL Deallocate_CVs(Ctrl_Cor_with_4_inner)
  CALL Deallocate_CVs(Ctrl_Cor_with_5_inner)

END IF






IF (Diagswitch(3)) THEN
  PRINT *, '========================================================================='
  PRINT *, 'Performing correlations between actual control variables'
  PRINT *, '========================================================================='

  CALL Initialise_model_vars (ABC_mean, .FALSE.)
  CALL Initialise_model_vars (ABC_pert, .FALSE.)
  CALL Initialise_CVs(ControlVar, .FALSE.)
  CALL Initialise_CVs(StdDevs, .FALSE.)


  ! Initialise the correlations to zero
  CALL Initialise_CVs(Ctrl_Cor_with_1_outer, .FALSE.)
  CALL Initialise_CVs(Ctrl_Cor_with_2_outer, .FALSE.)
  CALL Initialise_CVs(Ctrl_Cor_with_3_outer, .FALSE.)
  CALL Initialise_CVs(Ctrl_Cor_with_4_outer, .FALSE.)
  CALL Initialise_CVs(Ctrl_Cor_with_5_outer, .FALSE.)


  PRINT *, 'There are ', Neffmems, ' effective ensemble members.'
  DO ens = 1, Nens

    ! Initialise the correlations to zero
    CALL Initialise_CVs(Ctrl_Cor_with_1_inner, .FALSE.)
    CALL Initialise_CVs(Ctrl_Cor_with_2_inner, .FALSE.)
    CALL Initialise_CVs(Ctrl_Cor_with_3_inner, .FALSE.)
    CALL Initialise_CVs(Ctrl_Cor_with_4_inner, .FALSE.)
    CALL Initialise_CVs(Ctrl_Cor_with_5_inner, .FALSE.)

    ! Read-in the ensemble mean for this ensemble
    WRITE (meanfile, '(A,A,I0.3,A)') TRIM(datadirABCperts), '/MeanABC', ens, '.nc'
    PRINT *, 'Reading mean state ', TRIM (meanfile)
    CALL Read_state_2d (meanfile, ABC_mean, dims, 1, .TRUE.)
    PRINT *, '-- done'

    ! Read-in the ABC perts and convert to parameter perts
    item = 0
    DO mem = 1, NEnsMems
      DO lat = 1, Nlats
        item = item + 1

        ! Read-in a pert found from stage 2
        WRITE (ABCfile, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirABCperts), '/PertABC_Ens', ens, '_Item', item, '.nc'
        PRINT *, 'Reading file ', TRIM(ABCfile)
        CALL Read_state_2d (ABCfile, ABC_pert, dims, 1, .TRUE.)
        PRINT *, '-- done'
        IF ((ens == 1) .AND. (item == 1)) THEN
          CALL Set_ht_dep_cons (dims)
        END IF

        ! Pass through the inverse parameter transform
        PRINT *, 'Performing inverse control variable transform'
        CALL U_trans_inv (ABC_mean,                         &
                          ControlVar,                       &
                          ABC_pert,                         &
                          CVT,                              &
                          dims)
        PRINT *, 'Done'


        ! ========================================================================
        PRINT *, 'Computing covariance contribution'
        ! Deal with the largest scale first
        ! ---------------------------------
        ! psi with psi
        Ctrl_Cor_with_1_inner % v1(1, 1:nlevs) = Ctrl_Cor_with_1_inner % v1(1, 1:nlevs) +  &
          ControlVar % v1(1, 1:nlevs) * ControlVar % v1(1, 1:nlevs)
        ! psi with chi
        Ctrl_Cor_with_1_inner % v2(1, 1:nlevs) = Ctrl_Cor_with_1_inner % v2(1, 1:nlevs) +  &
          ControlVar % v1(1, 1:nlevs) * ControlVar % v2(1, 1:nlevs)
        ! psi with (unbalanced) r
        Ctrl_Cor_with_1_inner % v3(1, 1:nlevs) = Ctrl_Cor_with_1_inner % v3(1, 1:nlevs) +  &
          ControlVar % v1(1, 1:nlevs) * ControlVar % v3(1, 1:nlevs)
        ! psi with (unbalanced) b
        Ctrl_Cor_with_1_inner % v4(1, 1:nlevs) = Ctrl_Cor_with_1_inner % v4(1, 1:nlevs) +  &
          ControlVar % v1(1, 1:nlevs) * ControlVar % v4(1, 1:nlevs)
        ! psi with (unbalanced) w
        Ctrl_Cor_with_1_inner % v5(1, 1:nlevs) = Ctrl_Cor_with_1_inner % v5(1, 1:nlevs) +  &
          ControlVar % v1(1, 1:nlevs) * ControlVar % v5(1, 1:nlevs)

        ! chi with psi
        Ctrl_Cor_with_2_inner % v1(1, 1:nlevs) = Ctrl_Cor_with_2_inner % v1(1, 1:nlevs) +  &
          ControlVar % v2(1, 1:nlevs) * ControlVar % v1(1, 1:nlevs)
        ! chi with chi
        Ctrl_Cor_with_2_inner % v2(1, 1:nlevs) = Ctrl_Cor_with_2_inner % v2(1, 1:nlevs) +  &
          ControlVar % v2(1, 1:nlevs) * ControlVar % v2(1, 1:nlevs)
        ! chi with (unbalanced) r
        Ctrl_Cor_with_2_inner % v3(1, 1:nlevs) = Ctrl_Cor_with_2_inner % v3(1, 1:nlevs) +  &
          ControlVar % v2(1, 1:nlevs) * ControlVar % v3(1, 1:nlevs)
        ! chi with (unbalanced) b
        Ctrl_Cor_with_2_inner % v4(1, 1:nlevs) = Ctrl_Cor_with_2_inner % v4(1, 1:nlevs) +  &
          ControlVar % v2(1, 1:nlevs) * ControlVar % v4(1, 1:nlevs)
        ! chi with (unbalanced) w
        Ctrl_Cor_with_2_inner % v5(1, 1:nlevs) = Ctrl_Cor_with_2_inner % v5(1, 1:nlevs) +  &
          ControlVar % v2(1, 1:nlevs) * ControlVar % v5(1, 1:nlevs)

        ! (unbalanced) r with psi
        Ctrl_Cor_with_3_inner % v1(1, 1:nlevs) = Ctrl_Cor_with_3_inner % v1(1, 1:nlevs) +  &
          ControlVar % v3(1, 1:nlevs) * ControlVar % v1(1, 1:nlevs)
        ! (unbalanced) r with chi
        Ctrl_Cor_with_3_inner % v2(1, 1:nlevs) = Ctrl_Cor_with_3_inner % v2(1, 1:nlevs) +  &
          ControlVar % v3(1, 1:nlevs) * ControlVar % v2(1, 1:nlevs)
        ! (unbalanced) r with (unbalanced) r
        Ctrl_Cor_with_3_inner % v3(1, 1:nlevs) = Ctrl_Cor_with_3_inner % v3(1, 1:nlevs) +  &
          ControlVar % v3(1, 1:nlevs) * ControlVar % v3(1, 1:nlevs)
        ! (unbalanced) r with (unbalanced) b
        Ctrl_Cor_with_3_inner % v4(1, 1:nlevs) = Ctrl_Cor_with_3_inner % v4(1, 1:nlevs) +  &
          ControlVar % v3(1, 1:nlevs) * ControlVar % v4(1, 1:nlevs)
        ! (unbalanced) r with (unbalanced) w
        Ctrl_Cor_with_3_inner % v5(1, 1:nlevs) = Ctrl_Cor_with_3_inner % v5(1, 1:nlevs) +  &
          ControlVar % v3(1, 1:nlevs) * ControlVar % v5(1, 1:nlevs)

        ! (unbalanced) b with psi
        Ctrl_Cor_with_4_inner % v1(1, 1:nlevs) = Ctrl_Cor_with_4_inner % v1(1, 1:nlevs) +  &
          ControlVar % v4(1, 1:nlevs) * ControlVar % v1(1, 1:nlevs)
        ! (unbalanced) b with chi
        Ctrl_Cor_with_4_inner % v2(1, 1:nlevs) = Ctrl_Cor_with_4_inner % v2(1, 1:nlevs) +  &
          ControlVar % v4(1, 1:nlevs) * ControlVar % v2(1, 1:nlevs)
        ! (unbalanced) b with (unbalanced) r
        Ctrl_Cor_with_4_inner % v3(1, 1:nlevs) = Ctrl_Cor_with_4_inner % v3(1, 1:nlevs) +  &
          ControlVar % v4(1, 1:nlevs) * ControlVar % v3(1, 1:nlevs)
        ! (unbalanced) b with (unbalanced) b
        Ctrl_Cor_with_4_inner % v4(1, 1:nlevs) = Ctrl_Cor_with_4_inner % v4(1, 1:nlevs) +  &
          ControlVar % v4(1, 1:nlevs) * ControlVar % v4(1, 1:nlevs)
        ! (unbalanced) b with (unbalanced) w
        Ctrl_Cor_with_4_inner % v5(1, 1:nlevs) = Ctrl_Cor_with_4_inner % v5(1, 1:nlevs) +  &
          ControlVar % v4(1, 1:nlevs) * ControlVar % v5(1, 1:nlevs)

        ! (unbalanced) w with psi
        Ctrl_Cor_with_5_inner % v1(1, 1:nlevs) = Ctrl_Cor_with_5_inner % v1(1, 1:nlevs) +  &
          ControlVar % v5(1, 1:nlevs) * ControlVar % v1(1, 1:nlevs)
        ! (unbalanced) w with chi
        Ctrl_Cor_with_5_inner % v2(1, 1:nlevs) = Ctrl_Cor_with_5_inner % v2(1, 1:nlevs) +  &
          ControlVar % v5(1, 1:nlevs) * ControlVar % v2(1, 1:nlevs)
        ! (unbalanced) w with (unbalanced) r
        Ctrl_Cor_with_5_inner % v3(1, 1:nlevs) = Ctrl_Cor_with_5_inner % v3(1, 1:nlevs) +  &
          ControlVar % v5(1, 1:nlevs) * ControlVar % v3(1, 1:nlevs)
        ! (unbalanced) w with (unbalanced) b
        Ctrl_Cor_with_5_inner % v4(1, 1:nlevs) = Ctrl_Cor_with_5_inner % v4(1, 1:nlevs) +  &
          ControlVar % v5(1, 1:nlevs) * ControlVar % v4(1, 1:nlevs)
        ! (unbalanced) w with (unbalanced) w
        Ctrl_Cor_with_5_inner % v5(1, 1:nlevs) = Ctrl_Cor_with_5_inner % v5(1, 1:nlevs) +  &
          ControlVar % v5(1, 1:nlevs) * ControlVar % v5(1, 1:nlevs)

        ! Deal with the bulk of the scales
        ! --------------------------------
        DO k = 2, nlongs/2
          real_index = 2*k-2
          imag_index = 2*k-1

          ! psi with psi
          Ctrl_Cor_with_1_inner % v1(k, 1:nlevs) = Ctrl_Cor_with_1_inner % v1(k, 1:nlevs) +  &
            ControlVar % v1(real_index, 1:nlevs) * ControlVar % v1(real_index, 1:nlevs) +  &
            ControlVar % v1(imag_index, 1:nlevs) * ControlVar % v1(imag_index, 1:nlevs)
          ! psi with chi
          Ctrl_Cor_with_1_inner % v2(k, 1:nlevs) = Ctrl_Cor_with_1_inner % v2(k, 1:nlevs) +  &
            ControlVar % v1(real_index, 1:nlevs) * ControlVar % v2(real_index, 1:nlevs) +  &
            ControlVar % v1(imag_index, 1:nlevs) * ControlVar % v2(imag_index, 1:nlevs)
          ! psi with (unbalanced) r
          Ctrl_Cor_with_1_inner % v3(k, 1:nlevs) = Ctrl_Cor_with_1_inner % v3(k, 1:nlevs) +  &
            ControlVar % v1(real_index, 1:nlevs) * ControlVar % v3(real_index, 1:nlevs) +  &
            ControlVar % v1(imag_index, 1:nlevs) * ControlVar % v3(imag_index, 1:nlevs)
          ! psi with (unbalanced) b
          Ctrl_Cor_with_1_inner % v4(k, 1:nlevs) = Ctrl_Cor_with_1_inner % v4(k, 1:nlevs) +  &
            ControlVar % v1(real_index, 1:nlevs) * ControlVar % v4(real_index, 1:nlevs) +  &
            ControlVar % v1(imag_index, 1:nlevs) * ControlVar % v4(imag_index, 1:nlevs)
          ! psi with (unbalanced) w
          Ctrl_Cor_with_1_inner % v5(k, 1:nlevs) = Ctrl_Cor_with_1_inner % v5(k, 1:nlevs) +  &
            ControlVar % v1(real_index, 1:nlevs) * ControlVar % v5(real_index, 1:nlevs) +  &
            ControlVar % v1(imag_index, 1:nlevs) * ControlVar % v5(imag_index, 1:nlevs)

          ! chi with psi
          Ctrl_Cor_with_2_inner % v1(k, 1:nlevs) = Ctrl_Cor_with_2_inner % v1(k, 1:nlevs) +  &
            ControlVar % v2(real_index, 1:nlevs) * ControlVar % v1(real_index, 1:nlevs) +  &
            ControlVar % v2(imag_index, 1:nlevs) * ControlVar % v1(imag_index, 1:nlevs)
          ! chi with chi
          Ctrl_Cor_with_2_inner % v2(k, 1:nlevs) = Ctrl_Cor_with_2_inner % v2(k, 1:nlevs) +  &
            ControlVar % v2(real_index, 1:nlevs) * ControlVar % v2(real_index, 1:nlevs) +  &
            ControlVar % v2(imag_index, 1:nlevs) * ControlVar % v2(imag_index, 1:nlevs)
          ! chi with (unbalanced) r
          Ctrl_Cor_with_2_inner % v3(k, 1:nlevs) = Ctrl_Cor_with_2_inner % v3(k, 1:nlevs) +  &
            ControlVar % v2(real_index, 1:nlevs) * ControlVar % v3(real_index, 1:nlevs) +  &
            ControlVar % v2(imag_index, 1:nlevs) * ControlVar % v3(imag_index, 1:nlevs)
          ! chi with (unbalanced) b
          Ctrl_Cor_with_2_inner % v4(k, 1:nlevs) = Ctrl_Cor_with_2_inner % v4(k, 1:nlevs) +  &
            ControlVar % v2(real_index, 1:nlevs) * ControlVar % v4(real_index, 1:nlevs) +  &
            ControlVar % v2(imag_index, 1:nlevs) * ControlVar % v4(imag_index, 1:nlevs)
          ! chi with (unbalanced) w
          Ctrl_Cor_with_2_inner % v5(k, 1:nlevs) = Ctrl_Cor_with_2_inner % v5(k, 1:nlevs) +  &
            ControlVar % v2(real_index, 1:nlevs) * ControlVar % v5(real_index, 1:nlevs) +  &
            ControlVar % v2(imag_index, 1:nlevs) * ControlVar % v5(imag_index, 1:nlevs)

          ! (unbalanced) r with psi
          Ctrl_Cor_with_3_inner % v1(k, 1:nlevs) = Ctrl_Cor_with_3_inner % v1(k, 1:nlevs) +  &
            ControlVar % v3(real_index, 1:nlevs) * ControlVar % v1(real_index, 1:nlevs) +  &
            ControlVar % v3(imag_index, 1:nlevs) * ControlVar % v1(imag_index, 1:nlevs)
          ! (unbalanced) r with chi
          Ctrl_Cor_with_3_inner % v2(k, 1:nlevs) = Ctrl_Cor_with_3_inner % v2(k, 1:nlevs) +  &
            ControlVar % v3(real_index, 1:nlevs) * ControlVar % v2(real_index, 1:nlevs) +  &
            ControlVar % v3(imag_index, 1:nlevs) * ControlVar % v2(imag_index, 1:nlevs)
          ! (unbalanced) r with (unbalanced) r
          Ctrl_Cor_with_3_inner % v3(k, 1:nlevs) = Ctrl_Cor_with_3_inner % v3(k, 1:nlevs) +  &
            ControlVar % v3(real_index, 1:nlevs) * ControlVar % v3(real_index, 1:nlevs) +  &
            ControlVar % v3(imag_index, 1:nlevs) * ControlVar % v3(imag_index, 1:nlevs)
          ! (unbalanced) r with (unbalanced) b
          Ctrl_Cor_with_3_inner % v4(k, 1:nlevs) = Ctrl_Cor_with_3_inner % v4(k, 1:nlevs) +  &
            ControlVar % v3(real_index, 1:nlevs) * ControlVar % v4(real_index, 1:nlevs) +  &
            ControlVar % v3(imag_index, 1:nlevs) * ControlVar % v4(imag_index, 1:nlevs)
          ! (unbalanced) r with (unbalanced) w
          Ctrl_Cor_with_3_inner % v5(k, 1:nlevs) = Ctrl_Cor_with_3_inner % v5(k, 1:nlevs) +  &
            ControlVar % v3(real_index, 1:nlevs) * ControlVar % v5(real_index, 1:nlevs) +  &
            ControlVar % v3(imag_index, 1:nlevs) * ControlVar % v5(imag_index, 1:nlevs)

          ! (unbalanced) b with psi
          Ctrl_Cor_with_4_inner % v1(k, 1:nlevs) = Ctrl_Cor_with_4_inner % v1(k, 1:nlevs) +  &
            ControlVar % v4(real_index, 1:nlevs) * ControlVar % v1(real_index, 1:nlevs) +  &
            ControlVar % v4(imag_index, 1:nlevs) * ControlVar % v1(imag_index, 1:nlevs)
          ! (unbalanced) b with chi
          Ctrl_Cor_with_4_inner % v2(k, 1:nlevs) = Ctrl_Cor_with_4_inner % v2(k, 1:nlevs) +  &
            ControlVar % v4(real_index, 1:nlevs) * ControlVar % v2(real_index, 1:nlevs) +  &
            ControlVar % v4(imag_index, 1:nlevs) * ControlVar % v2(imag_index, 1:nlevs)
          ! (unbalanced) b with (unbalanced) r
          Ctrl_Cor_with_4_inner % v3(k, 1:nlevs) = Ctrl_Cor_with_4_inner % v3(k, 1:nlevs) +  &
            ControlVar % v4(real_index, 1:nlevs) * ControlVar % v3(real_index, 1:nlevs) +  &
            ControlVar % v4(imag_index, 1:nlevs) * ControlVar % v3(imag_index, 1:nlevs)
          ! (unbalanced) b with (unbalanced) b
          Ctrl_Cor_with_4_inner % v4(k, 1:nlevs) = Ctrl_Cor_with_4_inner % v4(k, 1:nlevs) +  &
            ControlVar % v4(real_index, 1:nlevs) * ControlVar % v4(real_index, 1:nlevs) +  &
            ControlVar % v4(imag_index, 1:nlevs) * ControlVar % v4(imag_index, 1:nlevs)
          ! (unbalanced) b with (unbalanced) w
          Ctrl_Cor_with_4_inner % v5(k, 1:nlevs) = Ctrl_Cor_with_4_inner % v5(k, 1:nlevs) +  &
            ControlVar % v4(real_index, 1:nlevs) * ControlVar % v5(real_index, 1:nlevs) +  &
            ControlVar % v4(imag_index, 1:nlevs) * ControlVar % v5(imag_index, 1:nlevs)

          ! (unbalanced) w with psi
          Ctrl_Cor_with_5_inner % v1(k, 1:nlevs) = Ctrl_Cor_with_5_inner % v1(k, 1:nlevs) +  &
            ControlVar % v5(real_index, 1:nlevs) * ControlVar % v1(real_index, 1:nlevs) +  &
            ControlVar % v5(imag_index, 1:nlevs) * ControlVar % v1(imag_index, 1:nlevs)
          ! (unbalanced) w with chi
          Ctrl_Cor_with_5_inner % v2(k, 1:nlevs) = Ctrl_Cor_with_5_inner % v2(k, 1:nlevs) +  &
            ControlVar % v5(real_index, 1:nlevs) * ControlVar % v2(real_index, 1:nlevs) +  &
            ControlVar % v5(imag_index, 1:nlevs) * ControlVar % v2(imag_index, 1:nlevs)
          ! (unbalanced) w with (unbalanced) r
          Ctrl_Cor_with_5_inner % v3(k, 1:nlevs) = Ctrl_Cor_with_5_inner % v3(k, 1:nlevs) +  &
            ControlVar % v5(real_index, 1:nlevs) * ControlVar % v3(real_index, 1:nlevs) +  &
            ControlVar % v5(imag_index, 1:nlevs) * ControlVar % v3(imag_index, 1:nlevs)
          ! (unbalanced) w with (unbalanced) b
          Ctrl_Cor_with_5_inner % v4(k, 1:nlevs) = Ctrl_Cor_with_5_inner % v4(k, 1:nlevs) +  &
            ControlVar % v5(real_index, 1:nlevs) * ControlVar % v4(real_index, 1:nlevs) +  &
            ControlVar % v5(imag_index, 1:nlevs) * ControlVar % v4(imag_index, 1:nlevs)
          ! (unbalanced) w with (unbalanced) w
          Ctrl_Cor_with_5_inner % v5(k, 1:nlevs) = Ctrl_Cor_with_5_inner % v5(k, 1:nlevs) +  &
            ControlVar % v5(real_index, 1:nlevs) * ControlVar % v5(real_index, 1:nlevs) +  &
            ControlVar % v5(imag_index, 1:nlevs) * ControlVar % v5(imag_index, 1:nlevs)
        END DO

        ! Deal with the smallest scale
        ! ----------------------------
        last_index = nlongs/2+1

        ! psi with psi
        Ctrl_Cor_with_1_inner % v1(last_index, 1:nlevs) = Ctrl_Cor_with_1_inner % v1(last_index, 1:nlevs) +  &
          ControlVar % v1(nlongs, 1:nlevs) * ControlVar % v1(nlongs, 1:nlevs)
        ! psi with chi
        Ctrl_Cor_with_1_inner % v2(last_index, 1:nlevs) = Ctrl_Cor_with_1_inner % v2(last_index, 1:nlevs) +  &
          ControlVar % v1(nlongs, 1:nlevs) * ControlVar % v2(nlongs, 1:nlevs)
        ! psi with (unbalanced) r
        Ctrl_Cor_with_1_inner % v3(last_index, 1:nlevs) = Ctrl_Cor_with_1_inner % v3(last_index, 1:nlevs) +  &
          ControlVar % v1(nlongs, 1:nlevs) * ControlVar % v3(nlongs, 1:nlevs)
        ! psi with (unbalanced) b
        Ctrl_Cor_with_1_inner % v4(last_index, 1:nlevs) = Ctrl_Cor_with_1_inner % v4(last_index, 1:nlevs) +  &
          ControlVar % v1(nlongs, 1:nlevs) * ControlVar % v4(nlongs, 1:nlevs)
        ! psi with (unbalanced) w
        Ctrl_Cor_with_1_inner % v5(last_index, 1:nlevs) = Ctrl_Cor_with_1_inner % v5(last_index, 1:nlevs) +  &
          ControlVar % v1(nlongs, 1:nlevs) * ControlVar % v5(nlongs, 1:nlevs)

        ! chi with psi
        Ctrl_Cor_with_2_inner % v1(last_index, 1:nlevs) = Ctrl_Cor_with_2_inner % v1(last_index, 1:nlevs) +  &
          ControlVar % v2(nlongs, 1:nlevs) * ControlVar % v1(nlongs, 1:nlevs)
        ! chi with chi
        Ctrl_Cor_with_2_inner % v2(last_index, 1:nlevs) = Ctrl_Cor_with_2_inner % v2(last_index, 1:nlevs) +  &
          ControlVar % v2(nlongs, 1:nlevs) * ControlVar % v2(nlongs, 1:nlevs)
        ! chi with (unbalanced) r
        Ctrl_Cor_with_2_inner % v3(last_index, 1:nlevs) = Ctrl_Cor_with_2_inner % v3(last_index, 1:nlevs) +  &
          ControlVar % v2(nlongs, 1:nlevs) * ControlVar % v3(nlongs, 1:nlevs)
        ! chi with (unbalanced) b
        Ctrl_Cor_with_2_inner % v4(last_index, 1:nlevs) = Ctrl_Cor_with_2_inner % v4(last_index, 1:nlevs) +  &
          ControlVar % v2(nlongs, 1:nlevs) * ControlVar % v4(nlongs, 1:nlevs)
        ! chi with (unbalanced) w
        Ctrl_Cor_with_2_inner % v5(last_index, 1:nlevs) = Ctrl_Cor_with_2_inner % v5(last_index, 1:nlevs) +  &
          ControlVar % v2(nlongs, 1:nlevs) * ControlVar % v5(nlongs, 1:nlevs)

        ! (unbalanced) r with psi
        Ctrl_Cor_with_3_inner % v1(last_index, 1:nlevs) = Ctrl_Cor_with_3_inner % v1(last_index, 1:nlevs) +  &
          ControlVar % v3(nlongs, 1:nlevs) * ControlVar % v1(nlongs, 1:nlevs)
        ! (unbalanced) r with chi
        Ctrl_Cor_with_3_inner % v2(last_index, 1:nlevs) = Ctrl_Cor_with_3_inner % v2(last_index, 1:nlevs) +  &
          ControlVar % v3(nlongs, 1:nlevs) * ControlVar % v2(nlongs, 1:nlevs)
        ! (unbalanced) r with (unbalanced) r
        Ctrl_Cor_with_3_inner % v3(last_index, 1:nlevs) = Ctrl_Cor_with_3_inner % v3(last_index, 1:nlevs) +  &
          ControlVar % v3(nlongs, 1:nlevs) * ControlVar % v3(nlongs, 1:nlevs)
        ! (unbalanced) r with (unbalanced) b
        Ctrl_Cor_with_3_inner % v4(last_index, 1:nlevs) = Ctrl_Cor_with_3_inner % v4(last_index, 1:nlevs) +  &
          ControlVar % v3(nlongs, 1:nlevs) * ControlVar % v4(nlongs, 1:nlevs)
        ! (unbalanced) r with (unbalanced) w
        Ctrl_Cor_with_3_inner % v5(last_index, 1:nlevs) = Ctrl_Cor_with_3_inner % v5(last_index, 1:nlevs) +  &
          ControlVar % v3(nlongs, 1:nlevs) * ControlVar % v5(nlongs, 1:nlevs)

        ! (unbalanced) b with psi
        Ctrl_Cor_with_4_inner % v1(last_index, 1:nlevs) = Ctrl_Cor_with_4_inner % v1(last_index, 1:nlevs) +  &
          ControlVar % v4(nlongs, 1:nlevs) * ControlVar % v1(nlongs, 1:nlevs)
        ! (unbalanced) b with chi
        Ctrl_Cor_with_4_inner % v2(last_index, 1:nlevs) = Ctrl_Cor_with_4_inner % v2(last_index, 1:nlevs) +  &
          ControlVar % v4(nlongs, 1:nlevs) * ControlVar % v2(nlongs, 1:nlevs)
        ! (unbalanced) b with (unbalanced) r
        Ctrl_Cor_with_4_inner % v3(last_index, 1:nlevs) = Ctrl_Cor_with_4_inner % v3(last_index, 1:nlevs) +  &
          ControlVar % v4(nlongs, 1:nlevs) * ControlVar % v3(nlongs, 1:nlevs)
        ! (unbalanced) b with (unbalanced) b
        Ctrl_Cor_with_4_inner % v4(last_index, 1:nlevs) = Ctrl_Cor_with_4_inner % v4(last_index, 1:nlevs) +  &
          ControlVar % v4(nlongs, 1:nlevs) * ControlVar % v4(nlongs, 1:nlevs)
        ! (unbalanced) b with (unbalanced) w
        Ctrl_Cor_with_4_inner % v5(last_index, 1:nlevs) = Ctrl_Cor_with_4_inner % v5(last_index, 1:nlevs) +  &
          ControlVar % v4(nlongs, 1:nlevs) * ControlVar % v5(nlongs, 1:nlevs)

        ! (unbalanced) w with psi
        Ctrl_Cor_with_5_inner % v1(last_index, 1:nlevs) = Ctrl_Cor_with_5_inner % v1(last_index, 1:nlevs) +  &
          ControlVar % v5(nlongs, 1:nlevs) * ControlVar % v1(nlongs, 1:nlevs)
        ! (unbalanced) w with chi
        Ctrl_Cor_with_5_inner % v2(last_index, 1:nlevs) = Ctrl_Cor_with_5_inner % v2(last_index, 1:nlevs) +  &
          ControlVar % v5(nlongs, 1:nlevs) * ControlVar % v2(nlongs, 1:nlevs)
        ! (unbalanced) w with (unbalanced) r
        Ctrl_Cor_with_5_inner % v3(last_index, 1:nlevs) = Ctrl_Cor_with_5_inner % v3(last_index, 1:nlevs) +  &
          ControlVar % v5(nlongs, 1:nlevs) * ControlVar % v3(nlongs, 1:nlevs)
        ! (unbalanced) w with (unbalanced) b
        Ctrl_Cor_with_5_inner % v4(last_index, 1:nlevs) = Ctrl_Cor_with_5_inner % v4(last_index, 1:nlevs) +  &
          ControlVar % v5(nlongs, 1:nlevs) * ControlVar % v4(nlongs, 1:nlevs)
        ! (unbalanced) w with (unbalanced) w
        Ctrl_Cor_with_5_inner % v5(last_index, 1:nlevs) = Ctrl_Cor_with_5_inner % v5(last_index, 1:nlevs) +  &
          ControlVar % v5(nlongs, 1:nlevs) * ControlVar % v5(nlongs, 1:nlevs)
        PRINT *, 'Done'
        ! ========================================================================

      END DO
    END DO

    ! Normalise to give covariances
    PRINT *, 'Normalising'
    CALL Div_CV_cons(Ctrl_Cor_with_1_inner, Neffmems_r)
    CALL Div_CV_cons(Ctrl_Cor_with_2_inner, Neffmems_r)
    CALL Div_CV_cons(Ctrl_Cor_with_3_inner, Neffmems_r)
    CALL Div_CV_cons(Ctrl_Cor_with_4_inner, Neffmems_r)
    CALL Div_CV_cons(Ctrl_Cor_with_5_inner, Neffmems_r)
    PRINT *, 'Done'

    ! Compute the standard deviations
    PRINT *, 'Computing standard deviations'
    StdDevs % v1(1:nlongs, 1:nlevs) = SQRT(Ctrl_Cor_with_1_inner % v1(1:nlongs, 1:nlevs))
    StdDevs % v2(1:nlongs, 1:nlevs) = SQRT(Ctrl_Cor_with_2_inner % v2(1:nlongs, 1:nlevs))
    StdDevs % v3(1:nlongs, 1:nlevs) = SQRT(Ctrl_Cor_with_3_inner % v3(1:nlongs, 1:nlevs))
    StdDevs % v4(1:nlongs, 1:nlevs) = SQRT(Ctrl_Cor_with_4_inner % v4(1:nlongs, 1:nlevs))
    StdDevs % v5(1:nlongs, 1:nlevs) = SQRT(Ctrl_Cor_with_5_inner % v5(1:nlongs, 1:nlevs))
    PRINT *, 'Done'

    ! Compute correlations
    PRINT *, 'Computing correlations'
    Ctrl_Cor_with_1_inner % v1(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_1_inner % v1(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v1(1:nlongs/2+1, 1:nlevs) * StdDevs % v1(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_1_inner % v2(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_1_inner % v2(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v1(1:nlongs/2+1, 1:nlevs) * StdDevs % v2(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_1_inner % v3(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_1_inner % v3(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v1(1:nlongs/2+1, 1:nlevs) * StdDevs % v3(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_1_inner % v4(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_1_inner % v4(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v1(1:nlongs/2+1, 1:nlevs) * StdDevs % v4(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_1_inner % v5(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_1_inner % v5(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v1(1:nlongs/2+1, 1:nlevs) * StdDevs % v5(1:nlongs/2+1, 1:nlevs) )

    Ctrl_Cor_with_2_inner % v1(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_2_inner % v1(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v2(1:nlongs/2+1, 1:nlevs) * StdDevs % v1(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_2_inner % v2(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_2_inner % v2(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v2(1:nlongs/2+1, 1:nlevs) * StdDevs % v2(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_2_inner % v3(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_2_inner % v3(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v2(1:nlongs/2+1, 1:nlevs) * StdDevs % v3(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_2_inner % v4(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_2_inner % v4(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v2(1:nlongs/2+1, 1:nlevs) * StdDevs % v4(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_2_inner % v5(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_2_inner % v5(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v2(1:nlongs/2+1, 1:nlevs) * StdDevs % v5(1:nlongs/2+1, 1:nlevs) )

    Ctrl_Cor_with_3_inner % v1(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_3_inner % v1(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v3(1:nlongs/2+1, 1:nlevs) * StdDevs % v1(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_3_inner % v2(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_3_inner % v2(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v3(1:nlongs/2+1, 1:nlevs) * StdDevs % v2(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_3_inner % v3(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_3_inner % v3(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v3(1:nlongs/2+1, 1:nlevs) * StdDevs % v3(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_3_inner % v4(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_3_inner % v4(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v3(1:nlongs/2+1, 1:nlevs) * StdDevs % v4(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_3_inner % v5(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_3_inner % v5(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v3(1:nlongs/2+1, 1:nlevs) * StdDevs % v5(1:nlongs/2+1, 1:nlevs) )

    Ctrl_Cor_with_4_inner % v1(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_4_inner % v1(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v4(1:nlongs/2+1, 1:nlevs) * StdDevs % v1(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_4_inner % v2(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_4_inner % v2(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v4(1:nlongs/2+1, 1:nlevs) * StdDevs % v2(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_4_inner % v3(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_4_inner % v3(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v4(1:nlongs/2+1, 1:nlevs) * StdDevs % v3(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_4_inner % v4(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_4_inner % v4(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v4(1:nlongs/2+1, 1:nlevs) * StdDevs % v4(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_4_inner % v5(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_4_inner % v5(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v4(1:nlongs/2+1, 1:nlevs) * StdDevs % v5(1:nlongs/2+1, 1:nlevs) )

    Ctrl_Cor_with_5_inner % v1(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_5_inner % v1(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v5(1:nlongs/2+1, 1:nlevs) * StdDevs % v1(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_5_inner % v2(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_5_inner % v2(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v5(1:nlongs/2+1, 1:nlevs) * StdDevs % v2(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_5_inner % v3(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_5_inner % v3(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v5(1:nlongs/2+1, 1:nlevs) * StdDevs % v3(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_5_inner % v4(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_5_inner % v4(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v5(1:nlongs/2+1, 1:nlevs) * StdDevs % v4(1:nlongs/2+1, 1:nlevs) )
    Ctrl_Cor_with_5_inner % v5(1:nlongs/2+1, 1:nlevs) = Ctrl_Cor_with_5_inner % v5(1:nlongs/2+1, 1:nlevs) / &
      ( StdDevs % v5(1:nlongs/2+1, 1:nlevs) * StdDevs % v5(1:nlongs/2+1, 1:nlevs) )
    PRINT *, 'Done'


    ! Add on to the outer for averages
    PRINT *, 'Averaging over ensembles'
    CALL Add_CVs (Ctrl_Cor_with_1_outer, Ctrl_Cor_with_1_inner)
    CALL Add_CVs (Ctrl_Cor_with_2_outer, Ctrl_Cor_with_2_inner)
    CALL Add_CVs (Ctrl_Cor_with_3_outer, Ctrl_Cor_with_3_inner)
    CALL Add_CVs (Ctrl_Cor_with_4_outer, Ctrl_Cor_with_4_inner)
    CALL Add_CVs (Ctrl_Cor_with_5_outer, Ctrl_Cor_with_5_inner)
    PRINT *, 'Done'


  END DO


  ! Normalise correlations
  PRINT *, 'Normalising'
  CALL Div_CV_cons(Ctrl_Cor_with_1_outer, Nens_r)
  CALL Div_CV_cons(Ctrl_Cor_with_2_outer, Nens_r)
  CALL Div_CV_cons(Ctrl_Cor_with_3_outer, Nens_r)
  CALL Div_CV_cons(Ctrl_Cor_with_4_outer, Nens_r)
  CALL Div_CV_cons(Ctrl_Cor_with_5_outer, Nens_r)
  PRINT *, 'Done'


  ! Write out the correlations
  PRINT *, 'Outputting the correlations'
  CALL Write_one_field ('ConVar_psi_psi.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_1_outer % v1(1:nlongs/2+1,1:nlevs), &
                        'ConVar_psi_psi')
  CALL Write_one_field ('ConVar_psi_chi.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_1_outer % v2(1:nlongs/2+1,1:nlevs), &
                        'ConVar_psi_chi')
  CALL Write_one_field ('ConVar_psi_ru.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_1_outer % v3(1:nlongs/2+1,1:nlevs), &
                        'ConVar_psi_ru')
  CALL Write_one_field ('ConVar_psi_bu.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_1_outer % v4(1:nlongs/2+1,1:nlevs), &
                        'ConVar_psi_bu')
  CALL Write_one_field ('ConVar_psi_w.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_1_outer % v5(1:nlongs/2+1,1:nlevs), &
                        'ConVar_psi_w')

  CALL Write_one_field ('ConVar_chi_psi.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_2_outer % v1(1:nlongs/2+1,1:nlevs), &
                        'ConVar_chi_psi')
  CALL Write_one_field ('ConVar_chi_chi.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_2_outer % v2(1:nlongs/2+1,1:nlevs), &
                        'ConVar_chi_chi')
  CALL Write_one_field ('ConVar_chi_ru.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_2_outer % v3(1:nlongs/2+1,1:nlevs), &
                        'ConVar_chi_ru')
  CALL Write_one_field ('ConVar_chi_bu.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_2_outer % v4(1:nlongs/2+1,1:nlevs), &
                        'ConVar_chi_bu')
  CALL Write_one_field ('ConVar_chi_w.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_2_outer % v5(1:nlongs/2+1,1:nlevs), &
                        'ConVar_chi_w')

  CALL Write_one_field ('ConVar_ru_psi.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_3_outer % v1(1:nlongs/2+1,1:nlevs), &
                        'ConVar_ru_psi')
  CALL Write_one_field ('ConVar_ru_chi.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_3_outer % v2(1:nlongs/2+1,1:nlevs), &
                        'ConVar_ru_chi')
  CALL Write_one_field ('ConVar_ru_ru.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_3_outer % v3(1:nlongs/2+1,1:nlevs), &
                        'ConVar_ru_ru')
  CALL Write_one_field ('ConVar_ru_bu.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_3_outer % v4(1:nlongs/2+1,1:nlevs), &
                        'ConVar_ru_bu')
  CALL Write_one_field ('ConVar_ru_w.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_3_outer % v5(1:nlongs/2+1,1:nlevs), &
                        'ConVar_ru_w')

  CALL Write_one_field ('ConVar_bu_psi.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_4_outer % v1(1:nlongs/2+1,1:nlevs), &
                        'ConVar_bu_psi')
  CALL Write_one_field ('ConVar_bu_chi.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_4_outer % v2(1:nlongs/2+1,1:nlevs), &
                        'ConVar_bu_chi')
  CALL Write_one_field ('ConVar_bu_ru.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_4_outer % v3(1:nlongs/2+1,1:nlevs), &
                        'ConVar_bu_ru')
  CALL Write_one_field ('ConVar_bu_bu.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_4_outer % v4(1:nlongs/2+1,1:nlevs), &
                        'ConVar_bu_bu')
  CALL Write_one_field ('ConVar_bu_w.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_4_outer % v5(1:nlongs/2+1,1:nlevs), &
                        'ConVar_bu_w')

  CALL Write_one_field ('ConVar_w_psi.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_5_outer % v1(1:nlongs/2+1,1:nlevs), &
                        'ConVar_w_psi')
  CALL Write_one_field ('ConVar_w_chi.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_5_outer % v2(1:nlongs/2+1,1:nlevs), &
                        'ConVar_w_chi')
  CALL Write_one_field ('ConVar_w_ru.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_5_outer % v3(1:nlongs/2+1,1:nlevs), &
                        'ConVar_w_ru')
  CALL Write_one_field ('ConVar_w_bu.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_5_outer % v4(1:nlongs/2+1,1:nlevs), &
                        'ConVar_w_bu')
  CALL Write_one_field ('ConVar_w_w.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_5_outer % v5(1:nlongs/2+1,1:nlevs), &
                        'ConVar_w_w')


  ! Write out the standard deviations
  CALL Write_one_field ('ConVarstddev_psi.nc', nlongs/2+1, nlevs, StdDevs % v1(1:nlongs/2+1,1:nlevs), &
                        'ConVarstddev_psi')
  CALL Write_one_field ('ConVarstddev_chi.nc', nlongs/2+1, nlevs, StdDevs % v2(1:nlongs/2+1,1:nlevs), &
                        'ConVarstddev_chi')
  CALL Write_one_field ('ConVarstddev_ru.nc', nlongs/2+1, nlevs, StdDevs % v3(1:nlongs/2+1,1:nlevs), &
                        'ConVarstddev_ru')
  CALL Write_one_field ('ConVarstddev_bu.nc', nlongs/2+1, nlevs, StdDevs % v4(1:nlongs/2+1,1:nlevs), &
                        'ConVarstddev_bu')
  CALL Write_one_field ('ConVarstddev_w.nc', nlongs/2+1, nlevs, StdDevs % v5(1:nlongs/2+1,1:nlevs), &
                        'ConVarstddev_w')

  PRINT *, 'Done'

  ! Tidy up
  CALL Deallocate_model_vars (ABC_mean)
  CALL Deallocate_model_vars (ABC_pert)
  CALL Deallocate_CVs(ControlVar)
  CALL Deallocate_CVs(StdDevs)
  CALL Deallocate_CVs(Ctrl_Cor_with_1_outer)
  CALL Deallocate_CVs(Ctrl_Cor_with_2_outer)
  CALL Deallocate_CVs(Ctrl_Cor_with_3_outer)
  CALL Deallocate_CVs(Ctrl_Cor_with_4_outer)
  CALL Deallocate_CVs(Ctrl_Cor_with_5_outer)
  CALL Deallocate_CVs(Ctrl_Cor_with_1_inner)
  CALL Deallocate_CVs(Ctrl_Cor_with_2_inner)
  CALL Deallocate_CVs(Ctrl_Cor_with_3_inner)
  CALL Deallocate_CVs(Ctrl_Cor_with_4_inner)
  CALL Deallocate_CVs(Ctrl_Cor_with_5_inner)


END IF






IF (DiagSwitch(4)) THEN
  PRINT *, '========================================================================='
  PRINT *, 'Computing ensemble average KE spectrum'
  PRINT *, '========================================================================='
  ! This is based on code /home/ross/DataAssim/EnsDiags/LargeEnsembleDiags/Stefano/KEspec_RB.py

  CALL Initialise_model_vars (ABC_state, .FALSE.)
  ALLOCATE (EnsAvKE(1:nlongs/2+1, 1:nlevs))
  ALLOCATE (data_fl(1:nlevs))
  CALL Initialise_CVs(ControlVarS, .FALSE.)

  ! Initialise array (wn vs height)
  EnsAvKE(1:nlongs/2+1, 1:nlevs) = 0.0


  PRINT *, 'There are ', Neffmems, ' effective ensemble members.'
  DO ens = 1, Nens

    ! Read-in the ABC full states
    item = 0
    DO mem = 1, NEnsMems
      DO lat = 1, Nlats
        item = item + 1

        ! Read-in a state
        WRITE (ABCfile, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirABCfcs), '/FC_Ens', ens, '_Item', item, '.nc'
        PRINT *, 'Reading file ', TRIM(ABCfile)
        CALL Read_state_2d (ABCfile, ABC_state, dims, 1, .TRUE.)
        PRINT *, '-- done'
        IF ((ens == 1) .AND. (item == 1)) THEN
          CALL Set_ht_dep_cons (dims)
        END IF

        ! Replace the full rhotilde with its square-root
        ABC_state % rho(1:nlongs, 1:nlevs) = SQRT(ABC_state % rho(1:nlongs, 1:nlevs))

        ! Transfer u, v, and sqrt rho from half-levels to full-levels and place back in structure
        DO x = 1, nlongs
          ! *** Deal with u ***
          DO z = 1, nlevs
            data_fl(z) = INT_HF (ABC_state % u(x,z),       &
                                 ABC_state % u(x,z+1),     &
                                 z, dims)
          END DO
          ! Place back in the ABC_state structure
          ABC_state % u(x,1:nlevs) = data_fl(1:nlevs)
          ! *** Deal with v ***
          DO z = 1, nlevs
            data_fl(z) = INT_HF (ABC_state % v(x,z),       &
                                 ABC_state % v(x,z+1),     &
                                 z, dims)
          END DO
          ! Place back in the ABC_state structure
          ABC_state % v(x,1:nlevs) = data_fl(1:nlevs)
          ! *** Deal with sqrt rho ***
          DO z = 1, nlevs
            data_fl(z) = INT_HF (ABC_state % rho(x,z),     &
                                 ABC_state % rho(x,z+1),   &
                                 z, dims)
          END DO
          ! Place back in the ABC_state structure
          ABC_state % rho(x,1:nlevs) = data_fl(1:nlevs)
        END DO

        ! Multiply each of u, v, w by the square-root of density
        ABC_state % u(1:nlongs, 1:nlevs) = ABC_state % u(1:nlongs, 1:nlevs) *  &
                                           ABC_state % rho(1:nlongs, 1:nlevs)
        ABC_state % v(1:nlongs, 1:nlevs) = ABC_state % v(1:nlongs, 1:nlevs) *  &
                                           ABC_state % rho(1:nlongs, 1:nlevs)
        ABC_state % w(1:nlongs, 1:nlevs) = ABC_state % w(1:nlongs, 1:nlevs) *  &
                                           ABC_state % rho(1:nlongs, 1:nlevs)


        ! Fourier transform each (put inside ControlVarS for convenience)
        CALL fft_real2spec ( ABC_state % u(1:nlongs, 1:nlevs),                 &
                           ControlVarS % v1(1:nlongs, 1:nlevs) )
        CALL fft_real2spec ( ABC_state % v(1:nlongs, 1:nlevs),                 &
                           ControlVarS % v2(1:nlongs, 1:nlevs) )
        CALL fft_real2spec ( ABC_state % w(1:nlongs, 1:nlevs),                 &
                           ControlVarS % v3(1:nlongs, 1:nlevs) )

        ! Compute the KE in Fourier space
        ! Deal with the largest scale (real part only)
        EnsAvKE(1, 1:nlevs) = EnsAvKE(1, 1:nlevs) +                            &
                              0.5 * (ControlVarS % v1(1, 1:nlevs) *            &
                                     ControlVarS % v1(1, 1:nlevs) +            &
!
                                     ControlVarS % v2(1, 1:nlevs) *            &
                                     ControlVarS % v2(1, 1:nlevs) +            &
!
                                     ControlVarS % v3(1, 1:nlevs) *            &
                                     ControlVarS % v3(1, 1:nlevs))
        ! Deal with the bulk of the scales (real and imaginary components)
        DO k = 2, nlongs/2
          real_index = 2*k-2
          imag_index = 2*k-1
          EnsAvKE(k, 1:nlevs) = EnsAvKE(k, 1:nlevs) +                          &
                                0.5 * (ControlVarS % v1(real_index, 1:nlevs) * &
                                       ControlVarS % v1(real_index, 1:nlevs) + &
                                       ControlVarS % v1(imag_index, 1:nlevs) * &
                                       ControlVarS % v1(imag_index, 1:nlevs) + &
!
                                       ControlVarS % v2(real_index, 1:nlevs) * &
                                       ControlVarS % v2(real_index, 1:nlevs) + &
                                       ControlVarS % v2(imag_index, 1:nlevs) * &
                                       ControlVarS % v2(imag_index, 1:nlevs) + &
!
                                       ControlVarS % v3(real_index, 1:nlevs) * &
                                       ControlVarS % v3(real_index, 1:nlevs) + &
                                       ControlVarS % v3(imag_index, 1:nlevs) * &
                                       ControlVarS % v3(imag_index, 1:nlevs))
        END DO
        ! Deal with the smallest scale (real part only)
        last_index = nlongs/2+1
        EnsAvKE(last_index, 1:nlevs) = EnsAvKE(last_index, 1:nlevs) +          &
                                0.5 * (ControlVarS % v1(nlongs, 1:nlevs) *     &
                                       ControlVarS % v1(nlongs, 1:nlevs) +     &
!
                                       ControlVarS % v2(nlongs, 1:nlevs) *     &
                                       ControlVarS % v2(nlongs, 1:nlevs) +     &
!
                                       ControlVarS % v3(nlongs, 1:nlevs) *     &
                                       ControlVarS % v3(nlongs, 1:nlevs))
      END DO
    END DO
  END DO

  ! Normalise
  EnsAvKE(1:nlongs/2+1, 1:nlevs) = EnsAvKE(1:nlongs/2+1, 1:nlevs) /            &
                                   ( Neffmems_r * Nens_r )

  ! Write out the KE spectra
  OPEN (12, file='KEspec.dat')

  WRITE (12,*) 'KE spectra for'
  WRITE (12,*) nlongs/2+1
  WRITE (12,*) 'wavenumbers and'
  WRITE (12,*) nlevs
  WRITE (12,*) 'levels'
  WRITE (12,*) '---------------------------------------------------'
  DO z = 1, nlevs
    WRITE (12,*) EnsAvKE(1:nlongs/2+1,z)
  END DO

  WRITE (12,*) '---------------------------------------------------'

  CLOSE (12)

  ! Tidy up
  CALL Deallocate_model_vars (ABC_state)
  DEALLOCATE (EnsAvKE, data_fl)
  CALL Deallocate_CVs(ControlVarS)

END IF






IF (DiagSwitch(5)) THEN
  PRINT *, '========================================================================='
  PRINT *, 'Computing scale dependence of balance'
  PRINT *, '========================================================================='
  ! This is based on code /home/ross/DataAssim/Meso/Covs/Fortran_Ubuntu/Var_HRAA_Smooth4Bal.f90

  CALL Initialise_model_vars (ABC_mean, .FALSE.)
  CALL Initialise_model_vars (ABC_pert, .FALSE.)
  ALLOCATE (GeoCorr(0:nscales, 1:nlevs))
  ALLOCATE (LBCorr(0:nscales, 1:nlevs))
  ALLOCATE (HydCorr(0:nscales, 1:nlevs))
  ALLOCATE (GeoCorr_outer(1:nlongs, 1:nlevs))
  ALLOCATE (LBCorr_outer(1:nlongs, 1:nlevs))
  ALLOCATE (HydCorr_outer(1:nlongs, 1:nlevs))
  ALLOCATE (GeoMean(1:2, 1:nlongs, 1:nlevs))
  ALLOCATE (LBMean(1:2, 1:nlongs, 1:nlevs))
  ALLOCATE (HydMean(1:2, 1:nlongs, 1:nlevs))
  ALLOCATE (Geo(1:2, 1:nlongs, 1:nlevs))
  ALLOCATE (LB(1:2, 1:nlongs, 1:nlevs))
  ALLOCATE (Hyd(1:2, 1:nlongs, 1:nlevs))
  ALLOCATE (GeoVar(1:2, 1:nlongs, 1:nlevs))
  ALLOCATE (LBVar(1:2, 1:nlongs, 1:nlevs))
  ALLOCATE (HydVar(1:2, 1:nlongs, 1:nlevs))
  ALLOCATE (GeoCorr_inner(1:nlongs, 1:nlevs))
  ALLOCATE (LBCorr_inner(1:nlongs, 1:nlevs))
  ALLOCATE (HydCorr_inner(1:nlongs, 1:nlevs))
  ALLOCATE (Scales(0:nscales))

  ! Choose the smoothing scales (number of grid boxes either side of each grid point)
  ! Suggest that the zeroth value is 0
  ! These should be integers
  DO scale = 0, nscales
    Scales(scale) = scale * 50
  END DO


  GeoCorr(0:nscales, 1:nlevs) = 0.0  !Geostrophic
  LBCorr(0:nscales, 1:nlevs)  = 0.0  !Linear balance
  HydCorr(0:nscales, 1:nlevs) = 0.0  !Hydrostatic

  PRINT *, 'There are ', Neffmems, ' effective ensemble members.'

  ! Loop over different horizontal smoothing scales
  DO scale = 0, Nscales

    PRINT *, '==== Computations for scale = ', scale, ' which is ', Scales(scale), ' grid boxes'
    ! Initialise the outer correlations to zero for this scaling
    GeoCorr_outer(1:nlongs, 1:nlevs) = 0.0
    LBCorr_outer(1:nlongs, 1:nlevs)  = 0.0
    HydCorr_outer(1:nlongs, 1:nlevs) = 0.0

    ! Loop over different ensembles
    DO ens = 1, Nens

      PRINT *, 'Computing the ens mean balance quantities for ensemble no ', ens
      ! ==== Compute the ensemble mean of the balance quantities for this ensemble ====
      ! Initialise the means to zero
      GeoMean(1:2, 1:nlongs, 1:nlevs) = 0.0
      LBMean(1:2, 1:nlongs, 1:nlevs)  = 0.0
      HydMean(1:2, 1:nlongs, 1:nlevs) = 0.0

      ! Read-in the ensemble mean for this ensemble
      WRITE (meanfile, '(A,A,I0.3,A)') TRIM(datadirABCperts), '/MeanABC', ens, '.nc'
      !!PRINT *, 'Reading mean state ', TRIM (meanfile)
      CALL Read_state_2d (meanfile, ABC_mean, dims, 1, .TRUE.)
      !!PRINT *, '-- done'

      ! Read-in the ABC perts
      item = 0
      DO mem = 1, NEnsMems
        DO lat = 1, Nlats
          item = item + 1

          ! Read-in a pert found from stage 2
          WRITE (ABCfile, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirABCperts), '/PertABC_Ens', ens, '_Item', item, '.nc'
          !!PRINT *, 'Reading file ', TRIM(ABCfile)
          CALL Read_state_2d (ABCfile, ABC_pert, dims, 1, .TRUE.)
          !!PRINT *, '-- done'
          IF ((ens == 1) .AND. (item == 1)) THEN
            CALL Set_ht_dep_cons (dims)
          END IF

          ! ==== Smooth this field in the horizontal ====
          PRINT *, 'Smoothing fields for mean calculation ...'
          CALL HorizSmooth (ABC_pert, scale, .FALSE., .TRUE., .FALSE., .TRUE., .TRUE.)
          PRINT *, '... done'


          ! ==== Compute geostrophic, linear balance, and hydrostatic related terms ====
          !!PRINT *, 'Computing balance terms for mean calculation ...'
          CALL BalanceTerms (ABC_pert, dims, Geo(1:2, 1:nlongs, 1:nlevs), &
                                             LB (1:2, 1:nlongs, 1:nlevs), &
                                             Hyd(1:2, 1:nlongs, 1:nlevs))
          !!PRINT *, '... done'

          ! Increment the mean
          GeoMean(1:2, 1:nlongs, 1:nlevs) = GeoMean(1:2, 1:nlongs, 1:nlevs) +  &
                                            Geo(1:2, 1:nlongs, 1:nlevs)
          LBMean(1:2, 1:nlongs, 1:nlevs)  = LBMean(1:2, 1:nlongs, 1:nlevs) +   &
                                            LB(1:2, 1:nlongs, 1:nlevs)
          HydMean(1:2, 1:nlongs, 1:nlevs) = HydMean(1:2, 1:nlongs, 1:nlevs) +  &
                                            Hyd(1:2, 1:nlongs, 1:nlevs)

        END DO  ! lat
      END DO  ! mem

      ! Normalise
      GeoMean(1:2, 1:nlongs, 1:nlevs) = GeoMean(1:2, 1:nlongs, 1:nlevs) / Neffmems_r
      LBMean(1:2, 1:nlongs, 1:nlevs) = GeoMean(1:2, 1:nlongs, 1:nlevs)  / Neffmems_r
      HydMean(1:2, 1:nlongs, 1:nlevs) = HydMean(1:2, 1:nlongs, 1:nlevs) / Neffmems_r



      ! ==== Compute the ensemble variances and correlations between the balance quantities for this ensemble ====
      ! Initialise the variances and correlations to zero
      GeoVar(1:2, 1:nlongs, 1:nlevs)   = 0.0
      LBVar(1:2, 1:nlongs, 1:nlevs)    = 0.0
      HydVar(1:2, 1:nlongs, 1:nlevs)   = 0.0

      GeoCorr_inner(1:nlongs, 1:nlevs) = 0.0
      LBCorr_inner(1:nlongs, 1:nlevs)  = 0.0
      HydCorr_inner(1:nlongs, 1:nlevs) = 0.0

      ! Read-in the ABC perts
      item = 0
      DO mem = 1, NEnsMems
        DO lat = 1, Nlats
          item = item + 1

          ! Read-in a pert found from stage 2
          WRITE (ABCfile, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirABCperts), '/PertABC_Ens', ens, '_Item', item, '.nc'
          !!PRINT *, 'Reading file ', TRIM(ABCfile)
          CALL Read_state_2d (ABCfile, ABC_pert, dims, 1, .TRUE.)
          !!PRINT *, '-- done'
          IF ((ens == 1) .AND. (item == 1)) THEN
            CALL Set_ht_dep_cons (dims)
          END IF

          ! ==== Smooth this field in the horizontal ====
          PRINT *, 'Smoothing fields for correlation calculation ...'
          CALL HorizSmooth (ABC_pert, scale, .FALSE., .TRUE., .FALSE., .TRUE., .TRUE.)
          PRINT *, '... done'


          ! ==== Compute geostrophic and hydrostatic related terms ====
          !!PRINT *, 'Computing balance terms for correlation calculation ...'
          CALL BalanceTerms (ABC_pert, dims, Geo(1:2, 1:nlongs, 1:nlevs), &
                                             LB (1:2, 1:nlongs, 1:nlevs), &
                                             Hyd(1:2, 1:nlongs, 1:nlevs))
          !!PRINT *, '... done'


          ! Increment the variances
          GeoVar(1:2, 1:nlongs, 1:nlevs)   = GeoVar(1:2, 1:nlongs, 1:nlevs) +                                   &
                                             (Geo(1:2, 1:nlongs, 1:nlevs) - GeoMean(1:2, 1:nlongs, 1:nlevs)) *  &
                                             (Geo(1:2, 1:nlongs, 1:nlevs) - GeoMean(1:2, 1:nlongs, 1:nlevs))
          LBVar(1:2, 1:nlongs, 1:nlevs)    = LBVar(1:2, 1:nlongs, 1:nlevs) +                                    &
                                             (LB(1:2, 1:nlongs, 1:nlevs) - LBMean(1:2, 1:nlongs, 1:nlevs)) *    &
                                             (LB(1:2, 1:nlongs, 1:nlevs) - LBMean(1:2, 1:nlongs, 1:nlevs))
          HydVar(1:2, 1:nlongs, 1:nlevs)   = HydVar(1:2, 1:nlongs, 1:nlevs) +                                   &
                                             (Hyd(1:2, 1:nlongs, 1:nlevs) - HydMean(1:2, 1:nlongs, 1:nlevs)) *  &
                                             (Hyd(1:2, 1:nlongs, 1:nlevs) - HydMean(1:2, 1:nlongs, 1:nlevs))
          ! Increment the covariances
          GeoCorr_inner(1:nlongs, 1:nlevs) = GeoCorr_inner(1:nlongs, 1:nlevs) +                                 &
                                             (Geo(1, 1:nlongs, 1:nlevs) - GeoMean(1, 1:nlongs, 1:nlevs)) *      &
                                             (Geo(2, 1:nlongs, 1:nlevs) - GeoMean(2, 1:nlongs, 1:nlevs))
          LBCorr_inner(1:nlongs, 1:nlevs)  = LBCorr_inner(1:nlongs, 1:nlevs) +                                  &
                                             (LB(1, 1:nlongs, 1:nlevs) - LBMean(1, 1:nlongs, 1:nlevs)) *        &
                                             (LB(2, 1:nlongs, 1:nlevs) - LBMean(2, 1:nlongs, 1:nlevs))
          HydCorr_inner(1:nlongs, 1:nlevs) = HydCorr_inner(1:nlongs, 1:nlevs) +                                 &
                                             (Hyd(1, 1:nlongs, 1:nlevs) - HydMean(1, 1:nlongs, 1:nlevs)) *      &
                                             (Hyd(2, 1:nlongs, 1:nlevs) - HydMean(2, 1:nlongs, 1:nlevs))
        END DO  ! lat
      END DO  ! mem

      ! Normalise
      GeoVar(1:2, 1:nlongs, 1:nlevs)   = GeoVar(1:2, 1:nlongs, 1:nlevs)   / Neffmems_r
      LBVar(1:2, 1:nlongs, 1:nlevs)    = LBVar(1:2, 1:nlongs, 1:nlevs)    / Neffmems_r
      HydVar(1:2, 1:nlongs, 1:nlevs)   = HydVar(1:2, 1:nlongs, 1:nlevs)   / Neffmems_r
      GeoCorr_inner(1:nlongs, 1:nlevs) = GeoCorr_inner(1:nlongs, 1:nlevs) / Neffmems_r
      LBCorr_inner(1:nlongs, 1:nlevs)  = LBCorr_inner(1:nlongs, 1:nlevs)  / Neffmems_r
      HydCorr_inner(1:nlongs, 1:nlevs) = HydCorr_inner(1:nlongs, 1:nlevs) / Neffmems_r

      ! Compute correlations from covariances
      GeoCorr_inner(1:nlongs, 1:nlevs) = GeoCorr_inner(1:nlongs, 1:nlevs) /                                     &
                                         SQRT(GeoVar(1, 1:nlongs, 1:nlevs) * GeoVar(2, 1:nlongs, 1:nlevs))
      LBCorr_inner(1:nlongs, 1:nlevs)  = LBCorr_inner(1:nlongs, 1:nlevs) /                                      &
                                         SQRT(LBVar(1, 1:nlongs, 1:nlevs) * LBVar(2, 1:nlongs, 1:nlevs))
      HydCorr_inner(1:nlongs, 1:nlevs) = HydCorr_inner(1:nlongs, 1:nlevs) /                                     &
                                         SQRT(HydVar(1, 1:nlongs, 1:nlevs) * HydVar(2, 1:nlongs, 1:nlevs))

      ! Increment the average over ensembles
      GeoCorr_outer(1:nlongs, 1:nlevs) = GeoCorr_outer(1:nlongs, 1:nlevs) + GeoCorr_inner(1:nlongs, 1:nlevs)
      LBCorr_outer(1:nlongs, 1:nlevs)  = LBCorr_outer(1:nlongs, 1:nlevs)  + LBCorr_inner(1:nlongs, 1:nlevs)
      HydCorr_outer(1:nlongs, 1:nlevs) = HydCorr_outer(1:nlongs, 1:nlevs) + HydCorr_inner(1:nlongs, 1:nlevs)

    END DO ! ens

    ! Normalise
    GeoCorr_outer(1:nlongs, 1:nlevs) = GeoCorr_outer(1:nlongs, 1:nlevs) / Nens_r
    LBCorr_outer(1:nlongs, 1:nlevs)  = LBCorr_outer(1:nlongs, 1:nlevs)  / Nens_r
    HydCorr_outer(1:nlongs, 1:nlevs) = HydCorr_outer(1:nlongs, 1:nlevs) / Nens_r

    ! Average over longitude and store
    DO z = 1, nlevs
      GeoCorr(scale, z) = SUM(GeoCorr_outer(1:nlongs, z)) / Nlongs_r
      LBCorr(scale, z)  = SUM(LBCorr_outer(1:nlongs, z))  / Nlongs_r
      HydCorr(scale, z) = SUM(HydCorr_outer(1:nlongs, z)) / Nlongs_r
    END DO


  END DO  ! scale

  ! Output the results
  OPEN (12, file='GeoBalCorr.dat')
  WRITE (12,*) 'Smoothing scales (grid boxes)'
  WRITE (12,*) Scales(0:nscales)
  WRITE (12,*) '-----------------------------'
  DO z = 1, nlevs
    WRITE (12, *) z, GeoCorr(0:nscales,z)
  END DO
  CLOSE (12)

  OPEN (12, file='LinBalCorr.dat')
  WRITE (12,*) 'Smoothing scales (grid boxes)'
  WRITE (12,*) Scales(0:nscales)
  WRITE (12,*) '-----------------------------'
  DO z = 1, nlevs
    WRITE (12, *) z, LBCorr(0:nscales,z)
  END DO
  CLOSE (12)

  OPEN (12, file='HydBalCorr.dat')
  WRITE (12,*) 'Smoothing scales (grid boxes)'
  WRITE (12,*) Scales(0:nscales)
  WRITE (12,*) '-----------------------------'
  DO z = 1, nlevs
    WRITE (12, *) z, HydCorr(0:nscales,z)
  END DO
  CLOSE (12)

  CALL Deallocate_model_vars (ABC_mean)
  CALL Deallocate_model_vars (ABC_pert)
  DEALLOCATE (GeoCorr, LBCorr, HydCorr)
  DEALLOCATE (GeoCorr_outer, LBCorr_outer, HydCorr_outer)
  DEALLOCATE (GeoMean, LBMean, HydMean, Geo, LB, Hyd, GeoVar, LBVar, HydVar)
  DEALLOCATE (GeoCorr_inner, LBCorr_inner, HydCorr_inner)
  DEALLOCATE (Scales)

END IF






IF (DiagSwitch(6)) THEN
  PRINT *, '========================================================================='
  PRINT *, 'Variances of total, balanced, and unbalanced scaled densities'
  PRINT *, '========================================================================='

  CALL Initialise_model_vars (ABC_mean, .FALSE.)
  CALL Initialise_model_vars (ABC_pert, .FALSE.)
  CALL Initialise_CVs(ControlVar, .FALSE.)
  CALL Initialise_CVs(ControlVarS, .FALSE.)

  CALL Initialise_CVs(Ctrl_Cor_with_1_outer, .FALSE.)

  PRINT *, 'There are ', Neffmems, ' effective ensemble members.'
  DO ens = 1, Nens

    ! Read-in the ensemble mean for this ensemble
    WRITE (meanfile, '(A,A,I0.3,A)') TRIM(datadirABCperts), '/MeanABC', ens, '.nc'
    PRINT *, 'Reading mean state ', TRIM (meanfile)
    CALL Read_state_2d (meanfile, ABC_mean, dims, 1, .TRUE.)
    PRINT *, '-- done'

    CALL Initialise_CVs(Ctrl_Cor_with_1_inner, .FALSE.)


    ! Read-in the ABC perts and convert to parameter perts
    item = 0
    DO mem = 1, NEnsMems
      DO lat = 1, Nlats
        item = item + 1

        ! Read-in a pert found from stage 2
        WRITE (ABCfile, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirABCperts), '/PertABC_Ens', ens, '_Item', item, '.nc'
        PRINT *, 'Reading file ', TRIM(ABCfile)
        CALL Read_state_2d (ABCfile, ABC_pert, dims, 1, .TRUE.)
        PRINT *, '-- done'
        IF ((ens == 1) .AND. (item == 1)) THEN
          CALL Set_ht_dep_cons (dims)
        END IF

        ! Pass through the inverse parameter transform
        PRINT *, 'Performing inverse parameter transform'
        CALL U_p_inv (ABC_mean, ControlVar, ABC_pert, CVT % CVT_order,                        &
                      CVT % CVT_param_opt_gb, CVT % CVT_param_opt_hb, CVT % CVT_param_opt_ab,   &
                      CVT % CVT_param_opt_reg, LevMeanBalr,                                     &
                      CVT % Regression(1:nlevs,1:nlevs), dims,                                  &
                      ((ens == 1) .AND. (item ==1)),  &  ! Output diagnostics only if this is met
                      '.')                               ! For diagnostic location
        PRINT *, 'Done'

        ! ===== IMPORTANT =====
        ! At this stage we have:
        !   * total scaled density pert: ABC_pert % r
        !   * unbalanced scaled density pert: ControlVar % v3
        ! Use ControlVar % v1 as a space to place the total scaled density pert
        ControlVar % v1(1:nlongs, 1:nlevs) = ABC_pert % r(1:nlongs, 1:nlevs)

        ! Use ControlVar % v2 as a space to place the balanced scaled density pert
        ControlVar % v2(1:nlongs, 1:nlevs) = ControlVar % v1(1:nlongs, 1:nlevs) - ControlVar % v3(1:nlongs, 1:nlevs)

        ! The variances will be placed inside Ctrl_Cor_with_1_inner.  These will be a function of wavenumber.
        ! There are only n/2 + 1 unique scales, so some of this array will be unused.
        ! ===== IMPORTANT =====


        ! ========================================================================
        ! Fourier transform the scaled densities
        PRINT *, 'Performing FT'
        CALL fft_real2spec ( ControlVar % v1(1:nlongs, 1:nlevs), ControlVarS % v1(1:nlongs, 1:nlevs) ) ! total
        CALL fft_real2spec ( ControlVar % v2(1:nlongs, 1:nlevs), ControlVarS % v2(1:nlongs, 1:nlevs) ) ! bal
        CALL fft_real2spec ( ControlVar % v3(1:nlongs, 1:nlevs), ControlVarS % v3(1:nlongs, 1:nlevs) ) ! unbal
        PRINT *, 'Done'


        ! ========================================================================
        ! Contribute to the variances for total, balanced, and unbalanced scaled densities
        PRINT *, 'Computing variance contributions'

        ! Deal with the largest scale first
        ! total scaled density
        Ctrl_Cor_with_1_inner % v1(1, 1:nlevs) = Ctrl_Cor_with_1_inner % v1(1, 1:nlevs) + &
                                                 ControlVarS % v1(1, 1:nlevs) * ControlVarS % v1(1, 1:nlevs)
        ! balanced scaled density
        Ctrl_Cor_with_1_inner % v2(1, 1:nlevs) = Ctrl_Cor_with_1_inner % v2(1, 1:nlevs) + &
                                                 ControlVarS % v2(1, 1:nlevs) * ControlVarS % v2(1, 1:nlevs)
        ! unbalanced scaled density
        Ctrl_Cor_with_1_inner % v3(1, 1:nlevs) = Ctrl_Cor_with_1_inner % v3(1, 1:nlevs) + &
                                                 ControlVarS % v3(1, 1:nlevs) * ControlVarS % v3(1, 1:nlevs)

        ! Deal with the bulk of the scales
        DO k = 2, nlongs/2
          real_index = 2*k-2
          imag_index = 2*k-1
          ! total scaled density
          Ctrl_Cor_with_1_inner % v1(k, 1:nlevs) = Ctrl_Cor_with_1_inner % v1(k, 1:nlevs) + &
                                                   ControlVarS % v1(real_index, 1:nlevs) * ControlVarS % v1(real_index, 1:nlevs) + &
                                                   ControlVarS % v1(imag_index, 1:nlevs) * ControlVarS % v1(imag_index, 1:nlevs)
          ! balanced scaled density
          Ctrl_Cor_with_1_inner % v2(k, 1:nlevs) = Ctrl_Cor_with_1_inner % v2(k, 1:nlevs) + &
                                                   ControlVarS % v2(real_index, 1:nlevs) * ControlVarS % v2(real_index, 1:nlevs) + &
                                                   ControlVarS % v2(imag_index, 1:nlevs) * ControlVarS % v2(imag_index, 1:nlevs)
          ! unbalanced scaled density
          Ctrl_Cor_with_1_inner % v3(k, 1:nlevs) = Ctrl_Cor_with_1_inner % v3(k, 1:nlevs) + &
                                                   ControlVarS % v3(real_index, 1:nlevs) * ControlVarS % v3(real_index, 1:nlevs) + &
                                                   ControlVarS % v3(imag_index, 1:nlevs) * ControlVarS % v3(imag_index, 1:nlevs)
        END DO

        ! Deal with the smallest scale
        last_index = nlongs/2+1
        ! total scaled density
        Ctrl_Cor_with_1_inner % v1(last_index, 1:nlevs) = Ctrl_Cor_with_1_inner % v1(last_index, 1:nlevs) + &
                                                 ControlVarS % v1(nlongs, 1:nlevs) * ControlVarS % v1(nlongs, 1:nlevs)
        ! balanced scaled density
        Ctrl_Cor_with_1_inner % v2(last_index, 1:nlevs) = Ctrl_Cor_with_1_inner % v2(last_index, 1:nlevs) + &
                                                 ControlVarS % v2(nlongs, 1:nlevs) * ControlVarS % v2(nlongs, 1:nlevs)
        ! unbalanced scaled density
        Ctrl_Cor_with_1_inner % v3(last_index, 1:nlevs) = Ctrl_Cor_with_1_inner % v3(last_index, 1:nlevs) + &
                                                 ControlVarS % v3(nlongs, 1:nlevs) * ControlVarS % v3(nlongs, 1:nlevs)

        PRINT *, 'Done'
        ! ========================================================================

      END DO
    END DO

    ! Normalise to give variances
    PRINT *, 'Normalising'
    CALL Div_CV_cons(Ctrl_Cor_with_1_inner, Neffmems_r)
    PRINT *, 'Done'


    ! Add on to the outer for averages
    PRINT *, 'Averaging over ensembles'
    CALL Add_CVs (Ctrl_Cor_with_1_outer, Ctrl_Cor_with_1_inner)
    PRINT *, 'Done'


  END DO


  ! Normalise correlations
  PRINT *, 'Normalising'
  CALL Div_CV_cons(Ctrl_Cor_with_1_outer, Nens_r)
  PRINT *, 'Done'


  ! Write out the variances as a function of wavenumber

  CALL Write_one_field ('TotrVarSpec.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_1_outer % v1(1:nlongs/2+1,1:nlevs), 'TotrVarSpec')
  CALL Write_one_field ('BalrVarSpec.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_1_outer % v2(1:nlongs/2+1,1:nlevs), 'BalrVarSpec')
  CALL Write_one_field ('UnbrVarSpec.nc', nlongs/2+1, nlevs, Ctrl_Cor_with_1_outer % v3(1:nlongs/2+1,1:nlevs), 'UnbrVarSpec')



  PRINT *, 'Done'

  ! Tidy up
  CALL Deallocate_model_vars (ABC_mean)
  CALL Deallocate_model_vars (ABC_pert)
  CALL Deallocate_CVs(ControlVar)
  CALL Deallocate_CVs(ControlVarS)
  CALL Deallocate_CVs(Ctrl_Cor_with_1_outer)
  CALL Deallocate_CVs(Ctrl_Cor_with_1_inner)

END IF








! Tidy up
IF (DiagSwitch(1) .OR. DiagSwitch(2) .OR. DiagSwitch(3)) CALL Deallocate_CVT (CVT)
CALL Deallocate_dims (dims)

IF (ALLOCATED(fft_wsave_x)) DEALLOCATE(fft_wsave_x)
IF (ALLOCATED(fft_work_x)) DEALLOCATE(fft_work_x)


PRINT *, 'Program finished'

END PROGRAM Util_SpecialDiags



! ==========================================================================
SUBROUTINE HorizSmooth1Field (Field, gridpoints)
! Subroutine to do horizontal smoothing of a single field

USE DefConsTypes, ONLY :         &
    ZREAL8,                      &
    nlongs,                      &
    nlevs

IMPLICIT NONE

! Subroutine parameters
REAL(ZREAL8), INTENT(INOUT) :: Field(1:nlongs, 1:nlevs)
INTEGER,      INTENT(IN)    :: gridpoints

! Local variables
REAL(ZREAL8), ALLOCATABLE   :: smoothed(:)
REAL(ZREAL8)                :: recip_2gridpoints, total
INTEGER                     :: x, xp, xpp, z

ALLOCATE (smoothed(1:nlongs))

recip_2gridpoints = 1.0 / REAL(2*gridpoints+1)

DO z = 1, nlevs
  DO x = 1, nlongs
    ! Do smoothing
    total = 0.0
    DO xp = x-gridpoints, x+gridpoints
      IF (xp < 1) THEN
        xpp = xp + nlongs
      ELSE IF (xp > nlongs) THEN
        xpp = xp - nlongs
      ELSE
        xpp = xp
      END IF
      total = total + Field(xpp, z)
    END DO
    smoothed(x) = total * recip_2gridpoints
  END DO

  ! Put smoothed result back into main array
  Field(1:nlongs, z) = smoothed(1:nlongs)
END DO

DEALLOCATE (smoothed)

END SUBROUTINE HorizSmooth1Field


! ==========================================================================
SUBROUTINE HorizSmooth (state, gridpoints,    &
                        smooth_u, smooth_v, smooth_w, smooth_rp, smooth_bp)
! Horizontally smooth the requested fields

USE DefConsTypes, ONLY :         &
    nlongs,                      &
    nlevs,                       &
    ABC_type

IMPLICIT NONE

! Subroutine parameters
TYPE(ABC_type), INTENT(INOUT) :: state
INTEGER,        INTENT(IN)    :: gridpoints
LOGICAL,        INTENT(IN)    :: smooth_u
LOGICAL,        INTENT(IN)    :: smooth_v
LOGICAL,        INTENT(IN)    :: smooth_w
LOGICAL,        INTENT(IN)    :: smooth_rp
LOGICAL,        INTENT(IN)    :: smooth_bp

IF (gridpoints > 0) THEN
  IF (smooth_u)  CALL HorizSmooth1Field (state % u(1:nlongs, 1:nlevs), gridpoints)
  IF (smooth_v)  CALL HorizSmooth1Field (state % v(1:nlongs, 1:nlevs), gridpoints)
  IF (smooth_w)  CALL HorizSmooth1Field (state % w(1:nlongs, 1:nlevs), gridpoints)
  IF (smooth_rp) CALL HorizSmooth1Field (state % r(1:nlongs, 1:nlevs), gridpoints)
  IF (smooth_bp) CALL HorizSmooth1Field (state % b(1:nlongs, 1:nlevs), gridpoints)
END IF

END SUBROUTINE HorizSmooth


! ==========================================================================
SUBROUTINE BalanceTerms (state, dims, Geo, LB, Hyd)
! Compute the following terms
! Geo part 1: C * d r_prime / dx
!     part 2: f * v
! LB  part 1: C * d2 r_prime / dx2
!     part 2: f * dv / dx
! Hyd part 1: C * d rho / dz
!     part 2: b_prime

USE DefConsTypes, ONLY :         &
    ZREAL8,                      &
    nlongs,                      &
    nlevs,                       &
    ABC_type,                    &
    dims_type,                   &
    C,                           &
    f,                           &
    recip2dx2,                   &
    recip2dx,                    &
    recipdx,                     &
    half

IMPLICIT NONE

! Subroutine parameters
TYPE(ABC_type),  INTENT(IN)    :: state
TYPE(dims_type), INTENT(IN)    :: dims
REAL(ZREAL8),    INTENT(INOUT) :: Geo(1:2, 1:nlongs, 1:nlevs)
REAL(ZREAL8),    INTENT(INOUT) :: LB(1:2, 1:nlongs, 1:nlevs)
REAL(ZREAL8),    INTENT(INOUT) :: Hyd(1:2, 1:nlongs, 1:nlevs)

! Local variables
INTEGER                       :: x, z

DO z = 1, nlevs
  DO x = 1, nlongs

    ! Geostrophic terms will be at a u point
    Geo(1, x, z) = C * (state % r(x+1,z) - state % r(x,z)) * recipdx
    Geo(2, x, z) = f * (state % v(x,z) + state % v(x+1,z)) * half

    ! Linear balance terms will be at a rho point
    LB(1, x, z) = C * (state % r(x+1,z) + state % r(x-1,z) - 2.0 * state % r(x,z)) * recip2dx2
    LB(2, x, z) = f * (state % v(x+1,z) - state % v(x-1,z)) * recip2dx

    ! Hydrostatic terms will be at b point
    Hyd(1, x, z) = C * (state % r(x,z+1) - state % r(x,z)) / &
                       (dims % half_levs(z+1) - dims % half_levs(z))
    Hyd(2, x, z) = state % b(x,z)

  END DO
END DO

END SUBROUTINE BalanceTerms
