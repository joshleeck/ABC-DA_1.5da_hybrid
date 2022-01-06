PROGRAM Util_TestConvScaleDens

!*****************************************************
!*   Code to test the routine ConvScaleDens          *
!*   Ross Bannister, June 2020                       *
!*                                                   *
!*****************************************************


! Use Statements
!===============

USE DefConsTypes, ONLY :         &
    ZREAL8,                      &
    dims_type,                   &
    ABC_type,                    &
    nlongs, nlevs,               &
    datadirABC_in,               &
    LS_file1, LS_file2,          &
    datadirABC_out,              &
    output_ABC_file,             &
    fft_wsave_x,                 &
    fft_work_x

IMPLICIT NONE

! Declare variables
!==========================
TYPE(dims_type)           :: dims
TYPE(ABC_type)            :: ABC_ref, dABC
CHARACTER(LEN=320)        :: ABC_input_file1, ABC_input_file2, drpcs_outputfile, incfile
REAL(ZREAL8), ALLOCATABLE :: rpls(:,:)
REAL(ZREAL8), ALLOCATABLE :: rpcs(:,:)
REAL(ZREAL8), ALLOCATABLE :: rhols(:,:)
REAL(ZREAL8), ALLOCATABLE :: bpcs(:,:)
REAL(ZREAL8), ALLOCATABLE :: drpls(:,:)
REAL(ZREAL8), ALLOCATABLE :: dbpcs(:,:)
REAL(ZREAL8), ALLOCATABLE :: drpcs(:,:)
INTEGER                   :: z
REAL(ZREAL8)              :: mean


PRINT*, '*************************************************************************'
PRINT*, 'Running Util_TestConvScaleDens'
PRINT*, '*************************************************************************'

! Read namelist
CALL SetOptions

ABC_input_file1  = TRIM(datadirABC_in) // '/' // TRIM(LS_file1)
ABC_input_file2  = TRIM(datadirABC_in) // '/' // TRIM(LS_file2)
incfile          = TRIM(datadirABC_out) // '/inc.nc'
drpcs_outputfile = TRIM(datadirABC_out) // '/drpcs.nc'

PRINT *, 'ABC_intput_file1 = ', TRIM(ABC_input_file1)
PRINT *, 'ABC_intput_file2 = ', TRIM(ABC_input_file2)
PRINT *, 'incfile          = ', TRIM(incfile)
PRINT *, 'drpcs_outputfile = ', TRIM(drpcs_outputfile)

CALL Initialise_dims (dims)
ALLOCATE (rpls(0:nlongs+1, 0:nlevs+1))
ALLOCATE (rpcs(0:nlongs+1, 0:nlevs+1))
ALLOCATE (rhols(0:nlongs+1, 0:nlevs+1))
ALLOCATE (bpcs(0:nlongs+1, 0:nlevs+1))
ALLOCATE (drpls(0:nlongs+1, 0:nlevs+1))
ALLOCATE (dbpcs(0:nlongs+1, 0:nlevs+1))
ALLOCATE (drpcs(0:nlongs+1, 0:nlevs+1))
CALL Initialise_model_vars (ABC_ref, .FALSE.)
CALL Initialise_model_vars (dABC, .FALSE.)

PRINT*, 'Reading in data file 1 ...'
CALL Read_state_2d (ABC_input_file1, ABC_ref, dims, -1, .TRUE.)
PRINT*, '-- done'
PRINT*, 'Reading in data file 2 ...'
CALL Read_state_2d (ABC_input_file2, dABC, dims, -1, .TRUE.)
PRINT*, '-- done'

! Set some commonly-used constants
CALL Set_ht_dep_cons (dims)


PRINT *, 'Subtracting for incremental state'
! Let dABC = dABC - ABC_ref
CALL Subtract_model_vars(dABC, ABC_ref, .TRUE.)
PRINT*, '-- done'


PRINT *, 'Writing out the increment'
CALL Write_state_2d (incfile, dABC, dims, 1, 0, 1, .TRUE.)
PRINT*, '-- done'


! Calculate rpls (the global level mean in this test) and rpcs (the field minus the level mean)
PRINT *, 'Computing rpls and rpcs'
DO z = 0, nlevs+1
  mean = SUM(ABC_ref % r(1:nlongs,z)) / REAL(nlongs)
  rpls(0:nlongs+1,z) = mean
  rpcs(0:nlongs+1,z) = ABC_ref % r(0:nlongs+1,z) - mean
END DO
PRINT*, '-- done'

PRINT *, 'Computing rhols'
! Calculate rhols (the global level mean in this test)
DO z = 0, nlevs+1
  mean = SUM(ABC_ref % rho(1:nlongs,z)) / REAL(nlongs)
  rhols(0:nlongs+1,z) = mean
END DO
PRINT*, '-- done'

PRINT *, 'Computing bpcs'
! Calculate bpcs (the field minus the level mean)
DO z = 0, nlevs+1
  mean = SUM(ABC_ref % b(1:nlongs,z)) / REAL(nlongs)
  bpcs(0:nlongs+1,z) = ABC_ref % b(0:nlongs+1,z) - mean
END DO
PRINT*, '-- done'

PRINT *, 'Computing drpls'
! Calculate drpls (the global level mean in this test)
DO z = 0, nlevs+1
  mean = SUM(dABC % r(1:nlongs,z)) / REAL(nlongs)
  drpls(0:nlongs+1,z) = mean
END DO
PRINT*, '-- done'

PRINT *, 'Computing dbpcs'
! Calculate dbpls (the global level mean in this test)
DO z = 0, nlevs+1
  mean = SUM(dABC % b(1:nlongs,z)) / REAL(nlongs)
  dbpcs(0:nlongs+1,z) = dABC % b(0:nlongs+1,z) - mean
END DO
PRINT*, '-- done'


PRINT *, '****** Comuting drpcs *******'
! Call the routine to compute the convective-scale density pert
CALL ConvScaleDens (ABC_ref % u(0:nlongs+1, 0:nlevs+1),         &   ! Reference state
                    ABC_ref % v(0:nlongs+1, 0:nlevs+1),         &   ! Reference state
                    ABC_ref % w(0:nlongs+1, 0:nlevs+1),         &   ! Reference state
                    rpls(0:nlongs+1, 0:nlevs+1),                &   ! Reference state
                    rpcs(0:nlongs+1, 0:nlevs+1),                &   ! Reference state
                    rhols(0:nlongs+1, 0:nlevs+1),               &   ! Reference state
                    bpcs(0:nlongs+1, 0:nlevs+1),                &   ! Reference state
                    dABC % u(0:nlongs+1, 0:nlevs+1),            &   ! Input inc
                    dABC % v(0:nlongs+1, 0:nlevs+1),            &   ! Input inc
                    dABC % w(0:nlongs+1, 0:nlevs+1),            &   ! Input inc
                    drpls(0:nlongs+1, 0:nlevs+1),               &   ! Input inc
                    dbpcs(0:nlongs+1, 0:nlevs+1),               &   ! Input inc
                    drpcs(0:nlongs+1, 0:nlevs+1),               &   ! Output inc
                    dims )                                          ! Meta data
PRINT*, '-- done'


! Output the result
PRINT *, 'Outputting drpcs'
CALL Write_one_field (drpcs_outputfile,                         &
                      nlongs, nlevs,                            &
                      drpcs(1:nlongs, 1:nlevs),                 &
                      'drpcs')
PRINT*, '-- done'


! Tidy up
CALL Deallocate_dims (dims)
DEALLOCATE (rpls, rpcs, rhols, bpcs, drpls, dbpcs, drpcs)
CALL Deallocate_model_vars (ABC_ref)
CALL Deallocate_model_vars (dABC)
IF (ALLOCATED(fft_wsave_x)) DEALLOCATE(fft_wsave_x)
IF (ALLOCATED(fft_work_x)) DEALLOCATE(fft_work_x)

END PROGRAM Util_TestConvScaleDens
