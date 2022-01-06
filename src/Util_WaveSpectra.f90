PROGRAM Util_WaveSpectra

!*****************************************************
!*   Code to run the non-linear ABC model            *
!*   and to compute wave spectra                     *
!*                                                   *
!*   R. Petrie,    2.0:  10-06-2011                  *
!*   R. Petrie,    3.0:  30-07-2013                  *
!*   R. Bannister, 3.1:  30-07-2013                  *
!*   R. Bannister, 1.4da 22-10-2017                  *
!*   R. Bannister, 1.5da 17-06/2020                  *
!*                                                   *
!*****************************************************


! Use Statements
!===============

USE DefConsTypes, ONLY :         &
    ZREAL8,                      &
    nlongs,                      &
    dims_type,                   &
    ABC_type,                    &
    datadirABC_in,               &
    datadirABC_out,              &
    init_ABC_file,               &
    output_ABC_file,             &
    ntimesteps,                  &
    ndumps,                      &
    diagnostics_file


IMPLICIT NONE

! Declare variables
!==========================
TYPE(dims_type)           :: dims
TYPE(ABC_type)            :: ABC_data
REAL(ZREAL8), ALLOCATABLE :: spacetime(:,:)
INTEGER                   :: store_level = 30

CHARACTER(LEN=320)        :: ABC_init_file, ABC_output_file, ABC_diags_file
INTEGER                   :: x, t
REAL(ZREAL8), ALLOCATABLE :: fftx(:,:)            ! After doing fft in longitude
REAL(ZREAL8), ALLOCATABLE :: fftxt(:,:)           ! After doing fft in longitude and time
INTEGER                   :: ierr, kx, kt, qRex, qImx, qRet, qImt
REAL(ZREAL8)              :: RealPart, ImagPart
REAL(ZREAL8), ALLOCATABLE :: SpecDens(:,:)

! Variables to do with fft in longitudinal direction
INTEGER                   :: fft_worklen_x
REAL(ZREAL8), ALLOCATABLE :: fft_wsave_x(:)
REAL(ZREAL8), ALLOCATABLE :: fft_work_x(:)
! Variables to do with fft in time
INTEGER                   :: fft_worklen_t
REAL(ZREAL8), ALLOCATABLE :: fft_wsave_t(:)
REAL(ZREAL8), ALLOCATABLE :: fft_work_t(:)


! Read namelist
CALL SetOptions

CALL Initialise_dims (dims)
CALL Initialise_model_vars (ABC_data, .FALSE.)


PRINT*, '*************************************************************************'
PRINT*, 'Running Master_RunNLModel'
PRINT*, '*************************************************************************'

ABC_init_file   = TRIM(datadirABC_in)  // '/' // TRIM(init_ABC_file)
ABC_output_file = TRIM(datadirABC_out) // '/' // TRIM(output_ABC_file)
ABC_diags_file  = TRIM(datadirABC_out) // '/' // TRIM(diagnostics_file)

! Set state to zero
CALL Initialise_model_vars (ABC_data, .FALSE.)
CALL Initialise_dims (dims)

! Read in preprocessed UM data store in ABC_data
PRINT*, 'Reading in processed data ...'
CALL Read_state_2d (ABC_init_file, ABC_data, dims, -1, .TRUE.)
PRINT*, '-- done'

! Set some commonly-used constants
CALL Set_ht_dep_cons (dims)



ALLOCATE (spacetime(1:nlongs, 0:ntimesteps))

CALL ABC_NL_ModelDriverWave ( ABC_data, dims, ntimesteps, ndumps,      &
                              ABC_output_file, ABC_diags_file,         &
                              store_level, spacetime(1:nlongs, 0:ntimesteps))


! Output the space/time plot
CALL Write_one_field ('SpaceTime.nc', nlongs, ntimesteps+1, &
                      spacetime(1:nlongs, 0:ntimesteps), 'w')


! Set-up stuff for ffts
PRINT *, 'Setting-up ffts'
fft_worklen_x = 2*nlongs
ALLOCATE (fft_wsave_x(1:fft_worklen_x))
ALLOCATE (fft_work_x(1:nlongs))
fft_worklen_t = 2*(ntimesteps+1)
ALLOCATE (fft_wsave_t(1:fft_worklen_t))
ALLOCATE (fft_work_t(1:ntimesteps+1))

CALL rfft1i (nlongs, fft_wsave_x, fft_worklen_x, ierr)
CALL rfft1i (ntimesteps+1, fft_wsave_t, fft_worklen_t, ierr)

! Do FFT in x-direction
PRINT *, 'Doing fft in x-direction'
ALLOCATE (fftx(1:nlongs, 0:ntimesteps))
fftx(1:nlongs, 0:ntimesteps) = spacetime(1:nlongs, 0:ntimesteps)
DEALLOCATE (spacetime)
DO t = 0, ntimesteps
  CALL rfft1f (nlongs, 1,                                           &
               fftx(1:nlongs,t), nlongs,                            &
               fft_wsave_x, fft_worklen_x, fft_work_x, nlongs, ierr)
  IF (ierr /= 0) THEN
    PRINT *, 'Error with rfft1f in x t = ', t
    STOP
  END IF
END DO
DEALLOCATE (fft_wsave_x, fft_work_x)


! Do FFT in t-direction
PRINT *, 'Doing fft in t-direction'
ALLOCATE (fftxt(1:nlongs, 0:ntimesteps))
fftxt(1:nlongs, 0:ntimesteps) = fftx(1:nlongs, 0:ntimesteps)
DEALLOCATE (fftx)
DO x = 1, nlongs
  CALL rfft1f (ntimesteps+1, 1,                                                &
               fftxt(x,0:ntimesteps), ntimesteps+1,                            &
               fft_wsave_t, fft_worklen_t, fft_work_t, ntimesteps+1, ierr)
  IF (ierr /= 0) THEN
    PRINT *, 'Error with rfft1f in t x = ', x
    STOP
  END IF
END DO
DEALLOCATE (fft_wsave_t, fft_work_t)

! Output the spectrum
PRINT *, 'Outputting spectrum'
CALL Write_one_field ('Spec.nc', nlongs, ntimesteps+1, &
                      fftxt(1:nlongs, 0:ntimesteps), 'Spec')



! Calculate the spectral density
PRINT *, 'Calculating spectral density'
ALLOCATE (SpecDens(0:nlongs/2, 0:(ntimesteps+1)/2))
SpecDens(0:nlongs/2, 0:(ntimesteps+1)/2) = 0.0

DO kx = 1, nlongs/2 - 1
  qRex = 2*kx - 1  ! array index of real part with wn kx
  qImx = 2*kx      ! array index of imag part with wn kx
  DO kt = 1, (ntimesteps+1)/2 - 1
    qRet = 2*kt - 2  ! array index of real part with wn kt (minus 1 because index starts at 0)
    qImt = 2*kt - 1  ! array index of imag part with wn kt (minus 1 because index starts at 0)

    RealPart = (fftxt(qRex,qRet) - fftxt(qImx,qImt)) / 2.0
    ImagPart = (fftxt(qRex,qImt) + fftxt(qImx,qRet)) / 2.0

    SpecDens(kx, kt) = SQRT(RealPart * RealPart + ImagPart * ImagPart)
  END DO
END DO

DEALLOCATE (fftxt)

! Output the spectral density
PRINT *, 'Outputting spectral density'
CALL Write_one_field ('SpecDens.nc', nlongs/2+1, (ntimesteps+1)/2+1, &
                      SpecDens(0:nlongs/2, 0:(ntimesteps+1)/2), 'SpecDens')


PRINT *, 'Tidying up'
DEALLOCATE (SpecDens)

CALL Deallocate_dims (dims)
CALL Deallocate_model_vars (ABC_data)

END PROGRAM Util_WaveSpectra
