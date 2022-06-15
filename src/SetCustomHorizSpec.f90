SUBROUTINE SetCustomHorizSpec (hScale, dims, CVT)

! Code to use custom horizontal spectrum in the static covariance

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  dims_type,              &
  CVT_type,               &
  CV_type,                &
  nlongs,                 &
  nlevs, dx               


IMPLICIT NONE

REAL(ZREAL8),    INTENT(IN)     :: hScale
TYPE(dims_type), INTENT(IN)     :: dims
TYPE(CVT_type),  INTENT(INOUT)  :: CVT

INTEGER                      :: k, x1, x2, lev, nover2p1, real_index, imag_index
REAL(ZREAL8)                 :: xdiff1, xdiff2, xdiff
REAL(ZREAL8), ALLOCATABLE    :: Lh(:,:), fftIn(:,:)
REAL(ZREAL8), EXTERNAL       :: Fn_GC_corr
TYPE(CV_type)                :: Interim

ALLOCATE (Lh(1:nlongs, 1:nlongs))
ALLOCATE (fftIn(1:nlongs,1:nlevs))
CALL Initialise_CVs (Interim, .FALSE.)

! Compute Gaspari Cohn/Gaussian function
DO x1 = 1, nlongs
  DO x2 = 1, nlongs
    ! Need to account for 'two distances' between any two points because of periodic BC and construct circulant Lh 
    xdiff1 = ABS(x2 - x1)*dx
    IF (xdiff1 <= 0.5*dims % longs_u(nlongs)) THEN
      !Lh(x1,x2) = Fn_GC_corr(xdiff1, hScale)
      Lh(x1,x2) = EXP(-0.5*(xdiff1/hScale)**2)
    ELSE
      IF (x2 > x1) THEN
        xdiff2 = (x1 + (nlongs-x2))*dx
      ELSE
        xdiff2 = (x2 + (nlongs-x1))*dx
      END IF
    !Lh(x1,x2) = Fn_GC_corr(xdiff2, hScale)
    Lh(x1,x2) = EXP(-0.5*(xdiff2/hScale)**2)
    END IF
  END DO
END DO

nover2p1 = nlongs/2 + 1


! Do transform from real space to spectral space, apply for each vertical mode
DO lev = 1, nlevs
  fftIn(:,lev) = Lh(nlongs/2,1:nlongs)
END DO
CALL fft_real2spec (fftIn(1:nlongs,1:nlevs), Interim % v1(1:nlongs,1:nlevs))

! Deal with the largest scale
CVT % HorizEV1(1,1:nlevs) = Interim % v1(1,1:nlevs) * Interim % v1(1,1:nlevs)

! Deal with the intermediate scales
DO k = 2, nlongs/2
  real_index = 2*k-2
  imag_index = 2*k-1
  CVT % HorizEV1(k,1:nlevs) = Interim % v1(real_index,1:nlevs) * Interim % v1(real_index,1:nlevs) + &
                              Interim % v1(imag_index,1:nlevs) * Interim % v1(imag_index,1:nlevs)
END DO 

! Deal with the smallest scale
CVT % HorizEV1(nover2p1,1:nlevs) = Interim % v1(nlongs,1:nlevs) * Interim % v1(nlongs,1:nlevs)

! Square root
CVT % HorizEV1(1:nlongs/2+1,1:nlevs) = SQRT(CVT % HorizEV1(1:nlongs/2+1,1:nlevs))


DEALLOCATE (Lh)
DEALLOCATE (fftIn)

CALL Deallocate_CVs (Interim)
 

END SUBROUTINE SetCustomHorizSpec
