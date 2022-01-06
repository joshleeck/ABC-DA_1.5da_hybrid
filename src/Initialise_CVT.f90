SUBROUTINE Initialise_CVT (CVT)

! Initialise control variable transform structure
! Does not alter the ABC parameters and options though

USE DefConsTypes, ONLY :     &
    CVT_type,                &
    nlongs, nlevs

IMPLICIT NONE

TYPE(CVT_type), INTENT(INOUT)   :: CVT

  ! Standard deviations of the 6 control parameters
  ! -----------------------------------------------
  IF (.NOT.ALLOCATED(CVT % sigma1)) ALLOCATE (CVT % sigma1(1:nlongs, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % sigma2)) ALLOCATE (CVT % sigma2(1:nlongs, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % sigma3)) ALLOCATE (CVT % sigma3(1:nlongs, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % sigma4)) ALLOCATE (CVT % sigma4(1:nlongs, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % sigma5)) ALLOCATE (CVT % sigma5(1:nlongs, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % sigma6)) ALLOCATE (CVT % sigma6(1:nlongs, 1:nlevs))
  CVT % sigma1(1:nlongs, 1:nlevs) = 0.0
  CVT % sigma2(1:nlongs, 1:nlevs) = 0.0
  CVT % sigma3(1:nlongs, 1:nlevs) = 0.0
  CVT % sigma4(1:nlongs, 1:nlevs) = 0.0
  CVT % sigma5(1:nlongs, 1:nlevs) = 0.0
  CVT % sigma6(1:nlongs, 1:nlevs) = 0.0


  ! Vertical modes of the 6 control parameters
  ! ------------------------------------------

  IF (.NOT.ALLOCATED(CVT % VertMode1)) ALLOCATE (CVT % VertMode1(1:nlevs, 1:nlevs, 1:nlongs/2+1))
  IF (.NOT.ALLOCATED(CVT % VertMode2)) ALLOCATE (CVT % VertMode2(1:nlevs, 1:nlevs, 1:nlongs/2+1))
  IF (.NOT.ALLOCATED(CVT % VertMode3)) ALLOCATE (CVT % VertMode3(1:nlevs, 1:nlevs, 1:nlongs/2+1))
  IF (.NOT.ALLOCATED(CVT % VertMode4)) ALLOCATE (CVT % VertMode4(1:nlevs, 1:nlevs, 1:nlongs/2+1))
  IF (.NOT.ALLOCATED(CVT % VertMode5)) ALLOCATE (CVT % VertMode5(1:nlevs, 1:nlevs, 1:nlongs/2+1))
  IF (.NOT.ALLOCATED(CVT % VertMode6)) ALLOCATE (CVT % VertMode6(1:nlevs, 1:nlevs, 1:nlongs/2+1))
  CVT % VertMode1(1:nlevs, 1:nlevs, 1:nlongs/2+1) = 0.0
  CVT % VertMode2(1:nlevs, 1:nlevs, 1:nlongs/2+1) = 0.0
  CVT % VertMode3(1:nlevs, 1:nlevs, 1:nlongs/2+1) = 0.0
  CVT % VertMode4(1:nlevs, 1:nlevs, 1:nlongs/2+1) = 0.0
  CVT % VertMode5(1:nlevs, 1:nlevs, 1:nlongs/2+1) = 0.0
  CVT % VertMode6(1:nlevs, 1:nlevs, 1:nlongs/2+1) = 0.0

  ! Vertical eigenvalues of the 6 control parameters (these are actually the square-roots)
  ! --------------------------------------------------------------------------------------
  IF (.NOT.ALLOCATED(CVT % VertEV1)) ALLOCATE (CVT % VertEV1(1:nlevs, 1:nlongs/2+1))
  IF (.NOT.ALLOCATED(CVT % VertEV2)) ALLOCATE (CVT % VertEV2(1:nlevs, 1:nlongs/2+1))
  IF (.NOT.ALLOCATED(CVT % VertEV3)) ALLOCATE (CVT % VertEV3(1:nlevs, 1:nlongs/2+1))
  IF (.NOT.ALLOCATED(CVT % VertEV4)) ALLOCATE (CVT % VertEV4(1:nlevs, 1:nlongs/2+1))
  IF (.NOT.ALLOCATED(CVT % VertEV5)) ALLOCATE (CVT % VertEV5(1:nlevs, 1:nlongs/2+1))
  IF (.NOT.ALLOCATED(CVT % VertEV6)) ALLOCATE (CVT % VertEV6(1:nlevs, 1:nlongs/2+1))
  CVT % VertEV1(1:nlevs, 1:nlongs/2+1) = 0.0
  CVT % VertEV2(1:nlevs, 1:nlongs/2+1) = 0.0
  CVT % VertEV3(1:nlevs, 1:nlongs/2+1) = 0.0
  CVT % VertEV4(1:nlevs, 1:nlongs/2+1) = 0.0
  CVT % VertEV5(1:nlevs, 1:nlongs/2+1) = 0.0
  CVT % VertEV6(1:nlevs, 1:nlongs/2+1) = 0.0

  ! Horizontal eigenvalues of the 6 control parameters (these are actually the square-roots)
  ! ----------------------------------------------------------------------------------------
  IF (.NOT.ALLOCATED(CVT % HorizEV1)) ALLOCATE (CVT % HorizEV1(1:nlongs/2+1, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % HorizEV2)) ALLOCATE (CVT % HorizEV2(1:nlongs/2+1, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % HorizEV3)) ALLOCATE (CVT % HorizEV3(1:nlongs/2+1, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % HorizEV4)) ALLOCATE (CVT % HorizEV4(1:nlongs/2+1, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % HorizEV5)) ALLOCATE (CVT % HorizEV5(1:nlongs/2+1, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % HorizEV6)) ALLOCATE (CVT % HorizEV6(1:nlongs/2+1, 1:nlevs))
  CVT % HorizEV1(1:nlongs/2+1, 1:nlevs) = 0.0
  CVT % HorizEV2(1:nlongs/2+1, 1:nlevs) = 0.0
  CVT % HorizEV3(1:nlongs/2+1, 1:nlevs) = 0.0
  CVT % HorizEV4(1:nlongs/2+1, 1:nlevs) = 0.0
  CVT % HorizEV5(1:nlongs/2+1, 1:nlevs) = 0.0
  CVT % HorizEV6(1:nlongs/2+1, 1:nlevs) = 0.0

  ! Regression data for balanced density
  ! ------------------------------------
  IF (.NOT.ALLOCATED(CVT % Cov_rbalrbal)) ALLOCATE (CVT % Cov_rbalrbal(1:nlevs, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % Cov_rtotrbal)) ALLOCATE (CVT % Cov_rtotrbal(1:nlevs, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % Regression))   ALLOCATE (CVT % Regression(1:nlevs, 1:nlevs))
  CVT % Cov_rbalrbal(1:nlevs, 1:nlevs) = 0.0
  CVT % Cov_rtotrbal(1:nlevs, 1:nlevs) = 0.0
  CVT % Regression(1:nlevs, 1:nlevs)   = 0.0

END SUBROUTINE Initialise_CVT
