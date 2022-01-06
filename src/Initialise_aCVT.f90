SUBROUTINE Initialise_aCVT (aCVT)

! Initialise alpha control variable transform structure

USE DefConsTypes, ONLY :     &
    aCVT_type,                &
    nlongs, nlevs

IMPLICIT NONE

TYPE(aCVT_type), INTENT(INOUT)   :: aCVT

  ! Vertical modes of the 5 variables in model space
  ! -----------------------------------------------
  IF (.NOT.ALLOCATED(aCVT % VertMode1)) ALLOCATE (aCVT % VertMode1(1:nlevs, 1:nlevs, 1:nlongs))
  IF (.NOT.ALLOCATED(aCVT % VertMode2)) ALLOCATE (aCVT % VertMode2(1:nlevs, 1:nlevs, 1:nlongs))
  IF (.NOT.ALLOCATED(aCVT % VertMode3)) ALLOCATE (aCVT % VertMode3(1:nlevs, 1:nlevs, 1:nlongs))
  IF (.NOT.ALLOCATED(aCVT % VertMode4)) ALLOCATE (aCVT % VertMode4(1:nlevs, 1:nlevs, 1:nlongs))
  IF (.NOT.ALLOCATED(aCVT % VertMode5)) ALLOCATE (aCVT % VertMode5(1:nlevs, 1:nlevs, 1:nlongs))
  aCVT % VertMode1(1:nlevs, 1:nlevs, 1:nlongs) = 0.0
  aCVT % VertMode2(1:nlevs, 1:nlevs, 1:nlongs) = 0.0
  aCVT % VertMode3(1:nlevs, 1:nlevs, 1:nlongs) = 0.0
  aCVT % VertMode4(1:nlevs, 1:nlevs, 1:nlongs) = 0.0
  aCVT % VertMode5(1:nlevs, 1:nlevs, 1:nlongs) = 0.0

  ! Vertical eigenvalues of the 5 variables in model space (these are actually the square-roots)
  ! -----------------------------------------------
  IF (.NOT.ALLOCATED(aCVT % VertEV1)) ALLOCATE (aCVT % VertEV1(1:nlevs, 1:nlongs))
  IF (.NOT.ALLOCATED(aCVT % VertEV2)) ALLOCATE (aCVT % VertEV2(1:nlevs, 1:nlongs))
  IF (.NOT.ALLOCATED(aCVT % VertEV3)) ALLOCATE (aCVT % VertEV3(1:nlevs, 1:nlongs))
  IF (.NOT.ALLOCATED(aCVT % VertEV4)) ALLOCATE (aCVT % VertEV4(1:nlevs, 1:nlongs))
  IF (.NOT.ALLOCATED(aCVT % VertEV5)) ALLOCATE (aCVT % VertEV5(1:nlevs, 1:nlongs))
  aCVT % VertEV1(1:nlevs, 1:nlongs)            = 0.0
  aCVT % VertEV2(1:nlevs, 1:nlongs)            = 0.0
  aCVT % VertEV3(1:nlevs, 1:nlongs)            = 0.0
  aCVT % VertEV4(1:nlevs, 1:nlongs)            = 0.0
  aCVT % VertEV5(1:nlevs, 1:nlongs)            = 0.0

  ! Horizontal modes of the 5 variables in model space
  ! -----------------------------------------------
  IF (.NOT.ALLOCATED(aCVT % HorizMode1)) ALLOCATE (aCVT % HorizMode1(1:nlongs, 1:nlongs, 1:nlevs))
  IF (.NOT.ALLOCATED(aCVT % HorizMode2)) ALLOCATE (aCVT % HorizMode2(1:nlongs, 1:nlongs, 1:nlevs))
  IF (.NOT.ALLOCATED(aCVT % HorizMode3)) ALLOCATE (aCVT % HorizMode3(1:nlongs, 1:nlongs, 1:nlevs))
  IF (.NOT.ALLOCATED(aCVT % HorizMode4)) ALLOCATE (aCVT % HorizMode4(1:nlongs, 1:nlongs, 1:nlevs))
  IF (.NOT.ALLOCATED(aCVT % HorizMode5)) ALLOCATE (aCVT % HorizMode5(1:nlongs, 1:nlongs, 1:nlevs))
  aCVT % HorizMode1(1:nlongs, 1:nlongs, 1:nlevs) = 0.0
  aCVT % HorizMode2(1:nlongs, 1:nlongs, 1:nlevs) = 0.0
  aCVT % HorizMode3(1:nlongs, 1:nlongs, 1:nlevs) = 0.0
  aCVT % HorizMode4(1:nlongs, 1:nlongs, 1:nlevs) = 0.0
  aCVT % HorizMode5(1:nlongs, 1:nlongs, 1:nlevs) = 0.0

  ! Horizontal eigenvalues of the 5 variables in model space (these are actually the square-roots)
  ! -----------------------------------------------
  IF (.NOT.ALLOCATED(aCVT % HorizEV1)) ALLOCATE (aCVT % HorizEV1(1:nlongs, 1:nlevs))
  IF (.NOT.ALLOCATED(aCVT % HorizEV2)) ALLOCATE (aCVT % HorizEV2(1:nlongs, 1:nlevs))
  IF (.NOT.ALLOCATED(aCVT % HorizEV3)) ALLOCATE (aCVT % HorizEV3(1:nlongs, 1:nlevs))
  IF (.NOT.ALLOCATED(aCVT % HorizEV4)) ALLOCATE (aCVT % HorizEV4(1:nlongs, 1:nlevs))
  IF (.NOT.ALLOCATED(aCVT % HorizEV5)) ALLOCATE (aCVT % HorizEV5(1:nlongs, 1:nlevs))
  aCVT % HorizEV1(1:nlongs, 1:nlevs)           = 0.0
  aCVT % HorizEV2(1:nlongs, 1:nlevs)           = 0.0
  aCVT % HorizEV3(1:nlongs, 1:nlevs)           = 0.0
  aCVT % HorizEV4(1:nlongs, 1:nlevs)           = 0.0
  aCVT % HorizEV5(1:nlongs, 1:nlevs)           = 0.0

END SUBROUTINE Initialise_aCVT
