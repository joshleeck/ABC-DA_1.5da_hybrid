SUBROUTINE Deallocate_aCVT (aCVT)

USE DefConsTypes, ONLY :     &
    aCVT_type

IMPLICIT NONE

TYPE(aCVT_type), INTENT(INOUT)   :: aCVT

DEALLOCATE (aCVT % VertMode1, aCVT % VertMode2, aCVT % VertMode3)
DEALLOCATE (aCVT % VertMode4, aCVT % VertMode5)

DEALLOCATE (aCVT % HorizMode1, aCVT % HorizMode2, aCVT % HorizMode3)
DEALLOCATE (aCVT % HorizMode4, aCVT % HorizMode5)

DEALLOCATE (aCVT % VertEV1, aCVT % VertEV2, aCVT % VertEV3)
DEALLOCATE (aCVT % VertEV4, aCVT % VertEV5)

DEALLOCATE (aCVT % HorizEV1, aCVT % HorizEV2, aCVT % HorizEV3)
DEALLOCATE (aCVT % HorizEV4, aCVT % HorizEV5)

END SUBROUTINE Deallocate_aCVT
