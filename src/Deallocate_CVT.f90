SUBROUTINE Deallocate_CVT (CVT)

USE DefConsTypes, ONLY :     &
    CVT_type

IMPLICIT NONE

TYPE(CVT_type), INTENT(INOUT)   :: CVT

DEALLOCATE (CVT % sigma1, CVT % sigma2, CVT % sigma3)
DEALLOCATE (CVT % sigma4, CVT % sigma5, CVT % sigma6)

DEALLOCATE (CVT % VertMode1, CVT % VertMode2, CVT % VertMode3)
DEALLOCATE (CVT % VertMode4, CVT % VertMode5, CVT % VertMode6)

DEALLOCATE (CVT % VertEV1, CVT % VertEV2, CVT % VertEV3)
DEALLOCATE (CVT % VertEV4, CVT % VertEV5, CVT % VertEV6)

DEALLOCATE (CVT % HorizEV1, CVT % HorizEV2, CVT % HorizEV3)
DEALLOCATE (CVT % HorizEV4, CVT % HorizEV5, CVT % HorizEV6)

DEALLOCATE (CVT % Cov_rbalrbal, CVT % Cov_rtotrbal, CVT % Regression)

END SUBROUTINE Deallocate_CVT
