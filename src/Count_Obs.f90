SUBROUTINE Count_Obs (Obs, count)

!*********************************************************************************
!*                                                                               *
!*  Run through observations structure and count number of observations          *
!*                                                                               *
!*  Obs                    - pointer to the start of the observation linked list *
!*  count                  - (output) the number of observations                 *
!*                                                                               *
!*   J. Lee, 1.4da 02-11-2021                                                    *
!*                                                                               *
!*********************************************************************************

USE DefConsTypes, ONLY :  &
  Obs_type


IMPLICIT NONE

! Parameters
!-----------
TYPE(Obs_type), POINTER, INTENT(IN) :: Obs
INTEGER,                 INTENT(OUT):: count

! Local Variables
!----------------
TYPE(Obs_type), POINTER             :: thisob
INTEGER                             :: obc


! Loop through observations and count
thisob => Obs

obc = 0
DO
  IF (ASSOCIATED(thisob)) THEN
    obc = obc + 1
    thisob => thisob % next
  ELSE
    EXIT
  END IF
END DO

count = obc

END SUBROUTINE Count_Obs
