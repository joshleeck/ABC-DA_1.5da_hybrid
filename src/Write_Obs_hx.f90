SUBROUTINE Write_Obs_hx ( Filename, Obs, EnsNum, NEnsMems, count, hx, R )


!****************************************************************************************
!*                                                                                      *
!*  Output the ensemble background in observation space                                 *
!*                                                                                      *
!*  Filename          - Output filename                                                 *
!*  Obs               - Pointer to the first observation (linked list)                  *
!*  EnsNum            - Ensemble member index                                           *
!*  NEnsMems          - Number of ensemble members                                      *
!*  count             - Number of observations                                          *
!*  hx                - Array containing values of h(x)                                 *
!*  R                 - Observation error covariance matrix                             *
!*                                                                                      *
!*  J. Lee, 1.5da_hybrid 04-10-2021                                                     *
!*                                                                                      *
!****************************************************************************************

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  Obs_type

IMPLICIT NONE

! Parameters
!-----------
CHARACTER(LEN=320),      INTENT(IN)   :: Filename
TYPE(Obs_type), POINTER, INTENT(IN)   :: Obs
INTEGER,                 INTENT(IN)   :: EnsNum
INTEGER,                 INTENT(IN)   :: NEnsMems
INTEGER,                 INTENT(IN)   :: count
REAL(ZREAL8),            INTENT(INOUT):: hx(1:count,1:NEnsMems)
REAL(ZREAL8),            INTENT(INOUT):: R(1:count,1:count)

! Local Variables
!----------------
TYPE(Obs_type), POINTER             :: thisob
INTEGER                             :: obc

  ! Open file
  OPEN (13, file=Filename)

  thisob  => Obs
  obc = 1
  DO
    IF (ASSOCIATED(thisob)) THEN
      WRITE (13,'(A)') '-----------------------------------------'
      WRITE (13,'(A18,I10)')   'Observation No  : ', thisob % obnumber_thisfile

      ! y_ref is the observation operator applied on a model state, i.e. hx
      WRITE (13,'(A18,E14.5)') 'y_ref           : ', thisob % y_ref

      ! Store value for respective ensemble member and observation number in hx and R
      hx(obc,EnsNum) = thisob % y_ref
      R(obc,obc)     = thisob % stddev * thisob % stddev 
      obc = obc + 1

      thisob => thisob % next
    ELSE
      EXIT
    END IF

  END DO

  CLOSE (13)

END SUBROUTINE Write_Obs_hx
