FUNCTION Fn_GC_corr (zdiff, vScale_alpha) RESULT (corr)

! Code to compute the Gaspari-Cohn correlation as a function of difference in z (model height)
! Correlation function given by equation (4.10) in Gaspari and Cohn,
! Q. J. R. Meteorol. Soc., 125, 723-757. vScale_alpha is the separation at
! which the correlation reaches zero, when 2c <= |zdiff| as shown in equation.
! 
! This is also now used for horizontal correlation, i.e. function of difference in horizontal distance
! J Lee. 22-03-2021 

USE DefConsTypes, ONLY :  &
  ZREAL8

IMPLICIT NONE

! Function arguments

REAL(ZREAL8),    INTENT(IN)     :: zdiff
REAL(ZREAL8),    INTENT(IN)     :: vScale_alpha

REAL(ZREAL8)                   :: corr         ! the result

! Local variables
REAL(ZREAL8)                    :: c
REAL(ZREAL8)                    :: temp_5
REAL(ZREAL8)                    :: temp_4
REAL(ZREAL8)                    :: temp_3
REAL(ZREAL8)                    :: temp_2
REAL(ZREAL8)                    :: temp_1
REAL(ZREAL8)                    :: temp_0
REAL(ZREAL8)                    :: temp_inv1
REAL(ZREAL8)                    :: absz_divc



! Here we pick a=1/2 in the paper in the equation so that we get equation 4.10
c = vScale_alpha

IF (zdiff == 0.0) THEN

  corr = 1.0

ELSE IF (ABS(zdiff) <= 2.0*c) THEN

  ! Set up coefficients for polynomial in piecewise rational function:
  IF (ABS(zdiff) <= c) THEN
    temp_5    = -0.25
    temp_4    =  0.5
    temp_3    =  0.625
    temp_2    = -5.0/3.0
    temp_1    =  0.0
    temp_0    =  1.0
    temp_inv1 =  0.0
  ELSE
    temp_5    =  1.0/12.0
    temp_4    = -0.5
    temp_3    =  0.625
    temp_2    =  5.0/3.0
    temp_1    = -5.0
    temp_0    =  4.0
    temp_inv1 = -2.0/3.0
  END IF

  absz_divc = ABS(zdiff/c)

  corr = temp_5    * absz_divc**5  &
       + temp_4    * absz_divc**4  &
       + temp_3    * absz_divc**3  &
       + temp_2    * absz_divc**2  &
       + temp_1    * absz_divc     &
       + temp_0    * 1.0           &
       + temp_inv1 * 1.0/absz_divc

ELSE

  corr = 0.0

END IF


END FUNCTION Fn_GC_corr

