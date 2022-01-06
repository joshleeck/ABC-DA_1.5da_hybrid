SUBROUTINE U_p_inv (LS, ControlVar, ModelVar, order,             &
                    option_gb, option_hb, option_ab, option_reg, &
                    LevMeanBalr,                                 &
                    Regression, dims,                            &
                    diags_flag, outputdir)

! Code to perform the parameter cvt: ControlVar = U_p^-1 ModelVar

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  ABC_type,               &
  CV_type,                &
  dims_type,              &
  nlongs,                 &
  nlevs,                  &
  f, C

IMPLICIT NONE

INCLUDE "Boundaries_CV.interface"


TYPE(ABC_type),   INTENT(IN)    :: LS
TYPE(CV_type),    INTENT(INOUT) :: ControlVar
TYPE(ABC_type),   INTENT(IN)    :: ModelVar
INTEGER,          INTENT(IN)    :: order
INTEGER,          INTENT(IN)    :: option_gb
INTEGER,          INTENT(IN)    :: option_hb
INTEGER,          INTENT(IN)    :: option_ab
INTEGER,          INTENT(IN)    :: option_reg
LOGICAL,          INTENT(IN)    :: LevMeanBalr  ! Set to make level mean r balanced
REAL(ZREAL8),     INTENT(IN)    :: Regression(1:nlevs, 1:nlevs)
TYPE(dims_type),  INTENT(IN)    :: dims
LOGICAL,          INTENT(IN)    :: diags_flag
CHARACTER(LEN=*), INTENT(IN)    :: outputdir  ! For diags if needed


INTEGER                         :: x, z
REAL(ZREAL8), ALLOCATABLE       :: r_b(:,:)
REAL(ZREAL8), ALLOCATABLE       :: b_b(:,:)
REAL(ZREAL8), ALLOCATABLE       :: w_b(:,:)
CHARACTER(LEN=320)              :: state_file
REAL(ZREAL8)                    :: LevMeanr, LevMeanrb, deltar, deltapsi



! *******************************************************************
! *** The same numbering system in the comments is used as in U_p ***
! *******************************************************************


ALLOCATE(r_b(1:nlongs, 1:nlevs))
ALLOCATE(b_b(1:nlongs, 1:nlevs))
ALLOCATE(w_b(1:nlongs, 1:nlevs))

IF ((order == 1).OR.(order == 2)) THEN
  ! Traditional kind of control variable transform
  ! ----------------------------------------------

  ! ----------------------------------------------
  ! 9. Compute the tracer
  ControlVar % v6(1:nlongs,1:nlevs) = ModelVar % tracer(1:nlongs,1:nlevs)
  CALL Boundaries_CV (ControlVar, set_6=.TRUE.)
  ! ----------------------------------------------


  ! ----------------------------------------------
  IF (option_ab == 1) THEN
    ! 7. Compute the balanced w
    CALL Anbalw (LS,                                   &
                 ModelVar % u(0:nlongs+1,0:nlevs+1),   &
                 ModelVar % r(0:nlongs+1,0:nlevs+1),   &
                 w_b(1:nlongs,1:nlevs),                &
                 dims)

    IF (diags_flag) THEN
      state_file = TRIM(outputdir) // '/w_b.nc'
      CALL Write_one_field (state_file, nlongs, nlevs, &
                            w_b(1:nlongs,1:nlevs), 'w_b')
    END IF

    ! 8. Compute the unbalanced w
    ControlVar % v5(1:nlongs,1:nlevs) = ModelVar % w(1:nlongs,1:nlevs) - w_b(1:nlongs,1:nlevs)
  ELSE
    ! No anelastic balance relations used
    ControlVar % v5(1:nlongs,1:nlevs) = ModelVar % w(1:nlongs,1:nlevs)
  END IF
  CALL Boundaries_CV (ControlVar, set_5=.TRUE.)
  ! ----------------------------------------------


  ! ----------------------------------------------
  IF ((option_hb == 1) .OR. (option_hb == 2)) THEN
    IF (option_hb == 1) THEN
      ! 5. Compute the balanced b - see Eq (19) of model paper
      CALL HydroBal_b (ModelVar % r(1:nlongs,0:nlevs+1),   &
                       b_b(1:nlongs,1:nlevs),              &
                       dims)
    ELSE
      ! Statistical balance
    END IF

      IF (diags_flag) THEN
        state_file = TRIM(outputdir) // '/b_b.nc'
        CALL Write_one_field (state_file, nlongs, nlevs, &
                              b_b(1:nlongs,1:nlevs), 'b_b')
      END IF

    ! 6. Compute the unbalanced b
    ControlVar % v4(1:nlongs,1:nlevs) = ModelVar % b(1:nlongs,1:nlevs) - b_b(1:nlongs,1:nlevs)
  ELSE
    ! No hydrostatic balance relations used
    ControlVar % v4(1:nlongs,1:nlevs) = ModelVar % b(1:nlongs,1:nlevs)
  END IF
  CALL Boundaries_CV (ControlVar, set_4=.TRUE.)
  ! ----------------------------------------------


  ! ----------------------------------------------
  ! 1. Compute psi and chi from u and v
  CALL Helmholtz_inv (ControlVar % v1(1:nlongs,1:nlevs),     &
                      ControlVar % v2(1:nlongs,1:nlevs),     &
                      ModelVar % u(0:nlongs+1,1:nlevs),      &
                      ModelVar % v(0:nlongs+1,1:nlevs))
  CALL Boundaries_CV (ControlVar, set_1=.TRUE., set_2=.TRUE.)
  ! ----------------------------------------------


  ! ----------------------------------------------
  IF ((option_gb == 1) .OR. (option_gb == 2)) THEN
    ! 2. Compute the balanced r (r_b) from psi (not known up to a constant on each level)
    IF (option_gb ==1) THEN
      ! Analytical balance
      CALL LinearBal_r (ControlVar % v1(0:nlongs+1,1:nlevs),   &
                        r_b(1:nlongs, 1:nlevs))

      IF (LevMeanBalr) THEN
        DO z = 1, nlevs
          LevMeanr        = SUM(ModelVar % r(1:nlongs,z)) / REAL(nlongs)
          LevMeanrb       = SUM(r_b(1:nlongs,z)) / REAL(nlongs)
          deltar          = LevMeanr - LevMeanrb
          ! Transfer the level mean of total scaled density to balanced scaled density
          r_b(1:nlongs,z) = r_b(1:nlongs,z) + deltar
          ! Modify the streamfunction to be in balance with this
          deltapsi        = C * deltar / f
          ControlVar % v1(0:nlongs+1,z) = ControlVar % v1(0:nlongs+1,z) + deltapsi
        END DO
      END IF

    ELSE
      ! Statistical balance
    END IF

    IF (diags_flag) THEN
      state_file = TRIM(outputdir) // '/r_b_preregress.nc'
      CALL Write_one_field (state_file, nlongs, nlevs, &
                            r_b(1:nlongs,1:nlevs), 'r_b_preregress')
    END IF

    IF (option_reg == 1) THEN
      ! 3. Perform vertical regression
      DO x = 1, nlongs
        r_b(x,1:nlevs) = MATMUL(Regression(1:nlevs, 1:nlevs), r_b(x,1:nlevs))
      END DO
      IF (diags_flag) THEN
        state_file = TRIM(outputdir) // '/r_b_postregress.nc'
        CALL Write_one_field (state_file, nlongs, nlevs, &
                              r_b(1:nlongs,1:nlevs), 'r_b_postregress')
      END IF
    END IF

    ! 4. Compute the unbalanced r
    ControlVar % v3(1:nlongs,1:nlevs) = ModelVar % r(1:nlongs,1:nlevs) - r_b(1:nlongs,1:nlevs)

  ELSE
    ! No geostrophic balance relations used
    ControlVar % v3(1:nlongs,1:nlevs) = ModelVar % r(1:nlongs,1:nlevs)
  END IF
  CALL Boundaries_CV (ControlVar, set_3=.TRUE.)
  ! ----------------------------------------------



  ! Control variables are psi, chi, (unbalanced) r, (unbalanced) b, (unbalanced) w, tracer
ELSE
  ! Control variable transform as REP's thesis
  ! ----------------------------------------------

END IF

DEALLOCATE(r_b, b_b, w_b)


END SUBROUTINE U_p_inv
