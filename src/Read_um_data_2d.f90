SUBROUTINE Read_um_data_2d (umdata, filename, latitude)

!**********************************
!* Subroutine to read a latitude  *
!* slice of um netcdf output into *
!* a compound data type um_data   *
!*                                *
!* R. Petrie                      *
!* version 2                      *
!* 06/06/2011                     *
!* RNB - allow for different res  *
!* 04/06/2020                     *
!**********************************

USE DefConsTypes, ONLY :  &
    ZREAL8,               &
    UM_type,              &
    nlongs, nlevs

IMPLICIT NONE

INCLUDE '/scratch/singadm/pkg_hera/netcdf/4.3.3/gnu/4.9.2/include/netcdf.inc'

!Declare parameters
!------------------
TYPE(UM_type),     INTENT(INOUT) :: umdata
CHARACTER (LEN=*), INTENT(IN)    :: filename
INTEGER,           INTENT(IN)    :: latitude

! Declare local variables
!-------------------------
INTEGER                          :: ncid, ierr
INTEGER                          :: dimidLongs_u, dimidLongs_v, dimidhalf_levs, dimidfull_levs
INTEGER                          :: varidLongs_u, varidLongs_v, varidhalf_levs, varidfull_levs
INTEGER                          :: varidu, varidv, varidw, variddensity, varidtheta
INTEGER                          :: varidorog, varidexpres
INTEGER                          :: startA(1), countA(1), startB(4), countB(4), z, x
INTEGER                          :: ierr1, ierr2, ierr3, ierr4, ierr5, ierr6, ierr7
REAL(ZREAL8)                     :: temp_x(nlongs)
REAL(ZREAL8), ALLOCATABLE        :: fieldx_um(:), fieldz_um(:), fieldxz_um(:,:)
REAL(ZREAL8)                     :: axisx(1:2), axisz(1:2)
REAL(ZREAL8)                     :: fracx, fracz
INTEGER                          :: nlongs_um, nlevs_um, lowerx, upperx, lowerz, upperz
LOGICAL                          :: same_grid

! Function
REAL(ZREAL8)                     :: Interpolate1D, Interpolate2D

! Open the netCDF file
!----------------------
ierr = NF_OPEN(TRIM(filename), NF_NOWRITE, ncid)
IF ( ierr .NE. 0 ) THEN
  PRINT*, ' *** Error opening file ***'
  PRINT*, 'Filename: ', TRIM(filename)
  PRINT*, ierr, NF_STRERROR(ierr)
  STOP
ENDIF

!Get the dimension ids
!-----------------------
ierr1 = NF_INQ_DIMID(ncid, 'x', dimidLongs_u)
ierr2 = NF_INQ_DIMID(ncid, 'x_1', dimidLongs_v)
ierr3 = NF_INQ_DIMID(ncid, 'hybrid_ht_3', dimidhalf_levs)
ierr4 = NF_INQ_DIMID(ncid, 'hybrid_ht_2', dimidfull_levs)
IF ( (ierr1 .NE. 0) .OR. (ierr2 .NE. 0) .OR. (ierr3 .NE. 0) .OR. (ierr4 .NE. 0) ) THEN
  PRINT*, '***Error getting dimension ids ***'
  PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
  PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
  PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
  PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
  STOP
!ELSE
! PRINT*, ' Dimension ids ok'
ENDIF


! Find out how many longitudes and levels are in the file
! -------------------------------------------------------

ierr1 = NF_INQ_DIMLEN (ncid, dimidLongs_u, nlongs_um)
ierr2 = NF_INQ_DIMLEN (ncid, dimidhalf_levs, nlevs_um)
IF ( (ierr1 .NE. 0) .OR. (ierr2 .NE. 0) ) THEN
  PRINT*, '***Error getting dimension lengths ***'
  PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
  PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
  STOP
!ELSE
! PRINT*, ' Dimension lengths read ok'
ENDIF


PRINT *, 'Reporting from UM file:'
PRINT *, 'nlongs_um = ', nlongs_um
PRINT *, 'nlevs_um  = ', nlevs_um

! Hard-wire to use 60 levels even though there may be more than 60 levels
nlevs_um = 60

! Is the UM grid the same as the new grid?
same_grid = ((nlongs_um .EQ. nlongs) .AND. (nlevs_um .EQ. nlevs))


! Get the dimension variable ids
!---------------------------------
ierr1 = NF_INQ_VARID(ncid, 'x', varidLongs_u)
ierr2 = NF_INQ_VARID(ncid, 'x_1', varidLongs_v)
ierr3 = NF_INQ_VARID(ncid, 'hybrid_ht_3', varidhalf_levs)
ierr4 = NF_INQ_VARID(ncid, 'hybrid_ht_2', varidfull_levs)
IF ( (ierr1 .NE. 0) .OR. (ierr2 .NE. 0) .OR. (ierr3 .NE. 0) .OR. (ierr4 .NE. 0) ) THEN
  PRINT*, '***Error getting dimension variable ids ***'
  PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
  PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
  PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
  PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
  PRINT*,'ierr5', ierr5, NF_STRERROR(ierr5)
  STOP
!ELSE
!  PRINT*, ' Dimension variable ids ok'
ENDIF


! Get the variable ids
!----------------------
ierr1 = NF_INQ_VARID(ncid, 'u', varidu)
ierr2 = NF_INQ_VARID(ncid, 'v', varidv)
ierr3 = NF_INQ_VARID(ncid, 'dz_dt', varidw)
ierr4 = NF_INQ_VARID(ncid, 'unspecified', variddensity)
ierr5 = NF_INQ_VARID(ncid, 'theta', varidtheta)
ierr6 = NF_INQ_VARID(ncid, 'ht', varidorog)
ierr7 = NF_INQ_VARID(ncid, 'field7', varidexpres)
IF ( (ierr1 .NE. 0) .OR. (ierr2 .NE. 0) .OR. (ierr3 .NE. 0) .OR. (ierr4 .NE. 0) .OR. (ierr5 .NE. 0)&
      .OR. (ierr6 .NE. 0) .OR. (ierr7 .NE. 0) ) THEN
  PRINT*, '***Error getting variable ids ***'
  PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
  PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
  PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
  PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
  PRINT*,'ierr5', ierr5, NF_STRERROR(ierr5)
  PRINT*,'ierr6', ierr6, NF_STRERROR(ierr6)
  PRINT*,'ierr7', ierr6, NF_STRERROR(ierr6)
  STOP
!ELSE
!  PRINT*, ' Variable ids ok'
ENDIF





IF (same_grid) THEN

  ! THE UM GRID IS THE SAME AS THE ABC MODEL GRID SPECIFIED
  ! =======================================================
  PRINT *, 'Using the same grid as the UM data'

  ! Get the dimension data from the file
  !-------------------------------------
  ! Longitudinal distances
  !-----------------------
  startA(1) = 1
  countA(1) = nlongs

  ierr2 = NF_GET_VARA_DOUBLE (ncid, varidLongs_u, startA, countA, umdata % longs_u(1:nlongs))
  ierr3 = NF_GET_VARA_DOUBLE (ncid, varidLongs_v, startA, countA, umdata % longs_v(1:nlongs))

  ! Level Heights
  !--------------
  startA(1) = 1
  countA(1) = nlevs+1
  ierr4 = NF_GET_VARA_DOUBLE (ncid, varidfull_levs, startA, countA, umdata % full_levs(0:nlevs))
  ierr5 = NF_GET_VARA_DOUBLE (ncid, varidhalf_levs, startA, countA, umdata % half_levs(1:nlevs+1))

  IF ( (ierr2 .NE. 0) .OR. (ierr3 .NE. 0) .OR. &
       (ierr4 .NE. 0) .OR. (ierr5 .NE. 0) ) THEN
    PRINT*, '***Error getting dimension data ***'
    PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
    PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
    PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
    PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
    PRINT*,'ierr5', ierr5, NF_STRERROR(ierr5)
    STOP
  ! ELSE
  !  PRINT*, ' Dimension data ok'
  ENDIF

  ! Get the main data from the file
  !--------------------------------
  startB(1) = 1
  countB(1) = nlongs    ! All longitude points
  startB(2) = latitude  ! Selected latitude slice
  countB(2) = 1         !lats
  startB(3) = 1
  countB(3) = 1         !levs
  startB(4) = 1
  countB(4) = 1         !time

  ! u
  !----
  DO z=1, nlevs
    startB(3)=z
    ierr = NF_GET_VARA_DOUBLE (ncid, varidu, startB, countB,  temp_x)
    umdata % u(1:nlongs,z) = temp_x(1:nlongs)
  END DO
  ierr1 = ierr
  !OPEN (15, FILE='u.dat')
  !DO x = 1, nlongs
  !  WRITE (15,*) umdata % longs_u(x), umdata % u(x,30)
  !END DO
  !CLOSE (15)

  ! v
  !----
  DO z=1, nlevs
    startB(3)=z
    ierr2 = NF_GET_VARA_DOUBLE (ncid, varidv, startB, countB, temp_x)
    umdata % v(1:nlongs, z) = temp_x(1:nlongs)
  END DO
  !OPEN (15, FILE='v.dat')
  !DO x = 1, nlongs
  !  WRITE (15,*) umdata % longs_v(x), umdata % v(x,30)
  !END DO
  !CLOSE (15)

  ! w
  !----
  DO z=1, nlevs+1
    startB(3) = z
    ierr3 = NF_GET_VARA_DOUBLE (ncid, varidw, startB, countB, temp_x)
    umdata % w(1:nlongs, z-1) = temp_x(1:nlongs)
  ENDDO

  ! density
  !---------
  DO z=1, nlevs
    startB(3)=z
    ierr4 = NF_GET_VARA_DOUBLE (ncid, variddensity, startB, countB, temp_x)
    umdata % density(1:nlongs, z) = temp_x(1:nlongs)
  ENDDO

  ! theta
  !-------
  DO z=1, nlevs
    startB(3)= z
    ierr5 = NF_GET_VARA_DOUBLE (ncid, varidtheta, startB, countB, temp_x)
    umdata % theta(1:nlongs, z) = temp_x(1:nlongs)
  ENDDO

  ! exner presure
  !----------------
  DO z=1, nlevs + 1
    startB(3) = z
    ierr6 = NF_GET_VARA_DOUBLE (ncid, varidexpres, startB, countB, temp_x)
    umdata % exner_pressure(1:nlongs, z) = temp_x(1:nlongs)
  ENDDO

  ! orographic height
  !--------------------
  DO z=1,1
    startB(3)=z
    ierr7 = NF_GET_VARA_DOUBLE (ncid, varidorog, startB, countB, umdata % orog_height(1:nlongs))
  ENDDO

  IF ( (ierr1 .NE. 0) .OR. (ierr2 .NE. 0) .OR. (ierr3 .NE. 0) .OR. &
       (ierr4 .NE. 0) .OR. (ierr5 .NE. 0) .OR. (ierr6 .NE. 0) .OR. (ierr7 .NE. 0) ) THEN
    PRINT*, '***Error getting main variable data ***'
    PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
    PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
    PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
    PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
    PRINT*,'ierr5', ierr5, NF_STRERROR(ierr5)
    PRINT*,'ierr6', ierr6, NF_STRERROR(ierr6)
    PRINT*,'ierr7', ierr6, NF_STRERROR(ierr6)
    STOP
  ! ELSE
  !  PRINT*, ' Main variable data ok'
  ENDIF


ELSE

  ! THE UM GRID IS NOT THE SAME AS THE ABC MODEL GRID SPECIFIED
  ! ===========================================================
  PRINT *, 'Using a different grid as the UM data'

  ! Need to do some interpolation to the new grid

  ALLOCATE (fieldx_um(1:nlongs_um))
  ALLOCATE (fieldz_um(1:nlevs_um+1))
  ALLOCATE (fieldxz_um(1:nlongs_um, 0:nlevs_um+1))


  ! Get the dimension data from the file
  !-------------------------------------
  ! Longitudinal distances
  !-----------------------
  startA(1) = 1
  countA(1) = nlongs_um

  ! These are the um values
  ierr = NF_GET_VARA_DOUBLE (ncid, varidLongs_u, startA, countA, fieldx_um(1:nlongs_um))
  IF (ierr .NE. 0) THEN
    PRINT*, '***Error getting dimension data for Longs_u***'
    PRINT*,'ierr', ierr, NF_STRERROR(ierr)
    STOP
  ENDIF
  ! Interpolate to ABC grid
  DO x = 1, nlongs
    fracx    = REAL((x-1)*(nlongs_um-1)) / REAL(nlongs-1) + 1
    lowerx   = INT(fracx)
    IF (lowerx == nlongs_um) lowerx = nlongs_um - 1
    upperx   = lowerx + 1
    axisx(1) = REAL(lowerx)
    axisx(2) = REAL(upperx)
    umdata % longs_u(x) = Interpolate1D(fieldx_um(lowerx:upperx), axisx(1:2), fracx)
  END DO

  ierr = NF_GET_VARA_DOUBLE (ncid, varidLongs_v, startA, countA, fieldx_um(1:nlongs_um))
  IF (ierr .NE. 0) THEN
    PRINT*, '***Error getting dimension data for Longs_v***'
    PRINT*,'ierr', ierr, NF_STRERROR(ierr)
    STOP
  ENDIF
  ! Interpolate to ABC grid
  DO x = 1, nlongs
    fracx    = REAL((x-1)*(nlongs_um-1)) / REAL(nlongs-1) + 1
    lowerx   = INT(fracx)
    IF (lowerx == nlongs_um) lowerx = nlongs_um - 1
    upperx   = lowerx + 1
    axisx(1) = REAL(lowerx)
    axisx(2) = REAL(upperx)
    umdata % longs_v(x) = Interpolate1D(fieldx_um(lowerx:upperx), axisx(1:2), fracx)
  END DO

  ! Level Heights
  !--------------
  startA(1) = 1
  countA(1) = nlevs_um + 1

  ierr = NF_GET_VARA_DOUBLE (ncid, varidfull_levs, startA, countA, fieldz_um(1:nlevs_um+1))
  IF (ierr .NE. 0) THEN
    PRINT*, '***Error getting dimension data for full_levs***'
    PRINT*,'ierr', ierr, NF_STRERROR(ierr)
    STOP
  ENDIF
  ! Interpolate to ABC grid
  DO z = 1, nlevs+1
    fracz    = REAL((z-1)*nlevs_um) / REAL(nlevs+1) + 1.0
    lowerz   = INT(fracz)
    IF (lowerz == nlevs_um) lowerz = nlevs_um - 1
    upperz   = lowerz + 1
    axisz(1) = REAL(lowerz)
    axisz(2) = REAL(upperz)
    umdata % full_levs(z-1) = Interpolate1D(fieldz_um(lowerz:upperz), axisz(1:2), fracz)
  END DO
  umdata % full_levs(nlevs+1) = 2.0 * umdata % full_levs(nlevs) - umdata % full_levs(nlevs-1)
  IF (nlevs <= nlevs_um) THEN
    ! Read the half levels from the um file
    ierr = NF_GET_VARA_DOUBLE (ncid, varidhalf_levs, startA, countA, fieldz_um(1:nlevs_um+1))
    IF (ierr .NE. 0) THEN
      PRINT*, '***Error getting dimension data for half_levs***'
      PRINT*,'ierr', ierr, NF_STRERROR(ierr)
      STOP
    ENDIF
    ! Interpolate to ABC grid
    DO z = 1, nlevs+1
      fracz    = REAL((z-1)*nlevs_um) / REAL(nlevs) + 1.0
      lowerz   = INT(fracz)
      IF (lowerz == nlevs_um) lowerz = nlevs_um - 1
      upperz   = lowerz + 1
      axisz(1) = REAL(lowerz)
      axisz(2) = REAL(upperz)
      umdata % half_levs(z) = Interpolate1D(fieldz_um(lowerz:upperz), axisz(1:2), fracz)
    END DO
  ELSE
    ! When nlevs > nlevs_um it is easier to select the half levels as half-way
    DO z = 1, nlevs
      umdata % half_levs(z) = (umdata % full_levs(z-1) + umdata % full_levs(z)) / 2.0
    END DO
    umdata % half_levs(nlevs+1) = 2.0 * umdata % half_levs(nlevs) - umdata % half_levs(nlevs-1)
  END IF


  ! Get the main data from the file
  !--------------------------------
  startB(1) = 1
  countB(1) = nlongs_um ! All longitude points
  startB(2) = latitude  ! Selected latitude slice
  countB(2) = 1         !lats
  startB(3) = 1
  countB(3) = 1         !levs
  startB(4) = 1
  countB(4) = 1         !time

  ! u
  !--
  DO z = 1, nlevs_um
    startB(3)                 = z
    ierr                      = NF_GET_VARA_DOUBLE (ncid, varidu, startB, countB, fieldx_um(1:nlongs_um))
    fieldxz_um(1:nlongs_um,z) = fieldx_um(1:nlongs_um)
  END DO
  IF (ierr .NE. 0) THEN
    PRINT*, '***Error getting main variable data u***'
    PRINT*,'ierr', ierr, NF_STRERROR(ierr)
    STOP
  ENDIF
  ! Interpolate to ABC grid
  DO x = 1, nlongs
    fracx    = REAL((x-1)*(nlongs_um-1)) / REAL(nlongs-1) + 1
    lowerx   = INT(fracx)
    IF (lowerx == nlongs_um) lowerx = nlongs_um - 1
    upperx   = lowerx + 1
    axisx(1) = REAL(lowerx)
    axisx(2) = REAL(upperx)
    DO z = 1, nlevs
      fracz    = REAL((z-1)*(nlevs_um-1)) / REAL(nlevs-1) + 1
      lowerz   = INT(fracz)
      IF (lowerz == nlevs_um) lowerz = nlevs_um - 1
      upperz   = lowerz + 1
      axisz(1) = REAL(lowerz)
      axisz(2) = REAL(upperz)
      umdata % u(x,z) = Interpolate2D (fieldxz_um(lowerx:upperx, lowerz:upperz), axisx(1:2), axisz(1:2), fracx, fracz)
    END DO
  END DO

  !OPEN (15, FILE='u.dat')
  !DO x = 1, nlongs
  !  WRITE (15,*) umdata % longs_u(x), umdata % u(x,30)
  !END DO
  !CLOSE (15)



  ! v
  !--
  DO z = 1, nlevs_um
    startB(3)                  = z
    ierr                       = NF_GET_VARA_DOUBLE (ncid, varidv, startB, countB, fieldx_um(1:nlongs_um))
    fieldxz_um(1:nlongs_um, z) = fieldx_um(1:nlongs_um)
  END DO
  IF (ierr .NE. 0) THEN
    PRINT*, '***Error getting main variable data v***'
    PRINT*,'ierr', ierr, NF_STRERROR(ierr)
    STOP
  ENDIF
  ! Interpolate to ABC grid
  DO x = 1, nlongs
    fracx    = REAL((x-1)*(nlongs_um-1)) / REAL(nlongs-1) + 1
    lowerx   = INT(fracx)
    IF (lowerx == nlongs_um) lowerx = nlongs_um - 1
    upperx   = lowerx + 1
    axisx(1) = REAL(lowerx)
    axisx(2) = REAL(upperx)
    DO z = 1, nlevs
      fracz    = REAL((z-1)*(nlevs_um-1)) / REAL(nlevs-1) + 1
      lowerz   = INT(fracz)
      IF (lowerz == nlevs_um) lowerz = nlevs_um - 1
      upperz   = lowerz + 1
      axisz(1) = REAL(lowerz)
      axisz(2) = REAL(upperz)
      umdata % v(x,z) = Interpolate2D (fieldxz_um(lowerx:upperx, lowerz:upperz), axisx(1:2), axisz(1:2), fracx, fracz)
    END DO
  END DO

  !OPEN (15, FILE='v.dat')
  !DO x = 1, nlongs
  !  WRITE (15,*) umdata % longs_v(x), umdata % v(x,30)
  !END DO
  !CLOSE (15)



  ! w
  !--
  DO z = 1, nlevs_um + 1
    startB(3)                  = z
    ierr                       = NF_GET_VARA_DOUBLE (ncid, varidw, startB, countB, fieldx_um(1:nlongs_um))
    fieldxz_um(1:nlongs_um, z) = fieldx_um(1:nlongs_um)
  ENDDO
  IF (ierr .NE. 0) THEN
    PRINT*, '***Error getting main variable data w***'
    PRINT*,'ierr', ierr, NF_STRERROR(ierr)
    STOP
  ENDIF
  ! Interpolate to ABC grid
  DO x = 1, nlongs
    fracx    = REAL((x-1)*(nlongs_um-1)) / REAL(nlongs-1) + 1
    lowerx   = INT(fracx)
    IF (lowerx == nlongs_um) lowerx = nlongs_um - 1
    upperx   = lowerx + 1
    axisx(1) = REAL(lowerx)
    axisx(2) = REAL(upperx)
    DO z = 1, nlevs+1
      fracz    = REAL((z-1)*(nlevs_um-1)) / REAL(nlevs) + 1
      lowerz   = INT(fracz)
      IF (lowerz == nlevs_um) lowerz = nlevs_um - 1
      upperz   = lowerz + 1
      axisz(1) = REAL(lowerz)
      axisz(2) = REAL(upperz)
      umdata % w(x,z-1) = Interpolate2D (fieldxz_um(lowerx:upperx, lowerz:upperz), axisx(1:2), axisz(1:2), fracx, fracz)
    END DO
  END DO

  ! density
  !--------
  DO z = 1, nlevs_um
    startB(3)                  = z
    ierr                       = NF_GET_VARA_DOUBLE (ncid, variddensity, startB, countB, fieldx_um(1:nlongs_um))
    fieldxz_um(1:nlongs_um, z) = fieldx_um(1:nlongs_um)
  ENDDO
  IF (ierr .NE. 0) THEN
    PRINT*, '***Error getting main variable data density***'
    PRINT*,'ierr', ierr, NF_STRERROR(ierr)
    STOP
  ENDIF
  ! Interpolate to ABC grid
  DO x = 1, nlongs
    fracx    = REAL((x-1)*(nlongs_um-1)) / REAL(nlongs-1) + 1
    lowerx   = INT(fracx)
    IF (lowerx == nlongs_um) lowerx = nlongs_um - 1
    upperx   = lowerx + 1
    axisx(1) = REAL(lowerx)
    axisx(2) = REAL(upperx)
    DO z = 1, nlevs
      fracz    = REAL((z-1)*(nlevs_um-1)) / REAL(nlevs-1) + 1
      lowerz   = INT(fracz)
      IF (lowerz == nlevs_um) lowerz = nlevs_um - 1
      upperz   = lowerz + 1
      axisz(1) = REAL(lowerz)
      axisz(2) = REAL(upperz)
      umdata % density(x,z) = Interpolate2D (fieldxz_um(lowerx:upperx, lowerz:upperz), axisx(1:2), axisz(1:2), fracx, fracz)
    END DO
  END DO

  ! theta
  !------
  DO z=1, nlevs_um
    startB(3)                  = z
    ierr                       = NF_GET_VARA_DOUBLE (ncid, varidtheta, startB, countB, fieldx_um(1:nlongs_um))
    fieldxz_um(1:nlongs_um, z) = fieldx_um(1:nlongs_um)
  ENDDO
  IF (ierr .NE. 0) THEN
    PRINT*, '***Error getting main variable data theta***'
    PRINT*,'ierr', ierr, NF_STRERROR(ierr)
    STOP
  ENDIF
  ! Interpolate to ABC grid
  DO x = 1, nlongs
    fracx    = REAL((x-1)*(nlongs_um-1)) / REAL(nlongs-1) + 1
    lowerx   = INT(fracx)
    IF (lowerx == nlongs_um) lowerx = nlongs_um - 1
    upperx   = lowerx + 1
    axisx(1) = REAL(lowerx)
    axisx(2) = REAL(upperx)
    DO z = 1, nlevs
      fracz    = REAL((z-1)*(nlevs_um-1)) / REAL(nlevs-1) + 1
      lowerz   = INT(fracz)
      IF (lowerz == nlevs_um) lowerz = nlevs_um - 1
      upperz   = lowerz + 1
      axisz(1) = REAL(lowerz)
      axisz(2) = REAL(upperz)
      umdata % theta(x,z) = Interpolate2D (fieldxz_um(lowerx:upperx, lowerz:upperz), axisx(1:2), axisz(1:2), fracx, fracz)
    END DO
  END DO

  ! exner presure
  !----------------
  DO z = 1, nlevs_um + 1
    startB(3)                  = z
    ierr                       = NF_GET_VARA_DOUBLE (ncid, varidexpres, startB, countB, fieldx_um(1:nlongs_um))
    fieldxz_um(1:nlongs_um, z) = fieldx_um(1:nlongs_um)
  ENDDO
  IF (ierr .NE. 0) THEN
    PRINT*, '***Error getting main variable data u***'
    PRINT*,'ierr', ierr, NF_STRERROR(ierr)
    STOP
  ENDIF
  ! Interpolate to ABC grid
  DO x = 1, nlongs
    fracx    = REAL((x-1)*(nlongs_um-1)) / REAL(nlongs-1) + 1
    lowerx   = INT(fracx)
    IF (lowerx == nlongs_um) lowerx = nlongs_um - 1
    upperx   = lowerx + 1
    axisx(1) = REAL(lowerx)
    axisx(2) = REAL(upperx)
    DO z = 1, nlevs+1
      fracz    = REAL((z-1)*(nlevs_um-1)) / REAL(nlevs) + 1
      lowerz   = INT(fracz)
      IF (lowerz == nlevs_um) lowerz = nlevs_um - 1
      upperz   = lowerz + 1
      axisz(1) = REAL(lowerz)
      axisz(2) = REAL(upperz)
      umdata % exner_pressure(x,z) = Interpolate2D (fieldxz_um(lowerx:upperx, lowerz:upperz), axisx(1:2), axisz(1:2), fracx, fracz)
    END DO
  END DO

  ! orographic height
  !--------------------
  DO z = 1, 1
    startB(3) = z
    ierr      = NF_GET_VARA_DOUBLE (ncid, varidorog, startB, countB, fieldx_um(1:nlongs_um))
  ENDDO
  IF (ierr .NE. 0) THEN
    PRINT*, '***Error getting data for orog_height***'
    PRINT*,'ierr', ierr, NF_STRERROR(ierr)
    STOP
  ENDIF
  ! Interpolate to ABC grid
  DO x = 1, nlongs
    fracx    = REAL((x-1)*(nlongs_um-1)) / REAL(nlongs-1) + 1
    lowerx   = INT(fracx)
    IF (lowerx == nlongs_um) lowerx = nlongs_um - 1
    upperx   = lowerx + 1
    axisx(1) = REAL(lowerx)
    axisx(2) = REAL(upperx)
    umdata % orog_height(x) = Interpolate1D(fieldx_um(lowerx:upperx), axisx(1:2), fracx)
  END DO

  DEALLOCATE (fieldx_um, fieldz_um, fieldxz_um)


END IF

!Close the netCDF file
ierr = NF_CLOSE(ncid)
IF ( ierr .NE. 0 ) THEN
  PRINT*, '***ERROR closing file***'
  STOP
ENDIF

END SUBROUTINE Read_um_data_2d
