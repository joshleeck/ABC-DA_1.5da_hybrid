SUBROUTINE Energy (state, dims)

! To calculate components of energy

USE DefConsTypes, ONLY :  &
    ABC_type,             &
    dims_type,            &
    ZREAL8,               &
    nlongs, nlevs,        &
    half,                 &
    rho0,                 &
    A, B, C, dx

IMPLICIT NONE

! Declare subroutine parameters
TYPE(ABC_type),  INTENT(INOUT)   :: state
TYPE(dims_type), INTENT(IN)      :: dims

! Declare local variables
INTEGER                          :: x, z
REAL(ZREAL8)                     :: u, v, w, rho, bp, r, vol
REAL(ZREAL8)                     :: E_k, E_b, E_e

! Declare functions
REAL(ZREAL8)                  :: INT_FH

! Initialise
E_k = 0.0
E_e = 0.0
E_b = 0.0

! Calculate energy
DO z = 1, nlevs
  vol  = dx * (dims % full_levs(z) - dims % full_levs(z-1))
  DO x = 1, nlongs
    ! Interpolate to rho points on the grid
    r     = state % r(x,z)
    u     = (state % u(x-1,z) + state % u(x,z)) * half
    v     = state % v(x,z)
    w     = INT_FH (state % w(x,z-1), state % w(x,z), z, dims)
    bp    = INT_FH (state % b(x,z-1), state % b(x,z), z, dims)
    rho   = state % rho(x,z)

    E_k  = E_k + vol * rho * (u*u + v*v + w*w) * half * rho0
    E_e  = E_e + vol * C * r*r * half * rho0 / B
    E_b  = E_b + vol * rho * bp*bp * half * rho0 / (A*A)
  END DO
END DO

state % Kinetic_Energy = E_k
state % Buoyant_Energy = E_b
state % Elastic_Energy = E_e
state % Total_Energy   = E_k + E_e + E_b

END SUBROUTINE Energy
