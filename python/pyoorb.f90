MODULE pyoorb

  USE Base_cl
  USE Orbit_cl
  USE Time_cl
  IMPLICIT NONE

CONTAINS

  SUBROUTINE propagation(elements, mjd_tt0, mjd_tt1)

    USE Base_cl
    USE Time_cl
    USE Orbit_cl
    IMPLICIT NONE
    REAL(bp), DIMENSION(6), INTENT(inout) :: elements
    REAL(bp), INTENT(in) :: mjd_tt0, mjd_tt1

    TYPE (Orbit) :: orb
    TYPE (Time) :: t0, t1

    CALL NEW(t0, mjd_tt0, "TT")
    CALL NEW(t1, mjd_tt1, "TT")
    CALL NEW(orb, elements, "cartesian", "ecliptic", t0)
    CALL propagate(orb, t1)
    elements = getElements(orb, "cartesian", "ecliptic")
    CALL NULLIFY(t0)
    CALL NULLIFY(t1)
    CALL NULLIFY(orb)

  END SUBROUTINE propagation

END MODULE pyoorb
