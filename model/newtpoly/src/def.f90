!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PARAMETER files for the a newtonian polytrope star       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE PARAMETER
USE DEFINITION
IMPLICIT NONE
SAVE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Polytope !
REAL*8, PARAMETER :: npoly = 1.0D0
REAL*8, PARAMETER :: kpoly = 100.0D0
REAL*8, PARAMETER :: rhocen = 1.28D-3

REAL*8, PARAMETER :: ggas = 1.0D0 + 1.0D0/npoly
REAL*8, PARAMETER :: alpha = SQRT((npoly + 1.0D0)*kpoly/4.0D0/pi)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Gravity !
REAL*8, ALLOCATABLE :: dphidx(:,:,:)

! Floors !
REAL*8, PARAMETER :: rho_floor = 1.0D-20
REAL*8, PARAMETER :: eps_floor = kpoly*rho_floor**(ggas - 1.0)/(ggas - 1.0)

END MODULE
