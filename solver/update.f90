!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! This subroutine prepare the data necessary for constructing
! the flux for the spatial discretization.
! It takes input/output of the U array and the 
! mode p which decides whether or not to update
! the gravitational potential
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE UPDATE (p_in)
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER, INTENT (IN) :: p_in

! Integer !
INTEGER :: i, j, k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Updates to hydrodynamic variables !

! Find pressure !
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC)
DO l = 1 - NGHOST, nz + NGHOST
  DO k = 1 - NGHOST, ny + NGHOST
    DO j = 1 - NGHOST, nx + NGHOST
      CALL EOS_PRESSURE (j, k, l, prim(irho,j,k,l), eps(j,k,l), prim(itau,j,k,l))
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Custom updates !

! Do custom updates !
CALL CUSTOM_UPDATE (p_in)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE