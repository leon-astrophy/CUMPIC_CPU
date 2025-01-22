!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom equations !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_EQN
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom arrays !
!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_HYDRO
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Populate custom arrays to GPU !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_POPULATE
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_GRID
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

!---------------------------------------------------------!

! Dummy variables
INTEGER :: i, j, k, l

! For radial grid !
REAL*8 :: DX_T, DY_T, DZ_T

!---------------------------------------------------------!
! Get interface coordiante !

! Assign grid size !
DX_T = DLOG(Rout/Rin)/DBLE(nxtot)
DY_T = (PI - 2.0d0*off_set)/DBLE(nytot)
DZ_T = 2.0d0*PI/DBLE(nztot)

! Log grid in the radial direction !
DO j = - NGHOST, nx + NGHOST
	xF(j) = DLOG(Rin) + (starts(1) + DBLE(j))*DX_T
  xF(j) = DEXP(xF(j))
END DO

! Compression along the theta, NOTE: might need boundary condition !
DO k = - NGHOST, ny + NGHOST
	yF(k) = off_set + (starts(2) + DBLE(k))*DY_T
  yF(k) = yF(k) + 0.5D0*(1.0D0 - s_comp)*DSIN(2.0d0*yF(k))
END DO

! Uniform along phi direction !
DO l = - NGHOST, nz + NGHOST
	zF(l) = 0.0D0 + (starts(3) + DBLE(l))*DZ_T
END DO

!---------------------------------------------------------!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_BOUNDARY
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC)
DO l = -2, nz + 3
  DO k = -2, ny + 3
    DO j = 1, 3
      prim(ivx,1-j,k,l) = MIN(prim(ivx,1-j,k,l), 0.0D0)
      prim(ivx,nx+j,k,l) = MAX(prim(ivx,nx+j,k,l), 0.0D0)
    END DO
  END DO               
ENDDO 
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!
! Custom variable floor !
!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_CHECKRHO
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k, l
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC)
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx

      ! Check density and internal energy density !
      IF(prim(irho,j,k,l) < rho_floor) THEN
        prim(irho,j,k,l) = rho_floor
      END IF
      IF(eps(j,k,l) < eps_floor) THEN
        eps(j,k,l) = eps_floor
      ENDIF

    END DO
  END DO
END DO
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom equations !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_SOURCE
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k, l

! Threshold for atmosphere density
REAL*8 :: dphidr
REAL*8 :: factor, diff

! Local !
REAL*8 :: x_loc, y_loc, z_loc
REAL*8 :: dx_loc, dy_loc, dz_loc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Add black hole gravity !
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE (dphidr, factor, diff, x_loc, y_loc, z_loc, dx_loc, dy_loc, dz_loc)
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      CALL GET_COORD(j,k,l,x_loc,y_loc,z_loc)
      CALL COORD_DX(j,k,l,dx_loc,dy_loc,dz_loc)
      diff = prim(irho,j,k,l) - rho_floor
      factor = MAX(SIGN(1.0D0, diff), 0.0D0)
      dphidr = 1.0d0/((x_loc - r_sh)*(x_loc - r_sh))
      sc(ivx,j,k,l) = sc(ivx,j,k,l) + (-factor*prim(irho,j,k,l)*dphidr)
      sc(itau,j,k,l) = sc(itau,j,k,l) + (-factor*prim(irho,j,k,l)*prim(ivx,j,k,l)*dphidr)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!
! Do custom updates !
!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_UPDATE (p_in)
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

! Integer !
INTEGER, INTENT (IN) :: p_in

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do custom operator split !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OPERATOR_SPLIT
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Print custom analysis file !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OPENFILE_CUSTOM
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

END SUBROUTINE
