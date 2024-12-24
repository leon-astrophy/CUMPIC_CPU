!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom equations !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_EQN
USE DEFINITION
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom arrays !
!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_HYDRO
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

! gravitational potential energy !
ALLOCATE (phi(-2:nx+3,-2:ny+3,-2:nz+3))
ALLOCATE (phi_old(-2:nx+3,-2:ny+3,-2:nz+3))

! for poisson equation !
ALLOCATE (ajp1(1:nx))
ALLOCATE (ajm1(1:nx))
ALLOCATE (bkp1(1:nx,1:ny))
ALLOCATE (bkm1(1:nx,1:ny))
ALLOCATE (clp1(1:nx,1:ny,1:nz))
ALLOCATE (clm1(1:nx,1:ny,1:nz))
ALLOCATE (epsc(1:nx,1:ny,1:nz))

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Populate custom arrays to GPU !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_POPULATE
USE PARAMETER
USE DEFINITION
IMPLICIT NONE

! Now populate all necessary, and reuseable arrays to the graphic cards !
!$ACC enter DATA COPYIN(phi, phi_old, ajp1, ajm1, bkp1, bkm1, clp1, clm1, epsc):q

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Clear custom arrays from GPU !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_CLEAR
USE PARAMETER
USE DEFINITION
IMPLICIT NONE

! Now we clear memory in the GPU device !
!$ACC exit DATA DELETE(phi, phi_old, ajp1, ajm1, bkp1, bkm1, clp1, clm1, epsc)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_GRID
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

! Integer !
INTEGER :: i, j, nlines

! Read the number of lines in the file !
nlines = 0 
OPEN (970, file = './hydro_x1_fgrid.dat') 
DO 
  READ (970,*, END=10) 
  nlines = nlines + 1 
END DO 
10 CLOSE (970) 

! Error message !
IF(nlines .ne. nx+7) THEN
  WRITE (*,*) 'number of grid faces from files', nlines
  WRITE (*,*) 'number of x grid faces in the program', nx+6
  STOP 'inconsistent number of grid faces, exit'
END IF

! Read !
OPEN(UNIT=970, FILE = './hydro_x1_fgrid.dat', ACTION='READ')
DO i = -3, nx+3
	READ(970,*) xF(i)
ENDDO
CLOSE(970)

xF = xF*lencgs2code

END SUBROUTINE
 
!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_BOUNDARY
USE DEFINITION
IMPLICIT NONE

INTEGER :: i, j, k, l

! x boundary !
DO l = - 2, nz + 3
  DO k = -2, ny + 3
    DO j = 1, 3
      prim(ivx,1-j,k,l) = 0.0D0
      prim(ivx,nx+j,k,l) = MAX(prim(ivx,nx+j,k,l), 0.0D0)
    END DO
  END DO
ENDDO

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

DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      prim(irho,j,k,l) = MAX(rho_floor, prim(irho,j,k,l))
      eps(j,k,l)       = MAX(eps_floor, eps(j,k,l))
    END DO
  END DO
END DO


END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom equations !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_SOURCE
USE PARAMETER
USE DEFINITION
IMPLICIT NONE

 ! Integer !
INTEGER :: i, j, k, l

! Derivatives of gravitational potential
REAL*8 :: dphidx, dphidy, dphidz

! Threshold for atmosphere density
REAL*8 :: factor, diff

! Local, plus 1, minus 1 !
REAL*8 :: xl, yl, zl
REAL*8 :: xp, yp, zp
REAL*8 :: xm, ym, zm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      CALL GET_COORD(j-1,k-1,l-1,xm,ym,zm)
      CALL GET_COORD(j  ,k  ,l  ,xl,yl,zl)
      CALL GET_COORD(j+1,k+1,l+1,xp,yp,zp)

      ! Include only non-atmosphere !
			diff = prim(irho,j,k,l) - rho_floor
      factor = MAX(SIGN(1.0D0, diff), 0.0D0)

      ! Gravitational potential of the matter !
      dphidx = first_derivative (xm, xl, xp, phi(j-1,k,l), phi(j,k,l), phi(j+1,k,l))
      dphidy = first_derivative (ym, yl, yp, phi(j,k-1,l), phi(j,k,l), phi(j,k+1,l))

      ! Add them to the source term !
      sc(ivx,j,k,l) = sc(ivx,j,k,l) - factor*prim(irho,j,k,l)*dphidx
      sc(ivy,j,k,l) = sc(ivy,j,k,l) - factor*prim(irho,j,k,l)*dphidy/xl
      ! sc(itau,j,k,l) = sc(itau,j,k,l) - factor*prim(irho,j,k,l)* &
                        ! (prim(ivx,j,k,l)*dphidx + prim(ivy,j,k,l)*dphidy/xl)
      sc(itau,j,k,l) = sc(itau,j,k,l) - factor*prim(irho,j,k,l)* &
                        (prim(ivx,j,k,l)*dphidx + prim(ivy,j,k,l)*dphidy/xl)
    END DO
  END DO
END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains

	REAL*8 function first_derivative (xm1, xc, xp1, fm1, fc, fp1)
	!$acc routine seq
	implicit none
	REAL*8 :: xm1, xc, xp1, fm1, fc, fp1, h1, h2
  h2 = xp1 - xc
  h1 = xc - xm1
	first_derivative = ((fp1-fc)*h1*h1+(fc-fm1)*h2*h2)/(h1*h2*(h1+h2))
	end function

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

! Integer !
INTEGER :: j, k, l, n

! For poisson solver !
REAL*8 :: abserror, rhs

! Density threshold !
REAL*8 :: rho_in, factor, diff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update gravitational potentials !

DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      CALL EOS_EPSILON(prim(irho,j,k,l),prim(itau,j,k,l),eps(j,k,l))
    END DO
  END DO
END DO



IF (p_in == 0 .OR. MOD(n_step, n_pot) == 0) THEN

  ! special treatment for initial model !
  IF(p_in == 0) THEN
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! First, give a guessing potential !
    !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)
    DO l = 0, nz + 1
      DO k = 0, ny + 1
        DO j = 0, nx + 1
          phi(j,k,l) = 0.0d0
        END DO
      END DO
    END DO
    !$OMP END PARALLEL DO
    
  END IF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Calucaltes potential by RBSOR
  DO n = 1, relax_max
    ! Back up potential !
    !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)
    DO l = 1, nz
      DO k = 1, ny
        DO j = 1, nx
          phi_old(j,k,l) = phi(j,k,l)
        END DO
      END DO
    END DO
    !$ACC END PARALLEL

    ! Set error !
    !$ACC SERIAL
	  abserror = 1.0D-50
	  !$ACC END SERIAL
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Red chess !
    !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE(diff,rhs,rho_in,factor)
    DO l = 1, nz
      DO k = 1, ny
        DO j = 1, nx
          IF ((-1)**(j+k+l)>0) THEN	
            diff = prim(irho,j,k,l) - rho_floor
            factor = MERGE(1.0d0, 0.0d0, diff > 0.0d0)
            rho_in = factor*prim(irho,j,k,l)
            rhs = (4.0d0*pi*rho_in - & 
								(ajp1(j)*phi(j+1,k,l) + ajm1(j)*phi(j-1,k,l) + & 
								 bkp1(j,k)*phi(j,k+1,l) + bkm1(j,k)*phi(j,k-1,l) + & 
								 clp1(j,k,l)*phi(j,k,l+1) + clm1(j,k,l)*phi(j,k,l-1)))/epsc(j,k,l)
					  phi(j,k,l) = (1.0d0 - omega_weight)*phi(j,k,l) + omega_weight*rhs
          ELSE 
            CYCLE
          END IF
        END DO
      END DO
    END DO
    !$ACC END PARALLEL

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Black chess !
    !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE(diff,rhs,rho_in,factor)
    DO l = 1, nz
      DO k = 1, ny
        DO j = 1, nx
          IF ((-1)**(j+k+l)<0) THEN
            diff = prim(irho,j,k,l) - rho_floor
            factor = MERGE(1.0d0, 0.0d0, diff > 0.0d0)
            rho_in = factor*prim(irho,j,k,l)
            rhs = (4.0d0*pi*rho_in - & 
								(ajp1(j)*phi(j+1,k,l) + ajm1(j)*phi(j-1,k,l) + & 
								 bkp1(j,k)*phi(j,k+1,l) + bkm1(j,k)*phi(j,k-1,l) + & 
								 clp1(j,k,l)*phi(j,k,l+1) + clm1(j,k,l)*phi(j,k,l-1)))/epsc(j,k,l)
					  phi(j,k,l) = (1.0d0 - omega_weight)*phi(j,k,l) + omega_weight*rhs
          ELSE 
            CYCLE
          END IF
        END DO
      END DO
    END DO
    !$ACC END PARALLEL

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Look for maximum abserror !
    !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) Reduction(MAX:abserror)
    DO l = 1, nz
      DO k = 1, ny
        DO j = 1, nx
          IF (phi_old(j,k,l) /= 0.0d0) THEN
            abserror = max(abserror, abs((phi(j,k,l) - phi_old(j,k,l)) / phi_old(j,k,l)))
          ELSE
            abserror = max(abserror, abs(phi(j,k,l)))
          END IF
        END DO
      END DO
    END DO
    !$ACC END PARALLEL

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Boundary conditions !
    !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(2) DEFAULT(PRESENT)
    DO l = 1, nz
      DO k = 1, ny
        phi(0,k,l) = phi(1,k,l)
        phi(nx+1,k,l) = 0.0d0
      END DO
    END DO
    !$ACC END PARALLEL
    !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(2) DEFAULT(PRESENT)
    DO l = 1, nz
      DO j = 1, nx
        phi(j,0,l) = phi(j,1,l)
        phi(j,ny+1,l) = phi(j,ny,l)
      END DO
    END DO
    !$ACC END PARALLEL
    !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(2) DEFAULT(PRESENT)
    DO k = 1, ny
      DO j = 1, nx
        phi(j,k,0) = phi(j,k,1)
        phi(j,k,nz+1) = phi(j,k,nz)
      END DO
    END DO
    !$ACC END PARALLEL
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Debug and exit !
	  !WRITE (*,*) n, abserror
    IF(abserror <= tolerance) EXIT 
  END DO

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Stop condition !
  IF(n == relax_max) THEN
    WRITE (*,*) n, relax_max
    STOP 'Convergence error in poisson solver'
  END IF
ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom equations !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OPERATOR_SPLIT
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OPENFILE_CUSTOM
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_ANALYSIS
IMPLICIT NONE

END SUBROUTINE
