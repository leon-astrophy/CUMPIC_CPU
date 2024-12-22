!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Builiding initial model                                                           !
! Simulation for rotating magnetised white dwarf AIC                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GET_MODEL
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k, l

! Coordinate !
REAL*8 :: xl, yl, zl
REAL*8 :: dxl, dyl, dzl

! Magnetic field !
REAL*8 :: maxdb
REAL*8 :: div_b
REAL*8 :: axp,axm,ayp,aym,azp,azm

REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: a_phi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Preperation !
! Allocate
Allocate(a_phi(-3:nx+3,-3:ny+3,-3:nz+3))

! Poisson interpolation coefficient !
call get_poisson

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
prim = 0.0d0
bcell = 0.0d0
eps = 0.0d0
a_phi = 0.0d0
! Read and assign density !
! OPEN(UNIT=970, FILE = './profile/hydro_rho.dat', ACTION='READ')
! READ(970,*) ((prim(irho,j,k,1), j = 1, nx), k = 1, ny)
! CLOSE(970)

OPEN(UNIT=970, FILE = './hydro_rho.dat', ACTION='READ')
READ(970,*) (prim(irho,j,1,1), j = 1, nx)
CLOSE(970)

DO k = 1, ny
  prim(irho,:,k,1) = prim(irho,:,1,1)
END DO
prim(irho,:,:,:) = prim(irho,:,:,:)*rhocgs2code
PRINT *, "Finished reading rho"

! Assign velocity !
prim(ivx,:,:,1) = 0.0d0
prim(ivy,:,:,1) = 0.0d0

! Read and assign velocity_z !
! OPEN(UNIT=970, FILE = './profile/hydro_vphi.dat', ACTION='READ')
! READ(970,*) ((prim(ivz,j,k,1), j = 1, nx), k = 1, ny)
! CLOSE(970)

! prim(ivz,:,:,:) = prim(ivz,:,:,:)*lencgs2code/tcgs2code
! PRINT *, "Finished reading vphi"

! Read for magnetic vector potential !
! OPEN(UNIT=970, FILE = './profile/hydro_Aphi.dat', ACTION='READ')
! READ(970,*) ((a_phi(j,k,1), j = 0, nx), k = 0, ny)
! CLOSE(970)

! PRINT *, "Finished reading Aphi"

! Calculate magnetic field !
! Coordinate here are in code unit but aphi is in gaussian unit !
! Unit conversion below !

! Get face-magnetic fields by cross product !
DO l = 0, nz
  DO k = 0, ny
    DO j = 0, nx
      CALL GET_COORD(j,k,l,xl,yl,zl)
      CALL COORD_DX(j,k,l,dxl,dyl,dzl)
      prim(ibx,j,k,l) = (DSIN(yF(k))*a_phi(j,k,l) - SIN(yF(k-1))*a_phi(j,k-1,l))/(xF(j)*(DCOS(yF(k-1)) - DCOS(yF(k))))
      prim(iby,j,k,l) = - (xF(j)*a_phi(j,k,l) - xF(j-1)*a_phi(j-1,k,l))/(xl*dxl)
    END DO
  END DO
END DO

! Get cell-centered magnetic fields by averaging !
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      bcell(ibx,j,k,l) = 0.5D0*(prim(ibx,j,k,l) + prim(ibx,j-1,k,l))
      bcell(iby,j,k,l) = 0.5D0*(prim(iby,j,k,l) + prim(iby,j,k-1,l))
    END DO
  END DO
END DO

prim(ibx:iby,:,:,:) = prim(ibx:iby,:,:,:)*gauss2code*lencgs2code  ! length conversion for curl !

PRINT *, "Finished calculating poloidal field"

! OPEN(UNIT=970, FILE = './profile/hydro_bphi.dat', ACTION='READ')
! READ(970,*) ((prim(ibz,j,k,1), j = 1, nx), k = 1, ny)
! CLOSE(970)
! prim(ibz,:,:,:) = prim(ibz,:,:,:)*gauss2code
! bcell(ibz,:,:,:) = prim(ibz,:,:,:)

PRINT *, "Finished reading torodial field "

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Deallocate
Deallocate(a_phi)

prim(ivx:ivz,:,:,:) = 0.0d0
prim(ibx:ibz,:,:,:) = 0.0d0
bcell(ibx:ibz,:,:,:) = 0.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      CALL EOS_PRESSURE(j,k,l,prim(irho,j,k,l),eps(j,k,l),prim(itau,j,k,l))
      CALL EOS_EPSILON(prim(irho,j,k,l),prim(itau,j,k,l),eps(j,k,l))
    END DO
  END DO
END DO
PRINT *, "Finished calculating pressure and eps"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Check divergence-B = 0 constraint !
maxdb = 0.0d0       
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      CALL GEOM_AREA(j,k,l,axp,axm,ayp,aym,azp,azm)
      div_b = (axp*prim(ibx,j,k,l) - axm*prim(ibx,j-1,k,l)) &
            + (ayp*prim(iby,j,k,l) - aym*prim(iby,j,k-1,l)) &
            + (azp*prim(ibz,j,k,l) - azm*prim(ibz,j,k,l-1))
      maxdb = MAX(maxdb, div_b)
    END DO
  END DO
END DO

! Communicate across all processor 
#ifdef MPI
CALL MPI_Allreduce(maxdb, maxdb, 1, MPI_DOUBLE, MPI_MAX, new_comm, ierror)
#endif

WRITE (*,*)
WRITE (*,*) 'Maximum initial divergence B', maxdb
WRITE (*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! set atmospheric primitive variables !
rho_floor = atmosphere
eps_floor = eps(nx,1,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign floor density !
DO j = 1, nx
  DO k = 1, ny
    IF(prim(irho,j,k,1) < rho_floor) THEN
      prim(irho,j,k,:)    = rho_floor
      eps(j,k,:)          = eps_floor
      prim(ivx:ivz,j,k,:) = 0.0d0
    END IF
  END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CALL BOUNDARY

END SUBROUTINE
