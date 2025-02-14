!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine prepares the data for spatial discretization,
! due ask the WENO_module to do the reconstruction
! and then combines the results for one Runge-Kutta sub-step
! Prototype developed by Wong Ka Wing in 2010 (or before?)
! Extended by Leung Shing Chi in 2016 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SPATIAL
USE OMP_LIB
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i,j,k,l

!$OMP PARALLEL
!---------------------------------------------------------------------------------------------!
! First, initialize !

! Initialize source term and rungekutta operator !
!$OMP DO COLLAPSE(4) SCHEDULE(STATIC)
DO l = 0, nz
	DO k = 0, ny
		DO j = 0, nx
			DO i = imin, imax
				sc(i,j,k,l) = 0.0D0
				l_rk(i,j,k,l) = 0.0D0
			END DO
		END DO
	END DO
END DO
!$OMP END DO

! Initialize electric field !
!$OMP DO COLLAPSE(4) SCHEDULE(STATIC)
DO l = 0, nz
	DO k = 0, ny
		DO j = 0, nx
			DO i = iex, iez
				ecorn(i,j,k,l) = 0.0D0
			END DO
		END DO
	END DO
END DO
!$OMP END DO

!---------------------------------------------------------------------------------------------!
!$OMP END PARALLEL

! Predefined source term !
CALL GET_SOURCE

! Custom source term !
CALL CUSTOM_SOURCE

!---------------------------------------------------------------------------------------------!
! Now compute the interface flux !

! Interface flux, sweep through all direction !
CALL GET_FLUXES

! Constrained transport !
CALL flux_ct

!---------------------------------------------------------------------------------------------!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Get geometric source terms
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GET_SOURCE
USE OMP_LIB
USE DEFINITION 
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k, l

! Geometric sources terms
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC)
!-----------------------------------------------------------------------------------!
DO l = 1, nz
	DO k = 1, ny
		DO j = 1, nx
				
			! Call warpper for computing geometric sources !
			CALL GEOM_SOURCE(j, k, l)

		END DO
	END DO
END DO
!-----------------------------------------------------------------------------------!
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! 1. Peform reconstruction to get primitive variables at cell interfaces 
! 2. Construct conservative variables and fluxes at cell interfaces 
! 3. Solve the Riemann problem
! 4. Solve for the flux difference
! 5. Update the rungekutta operator 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GET_FLUXES
USE OMP_LIB
USE DEFINITION 
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k, l

! Geometric factor !
REAL*8 :: geom_flux_p, geom_flux_c, geom_flux_m

!$OMP PARALLEL PRIVATE(geom_flux_p, geom_flux_c, geom_flux_m)
!==============================================================================================================!

! First loop through the x-direction
!--------------------------------------------------------------------------------------------------------------!
! Interpolate to get L/R state !
!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
DO l = 0, nz + 1 
	DO k = 0, ny + 1
		DO j = 0, nx + 1

			! Core hydrodynamic variables !
			DO i = imin, ibx - 1
				CALL INTERPOLATE (x_dir, j, k, l, prim(i,j-2,k,l), prim(i,j-1,k,l), prim(i,j,k,l), prim(i,j+1,k,l), prim(i,j+2,k,l), primR(i,j-1,k,l), primL(i,j,k,l))
			END DO

			! Extra scalar !

			! Cell center magnetic fields !
			CALL INTERPOLATE (x_dir, j, k, l, bcell(iby,j-2,k,l), bcell(iby,j-1,k,l), bcell(iby,j,k,l), bcell(iby,j+1,k,l), bcell(iby,j+2,k,l), primR(iby,j-1,k,l), primL(iby,j,k,l))
			CALL INTERPOLATE (x_dir, j, k, l, bcell(ibz,j-2,k,l), bcell(ibz,j-1,k,l), bcell(ibz,j,k,l), bcell(ibz,j+1,k,l), bcell(ibz,j+2,k,l), primR(ibz,j-1,k,l), primL(ibz,j,k,l))
				
			! Special treatment for normal field 
			primR(ibx,j-1,k,l) = prim(ibx,j-1,k,l)
			primL(ibx,j,k,l) = prim(ibx,j,k,l)

		END DO
	END DO
END DO
!$OMP END DO

! Then solve the Riemann Problem !
!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
DO l = 0, nz + 1 
	DO k = 0, ny + 1
		DO j = 0, nx

			! Core hydrodynamic variables !
			CALL RIEMANN (x_dir, j, k, l) 

		END DO
	END DO
END DO
!$OMP END DO

! Add the flux difference into the l-operator
!$OMP DO COLLAPSE(4) SCHEDULE(STATIC)
DO l = 1, nz
	DO k = 1, ny
		DO j = 1, nx
			DO i = imin, ibx - 1

				! Call the geometric factor !
				CALL GEOM_FLUX(x_dir,j,k,l,geom_flux_p,geom_flux_c,geom_flux_m)

				! Perform subtraction !
				l_rk(i,j,k,l) = l_rk(i,j,k,l) - (flux (i,j,k,l)*geom_flux_p - flux (i,j-1,k,l)*geom_flux_m) / geom_flux_c

			END DO
		END DO
	END DO
END DO
!$OMP END DO

! Now do the same for the electric field, using flux-CT scheme 
!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
DO l = 0, nz
  DO k = 0, ny
    DO j = 0, nx

      ! Add emf !
      ecorn(iey,j,k,l) = ecorn(iey,j,k,l) + 0.25d0*(flux(ibz,j,k,l) + flux(ibz,j,k,l+1))

      ecorn(iez,j,k,l) = ecorn(iez,j,k,l) + 0.25d0*(- flux(iby,j,k,l) - flux(iby,j,k+1,l))

    END DO
  END DO
END DO
!$OMP END DO

!==============================================================================================================!

! Then loop through the y-direction
!--------------------------------------------------------------------------------------------------------------!
! Interpolate to get L/R state !
!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
DO l = 0, nz + 1 
	DO k = 0, ny + 1
		DO j = 0, nx + 1

			! Core hydrodynamic variables !
			DO i = imin, ibx - 1
				CALL INTERPOLATE (y_dir, j, k, l, prim(i,j,k-2,l), prim(i,j,k-1,l), prim(i,j,k,l), prim(i,j,k+1,l), prim(i,j,k+2,l), primR(i,j,k-1,l), primL(i,j,k,l))
			END DO

			! Extra scalar !

			! Cell center magnetic fields !
			CALL INTERPOLATE (y_dir, j, k, l, bcell(ibx,j,k-2,l), bcell(ibx,j,k-1,l), bcell(ibx,j,k,l), bcell(ibx,j,k+1,l), bcell(ibx,j,k+2,l), primR(ibx,j,k-1,l), primL(ibx,j,k,l))
			CALL INTERPOLATE (y_dir, j, k, l, bcell(ibz,j,k-2,l), bcell(ibz,j,k-1,l), bcell(ibz,j,k,l), bcell(ibz,j,k+1,l), bcell(ibz,j,k+2,l), primR(ibz,j,k-1,l), primL(ibz,j,k,l))

			! Special treatment for normal field 
			primR(iby,j,k-1,l) = prim(iby,j,k-1,l)
			primL(iby,j,k,l) = prim(iby,j,k,l)

		END DO
	END DO
END DO
!$OMP END DO

! Then solve the Riemann Problem !
!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
DO l = 0, nz + 1 
	DO k = 0, ny 
		DO j = 0, nx + 1

			! Core hydrodynamic variables !
			CALL RIEMANN (y_dir, j, k, l)

		END DO
	END DO
END DO
!$OMP END DO

! Add the flux difference into the l-operator
!$OMP DO COLLAPSE(4) SCHEDULE(STATIC)
DO l = 1, nz
	DO k = 1, ny
		DO j = 1, nx
			DO i = imin, ibx - 1

				! Call the geometric factor !
				CALL GEOM_FLUX(y_dir,j,k,l,geom_flux_p,geom_flux_c,geom_flux_m)

				! Perform subtraction !
				l_rk(i,j,k,l) = l_rk(i,j,k,l) - (flux (i,j,k,l)*geom_flux_p - flux (i,j,k-1,l)*geom_flux_m) / geom_flux_c

			END DO
		END DO
	END DO
END DO
!$OMP END DO

! Now do the same for the electric field, using flux-CT scheme 
!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
DO l = 0, nz
  DO k = 0, ny
    DO j = 0, nx

      ! Add emf !
      ecorn(iex,j,k,l) = ecorn(iex,j,k,l) + 0.25d0*(- flux(ibz,j,k,l) - flux(ibz,j,k,l+1))

      ecorn(iez,j,k,l) = ecorn(iez,j,k,l) + 0.25d0*(flux(ibx,j,k,l) + flux(ibx,j+1,k,l))
			
    END DO
  END DO
END DO
!$OMP END DO

!==============================================================================================================!

! Finally loop through the z-direction
!--------------------------------------------------------------------------------------------------------------!
! Interpolate to get L/R state !
!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
DO l = 0, nz + 1 
	DO k = 0, ny + 1
		DO j = 0, nx + 1
			
			! Core hydrodynamic variables !
			DO i = imin, ibx - 1
				CALL INTERPOLATE (z_dir, j, k, l, prim(i,j,k,l-2), prim(i,j,k,l-1), prim(i,j,k,l), prim(i,j,k,l+1), prim(i,j,k,l+2), primR(i,j,k,l-1), primL(i,j,k,l))
			END DO

			! Extra scalar !

			! Cell center magnetic fields !
			CALL INTERPOLATE (z_dir, j, k, l, bcell(ibx,j,k,l-2), bcell(ibx,j,k,l-1), bcell(ibx,j,k,l), bcell(ibx,j,k,l+1), bcell(ibx,j,k,l+2), primR(ibx,j,k,l-1), primL(ibx,j,k,l))
			CALL INTERPOLATE (z_dir, j, k, l, bcell(iby,j,k,l-2), bcell(iby,j,k,l-1), bcell(iby,j,k,l), bcell(iby,j,k,l+1), bcell(iby,j,k,l+2), primR(iby,j,k,l-1), primL(iby,j,k,l))

			! Special treatment for normal field 
			primR(ibz,j,k,l-1) = prim(ibz,j,k,l-1)
			primL(ibz,j,k,l) = prim(ibz,j,k,l)

		END DO
	END DO
END DO
!$OMP END DO

! Then solve the Riemann Problem !
!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
DO l = 0, nz 
	DO k = 0, ny + 1
		DO j = 0, nx + 1

			! Core hydrodynamic variables !
			CALL RIEMANN (z_dir, j, k, l)

		END DO
	END DO
END DO
!$OMP END DO

! Add the flux difference into the l-operator
!$OMP DO COLLAPSE(4) SCHEDULE(STATIC)
DO l = 1, nz
	DO k = 1, ny
		DO j = 1, nx
			DO i = imin, ibx - 1

				! Call the geometric factor !
				CALL GEOM_FLUX(z_dir,j,k,l,geom_flux_p,geom_flux_c,geom_flux_m)

				! Perform subtraction !
				l_rk(i,j,k,l) = l_rk(i,j,k,l) - (flux (i,j,k,l)*geom_flux_p - flux (i,j,k,l-1)*geom_flux_m) / geom_flux_c

			END DO
		END DO
	END DO
END DO
!$OMP END DO

! Now do the same for the electric field, using flux-CT scheme 
!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
DO l = 0, nz
  DO k = 0, ny
    DO j = 0, nx

      ! Add emf !
      ecorn(iex,j,k,l) = ecorn(iex,j,k,l) + 0.25d0*(flux(iby,j,k,l) + flux(iby,j,k+1,l))

      ecorn(iey,j,k,l) = ecorn(iey,j,k,l) + 0.25d0*(- flux(ibx,j,k,l) - flux(ibx,j+1,k,l))

    END DO
  END DO
END DO
!$OMP END DO

!--------------------------------------------------------------------------------------------------------------!

! At the end, add source terms !
!--------------------------------------------------------------------------------------------------------------!
! Final step, get rungekutta operator, LHS of the hydro equation !
!$OMP DO COLLAPSE(4) SCHEDULE(STATIC)
DO l = 1, nz
	DO k = 1, ny
		DO j = 1, nx
			DO i = imin, ibx - 1
				l_rk(i,j,k,l) = l_rk(i,j,k,l) + sc(i,j,k,l)
			END DO
		END DO
	END DO
END DO
!$OMP END DO
!==============================================================================================================!
!$OMP END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE