!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Definition files for specific simulation models
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PARAMETER
IMPLICIT NONE
SAVE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Grid parameters !

! Inner radial boundary !
REAL*8, PARAMETER :: Rin = 3.0d0

! Outer radial boundary !
REAL*8, PARAMETER :: Rout = 23.0d0

! Off-set along the pole to avoid singularity !
REAL*8, PARAMETER :: off_set = 1.0e-6

! Compression along the theta direction !
REAL*8, PARAMETER :: s_comp = 0.5D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Adiabatic index !
REAL*8, PARAMETER :: ggas = 5.0d0/3.0d0

! Normalization constant !
REAL*8, PARAMETER :: beta_norm = 100.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Problem setup !

! schwarzschild radius !
REAL*8, PARAMETER :: r_sh = 2.0

! angular velocity gradient !
REAL*8, PARAMETER :: q_grad = 2.0d0

! inner (equatorial) radius of the torus !
REAL*8, PARAMETER :: s_in = 6.0d0

! (equatorial) radius where the density is at maximum !
REAL*8, PARAMETER :: s_max = 9.4d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For the magnetic field setup 

! minimum density (fraction) to define the last contour of the vector potential !
REAL*8, PARAMETER :: rho_cut = 0.2d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! atmospheric density !
REAL*8, PARAMETER :: rho_fac = 1.0d-4

! maximum density !
REAL*8, PARAMETER :: rho_max = 1.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Atmospheric value !
REAL*8 :: rho_floor
REAL*8 :: eps_floor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE PARAMETER
