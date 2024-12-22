!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PARAMETER files for the rotating magnetised white dwarf  !
! The cold fermi gas is used as the equation of state      !
! Normal Newtonian gravity is assumed                      !
! Electron fraction is assumed to be 0.5 always everywhere !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE PARAMETER
USE DEFINITION
IMPLICIT NONE
SAVE
! Unit constants !

! Physical constants to be as one !
REAL*8, PARAMETER :: gconst = 6.67430D-8
REAL*8, PARAMETER :: clight = 2.99792458D10
REAL*8, PARAMETER :: solar = 1.98847D33

! Here, mu_0 is not in cgs !
REAL*8, PARAMETER :: mu_0 = 1.25663706212D-6 ! unit in kg m s-2 A-2 or H/m !

! Solar Radius !
REAL*8, PARAMETER :: rsolar = 6.96342D10

! 1 GeV !
REAL*8, PARAMETER :: GeV2gram = 1.78266191D-24

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Unit conversions !
! Conversion between units !
REAL*8, PARAMETER :: lencgs2code = (clight**2)/(solar*gconst)
REAL*8, PARAMETER :: masscgs2code = (1.0D0/solar)
REAL*8, PARAMETER :: tcgs2code = (clight**3)/(solar*gconst)

! Derived conversion !
REAL*8, PARAMETER :: rhocgs2code = (masscgs2code/lencgs2code**3)
REAL*8, PARAMETER :: energycgs2code = (1.0D0/clight**2)
REAL*8, PARAMETER :: taucgs2code = (rhocgs2code*energycgs2code)        ! (masscgs2code*lencgs2code**2/tcgs2code**2) !
REAL*8, PARAMETER :: h_bar = (1.054571817D-27)*(lencgs2code**2*masscgs2code/tcgs2code)

! Current conversion !
REAL*8, PARAMETER :: amp2code = (mu_0*1.0D5*masscgs2code*lencgs2code)**(0.5D0)/tcgs2code

! Magnetic field !
REAL*8, PARAMETER :: gauss2code = 1.0D-1*masscgs2code/amp2code/tcgs2code**2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8, PARAMETER :: me2 = 9.1093837015D-28*masscgs2code
REAL*8, PARAMETER :: mb2 = 1.66053906660D-24*masscgs2code
REAL*8, PARAMETER :: ye = 0.5D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Parameters !
REAL*8, PARAMETER :: rhomax = 5.0D10*rhocgs2code
REAL*8, PARAMETER :: atmosphere = rhomax*1.0D-7
REAL*8, PARAMETER :: atmospheric = 1.0D-7

! Constant for fermi equation of state  !
! Note that the speed of light is unity !
REAL*8, PARAMETER :: amax = (me2**4)/(2.4D1*pi**2*h_bar**3)
REAL*8, PARAMETER :: bmax = (mb2*me2**3)/(3.0D0*pi**2*h_bar**3*ye)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for solving gravity !

! Solve the potential per how many steps
INTEGER, PARAMETER :: n_pot = 20

! maximum number of relaxation !
INTEGER, PARAMETER :: relax_max = 1000000

! Tolerance in relaxation of the potential			
REAL*8, PARAMETER :: tolerance = 1.0D-8

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Floors !
REAL*8 :: rho_floor
REAL*8 :: eps_floor

! gravitational potential !
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: phi, phi_old

! for poisson solver !
REAL*8, ALLOCATABLE, DIMENSION (:) :: ajp1, ajm1
REAL*8, ALLOCATABLE, DIMENSION (:,:) :: bkp1, bkm1
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: clp1, clm1, epsc

! for poisson solver !
REAL*8, PARAMETER :: omega_weight = 1.9d0


END MODULE
