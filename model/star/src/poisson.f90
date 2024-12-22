!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Compute the poisson equation coefficient for the relaxation method 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_poisson
USE PARAMETER 
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

REAL*8 :: xl, yl, zl
REAL*8 :: xp, yp, zp
REAL*8 :: xm, ym, zm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get poisson coefficient !
Do i = 1, nx
	Do j = 1, ny
		Do k = 1, nz
			CALL GET_COORD(i-1,j-1,k-1,xm,ym,zm)
			CALL GET_COORD(i  ,j  ,k  ,xl,yl,zl)
			CALL GET_COORD(i+1,j+1,k+1,xp,yp,zp)

			CALL poisson_coef(xm, xl, xp, &
							ym, yl, yp, &
							zm, zl, zp, &
							ajp1(i), ajm1(i), bkp1(i,j), bkm1(i,j), &
							clp1(i,j,k), clm1(i,j,k), epsc(i,j,k))
		end do
	end do
END DO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Get left hand side of the discrete poisson equation 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE poisson_coef(xm1, xc, xp1, ym1, yc, yp1, zm1, zc, zp1, &
  alphajp1, alphajm1, betakp1, betakm1, gammalp1, gammalm1, epsc)
implicit none

! input !
real*8, INTENT(IN) :: xm1, xc, xp1, ym1, yc, yp1, zm1, zc, zp1

! output !
real*8, INTENT(OUT) :: epsc
real*8, INTENT(OUT) :: alphajp1, alphajm1
real*8, INTENT(OUT) :: betakp1, betakm1
real*8, INTENT(OUT) :: gammalp1, gammalm1

! assign !
epsc = 2.0d0*(xp1+xm1-3.0d0*xc)/xc/(xp1-xc)/(xc-xm1) &
+ (yp1+ym1-2.0d0*yc-2.0d0*DTAN(yc))/xc**2/DTAN(yc)/(yp1-yc)/(yc-ym1) &
- 2.0d0/xc**2/DSIN(yc)**2/(zp1-zc)/(zc-zm1)

! assign !
alphajp1 = 2.0d0*(2.0d0*xc-xm1)/(xp1-xm1)/(xp1-xc)/xc
alphajm1 = 2.0d0*(2.0d0*xc-xp1)/(xp1-xm1)/(xc-xm1)/xc

! assign !
betakp1 = (2.0d0*DTAN(yc) + yc - ym1)/xc**2/DTAN(yc)/(yp1 - yc)/(yp1 - ym1)
betakm1 = (2.0d0*DTAN(yc) + yc - yp1)/xc**2/DTAN(yc)/(yp1 - ym1)/(yc - ym1)

! assign !
gammalp1 = 2.0d0/xc**2/DSIN(yc)**2/(zp1-zm1)/(zp1-zc)
gammalm1 = 2.0d0/xc**2/DSIN(yc)**2/(zp1-zm1)/(zc-zm1) 

END SUBROUTINE