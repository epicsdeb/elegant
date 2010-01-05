	subroutine brighte(i,alpha,cosphi,sinphi,
     1	           s0,s1,s2,s3,axr,ayr,axi,ayi)

C+
C 
C FUNCTIONAL DESCRIPTION:	
C 
C  Driver routine for calculation of the brightness function. This routine selects
C  the appropriate subroutine for either a plane device or an elliptical device.
C 
C AUTHORS: 
C 
C  Roger J. Dejus
C 
C  The Advanced Photon Source
C  Experimental Facilities Division
C  Argonne National Laboratory
C 
C CREATION DATE: 
C 
C   05-SEP-1994
C 
C FORMAL PARAMETERS:
C  
C  Input arguments:
C  i			Harmonic number; +1,+2,+3, ...
C  alpha		Gamma*Theta; Theta is the polar angle [rad]
C  cosphi		Cosine of azimuthal angle phi
C  sinphi		Sine   of azimuthal angle phi
C   
C  Output arguments:
C  s0                   Stokes parameter (s0)
C  s1                   Stokes parameter (s1)
C  s2                   Stokes parameter (s2)
C  s3                   Stokes parameter (s3)
C  axr   		Real  component of electric field in the x-direction (hor.)
C  ayr   		Real  component of electric field in the y-direction (vert.)
C  axi   		Imag. component of electric field in the x-direction (hor.)
C  ayi   		Imag. component of electric field in the y-direction (vert.)
C  
C COMMON BLOCKS:
C  
C  /calcr/              Real valued calculation parameters
C  
C DESIGN ISSUES:
C  
C  To be added.
C  
C MODIFICATION HISTORY:
C 
C	 Date     | Name  | Description
C ----------------+-------+-----------------------------------------------------
C [change_entry]
C-
C  Declarations of scalars:
	integer*4	i
	real*8		alpha,cosphi,sinphi
	real*8		s0,s1,s2,s3,axr,ayr,axi,ayi
  
C  Declarations of arrays:
  
C  Labeled constants:
	real*8		KMIN,ZERO
	parameter	(KMIN=1.0D-3,ZERO=0.0D0)

C  Common blocks:
	real*8		k,kx,ky,k3,npi,pe,gamma,er,ew,dew,len
	common		/calcr/	    k,kx,ky,k3,npi,pe,gamma,er,ew,dew,len

C+ Select call to appropriate Brightness routine
C  -----------------------------------------------------------------------------
	if (kx .lt. KMIN .and. ky .gt. KMIN) then      ! Regular plane device
	   call bright1(i,ky,alpha,cosphi,sinphi,
     1               s0,s1,s2,s3,axr,ayr)
	   axi = ZERO
	   ayi = ZERO
	else if (kx .gt. KMIN .and. ky .lt. KMIN) then ! Flipped plane device
	   call bright1(i,kx,alpha,sinphi,-cosphi,
     1               s0,s1,s2,s3,axr,ayr)
	   s1 = -s1
	   s2 = -s2
	   axi = ZERO
	   ayi = ZERO
	else                                           ! Elliptical device
	   call bright3(i,kx,ky,alpha,cosphi,sinphi,
     1               s0,s1,s2,s3,axr,ayr,axi,ayi)
	end if

	return
	end ! brighte

	subroutine bright1(i,k,alpha,cosphi,sinphi,
     1	           s0,s1,s2,s3,axr,ayr)
C+
C 
C FUNCTIONAL DESCRIPTION:	
C 
C  Routine to calculate the brightness using an expansion in Bessel
C  functions. This routine is for a plane ideal device.
C 
C AUTHORS: 
C 
C  Roger J. Dejus
C 
C  The Advanced Photon Source
C  Experimental Facilities Division
C  Argonne National Laboratory
C 
C CREATION DATE: 
C 
C   05-SEP-1994
C 
C FORMAL PARAMETERS:
C  
C  Input arguments:
C  i			Harmonic number; +1,+2,+3, ...
C  k		        Deflection parameter [dim. less]
C  alpha		Gamma*Theta; Theta is the polar angle [rad]
C  cosphi		Cosine of azimuthal angle phi
C  sinphi		Sine   of azimuthal angle phi
C   
C  Output arguments:
C  s0                   Stokes parameter (s0)
C  s1                   Stokes parameter (s1)
C  s2                   Stokes parameter (s2)
C  s3                   Stokes parameter (s3)
C  axr   		Real  component of electric field in the x-direction (hor.)
C  ayr   		Real  component of electric field in the y-direction (vert.)
C  
C COMMON BLOCKS:
C  
C  None
C  
C DESIGN ISSUES:
C  
C  The formalism is treated in "Handbook on Synchrotron Radiation, Vol. 1,
C  edited by E.E. Koch, North-Holland Publ. Co. (1983), ch. 2, "Characteristics
C  of Synchrotron Radiation and of its Sources", S. Krinsky, M.L. Perlman and
C  R.E. Watson. In particular section 6.4 and equations 269 - 272.
C  
C MODIFICATION HISTORY:
C 
C	 Date     | Name  | Description
C ----------------+-------+-----------------------------------------------------
C [change_entry]
C-
C  Size parameters:
	integer*4	B_SZ
	parameter	(B_SZ=1000)
  
C  Declarations of scalars:
	logical*4	lodd
	integer*4	i,n,nx,ny,np,nm,im
	real*8		k,alpha,cosphi,sinphi
	real*8		s0,s1,s2,s3,axr,ayr
	real*8		bxp,bxm,bx1,bx2,sum1,sum2,sum3,sign
	real*8		x,y,absx
	real*8		c,a,alpha2,k2
  
C  Declarations of arrays:
	real*8		bx(0:B_SZ),by(0:B_SZ)
  
C  Labeled constants:
	real*8		ZERO,ONE,TWO,FOUR
	real*8		HALF,QT
	parameter	(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,FOUR=4.0D0)
	parameter	(HALF=0.5D0,QT=0.25D0)
	real*8		A5,A6,A8,A10,B5,B6,B8,B10
	parameter	(A5=5.22D0,B5=1.40D0,A6 =6.20D0,B6 =1.41D0)
	parameter	(A8=8.03D0,B8=1.44D0,A10=9.74D0,B10=1.47D0)
	real*8		EPS,EPSX
	parameter	(EPS=1.0D-6) ! Must correspond to Ai and Bi below
	parameter	(EPSX=1.0D-5)! Provides approx. EPSX**2 sig. figures
				     ! Cannot make smaller than ~ 10**-8
C  Common blocks:

C+ Definitions
C  -----------------------------------------------------------------------------
	alpha2 = alpha*alpha
	k2     = k*k
	a      = ONE+HALF*k2+alpha2
	c      = -k*alpha*cosphi    ! c = (-a/i)*(x/2)
	x      = (TWO*i/a)*(-c)     ! x may be < 0.0, 0.0, > 0.0
	y      = (QT *i/a)*(k2)     ! y > 0.0
	absx   = abs(x)
C-
C  -----------------------------------------------------------------------------

C+ Calculate Bessel functions
C  -----------------------------------------------------------------------------
	nx = A6+B6*absx +ONE
	ny = A6+B6*y    +ONE
	if (nx .gt. B_SZ) then
	    nx = B_SZ
	    print 200,'&bright1-W-DIMERR, Dimension error'
            print 200,'Array of Bessel function values out ',
     1        'of bound for x'
	    print 210,'- upper bound reset to ',B_SZ
	end if ! nx
	if (ny .gt. B_SZ) then
	    ny = B_SZ
	    print 200,'&bright1-W-DIMERR, Dimension error'
	    print 200,'Array of Bessel function values out ',
     1        'of bound for y'
	    print 210,'- upper bound reset to ',B_SZ
	end if ! ny
200	format(' ',a)
210	format(' ',a,i4)

	if (absx .ge. EPSX) call bright_bessjn(x,nx,bx) ! J0(x) ... JNx(x)
	call bright_bessjn(y,ny,by)		        ! J0(y) ... JNy(y)

C ev. use abs value of by
	do while (by(ny) .lt. EPS)  ! Modify ny
	    ny = ny-1
	end do !
C-
C  -----------------------------------------------------------------------------

C+ Calculate sums
C  -----------------------------------------------------------------------------
	if (i/2*2 .ne. i) then ! Harmonic odd
	    lodd = .true.	
	else ! Harmonic even
	    lodd = .false.
	end if ! i

	if (absx .lt. EPSX) then ! Small arguments

	    if (lodd) then ! Harmonic odd
		sum1 = ZERO
		sum3 = ZERO
		np   = -(-i+1)/2
		nm   = -(-i-1)/2
		if (np .le. ny) then
		    if (np/2*2 .ne. np) then ! Index odd; change sign
			sum1 = -by(np)
 			sum3 = -by(np)
		    else
			sum1 = +by(np)
			sum3 = +by(np)
		    end if ! np
		end if ! np
		if (nm .le. ny) then
		    if (nm/2*2 .ne. nm) then ! Index odd; change sign
			sum1 = sum1+by(nm)
			sum3 = sum3-by(nm)
		    else
			sum1 = sum1-by(nm)
			sum3 = sum3+by(nm)
		    end if ! nm
		end if ! nm
		sum1 = HALF*x*sum1
	    else ! Harmonic even
		sum1 = ZERO
		nm   = +i/2
		if (nm .le. ny) then
		    if (nm/2*2 .ne. nm) then ! Index odd; change sign
			sum1 = -by(nm)
		    else
			sum1 = +by(nm)
		    end if ! nm
		end if ! nm

		sum3 = ZERO
		np   = -(-i+2)/2
		nm   = -(-i-2)/2
		if (np .le. ny) then
		    if (np/2*2 .ne. np) then ! Index odd; change sign
			sum3 = -by(np)
		    else
			sum3 = +by(np)
		    end if ! np
		end if ! np
		if (nm .le. ny) then
		    if (nm/2*2 .ne. nm) then ! Index odd; change sign
			sum3 = sum3+by(nm)
		    else
			sum3 = sum3-by(nm)
		    end if ! nm
		end if ! nm
		sum3 = HALF*x*sum3
	    end if ! lodd

	else ! abs(x) > EPSX; Large arguments

	    sum1 = ZERO
	    sum2 = ZERO
	    sign = ONE
	    if (i .le. nx) sum1 = by(0)*bx(i)
	    do n=1,ny
		sign = -sign
		np   =  2*n+i
		nm   = -2*n+i
		im   = abs(nm)

		if (np .le. nx) then
		    bxp = bx(np)
		else
		    bxp = ZERO
		end if ! np

		if (im .le. nx) then
		    if (lodd .and. nm .lt. 0) then
			bxm = -bx(im)
		    else
			bxm = +bx(im)
		    end if ! lodd
		else
		    bxm = ZERO
		end if ! im

		bx1  = bxp +sign*bxm
		bx2  = bxp -sign*bxm
		sum1 = sum1 +  by(n)*bx1
		sum2 = sum2 +n*by(n)*bx2
	    end do ! n
			
	    sum3 = TWO/x*(i*sum1 +TWO*sum2)

	end if ! abs(x)
C-
C  -----------------------------------------------------------------------------
	axr= (TWO*i/a)*(alpha*cosphi*sum1 -HALF*k*sum3)
	ayr= (TWO*i/a)*(alpha*sinphi*sum1)
c	ft = FOUR*i*i/a/a*(alpha2*sum1*sum1 +c*sum1*sum3 +QT*k2*sum3*sum3) ! old variable; same as s0
	s0 = axr*axr +ayr*ayr
	s1 = axr*axr -ayr*ayr
	s2 = TWO*axr*ayr
	s3 = ZERO

	return
	end ! bright1

	subroutine bright3(i,kx,ky,alpha,cosphi,sinphi,
     1	           s0,s1,s2,s3,axr,ayr,axi,ayi)
C+
C 
C FUNCTIONAL DESCRIPTION:	
C 
C  Routine to calculate the brightness using an expansion in Bessel
C  functions. This routine is for a elliptical ideal device.
C 
C AUTHORS: 
C 
C  Roger J. Dejus
C 
C  The Advanced Photon Source
C  Experimental Facilities Division
C  Argonne National Laboratory
C 
C CREATION DATE: 
C 
C   05-SEP-1994
C 
C FORMAL PARAMETERS:
C  
C  Input arguments:
C  i			Harmonic number; +1,+2,+3, ...
C  kx		        Deflection parameter in the x-direction (hor.)  [dim. less]
C  ky		        Deflection parameter in the y-direction (vert.) [dim. less]
C  alpha		Gamma*Theta; Theta is the polar angle [rad]
C  cosphi		Cosine of azimuthal angle phi
C  sinphi		Sine   of azimuthal angle phi
C   
C  Output arguments:
C  s0                   Stokes parameter (s0)
C  s1                   Stokes parameter (s1)
C  s2                   Stokes parameter (s2)
C  s3                   Stokes parameter (s3)
C  axr   		Real  component of electric field in the x-direction (hor.)
C  ayr   		Real  component of electric field in the y-direction (vert.)
C  axi   		Imag. component of electric field in the x-direction (hor.)
C  ayi   		Imag. component of electric field in the y-direction (vert.)
C  
C COMMON BLOCKS:
C  
C  None
C  
C DESIGN ISSUES:
C  
C  Proper expansion to first order in the variable x has been retained for small
C  arguments.  The expansion in the variable y for small arguments is only correct
C  to zero order.  This is OK because when y is small, we are interested in the 
C  case for the value of identical zero (Kx=Ky).
C  
C MODIFICATION HISTORY:
C 
C	 Date     | Name  | Description
C ----------------+-------+-----------------------------------------------------
C [change_entry]
C-
C  Size parameters:
	integer*4	B_SZ
	parameter	(B_SZ=1000)
  
C  Declarations of scalars:
	logical*4	lodd
	integer*4	i,n,nx,ny
	integer*4	np,np0,np1p,np1m
	integer*4       nm,nm0,nm1p,nm1m,im0,im1p,im1m
	real*8		kx,ky,alpha,cosphi,sinphi
	real*8		s0,s1,s2,s3,axr,ayr,axi,ayi
	real*8		sign
	real*8		x,y,absy,phi
	real*8		a,alpha2,k2
	real*8		co,so,cp0,cm0,sp0,sm0,cp1p,cm1p,sp1p,sm1p,
     1               cp1m,cm1m,sp1m,sm1m,cps,cp1s,cms,cm1s
	real*8          sum0r,sum0i,sum1pr,sum1pi,sum1mr,sum1mi
	real*8          bxp0r,bxp0i,bxp1pr,bxp1pi,bxp1mr,bxp1mi
	real*8          bxm0r,bxm0i,bxm1pr,bxm1pi,bxm1mr,bxm1mi
	real*8          bx10r,bx10i,bx11pr,bx11pi,bx11mr,bx11mi
	real*8          sum1p,sum1m
  
C  Declarations of arrays:
	real*8		bx(0:B_SZ),by(0:B_SZ)
  
C  Labeled constants:
	real*8		ZERO,ONE,TWO
	real*8		HALF,QT
	parameter	(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
	parameter	(HALF=0.5D0,QT=0.25D0)
	real*8		A5,A6,A8,A10,B5,B6,B8,B10
	parameter	(A5=5.22D0,B5=1.40D0,A6 =6.20D0,B6 =1.41D0)
	parameter	(A8=8.03D0,B8=1.44D0,A10=9.74D0,B10=1.47D0)
	real*8		EPS,EPSX,EPSY
	parameter	(EPS=1.0D-6) ! Must correspond to Ai and Bi below
	parameter	(EPSX=1.0D-5)! Provides approx. EPSX**2 sig. figures
				     ! Cannot make smaller than ~ 10**-8
	parameter	(EPSY=EPSX)  ! Make same as EPSX

C  Common blocks:

C+ Definitions
C  -----------------------------------------------------------------------------
	alpha2 = alpha*alpha
	k2     = kx*kx +ky*ky
	a      = ONE+HALF*k2+alpha2
	x      = (TWO*i/a)*alpha*
     1        sqrt(kx*kx*sinphi*sinphi +ky*ky*cosphi*cosphi) ! x >= 0.0
	y      = (QT *i/a)*(ky*ky -kx*kx)
	absy   = abs(y)
	phi    = atan2((kx*sinphi),(ky*cosphi))
	co     = cos(TWO*phi)
	so     = sin(TWO*phi)
	cp0    = cos(i*phi)
	cm0    = cp0
	sp0    = sin(i*phi)
	sm0    = sp0
	cp1p   = cos((i+1)*phi)
	cm1p   = cp1p
	sp1p   = sin((i+1)*phi)
	sm1p   = sp1p
	cp1m   = cos((i-1)*phi)
	cm1m   = cp1m
	sp1m   = sin((i-1)*phi)
	sm1m   = sp1m
	sum0r  = ZERO
	sum0i  = ZERO
	sum1pr = ZERO
	sum1pi = ZERO
	sum1mr = ZERO
	sum1mi = ZERO
C-
C  -----------------------------------------------------------------------------

C+ Calculate Bessel functions
C  -----------------------------------------------------------------------------
	nx = A6+B6*x    +ONE
	ny = A6+B6*absy +ONE
	if (nx .gt. B_SZ) then
	    nx = B_SZ
	    print 200,'&bright3-W-DIMERR, Dimension error'
	    print 200,'Array of Bessel function values out ',
     1         'of bound for x'
	    print 210,'- upper bound reset to ',B_SZ
	end if ! nx
	if (ny .gt. B_SZ) then
	    ny = B_SZ
	    print 200,'&bright3-W-DIMERR, Dimension error'
	    print 200,'Array of Bessel function values out ',
     1         'of bound for y'
	    print 210,'- upper bound reset to ',B_SZ
	end if ! ny
200	format(' ',a)
210	format(' ',a,i4)

	if (x .ge. EPSX) call bright_bessjn(x,nx,bx) ! J0(x) ... JNx(x)
	if (absy .ge. EPSY) then
	    call bright_bessjn(y,ny,by)              ! J0(y) ... JNy(y)
	    do while (by(ny) .lt. EPS) ! Modify ny
	        ny = ny-1
	    end do
	end if ! absy
C-
C  -----------------------------------------------------------------------------

C+ Calculate sums
C  -----------------------------------------------------------------------------
	if (i/2*2 .ne. i) then ! Harmonic odd
	    lodd = .true.	
	else ! Harmonic even
	    lodd = .false.
	end if ! i

	if (x .lt. EPSX .and. absy .lt. EPSY) then ! Two small arguments
	    if (i .eq. 1) then	    ! First harmonic
	        sum0r  =  cos(phi)*HALF*x
		sum0i  = -sin(phi)*HALF*x
	        sum1mr = ONE
	    else if (i .eq. 2) then ! Second harmonic
	        sum1mr =  cos(phi)*HALF*x
		sum1mi = -sin(phi)*HALF*x
	    end if ! i

	else if (x .lt. EPSX) then ! Small arguments for x
	    if (lodd) then ! Harmonic odd
		np   = -(-i+1)/2
		nm   = -(-i-1)/2
		if (np .le. ny) then
		    if (np/2*2 .ne. np) then ! Index odd; change sign
			sum0r  = -by(np)
			sum0i  = -by(np)
 			sum1mr = -by(np)
		    else
			sum0r  = +by(np)
			sum0i  = +by(np)
			sum1mr = +by(np)
		    end if ! np
		end if ! np
		if (nm .le. ny) then
		    if (nm/2*2 .ne. nm) then ! Index odd; change sign
			sum0r  =  sum0r +by(nm)
			sum0i  =  sum0i -by(nm)
			sum1pr = -by(nm)
		    else
			sum0r  =  sum0r -by(nm)
			sum0i  =  sum0i +by(nm)
			sum1pr = +by(nm)
		    end if ! nm
		end if ! nm
		sum0r =  cos(phi)*HALF*x*sum0r
		sum0i = -sin(phi)*HALF*x*sum0i
	    else ! Harmonic even
		nm   = +i/2
		if (nm .le. ny) then
		    if (nm/2*2 .ne. nm) then ! Index odd; change sign
			sum0r = -by(nm)
		    else
			sum0r = +by(nm)
		    end if ! nm
		end if ! nm

		np   = -(-i+2)/2
		nm   = -(-i-2)/2
		if (np .le. ny) then
		    if (np/2*2 .ne. np) then ! Index odd; change sign
			sum1m = -by(np)
		    else
			sum1m = +by(np)
		    end if ! np
		end if ! np
		if (nm .le. ny) then
		    if (nm/2*2 .ne. nm) then ! Index odd; change sign
			sum1p = -by(nm)
		    else
			sum1p = +by(nm)
		    end if ! nm
		end if ! nm
		
	        sum1pr =  cos(phi)*HALF*x*(sum0r -sum1p)
	        sum1pi = -sin(phi)*HALF*x*(sum0r +sum1p)
	        sum1mr =  cos(phi)*HALF*x*(sum1m -sum0r)
	        sum1mi = -sin(phi)*HALF*x*(sum1m +sum0r)
	    end if ! lodd

	else if (absy .lt. EPSY) then ! Small arguments for y
	    if (i .le. nx) then
	        sum0r  =  cp0 *bx(i)
		sum0i  = -sp0 *bx(i)
	    end if
	    if ((i+1) .le. nx) then
	        sum1pr =  cp1p*bx(i+1)
	        sum1pi = -sp1p*bx(i+1)
	    end if
	    if ((i-1) .le. nx) then
	        sum1mr =  cp1m*bx(i-1)
	        sum1mi = -sp1m*bx(i-1)
	    end if

	else ! x > EPSX and abs(y) > EPSY; Large arguments
	    if (i .le. nx) then
	        sum0r  =  cp0 *by(0)*bx(i)
		sum0i  = -sp0 *by(0)*bx(i)
	    end if
	    if ((i+1) .le. nx) then
	        sum1pr =  cp1p*by(0)*bx(i+1)
	        sum1pi = -sp1p*by(0)*bx(i+1)
	    end if
	    if ((i-1) .le. nx) then
	        sum1mr =  cp1m*by(0)*bx(i-1)
	        sum1mi = -sp1m*by(0)*bx(i-1)
	    end if

	    sign = ONE
	    do n=1,ny
		sign = -sign
		np0  =  2*n+i
		np1p =  2*n+i+1
		np1m =  2*n+i-1
		nm0  = -2*n+i
		nm1p = -2*n+i+1
		nm1m = -2*n+i-1
		im0  = abs(nm0)
		im1p = abs(nm1p)
		im1m = abs(nm1m)

c  Recurrence relations
		cps  = cp0
		cms  = cm0
		cp0  = cp0 *co -sp0 *so
		cm0  = cm0 *co +sm0 *so
		sp0  = sp0 *co +cps *so
		sm0  = sm0 *co -cms *so

		cp1s = cp1p
		cm1s = cm1p
		cp1p = cp1p*co -sp1p*so
		cm1p = cm1p*co +sm1p*so
		sp1p = sp1p*co +cp1s*so
		sm1p = sm1p*co -cm1s*so

		cp1s = cp1m
		cm1s = cm1m
		cp1m = cp1m*co -sp1m*so
		cm1m = cm1m*co +sm1m*so
		sp1m = sp1m*co +cp1s*so
		sm1m = sm1m*co -cm1s*so

		if (np0 .le. nx) then
		    bxp0r  = cp0 *bx(np0)
		    bxp0i  = sp0 *bx(np0)
		else
		    bxp0r  = ZERO
		    bxp0i  = ZERO
		end if ! np0
		if (np1p .le. nx) then
		    bxp1pr = cp1p*bx(np1p)
		    bxp1pi = sp1p*bx(np1p)
		else
		    bxp1pr = ZERO
		    bxp1pi = ZERO
		end if ! np1p
		if (np1m .le. nx) then
		    bxp1mr = cp1m*bx(np1m)
		    bxp1mi = sp1m*bx(np1m)
		else
		    bxp1mr = ZERO
		    bxp1mi = ZERO
		end if ! np1m

		if (im0 .le. nx) then
		    if (lodd .and. nm0 .lt. 0) then
			bxm0r  = -cm0 *bx(im0)
			bxm0i  = -sm0 *bx(im0)
		    else
			bxm0r  = +cm0 *bx(im0)
			bxm0i  = +sm0 *bx(im0)
		    end if ! lodd
		else
		    bxm0r  = ZERO
		    bxm0i  = ZERO
		end if ! im0
		if (im1p .le. nx) then
		    if (.not. lodd .and. nm1p .lt. 0) then
			bxm1pr = -cm1p*bx(im1p)
			bxm1pi = -sm1p*bx(im1p)
		    else
			bxm1pr = +cm1p*bx(im1p)
			bxm1pi = +sm1p*bx(im1p)
		    end if ! lodd
		else
		    bxm1pr = ZERO
		    bxm1pi = ZERO
		end if ! im1p
		if (im1m .le. nx) then
		    if (.not. lodd .and. nm1m .lt. 0) then
			bxm1mr = -cm1m*bx(im1m)
			bxm1mi = -sm1m*bx(im1m)
		    else
			bxm1mr = +cm1m*bx(im1m)
			bxm1mi = +sm1m*bx(im1m)
		    end if ! lodd
		else
		    bxm1mr = ZERO
		    bxm1mi = ZERO
		end if ! im1m

		bx10r  = bxp0r  +sign*bxm0r
		bx10i  = bxp0i  +sign*bxm0i
		bx11pr = bxp1pr +sign*bxm1pr
		bx11pi = bxp1pi +sign*bxm1pi
		bx11mr = bxp1mr +sign*bxm1mr
		bx11mi = bxp1mi +sign*bxm1mi
		sum0r  = sum0r  +by(n)*bx10r
		sum0i  = sum0i  -by(n)*bx10i
		sum1pr = sum1pr +by(n)*bx11pr
		sum1pi = sum1pi -by(n)*bx11pi
		sum1mr = sum1mr +by(n)*bx11mr
		sum1mi = sum1mi -by(n)*bx11mi
	    end do ! n
			
	end if
C-
C  -----------------------------------------------------------------------------

	axr= (i/a)*(TWO*alpha*cosphi*sum0r -ky*(sum1pr +sum1mr))
	axi= (i/a)*(TWO*alpha*cosphi*sum0i -ky*(sum1pi +sum1mi))
	ayr= (i/a)*(TWO*alpha*sinphi*sum0r +kx*(sum1pi -sum1mi))
	ayi= (i/a)*(TWO*alpha*sinphi*sum0i -kx*(sum1pr -sum1mr))

	s0 = axr*axr +axi*axi +ayr*ayr +ayi*ayi
	s1 = axr*axr +axi*axi -ayr*ayr -ayi*ayi
	s2 = TWO*(axr*ayr +axi*ayi)
	s3 = TWO*(axi*ayr -axr*ayi)

CC
c       print 220
c	print 220,'brighte:a,x,y,phi          ',a,x,y,phi
c	print 220,'brighte:sum0r,sum1pr,sum1mr',sum0r,sum1pr,sum1mr
c	print 220,'brighte:sum0i,sum1pi,sum1mi',sum0i,sum1pi,sum1mi
c	print 220,'brighte:ex**2,ey**2        ',(s0+s1)/TWO,(s0-s1)/TWO
c	print 220,'brighte:s0   ,s1    ,s2, s3',s0,s1,s2,s3
c220     format(' ',a,1p4e15.6)

	return
	end ! bright3

	subroutine bright_bessjn(x,nmax,bs)
C+
C 
C FUNCTIONAL DESCRIPTION:	
C 
C  Routine to calculate Bessel functions from order zero to order nmax.
C  Modified form of BESSJN in order to increase speed (speed increased more than
C  a factor of 2 in this version).
C 
C AUTHORS: 
C 
C  Roger J. Dejus
C 
C  The Advanced Photon Source
C  Experimental Facilities Division
C  Argonne National Laboratory
C 
C CREATION DATE: 
C 
C   3-FEB-1991
C 
C FORMAL PARAMETERS:
C  
C  Input arguments:
C  x  			Argument of Bessel function 
C                       (must not be zero; may be less than zero).
C  nmax  		Maximum order for which functions are calculated.
C                       (must be greater than or equal to zero; use symmetry
C			 relations to obtain values for negative order).
C   
C  Output arguments:
C  bs    		Array of Bessel function values. bs(0) = J0(x),
C			bs(1) = J1(x) ... bs(nmax) = Jnmax(x).
C  
C COMMON BLOCKS:
C  
C  None
C  
C DESIGN ISSUES:
C  
C  Adapted from "Numerical Recipes" by W.H. Press, B.P. Flannery, S.A. Teukolsky
C  and W.T. Vetterling, Cambridge University Press (1988), p. 175.
C  Note: Uses only Miller's downward recurrence relation.
C  
C MODIFICATION HISTORY:
C 
C	 Date     | Name  | Description
C ----------------+-------+-----------------------------------------------------
C 05-SEP-1994     | RJD   | Modified to allow for x < 0.0. x must not be 0.0
C ----------------+-------+-----------------------------------------------------
C [change_entry]
C-
C  Declarations of scalars:
	integer*4	nmax,m,i,j,jsum
	real*8		x,tox,bj,bjm,bjp,sum,sign

C  Declarations of arrays:
	real*8		bs(0:*)
 
C  Labeled constants:
	real*8		ZERO,ONE,TWO,BIGNO,SMLNO
	parameter	(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
	parameter	(BIGNO=1.0D10,SMLNO=1.0D-10)

	tox   = TWO/abs(x)
	m     = 2*((nmax+1)/2) ! Make M even
	jsum  = 0
	sum   = ZERO
	bjp   = ZERO
	bj    = ONE
	do j=m,1,-1
	    bjm = j*tox*bj-bjp
	    bjp = bj
	    bj  = bjm
	    if (abs(bj) .gt. BIGNO) then
		bj    = bj   *SMLNO
		bjp   = bjp  *SMLNO
		sum   = sum  *SMLNO
		do i=j,nmax
		    bs(i) = bs(i)*SMLNO
		end do ! i
	    end if
	    if (jsum .ne. 0) sum = sum+bj
	    jsum = 1-jsum
	    if (j .le. nmax) bs(j) = bjp
	end do ! j
	bs(0) = bj

	sum   = TWO*sum-bj
	if (x .lt. ZERO) then ! Change sign for odd indices
	    sign = -ONE
	    do j=0,nmax ! Normalize
	        sign  = -sign
	        bs(j) = sign*bs(j)/sum
	    end do ! j
	else ! x > 0.0
	    do j=0,nmax ! Normalize
	        bs(j) = bs(j)/sum
	    end do ! j
	end if ! x
	   
	return
	end ! bright_bessjn





