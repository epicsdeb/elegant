	program tc
C+
C PROGRAM DESCRIPTION: 
C  Program to calculate on-axis brilliance tuning curves for an ideal undulator
C  insertion device (regular planar device or a helical device). The program
C  uses the Bessel function approximation which is valid for an ideal device,
C  e.g., no magnetic field errors. The effect of the particle beam emittance and
C  the beam energy spread is taken into account.
C
C AUTHORS: 
C  Roger J. Dejus
C  The Advanced Photon Source
C  Experimental Facilities Division
C  Argonne National Laboratory
C 
C CREATION DATE: 
C  22-APR-1996
C 
C INPUT PARAMETERS:
C  The input parameters are divided into sections related to the storage ring,
C  the undulator device, and the brilliance calculation.
C Machine Parameters:
C  ENERGY		Storage ring energy              (GeV)
C  CUR		        Storage ring current             (mA)
C  SIGE		        Energy spread (sigma(E)/E)
C  SIGX		        RMS beam size (horizontal)       (mm)
C  SIGY		        RMS beam size (vertical)         (mm)
C  SIGX1	        RMS beam divergence (horizontal) (mrad)
C  SIGY1	        RMS beam divergence (vertical)   (mrad)
C Undulator Parameters:
C  PERIOD	        Undulator period length          (cm)
C  N			Number of periods        
C Scan Parameters:
C  EMIN			Lower energy limit for FIRST harmonic (eV)
C  EMAX			Upper energy limit for FIRST harmonic (eV)
C  NE	                Number of energy points
C  IHMIN		Minimum harmonic of interest
C  IHMAX		Maximum harmonic of interest
C  IHSTEP		Step size for the harmonics
C Brilliance Parameters:
C  IHEL			Type of device; ihel=0: regular planar, ihel=1: helical
C  METHOD               Method
C  			METHOD={0,1}Non-zero emittance;
C                                   infinite-N +convolution (Dejus' approach)
C  			METHOD=2    Non-zero emittance;
C                                   infinite-N +convolution (Walker's approach)
C  			METHOD=3    Non-zero emittance; finite-N (Walker's)
C
C  IK			Print K-values & powers; ik=0: quiet, ik=1: print
C  NEKS			Number of energy points for peak search: Default 100, or
C			use 0 for default.
C COMMON BLOCKS:
C  None.
C  
C DESIGN ISSUES:
C  The down shift in energy due to the beam emittance is most noticable for
C  small values of K, and therefore, for each harmonic, the energy shift is
C  calculated for K = Kmin. The shifted peak is then used to define the energy
C  interval over which the on-axis brilliance is calculated. The peak of this
C  function is stored and subsequently saved in a file vs. energy (eV). The
C  user enters the scanning range in energy (eV) for the first harmonic, even in
C  the case when higher harmonics are chosen (using IHMIN > 1). The beam energy
C  spread is included by using a straightforward convolution at the fixed energy
C  of the peak. Beam energy spreads typically in the range 1% to 0.01% can be
C  used. It is valid to set the beam energy spread to zero. The number of energy
C  points for peak search NEKS is used at each K-value and would typically be
C  set to the default value unless the down shift in energy is being sought with
C  high accuracy. For high accuracy typically use NEKS = 1000.
C  
C COPYRIGHT:
C  This routine must only be used at The Advanced Photon Source and must not
C  be tranferred or used at any other location without the consent of the
C  author.
C  
C FILES USED:
C Input files
C  tc.dat
C Output files
C  tc.plt
C
C KEYWORDS:
C  Undulator Tuning Curve, Undulator On-Axis Brilliance.
C  
C LINK/LIBRARY ISSUES:
C  Calls routine USB that calls routines BRIGHTE and HUNT.  
C  BRIGHTE calculates the brightness and HUNT searches an array of real
C  numbers (from Numerical Recipes).
C  
C PORTABILITY ISSUES:
C  Runs on DEC 3000/400 AXP alpha (Tru64Unix v5.0), SUN (Solaris: SunOS
C  Release v5.6), and Windows 95/98/NT (Pentium and higher).
C  
C TIMING:
C  Generally, the execution is fast. For example, the first three odd
C  harmonics (1, 3, 5) for Undulatar A at the APS over the full tuning range is 
C  calculated in about 10 s using the default parameters for the infinite-N
C  method with convolution and zero beam energy spread. The finite-N method is
C  about 10 times slower. Introduction of the beam energy spread increases the
C  execution time by typically 40%. The timing above is given for the default
C  value of NEKS (=100) and the number of points/harmonic ne = 20.
C
C VERSION:
C  1.91
C  
C MODIFICATION HISTORY:
C	 Date     | Name  | Description
C ----------------+-------+-----------------------------------------------------
C 10-MAY-1996     | RJD   | Tuning variable is K.
C ----------------+-------+-----------------------------------------------------
C 21-MAY-1996     | RJD   | Changed tuning variable to E1 and added descriptive
C                 |       | text describing the input/output. First official
C                 |       | release, v1.2.
C ----------------+-------+-----------------------------------------------------
C 13-MAR-1997     | RJD   | Added beam energy spread using a straightforward
C                 |       | convolution at the fixed energy of the peak.
C                 |       | Modified calculation of the variable nek. Added the
C                 |       | variable ik which controls printing of K-values.
C                 |       | Current version is v1.3.
C ----------------+-------+-----------------------------------------------------
C 18-MAR-1997     | RJD   | Modified upper limit of peak search from 1.0 to
C                 |       | fc2 = 1.002. Added parameter SPECLIM which defines
C                 |       | the minimum Brilliance to retain in the calculation.
C                 |       | Current version is v1.4.
C ----------------+-------+-----------------------------------------------------
C 21-MAR-1997     | RJD   | Modified lower limit of peak search from fc = 
C                 |       | 0.990d0*ep/eiz to fc = 0.985d0*ep/eiz so that the 
C                 |       | the peak of the higher odd harmonics will be found
C                 |       | (or is more likely to be found) when the beam energy
C                 |       | spread is taken into account.
C                 |       | Current version is v1.5.
C ----------------+-------+-----------------------------------------------------
C 15-JUL-1997     | RJD   | Added information about total emitted power and 
C                 |       | on-axis power density to the printout when the print
C                 |       | variable ik is set to 1. Current version is v1.6.
C ----------------+-------+-----------------------------------------------------
C 29-SEP-1997     | RJD   | Added printout of zero emittance energy (first
C                 |       | column in output file). Current version is v1.7.
C ----------------+-------+-----------------------------------------------------
C 06-OCT-1997     | RJD   | The parameter NEKS which determines the number of
C                 |       | energy points for the peak search at each K-value
C                 |       | was added to the input file. Default = 100 (or enter
C                 |       | 0). Min and max is 100 and 10000, respectively.
C ----------------+-------+-----------------------------------------------------
C 14-NOV-1997     | RJD   | The variable neks declared as integer*4.
C                 |       | Current version is v1.9.
C ----------------+-------+-----------------------------------------------------
C 16-JUL-2000     | RJD   | Minor change in the code to compile error-free on
C                 |       | Unix and Windows (no change in results vs. v1.9).
C                 |       | Current version is v1.91.
C ----------------+-------+-----------------------------------------------------
C 18-JAN-2001     | RJD   | Changed parameter NEKS1 from 200 to 500 and upper
C                 |       | range (ekmax = 1.00*eiz -> ekmax = 1.002*eiz) to
C                 |       | avoid peak searching error for second harmonic at
C                 |       | Kmin for small emittance and large N. Results for
C                 |       | for odd harmonics may also differ slightly from
C                 |       | previous version because peak shifts are different.
C                 |       | Current version is v1.92.
C ----------------+-------+-----------------------------------------------------
C-
C  Size parameters:
	integer*4	E_SZ,H_SZ,K_SZ
	parameter	(E_SZ=10001,H_SZ=20,K_SZ=1000)

C  Declarations of scalars:
	character*80	title
	character*8	tbuff
	character*9	dbuff

	logical*4	lodd

	integer*4	ierrorb,ierrorp,ierrorg
	integer*4	n,ne,nek,neks,ns,method
	integer*4	i,j,ie,ir,je
	integer*4	l,nch1
	integer*4	ih,ihmin,ihmax,ihstep,ihel,nharm,ik

	real*4		dtime,delta
	real*8		energy,cur,sige,sigx,sigy,sigx1,sigy1
	real*8		period,kx,ky,k,k2,ekmin,ekmax,emin,emax
	real*8		kmin,kmax
	real*8		gamma,g2,xg,yg,lamdar,er,e1z,e1zmin,e1zmax
	real*8		ek,eiz,ep,sp,fc,fc1,fc2,de,dek,dep1,dep2
	real*8		alpha2
	real*8		sigee,smax
	real*8		gk,gkp,gkh

C  Declarations of arrays:
	real*4		td(2)

	real*8		e(E_SZ),spec(E_SZ)
	real*8		ei(K_SZ,H_SZ),eb(K_SZ,H_SZ),sb(K_SZ,H_SZ)
	real*8		kxb(K_SZ),kyb(K_SZ),ptot(K_SZ),pd(K_SZ)
	real*8		dep(2)

C  Fundamental physical constants; Physics Today Aug. 1990:
	real*8		C,ME,MEE,EC,H,HBAR,MUZ,EPSZ
	parameter	(C    =2.99792458D8)	! Speed of light [m/s]
	parameter	(ME   =9.1093897D-31)	! Electron rest mass [kg]
	parameter	(MEE  =0.51099906D0)	! Electron rest mass [MeV]
	parameter	(EC   =1.60217733D-19)	! Elementary charge [C]
	parameter	(H    =6.6260755D-34)	! Planck's constant [Js]
	parameter	(HBAR =1.05457266D-34)	! Planck's constant/2Pi [Js]
	parameter	(MUZ  =1.2566370614D-6)	! Permeability of vacuum [NA-2]
	parameter	(EPSZ =8.854187817D-12)	! Permittivity of vacuum [Fm-1]

C  Conversion factors:
	real*8		C_EVANG,C_CM_ANG,C_MRAD_RAD,C_RAD_MRAD,C_MA_A
	real*8		C_MM_CM,C_CM_MM,C_CM_M
	parameter	(C_MM_CM=1.0D-1,C_CM_MM=1.0D+1,C_CM_M=1.0D-2)
	parameter	(C_EVANG=H*C/EC*1.0D10,C_CM_ANG=1.0D8)
	parameter	(C_MRAD_RAD=1.0D-3,C_RAD_MRAD=1.0D+3,C_MA_A=1.0D-3)

C  Labeled constants:
	character*(*)   VN
	parameter       (VN='v. 1.92')
	real*8		PI,PIHALF,TWOPI
	parameter	(PI    =3.1415 92653 58979 32384 62643D0)
	parameter	(PIHALF=1.5707 96326 79489 66192 31322D0)
	parameter	(TWOPI= 6.2831 85307 17958 64769 25287D0)
	real*8		ZERO,ONE,TWO,HALF
	parameter	(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,HALF=0.5D0)
	REAL*8		PTOT_FAC,PD_FAC
	PARAMETER	(PTOT_FAC=PI/3.0D0*EC/EPSZ/(MEE**2)*1.0D+6)   ! 0.07257
	PARAMETER	(PD_FAC  =21.0D0/(16.0D0*PI*MEE**2)*PTOT_FAC) ! 0.11611
	real*8		BW
	parameter	(BW=1.0D-3) ! 0.1%
	real*8		SPECLIM
	parameter	(SPECLIM=1.0D10)	! Min Brilliance to retain
	integer*4	NSIGMA,NPPSIGMA,NEKS1
	parameter	(NSIGMA=3,NPPSIGMA=6,NEKS1=500)
	integer*4	NEKS_MIN,NEKS_MAX
	parameter	(NEKS_MIN=100,NEKS_MAX=10000)

	data            td /0.0,0.0/

C  Common blocks:

C  Statement functions:
	gkp(k) = k*((k**6)   +(24.0d0*(k**4)/7.0D0) +	! planar device
	1        (4.0d0*k*k) +(16.0D0/7.0D0))/((1.0d0+(k*k))**3.5d0)
	gkh(k) = 32.0d0/7.0d0*k/((1.0d0+(k*k))**3.0d0)	! helical device

C+ First executable statement here
C  -----------------------------------------------------------------------------
	delta = dtime(td)	! start clock
	call date(dbuff)
	call time(tbuff)

C+ Read parameters from data file
C  -----------------------------------------------------------------------------
	open (unit=1,file='tc.dat',status='old',access='sequential',
	1     form='formatted',action='read')

	read (1,100) title
	read (1,*)   energy,cur,sige		! GeV, mA
	read (1,*)   sigx,sigy,sigx1,sigy1	! mm, mm, mrad, mrad
	read (1,*)   period,n			! cm
	read (1,*)   emin,emax,ne		! eV, eV (first harmonic)
	read (1,*)   ihmin,ihmax,ihstep
	read (1,*)   ihel,method,ik,neks	! 1:dejus,2:walker,3:walker f-N
	close (unit=1)				! 0:planar,1:helical
						! 0:quiet,1:print K-values & powers
						! # of energy points for peak seearch
C-
C  -----------------------------------------------------------------------------

C+ Default values
C  -----------------------------------------------------------------------------
	if (ihmin .eq. 0) ihmin = 1
	if (ihmax .eq. 0) ihmax = 5
	if (ihstep .eq. 0) ihstep = 2
	if (ne .eq. 0) ne = 20
	if (method .eq. 0) method = 1
	if (neks .eq. 0) neks = 100

C+ Check input
C  -----------------------------------------------------------------------------
	nharm = (ihmax -ihmin +ihstep)/ihstep
	if (ne .gt. K_SZ) goto 904
	if (nharm .gt. H_SZ) goto 908
	if (ihel .lt. 0 .or. ihel .gt. 1) goto 910
	if (ik .lt. 0 .or. ik .gt. 1) goto 911
	if (method .lt. 0 .or. method .gt. 3) goto 912
	if (neks .lt. NEKS_MIN .or. neks .gt. NEKS_MAX) goto 913

	l = 80
	do while(title(l:l) .eq. ' ') ! strip blanks
	    l = l-1
	end do ! while
	nch1 = l

C+ Definition of constants
C  -----------------------------------------------------------------------------
	gamma  = energy/MEE*1.0d3
	g2     = gamma*gamma
	lamdar = period*C_CM_ANG/(TWO*g2) ! reduced wavelength A
	er     = C_EVANG/lamdar		  ! reduced energy    eV
	xg     = gamma*sigx1*C_MRAD_RAD
	yg     = gamma*sigy1*C_MRAD_RAD
c	alpha2 = xg*xg +yg*yg
 	alpha2 = yg*yg

	kx = ZERO
	e1zmin= emin
	e1zmax= emax
	de = 0.0d0
	if (ne .gt. 1) de = (e1zmax -e1zmin)/(ne -1)

	if (e1zmax .ge. er) goto 915
	kmin = sqrt(TWO*(er/e1zmax -ONE))
	kmax = sqrt(TWO*(er/e1zmin -ONE))
	if (ihel .eq. 1) kmin = kmin/sqrt(TWO)
	if (ihel .eq. 1) kmax = kmax/sqrt(TWO)

C+ Open files
C  -----------------------------------------------------------------------------
	open (unit=2,file='tc.plt',status='unknown')

C+ Write input parameters
C  -----------------------------------------------------------------------------
	print 200,'****************** Start TC ',VN,
	1         ' at ',dbuff,' ',tbuff,' ******************'
	print 200,title(1:nch1)
	if (ihel .eq. 0) then
	  print 200,'Mode  : Regular Planar Undulator'
	else if (ihel .eq. 1) then
	  print 200,'Mode  : Helical Undulator'
	endif
	if (method .eq. 1) then
	  print 200,'Method: Infinite-N with convolution'
	else if (method .eq. 2) then
	  print 200,'Method: Infinite-N with convolution (Walker)'
	else if (method .eq. 3) then
	  print 200,'Method: Finite-N'
	endif
	print 257,'K     : ',kmin,' - ',kmax
	print 205,'Nharm :',nharm
	print *

	write (2,200) '****************** Start TC ',VN,
	1             ' at ',dbuff,' ',tbuff,' ******************'
	write (2,200)
	write (2,200) title(1:nch1)
	write (2,250) 'energy ',energy,' GeV','current',cur ,' mA',
	1             'sige',sige
	write (2,260) 'sigx   ',sigx  ,'  mm','sigy   ',sigy,' mm',
	1             'sigx1 ',sigx1,' mr','sigy1 ',sigy1,' mr'
	write (2,255) 'period ',period,'  cm','N      ',n
	write (2,256) 'emin   ',emin  ,'    ','emax',emax ,'   ',
	1             'ne     ',ne 
	write (2,258) 'hmin   ',ihmin ,'    ','hmax   ',ihmax,'   ',
	1             'hstep  ',ihstep 
	write (2,258) 'helical',ihel  ,'    ','method ',method,'   ',
	1	      'print_k',ik    ,'   ','neks   ',neks

	write (2,200)
	if (ihel .eq. 0) then
	  write (2,200) 'Mode  : Regular Planar Undulator'
	else if (ihel .eq. 1) then
	  write (2,200) 'Mode  : Helical Undulator'
	endif
	if (method .eq. 1) then
	  write (2,200) 'Method: Infinite-N with convolution'
	else if (method .eq. 2) then
	  write (2,200) 'Method: Infinite-N with convolution (Walker)'
	else if (method .eq. 3) then
	  write (2,200) 'Method: Finite-N'
	endif
	write (2,257) 'K     : ',kmin,' - ',kmax
	write (2,205) 'Nharm :',nharm
	write (2,200)
	write (2,200) 'On-axis Brilliance (ph/s/mrad^2/mm^2/0.1%bw):'
	if (ik .eq. 1) then
	  if (ihel .eq. 1) then
	    write (2,200) 'Energy (eV) (w/o w/ em.) Brilliance       Kx',
	1		  '      Ky  Ptot(W)  Pd(W/mr^2)'
	  else
	    write (2,200) 'Energy (eV) (w/o w/ em.) Brilliance       Ky',
	1		  '  Ptot(W)  Pd(W/mr^2)'
	  endif
	else
	  write (2,200) 'Energy (eV) (w/o w/ em.) Brilliance'
	endif
C-
C  -----------------------------------------------------------------------------

c  Determine peak shifts for first and second harmonics at Kmin
	nek = NEKS1			! # of energy points for peak search
	ek  = e1zmax
	ky  = sqrt(TWO*(er/ek -ONE))
	if (ihel .eq. 1) ky = ky/sqrt(TWO)
	if (ihel .eq. 1) kx = ky

	do i=1,2
	  e1z = ek
	  eiz = i*e1z

	  if (i .eq. 1) then		! first harmonic 
	    ekmin = 0.95d0*eiz
	    ekmax = 1.01d0*eiz
	  else				! second harmonic
	    ekmin = 0.820d0*eiz
	    ekmax = 1.002d0*eiz
	  endif
	  call usb(energy,cur,sigx,sigy,sigx1,sigy1,period,n,
	1          kx,ky,ekmin,ekmax,nek,method,e,spec,ns,ierrorb)
	  if (ierrorb .ne. 0) goto 916

	  call peak(e,spec,ns,ep,sp,ierrorp)
	  if (ierrorp .ne. 0) goto 919
	  dep(i) = eiz -ep
	enddo ! i
 	dep1 = dep(1)			! save peak shifts of harmonics
 	dep2 = dep(2)

C+ Main loop over harmonics and K-values
C  -----------------------------------------------------------------------------
	ih = 0
	do i=ihmin,ihmax,ihstep
	  ih = ih +1
	  if (i .ne. i/2*2) then
	    lodd = .true.
	  else
	    lodd = .false.
	  endif

c  Try Kmin initially to find fc for each harmonic
c  Omit Brilliances < SPECLIM (will be set to 0.0) and peak positions will be
c  calculated from the zero emittance formula.
	  nek  = neks			! # of energy points for peak search
	  smax = ZERO
	  je   = ne +1
	  do while (smax .lt. SPECLIM .and. je .gt. 1)
	    je = je -1
	    ek = e1zmin +(je-1)*de
	    ky = sqrt(TWO*(er/ek -ONE))
	    if (ihel .eq. 1) ky = ky/sqrt(TWO)
	    if (ihel .eq. 1) kx = ky
	    if (ik .eq. 1 .and. ih .eq. 1) then
	      kyb(je) = ky
	      if (ihel .eq. 1) kxb(je) = kx
	      k  = ky
	      k2 = kx*kx +ky*ky
	      ptot(je) = PTOT_FAC*n*k2*(energy**2)*cur*C_MA_A/
	1		 (period*C_CM_M) ! W
	      if (ihel .eq. 1) then
		gk = gkh(k)		! helical device
	      else
		gk = gkp(k)		! planar device
	      endif
	      pd(je) = PD_FAC*n*k*gk*(energy**4)*cur*C_MA_A/
	1	       (period*C_CM_M) ! W/mrad^2
	    endif

	    e1z = ek
	    eiz = i*e1z
	    if (lodd) then		! odd harmonics
	      ekmin = eiz -i*dep1
	      ekmax = eiz +i*dep1/TWO
	      if (i .eq. 1) ekmin = ekmin -dep1
	      ir = e1z/dep1
	      if (i .gt. ir) then
	        print 200,'&tc-W-HARMERR, Warning,', 
	1  		  ' Overlapping range for initial peak search'
  	        print 206,'- for harmonic number ',i
	        goto 800
	      endif
	    else 			! even harmonics
	      ekmin = eiz -4.0d0*dep2
	      ekmax = eiz
	    endif

	    call usb(energy,cur,sigx,sigy,sigx1,sigy1,period,n,
	1            kx,ky,ekmin,ekmax,nek,method,e,spec,ns,ierrorb)
	    if (ierrorb .ne. 0) goto 917

	    smax = ZERO
	    do ie=1,ns
	      if (spec(ie) .gt. smax) then
		smax = spec(ie)
	      endif
	    end do

	    ei(je,ih) = eiz
	    eb(je,ih) = eiz
	    sb(je,ih) = ZERO
	  enddo ! while

	  if (smax .lt. SPECLIM) then
	print 200,'&tc-W-HARMERR, Warning, Harmonic intensity too small'
	print 206,'- for harmonic number ',i
	    goto 800
	  endif

	  call peak(e,spec,ns,ep,sp,ierrorp)
	  if (ierrorp .ne. 0) goto 920

c  define fc
 	  fc = 0.985d0*ep/eiz
	  if (i .gt. int(1.0d0/(1.0d0 -fc))) then
	    print 200,'&tc-W-HARMERR, Warning,', 
	1	      ' Overlapping range for peak search'
	    print 206,'- for harmonic number ',i
	    goto 800
	  endif

	  do j=1,je
	    ek = e1zmin +(j-1)*de
	    ky = sqrt(TWO*(er/ek -ONE))
	    if (ihel .eq. 1) ky = ky/sqrt(TWO)
	    if (ihel .eq. 1) kx = ky
	    if (ik .eq. 1 .and. ih .eq. 1) then
	      kyb(j) = ky
	      if (ihel .eq. 1) kxb(j) = kx
	      k  = ky
	      k2 = kx*kx +ky*ky
	      ptot(j) = PTOT_FAC*n*k2*(energy**2)*cur*C_MA_A/
	1		(period*C_CM_M) ! W
	      if (ihel .eq. 1) then
		gk = gkh(k)		! helical device
	      else
		gk = gkp(k)		! planar device
	      endif
	      pd(j) = PD_FAC*n*k*gk*(energy**4)*cur*C_MA_A/
	1	      (period*C_CM_M) ! W/mrad^2
	    endif

	    fc1   = fc
	    if (lodd) then
	      fc2 = 1.002d0
	    else
	      fc2 = 1.000d0
	    endif
	    e1z   = ek
 	    eiz   = i*e1z
 	    ekmin = fc1*eiz
 	    ekmax = fc2*eiz

c  adjust ekmin, ekmax, and number of points if beam energy spread applied
	    if (sige .gt. 0.0d0) then
	      nek   = neks
	      dek   = (ekmax -ekmin)/nek
	      sigee = 2.0d0*sige*eiz		! estimated width (eV)
	      ekmin = ekmin -sigee*NSIGMA	! adjust for convolution
	      ekmax = ekmax +sigee*NSIGMA	! adjust for convolution
	      if (sigee/NPPSIGMA .lt. dek) dek = sigee/NPPSIGMA
  	      nek   = (ekmax-ekmin)/dek +1.0d0	! # of points
	      if (nek .ge. E_SZ) goto 914
	    endif

c  get undulator on-axis brilliance for given K value in energy range ekmin to ekmax
c  returns energy array e (eV) and spec (ph/s/mrad^2/mm^2/0.1%bw)
	    call usb(energy,cur,sigx,sigy,sigx1,sigy1,period,n,
 	1            kx,ky,ekmin,ekmax,nek,method,e,spec,ns,ierrorb)
	    if (ierrorb .ne. 0) goto 918

c for debugging to check spectrum
c	    if (i .eq. 13 .and. j .eq. 9) then
c	      do ie=1,ns
c	        write (4,*),e(ie),spec(ie)
c	      enddo
c	      call exit(0)
c	    endif

c  apply beam energy spread
c  Note: the size of the arrays e and spec will be adjusted (smaller by 
c        the half width of the Gaussian at each limit because of the
c        convolution).

	    if (sige .gt. 0.0d0) then
	      call gauss_convolve(e,spec,ns,sige,ierrorg)
	      if (ierrorg .ne. 0) goto 922
	    endif

c for debugging to check spectrum
c	    if (i .eq. 13 .and. j .eq. 9) then
c	      do ie=1,ns
c	        write (4,*),e(ie),spec(ie)
c	      enddo
c	      call exit(0)
c	    endif

	    call peak(e,spec,ns,ep,sp,ierrorp)
	    if (ierrorp .ne. 0) goto 925

	    ei(j,ih) = eiz
	    eb(j,ih) = ep
	    sb(j,ih) = sp

	  end do	! j
	  print 205, 'Harmonic ',i,' completed'
	end do 		! i
C- End main loop 
C  -----------------------------------------------------------------------------

C+ Save results
C  -----------------------------------------------------------------------------
800	continue
	ih = 0
	do i=ihmin,ihmax,ihstep
	  ih = ih +1
	  write (2,200)
	  write (2,205) 'Harmonic ',i
	  do j=1,ne
	    if (ik .eq. 1) then
	      if (ihel .eq. 1) then
 	        write (2,220) ei(j,ih),eb(j,ih),sb(j,ih),
	1		      kxb(j),kyb(j),ptot(j),pd(j)
	      else
 	        write (2,225) ei(j,ih),eb(j,ih),sb(j,ih),
	1		      kyb(j),ptot(j),pd(j)
	      endif
	    else
 	      write (2,220) ei(j,ih),eb(j,ih),sb(j,ih)
	    endif
	  end do	! j
	end do 		! i
	close (unit=2)
C-
C  -----------------------------------------------------------------------------

C+ Exit
C  -----------------------------------------------------------------------------
	delta = dtime(td)	! read clock
	print *
	print 210,'Time:',delta,' s.'
	print 200,'&tc-S-NORMAL, Successful completion'
	print *

	call exit(0)

100	format(a)
200	format(' ',8a)
205	format(' ',a,i2,2a,f8.3,i5,f8.3,f10.3)
206	format(' ',a,i3,a,i3)
207	format(' ',a,i2,a)
208	format(' ',a,a,i4,a,i4)
209	format(' ',a,f10.1,a)
210	format(' ',a,f10.3,a,f10.3,a)
212	format(' ',a,2(f10.1,a))
213	format(' ',a,f12.7,a)
214	format(' ',a,a,i5,a,i5)
220	format(' ',f10.3,tr2,f10.3,tr2,1pe12.5,0p2f8.3,tr2,f7.1,tr2,f10.1)
225	format(' ',f10.3,tr2,f10.3,tr2,1pe12.5,0p1f8.3,tr2,f7.1,tr2,f10.1)
250	format(' ',(a,f7.3,a,tr3),(a,f6.1,a,tr3),(a,f12.7))
255	format(' ',(a,f7.3,a,tr3),(a,i6,a,tr3),2(a,f6.3,a,tr2))
256	format(' ',(a,f7.1,a,tr3),(a,f9.1,a,tr3),(a,i6))
257	format(' ',a,f5.3,a,f5.3)
258	format(' ',(a,i7,a,tr3),(a,i6,a,tr3),2(a,i6,a,tr3))
260	format(' ',(a,f7.3,a,tr3),(a,f6.3,a,tr3),2(a,f7.4,a,tr3))

C+ Error returns
C  -----------------------------------------------------------------------------
900	continue
	print 200,'&tc-F-SUBERR, Future reference ...'
	goto 999

904	continue
	print 200,'&tc-F-BNDERR, Boundary error'
	print 208,'- energy dimension out of bounds;',
	1         ' number of points ',ne,' is greater than ',K_SZ
	goto 999

908	continue
	print 200,'&tc-F-BNDERR, Boundary error'
	print 208,'- harmonic dimension out of bounds;',
	1         ' number of points ',nharm,' is greater than ',H_SZ
	goto 999

910	continue
	print 200,'&tc-E-INVDAT, Invalid data'
	print 207,'- check input data file; hel ',ihel
	goto 999

911	continue
	print 200,'&tc-E-INVDAT, Invalid data'
	print 207,'- check input data file; print_k ',ik
	goto 999

912	continue
	print 200,'&tc-E-INVDAT, Invalid data'
	print 207,'- check input data file; method ',method 
	goto 999

913	continue
	print 200,'&tc-E-INVDAT, Invalid data'
	print 207,'- check input data file; neks ',neks 
	goto 999

914	continue
	print 200,'&tc-F-BNDERR, Boundary error'
	print 214,'- energy dimension out of bounds (sige too small);',
	1         ' number of points in energy ',nek,
	2         ' is greater than ',E_SZ-1
	print 213,'- check input data file; sige ',sige
	goto 999

915	continue
	print 200,'&tc-E-INVDAT, Invalid data'
	print 209,'- check input data file; emax ',emax,' eV.' 
	print 209,'- emax must be less than',er,' eV.'
	goto 999

916	continue
	print 200,'&tc-F-SUBERR, Error in subroutine USB'
  	print 206,'- for peak shift search for harmonic number ',i
	goto 999

917	continue
	print 200,'&tc-F-SUBERR, Error in subroutine USB'
  	print 206,'- for harmonic number ',i
	goto 999

918	continue
	print 200,'&tc-F-SUBERR, Error in subroutine USB'
  	print 206,'- for harmonic number ',i,' and point ',j
	goto 999

919	continue
	print 200,'&tc-F-PEAKERR, Peak error'
  	print 206,'- for peak shift search for harmonic number ',i
	print 212,'- no peak found in the energy range: ',
	1          ekmin,' -',ekmax,' eV.'
	goto 999

920	continue
	print 200,'&tc-F-PEAKERR, Peak error'
  	print 205,'- for harmonic number ',i
	print 212,'- no peak found in the energy range: ',
	1          ekmin,' -',ekmax,' eV.'
	goto 999

922	continue
	print 200,'&tc-F-CONVERR, Convolution error'
  	print 206,'- for harmonic number ',i,' and point ',j
	goto 999

925	continue
	print 200,'&tc-F-PEAKERR, Peak error'
  	print 206,'- for harmonic number ',i,' and point ',j
	print 212,'- no peak found in the energy range: ',
	1          ekmin,' -',ekmax,' eV.'
	goto 999

999	continue
	close (unit=2)
	print *
	print 200,'&tc-F-PRGERR, Program error'
	print 200,'- unsuccessful completion due to error status'
	print *
C-
C  -----------------------------------------------------------------------------
	call exit(0)
	end ! tc

	subroutine peak(e,spec,ns,ep,sp,ierror)
 
c FUNCTIONAL DESCRIPTION: 
c  Routine to find the abscissa and ordinate at the peak of the array spec.

C  Declarations of scalars:
	integer*4	ns,ierror,ie,ip
	real*8		ep,sp,smax
 
C  Declarations of arrays:
	real*8		e(ns),spec(ns)

C  Labeled constants:
	integer*4	IB
	parameter	(IB=3)
	real*8		ZERO
	parameter	(ZERO=0.0D0)

	ierror = 0
	smax = ZERO
	do ie=1,ns
	  if (spec(ie) .gt. smax) then
	    smax = spec(ie)
	    ip   = ie
	  endif
	end do

	if (ip .ge. (ns-IB) .or. ip .le. IB) then
	  ep = ZERO
	  sp = ZERO
	  ierror = -1
	else
	  ep = e(ip)
	  sp = smax
	endif

	return
	end 	! peak

	subroutine gauss_convolve(e,spec,ns,sige,ierror)
 
c FUNCTIONAL DESCRIPTION: 
c  Routine to make a Gaussian convolution of array spec.

C  Size parameters:
	integer*4	E_SZ
	parameter	(E_SZ=10001)

C  Declarations of scalars:
	integer*4	ns,ierror,ierrorp
	integer*4	ip,ie,nsigma2,np,ne1,ne2
	real*8		sige,ep,sp
	real*8		x,de,sigp,sum
 
C  Declarations of arrays:
	real*8		e(ns),spec(ns)
	real*8		gs(E_SZ),spec2(E_SZ)

C  Labeled constants:
	real*8		ZERO
	parameter	(ZERO=0.0D0)
	integer*4	NSIGMA,NPPSIGMA
	parameter	(NSIGMA=3,NPPSIGMA=6)

	ierror = 0
	nsigma2 = 2*NSIGMA

	call peak(e,spec,ns,ep,sp,ierrorp)
	if (ierrorp .eq. -1) goto 925

c  generate Gaussian with correct sigma in units of x-axis
	de   = e(2) -e(1)			! step-size in eV
	sigp = 2.0d0*sige*ep/de			! sigma in x-axis units
	if (sigp .lt. NPPSIGMA -1) goto 930
	np   = dble(nsigma2)*sigp +1.0d0
	if (np/2*2 .eq. np) np = np +1		! make odd
	sum  = ZERO
	do ip=1,np
	  x      = ((ip-1)/dble(np-1) -0.5d0)*dble(np-1)
	  gs(ip) = exp(-x**2/2.0d0/sigp**2)
	  sum    = sum +gs(ip)
	enddo

c  make convolution
	ne1 = np/2 +1
	ne2 = ns +1 -ne1
	do ie=ne1,ne2
	  spec2(ie) = ZERO
	  do ip=1,np
	    spec2(ie) = spec2(ie) +spec(ie+ne1-ip)*gs(ip)
	  enddo ! ip
	enddo   ! ie

c  return in original array and make adjustment of array sizes
	ns = ne2 -ne1 +1
	do ie=ne1,ne2
	  e   (ie-ne1+1) = e(ie)
	  spec(ie-ne1+1) = spec2(ie)/sum
	enddo   ! ie
	return

200	format(' ',8a)
213	format(' ',a,f12.7,a)

C+ Error returns
C  -----------------------------------------------------------------------------
925	continue
	print 200,'&gauss_convolve-F-PEAKERR, Peak error'
	ierror = -1
	goto 999

930	continue
	print 200,'&gauss_convolve-F-SIGERR, Sigma error'
	print 200,'- too few data points for Gaussian convolution'
	print 213,'- check input data file; sige ',sige
	ierror = -1
	goto 999

999	continue
	return
	end 	! gauss_convolve
