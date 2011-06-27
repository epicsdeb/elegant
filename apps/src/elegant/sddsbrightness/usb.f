	SUBROUTINE USB(ENERGY,CUR,SIGX,SIGY,SIGX1,SIGY1,PERIOD,N_F,
     1              KX_F,KY_F,EMINU,EMAXU,NEU_F,METHOD_F,
     2              E_F,SPEC0_F,NE_F,IERROR)
C+
C FUNCTIONAL DESCRIPTION:	
C  Subroutine to calculate the on-axis brilliance (ph/s/mrad^2/mm^2/0.1%bw)
C  using the Bessel function approximation for an ideal planar undulator or an
C  ideal helical undulator.
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
C FORMAL PARAMETERS:
C  The input parameters are divided into sections related to the storage ring,
C  the undulator device, and the brilliance calculation.
C  Input arguments:
C Machine Parameters:
C  ENERGY		Storage ring energy              (GeV)
C  CUR		        Storage ring current             (mA)
C  SIGX		        RMS beam size (horizontal)       (mm)
C  SIGY		        RMS beam size (vertical)         (mm)
C  SIGX1	        RMS beam divergence (horizontal) (mrad)
C  SIGY1	        RMS beam divergence (vertical)   (mrad)
C Undulator Parameters:
C  PERIOD	        Undulator period length          (cm)
C  N_F			Number of periods        
C  KX_F			Deflection parameter (hor.  field)
C  KY_F			Deflection parameter (vert. field)
C Scan Parameters:
C  EMINU		Lower energy limit               (eV)
C  EMAXU                Upper energy limit               (eV)
C  NEU_F                Number of energy points
C Brilliance Parameters:
C  METHOD_F             Method
C  			METHOD_F={0,1}Non-zero emittance;
C                                     infinite-N+convolution (Dejus' approach)
C  			METHOD_F=2    Non-zero emittance;
C                                     infinite-N+convolution (Walker's approach)
C  			METHOD_F=3    Non-zero emittance; finite-N (Walker's)
C  Output arguments:
C  E_F                  Energy array                     (eV)
C  SPEC0_F              On-axix brilliance array       (ph/s/mrad^2/mm^2/0.1%bw)
C  NE_F		        Number of points in the arrays
C                       METHOD_F={0,1} => NE_F=NEU_F+1
C                       METHOD_F=2     => NE_F~NEU_F
C                       METHOD_F=3     => NE_F=NEU_F+1
C  IERROR	        Error flag; normal return  => ierror= 0
C                                   error occurred => ierror=-1
C  
C COMMON BLOCKS:
C  Several labeled common blocks are used.
C  
C DESIGN ISSUES:
C  At the APS, the notation on-axis brilliance is the same as on-axis
C  brightness at many other synchrotron radiation sites.
C  The routine is based on the Bessel function approximation as programmed in 
C  the codes US and URGENT, and is valid in the far-field for an ideal sinus-
C  oidal magnetic field. Three methods (Dejus' and Walker's) can be chosen.
C  The Dejus-method uses an "internal" variable energy-mesh and returns 
C  an energy array that is exactly the array selected by the user.  
C  The Walker-method uses a constant energy-mesh that is based on the 
C  infinite-N convolution and which is not exactly the same as the array as
C  predicted from the user's selection. The third method is also by Walker, 
C  finite-N, and this method is considerably slower than the other two methods 
C  based on infinite-N. Generally, use the infinite-N method with convolution.
C  All harmonics are selected to be contribute at any given energy (IHARM=0).
C  Several "intrinsic parameters" are used, see the setting of the "default 
C  values" in the beginning of the code. The normalized Stokes parameters are
C  calculated including the unpolarized component but are not returned to the
C  caller.
C  
C COPYRIGHT:
C  This routine must only be used at The Advanced Photon Source and must not
C  be tranferred or used at any other location without the written consent
C  of the author.
C  
C FILES USED:
C  None.
C
C KEYWORDS:
C  Undulator Spectrum, Bessel Function Approximation, On-axis Brightness,
C  On-axis Brilliance.
C  
C LINK/LIBRARY ISSUES:
C  Calls routines BRIGHTE and HUNT.  BRIGHTE calculates the brightness and HUNT
C  searches an array of real numbers (from Numerical Recipes).
C  
C PORTABILITY ISSUES:
C  Runs on DEC 3000/400 AXP alpha (Unix v. 3.2c), SUNs (SUN-OS 4.1.3), 
C  HP 9000/735-series (HP-UX 9.03).
C  
C TIMING:
C  The execution is fast (less than a second) for the infinite-N methods, and
C  on the order of 10 times larger for the finite-N method.
C
C VERSION:
C  1.2
C  
C MODIFICATION HISTORY:
C	 Date     | Name  | Description
C ----------------+-------+-----------------------------------------------------
C 07-MAY-1996     | RJD   | First version.
C ----------------+-------+-----------------------------------------------------
C 17-MAY-1996     | RJD   | Added descriptive text for the input/output.
C ----------------+-------+-----------------------------------------------------
C-
C  Size parameters:
	INTEGER*4	E_SZ,A_SZ,B_SZ
	PARAMETER	(E_SZ=40001,A_SZ=400,B_SZ=A_SZ/4+1)
  
C  Declarations of scalars:
	INTEGER*4	IERROR
	INTEGER*4	N_F,NEU_F,METHOD_F,NE_F
	INTEGER*4	NOMEGA,NEI
	INTEGER*4	I,IE,ID,IW,ISIGN
	INTEGER*4	IE1,IE2,IDW
	INTEGER*4	IP,IP1,IP2

	REAL*8          ENERGY,CUR,PERIOD
	REAL*8		KX_F,KY_F
	REAL*8		SIGX,SIGY,SIGX1,SIGY1
	REAL*8		EMIN,EMAX
	REAL*8		EMINU,EMAXU
	REAL*8		COMEGA
	REAL*8		G2,LAMDAR,LAMDA1,E1Z
	REAL*8		ARG,D2,SL
	REAL*8		SIGU2,SIGV2,SIGU,SIGV
	REAL*8		XE,YE,DEU
	REAL*8		E1MIN,E1MAX,EP,EF,DEF
	REAL*8          K2

C  Declarations of arrays:
	REAL*8		E_F(*),SPEC0_F(*)

C  Fundamental physical constants; Physics Today Aug. 1990:
	REAL*8		C,ME,MEE,EC,H,HBAR,MUZ,EPSZ
	PARAMETER	(C    =2.99792458D8)	! Speed of light [m/s]
	PARAMETER	(ME   =9.1093897D-31)	! Electron rest mass [kg]
	PARAMETER	(MEE  =0.51099906D0)	! Electron rest mass [MeV]
	PARAMETER	(EC   =1.60217733D-19)	! Elementary charge [C]
	PARAMETER	(H    =6.6260755D-34)	! Planck's constant [Js]
	PARAMETER	(HBAR =1.05457266D-34)	! Planck's constant/2Pi [Js]
	PARAMETER	(MUZ  =1.2566370614D-6)	! Permeability of vacuum [NA-2]
	PARAMETER	(EPSZ =8.854187817D-12)	! Permittivity of vacuum [Fm-1]

C  Conversion factors:
	REAL*8		C_EVANG,C_MM_M,C_M_MM,C_CM_ANG,C_MRAD_RAD,C_MA_A
	REAL*8		C_CM_M,C_RAD_MRAD
	PARAMETER	(C_EVANG=H*C/EC*1.0D10,C_MM_M=1.0D-3,C_M_MM=1.0D3)
	PARAMETER	(C_CM_ANG=1.0D8,C_MRAD_RAD=1.0D-3,C_MA_A=1.0D-3)
	PARAMETER	(C_CM_M=1.0D-2,C_RAD_MRAD=1.0D+3)

C  Labeled constants:
	CHARACTER*(*)	VN
	PARAMETER       (VN='v. 1.2')
	INTEGER*4	IDWMAX
	PARAMETER	(IDWMAX=128)
	REAL*8		PI,PIHALF,TWOPI
	PARAMETER	(PI    =3.1415 92653 58979 32384 62643D0)
	PARAMETER	(PIHALF=1.5707 96326 79489 66192 31322D0)
	PARAMETER	(TWOPI= 6.2831 85307 17958 64769 25287D0)
	REAL*8		FINE_STRUCTURE_CONST
	PARAMETER	(FINE_STRUCTURE_CONST=1.0D19*EC*1.0D19*EC
     1		/(4.0D0*PI*EPSZ*1.0D38*HBAR*C))
	REAL*8		PTOT_FAC,PD_FAC,PDH_FAC
	PARAMETER	(PTOT_FAC=PI/3.0D0*EC/EPSZ/(MEE**2)*1.0D+6)   ! 0.07257
	PARAMETER	(PD_FAC  =21.0D0/(16.0D0*PI*MEE**2)*PTOT_FAC) ! 0.11611
	PARAMETER	(PDH_FAC =EC/EPSZ/(C_RAD_MRAD**2)) ! 1.80951e-14
	REAL*8		BW
	PARAMETER	(BW=1.0D-3) ! 0.1%
	REAL*8		ZERO,ONE,TWO,HALF,EPS,EPSE,EPSK
	PARAMETER	(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,HALF=0.5D0)
	PARAMETER	(EPS=1.0D-6,EPSE=1.0D-8,EPSK=1.0D-4)

C  Data:

C  Common blocks:
	INTEGER*4	MODE,METHOD,IHARM,N,NSIGMA
	INTEGER*4	NPHI4,NPHI,NPHI_BRIGHT
	INTEGER*4	NALPHA
	INTEGER*4	NE,NE1,NE2,NEU
	INTEGER*4	NW
	REAL*8		D
	REAL*8		K,KX,KY,K3,NPI,PE,GAMMA,ER,EW,DEW,LEN
	REAL*8		FAC,C1
	REAL*8		AP2MIN,AP2MAX,AP2CNT,ARGMAX,FU,FV,SIGX2,SIGY2
	REAL*8		DPHI
	REAL*8		CALPHA2
	REAL*8		DE
	REAL*8		DW
	INTEGER*4	INDEX_PHI(A_SZ)
	INTEGER*4	I1(E_SZ),I2(E_SZ)
	REAL*8		COSPHI(A_SZ),SINPHI(A_SZ),S2SIGN(A_SZ)
	REAL*8		E(E_SZ),EU(E_SZ)
	REAL*8		SPEC0(E_SZ),SPEC1(E_SZ),SPEC2(E_SZ),SPEC3(E_SZ)
	REAL*8		HE(E_SZ),HA(E_SZ)
	REAL*8		HW(E_SZ)
	REAL*8		BR0(B_SZ,B_SZ),BR1(B_SZ,B_SZ)
	REAL*8		BR2(B_SZ,B_SZ),BR3(B_SZ,B_SZ)
	REAL*8		RA0,RA1,RA2,RA3

	COMMON		/CALCI/     MODE,METHOD,IHARM,N,NSIGMA
	COMMON		/CALCR/	    K,KX,KY,K3,NPI,PE,GAMMA,ER,EW,DEW,LEN
	COMMON		/FACTOR/    FAC,C1
	COMMON		/BEAM/	    AP2MIN,AP2MAX,AP2CNT,ARGMAX,
     1			    FU,FV,SIGX2,SIGY2
	COMMON		/ANGLE_PHII/ INDEX_PHI,NPHI4,NPHI,NPHI_BRIGHT
	COMMON		/ANGLE_PHIR/ COSPHI,SINPHI,S2SIGN,DPHI
	COMMON		/ANGLE_ALPI/ NALPHA
	COMMON		/ANGLE_ALPR/ CALPHA2
	COMMON		/ENERGY1/    E,EU,SPEC0,SPEC1,SPEC2,SPEC3
	COMMON		/ENERGY2/    I1,I2,NE,NE1,NE2,NEU,DE
	COMMON		/STEP/	    HE,HA
	COMMON		/LINE_SHAPEI/ NW
	COMMON		/LINE_SHAPER/ HW,DW
	
	COMMON		/SPECTRA/   BR0,BR1,BR2,BR3,RA0,RA1,RA2,RA3

C+ First executable statement here
C  -----------------------------------------------------------------------------
	IERROR = 0

C  Copy formal parameters (that are in common blocks) to new variables
	N   = N_F	! number of periods
	KX  = KX_F	! deflection parameter
	KY  = KY_F	! deflection parameter
	NEU = NEU_F	! number of energy points
	
C+ Check for valid input
C  -----------------------------------------------------------------------------
	IF (METHOD_F .EQ. 0) THEN 	! default
	  METHOD = 4 			! Dejus'
	ELSE IF (METHOD_F .EQ. 1) THEN
	  METHOD = 4 			! Dejus'
	ELSE IF (METHOD_F .EQ. 2) THEN
	  METHOD = 14 			! Walker's
	ELSE IF (METHOD_F .EQ. 3) THEN
	  METHOD = 1 			! Walker's finite-N
	ELSE
	  GO TO 930
	END IF ! METHOD_F
	IF ((ABS(KX-KY) .GE. EPSK) 
     1    .AND. (KX .GE. EPSK .AND. KY .GE. EPSK)) GO TO 902
C-
C  -----------------------------------------------------------------------------

C+ Default values
C  -----------------------------------------------------------------------------
	NPHI    = 20
	NALPHA  = 15
	NOMEGA  = 64
	CALPHA2 = 2.0D0
	COMEGA  = 8.0D0
	NSIGMA  = 3
	IHARM   = 0
	D       = ONE
C-
C  -----------------------------------------------------------------------------

C+ Definition of constants
C  -----------------------------------------------------------------------------
	GAMMA  = ENERGY/MEE*1.0D3
	G2     = GAMMA*GAMMA
	LAMDAR = PERIOD*C_CM_ANG/(TWO*G2) ! Reduced wavelength [A]
	ER     = C_EVANG/LAMDAR		  ! Reduced energy    [eV]
	K2     = KX*KX +KY*KY
	K3     = ONE+K2/TWO
	LAMDA1 = LAMDAR*K3		  ! First harmonic on axis  [A]
	E1Z    = ER    /K3		  ! First harmonic on axis [eV]
	D2     = D*D			  ! Distance squared     [m**2]
	LEN    = N*PERIOD*C_CM_M	  ! Length of device        [m]
	NPI    = N*PI
C-
C  -----------------------------------------------------------------------------

C+ Beam emittance
C  -----------------------------------------------------------------------------
	FU    = ZERO
	FV    = ZERO
	SIGX2 = SIGX*SIGX
	SIGY2 = SIGY*SIGY
	SIGU2 = SIGX1*SIGX1*C_MRAD_RAD*C_MRAD_RAD
	SIGV2 = SIGY1*SIGY1*C_MRAD_RAD*C_MRAD_RAD
	SIGU  = SQRT(SIGU2) ! [rad]
	SIGV  = SQRT(SIGV2) ! [rad]
	IF (SIGU2 .NE. ZERO) FU = 0.5D0/SIGU2
	IF (SIGV2 .NE. ZERO) FV = 0.5D0/SIGV2
C-
C  -----------------------------------------------------------------------------

C+ Determine min and max emission angles and center
C  -----------------------------------------------------------------------------
	AP2MIN = ZERO
	AP2CNT = ZERO

	XE  = +NSIGMA*SIGU 		! Cartesian angle in x-dir.
	YE  = +NSIGMA*SIGV 		! Cartesian angle in y-dir.
	AP2MAX = G2*(XE*XE+YE*YE)
	ARGMAX = NSIGMA*NSIGMA/TWO
C-
C  -----------------------------------------------------------------------------

C+ Define energy scale and set up array for the line shape function
C  -----------------------------------------------------------------------------
	IF (METHOD .EQ. 4) THEN			! Dejus' method
	  E1MIN = E1Z    *K3/(K3+AP2MAX)
	  E1MAX = E1Z    *K3/(K3+AP2MIN)
	  DEW   = (E1Z/N)*K3/(K3+AP2CNT)
	  EW    = COMEGA*DEW
	  PE    = PI/DEW

C  Extend "internal" energy scale to make room for the line shape function
	  DW   = TWO*EW/NOMEGA ! Default step size
	  EMIN = EMINU-EW
	  EMAX = EMAXU+EW

C  Adjust EMIN and EMAX
	  I  = 1
	  EP = I*E1MAX
	  DO WHILE (EP .LT. EMIN)
	    I  = I+1
	    EP = I*E1MAX
	  END DO ! WHILE
	  IE1  = I
	  EMIN = MAX(IE1*E1MIN,EMIN)

	  EP = IE1*E1MIN
	  DO WHILE (EP .LT. EMAX)
	    I  = I+1
	    EP = I*E1MIN
	  END DO ! WHILE
	  IE2  = I-1
	  EMAX = MIN(IE2*E1MAX,EMAX)
	    
	  IF (EMAX .LE. EMIN) GO TO 910

	  I   = IE1
	  EP  = I*E1MAX
	  IE  = 1
	  E (IE) = EMIN

C  Generate energy scale with variable step size
	  DO WHILE (E(IE) .LT. EMAX)
	    IDW = 1
	    DEF = DW /IDW
	    EF  = DEW/IDW
	    DO WHILE (E(IE) .GT. EP-EF)
	      IDW = IDW*2
	      IF (IDW .GT. IDWMAX) THEN
	        EF = ZERO
	      ELSE
		DEF = DW /IDW
		EF  = DEW/IDW
	      END IF ! IDW
	    END DO ! WHILE

	    DO WHILE (E(IE) .LT. (EP-EPSE) .AND. E(IE) .LT. (EMAX-EPSE))
	      IF (E(IE) .GT. EP-EF) THEN
		IDW = IDW*2
		IF (IDW .GT. IDWMAX) THEN
		  EF = ZERO
		ELSE
		  DEF = DW /IDW
		  EF  = DEW/IDW
		END IF ! IDW
	      END IF ! E(IE)

	      IE = IE+1
	      IF (IE .GT. E_SZ) GO TO 900
	      HE(IE-1) = DEF
	      IF (IE .EQ. 2) THEN
	        HA(IE-1) = HE(IE-1)/TWO
	      ELSE
	        HA(IE-1) = (HE(IE-2)+HE(IE-1))/TWO ! Average
	      END IF ! IE
	      E(IE) = E(IE-1)+HE(IE-1)
	    END DO ! WHILE

C  Redefine the last point within range of harmonic # i
	    SL       = MIN(EP,EMAX)
	    E (IE)   = SL-EPSE
	    HE(IE-1) = E(IE)-E(IE-1)
	    HA(IE-1) = (HE(IE-2)+HE(IE-1))/TWO ! Average
	    IE       = IE+1
	    IF (IE .GT. E_SZ) GO TO 900
	    HE(IE-1) = TWO*EPSE
	    HA(IE-1) = (HE(IE-2)+HE(IE-1))/TWO ! Average
	    E (IE)   = E(IE-1)+HE(IE-1)

	    I   = I+1 ! Next harmonic
	    DEF = I*E1MIN-E(IE)
	    IF (DEF .GT. ZERO) THEN
	      IE     = IE+1
	      IF (IE .GT. E_SZ) GO TO 900
	      HE(IE-1) = DEF
	      HA(IE-1) = (HE(IE-2)+HE(IE-1))/TWO ! Average
	      E (IE)   = E(IE-1)+HE(IE-1)
	    END IF ! DEF

	    EP  = I*E1MAX
		
	  END DO ! WHILE

C  Adjust end points
	  IF (DEF .GT. ZERO) THEN
	    NE = IE-1
	  ELSE
	    NE = IE
	  END IF ! DEF
	  NE1 = 1
	  NE2 = NE
	  HE(NE) = ZERO
	  HA(NE) = HE(NE-1)/TWO
	
	  IF (NEU .EQ. 0) NEU = INT(NE2/100.0D0+ONE)*100.0D0 ! Default
C	  DEU = (EMAXU-EMINU)/NEU
	  DE  = (EMAXU-EMINU)/NEU
	  NEU = NEU+1
	  DO IE=1,NEU ! Energy scale (user's selection)
C	    EU(IE) = EMINU+(IE-1)*DEU
	    EU(IE) = EMINU+(IE-1)*DE
	  END DO ! IE

	ELSE IF (METHOD .EQ. 14) THEN 		! Walker's method
	  DEW = (E1Z/N)*K3/(K3+AP2CNT)
	  EW  = COMEGA*DEW ! [eV]
	  IF (NEU .GT. 0) THEN ! Redefine NOMEGA
	    SL     = (EMAXU-EMINU)/NEU
	    NOMEGA = TWO*EW/SL+ONE
	    NOMEGA = NOMEGA/2*2 ! Make even
	    IF (NOMEGA .LT. 16) THEN ! Reset to 16
	      NOMEGA = 16
	      PRINT *
	      PRINT *  ,'&USB-W-ISMALL, NOMEGA less than 16'
	      PRINT *  ,'- NOMEGA reset to 16'
	      PRINT *
	    END IF ! NOMEGA
	  END IF ! NEU
	    
	  DW   = TWO*COMEGA/NOMEGA
	  DE   = DW*DEW
	  EMIN = EMINU-EW		 ! Extend lower limit of energy scale
	  IF (EMIN .LT. EPS) EMIN = EPS
	  EMAX = EMAXU+EW		 ! Extend upper limit of energy scale
	  NE   = (EMAX-EMIN)/DE+ONE      ! Redefine the # of intervals
	  NE   = NE+1			 ! # of points
	  NE1  = NOMEGA/2+1		 ! Smallest index; odd
	  NE2  = NE+1-NE1		 ! Largest index
	  NEI  = NE2-NE1		 ! # of intervals

C  Generate line shape function (used in energy convolution)
	  NW = NOMEGA+1 ! # of points
	  DO IW=1,NW
	    ARG = (-COMEGA+(IW-1)*DW)*PI	    
	    IF (ABS(ARG) .GT. EPS) THEN
	      SL = SIN(ARG)/ARG
	      HW(IW) = SL*SL
	    ELSE
	      HW(IW) = ONE
	    END IF ! ARG
	  END DO ! IW
C  Generate energy scale
	  IF (NE .GT. E_SZ) GO TO 905
	  DO IE=1,NE
	    E(IE) = EMIN+(IE-1)*DE ! [eV]
	  END DO ! IE

	ELSE
	  NE  = NEU	    
	  IF (NE .GT. 0) DE  = (EMAXU-EMINU)/NE
	  NE  = NE+1		 	! # of points
	  NE1 = 1			! Smallest index
	  NE2 = NE+1-NE1		! Largest index
	  NEI = NE2-NE1		 	! # of intervals
C  Generate energy scale
	  IF (NE .GT. E_SZ) GO TO 905
	  DO IE=1,NE
	    E(IE) = EMINU+(IE-1)*DE ! [eV]
	  END DO ! IE
	END IF ! METHOD
C-
C  -----------------------------------------------------------------------------

C+ Set up arrays of cos(phi) and sin(phi), and index_phi
C  -----------------------------------------------------------------------------
	NPHI4 = NPHI   	    ! Use only 0 ... 90 deg.
	NPHI_BRIGHT = NPHI  ! # of indices for the brightness array (use sym.)
	DPHI   = PIHALF/NPHI
	ISIGN  = +1
	SL     = ISIGN*DPHI/TWO
	INDEX_PHI(1) = 1
	S2SIGN(1) = +ONE
	COSPHI(1) = COS(SL)
	SINPHI(1) = SIN(SL)
	DO ID=2,NPHI4
	  ARG           = SL+(ID-1)*DPHI
	  INDEX_PHI(ID) = INDEX_PHI(ID-1)+ISIGN
	  COSPHI(ID)    = COS(ARG)
	  SINPHI(ID)    = SIN(ARG)
	  IF (ID .EQ. NPHI) THEN
	    ISIGN = 0
	  ELSE IF (ID .EQ. NPHI+1) THEN
	    ISIGN = -1
	  ELSE IF (ID .EQ. 2*NPHI) THEN
	    ISIGN = 0
	  ELSE IF (ID .EQ. 2*NPHI+1) THEN
	    ISIGN = +1
	  ELSE IF (ID .EQ. 3*NPHI) THEN
	    ISIGN = 0
	  ELSE IF (ID .EQ. 3*NPHI+1) THEN
	    ISIGN = -1
	  END IF ! ID
	  IF (ISIGN .NE. 0) THEN
	    S2SIGN(ID) = ISIGN
	  ELSE ! ISIGN = 0
	    S2SIGN(ID) = +S2SIGN(ID-1)
	  END IF
	END DO ! ID
C-
C  -----------------------------------------------------------------------------

C+ Scale factors
C  -----------------------------------------------------------------------------
	FAC = 4.0D0			! use symmetry
	IF (SIGU .NE. ZERO .AND. SIGV .NE. ZERO)
     1   C1 =N*N*FINE_STRUCTURE_CONST*BW*CUR*C_MA_A/EC
     2        /(TWOPI*SIGU*SIGV*D2*C_M_MM*C_M_MM)
C-
C  -----------------------------------------------------------------------------

C+ Call spectral routine
C  -----------------------------------------------------------------------------
	CALL SPECTRAL_DISTRIBUTION(IERROR)
	IF (IERROR .NE. 0) RETURN
C-
C  -----------------------------------------------------------------------------

C+ Return results
C  -----------------------------------------------------------------------------
	IF (METHOD .EQ. 4) THEN
	  IP1 = 1
	  IP2 = NEU
	ELSE 				! METHOD={1,14}
	  IP1 = NE1
	  IP2 = NE2
	END IF ! METHOD

	NE_F = IP2-IP1+1
	DO IP=1,NE_F
	  IF (METHOD .EQ. 4) THEN
	    E_F(IP)     = EU(IP)
	    SPEC0_F(IP) = SPEC0(IP)
	  ELSE 				! METHOD={1,14}
	    E_F(IP)     = E(NE1-1+IP)
	    SPEC0_F(IP) = SPEC0(NE1-1+IP)
	  END IF ! METHOD
	END DO ! IP   
C-
C  -----------------------------------------------------------------------------
	RETURN

200	FORMAT(' ',8A)
205	FORMAT(' ',A,I6,A,I6)
207	FORMAT(' ',A,I2,A)
210	FORMAT(' ',A,2(F10.3,A))

C+ Error returns
C  -----------------------------------------------------------------------------
900	CONTINUE
	PRINT 200,'&USB-F-BNDERR, Boundary error'
	PRINT 205,'- energy array out of bounds; number of points',IE,
     1	  ' is greater than ',E_SZ
	GO TO 999

902	CONTINUE
	PRINT 200,'&USB-E-INVDAT, Deflection parameter errror'
	PRINT 200,'- deflection parameters Kx, Ky must be set equal or',
     1         ' one must be identically zero.'
	GO TO 999

905	CONTINUE
	PRINT 200,'&USB-F-BNDERR, Boundary error'
	PRINT 205,'- energy array out of bounds; number of points ',NE,
     1	  ' is greater than ',E_SZ
	GO TO 999

910	CONTINUE
	PRINT 200,'&USB-E-HARMERR, Harmonic errror'
	PRINT 210,'- no harmonics reachable in the range ',EMINU,
     1	  ' to ',EMAXU,' eV.'
	GO TO 999

930	CONTINUE
	PRINT 200,'&USB-E-INVDAT, Invalid data'
	PRINT 207,'- check input data file; method ',METHOD_F,
     1	  ' not valid for the on-axis brightness.'
	GO TO 999

999	CONTINUE
	PRINT *
	PRINT 200,'&USB-F-PRGERR, Program error'
	PRINT 200,'- unsuccessful completion due to error status'
	PRINT *

	WRITE (2,*)
	WRITE (2,200) '&USB-F-PRGERR, Program error'
	WRITE (2,200) '- unsuccessful completion due to error status'
	WRITE (2,*)
C-
C  -----------------------------------------------------------------------------

	IERROR = -1
	RETURN
	END ! USB

	SUBROUTINE SPECTRAL_DISTRIBUTION(IERROR)
 
C[subroutine_header_comments]

C  Size parameters:
	INTEGER*4	E_SZ,A_SZ,B_SZ
	PARAMETER	(E_SZ=40001,A_SZ=400,B_SZ=A_SZ/4+1)

C  Declarations of scalars:
	LOGICAL*4	LE1,LE2
	INTEGER*4	IERROR,ICOUNT
	INTEGER*4	I,IC,IE,IMIN,IMAX,IH1,IH2

	REAL*8		R,DA2,DALPHA,SL,CONST
	REAL*8		ALPHAI,ALPHA2I,ALPHAMIN,ALPHA2MIN,ALPHAMAX,ALPHA2MAX
	REAL*8		DEL,SIGR2

C  Declarations of arrays:
	REAL*8		ALPHA(B_SZ),THETA(B_SZ)

C  Fundamental physical constants; Physics Today Aug. 1990:
	REAL*8		C,ME,MEE,EC,H,HBAR,MUZ,EPSZ
	PARAMETER	(C    =2.99792458D8)	! Speed of light [m/s]
	PARAMETER	(ME   =9.1093897D-31)	! Electron rest mass [kg]
	PARAMETER	(MEE  =0.51099906D0)	! Electron rest mass [MeV]
	PARAMETER	(EC   =1.60217733D-19)	! Elementary charge [C]
	PARAMETER	(H    =6.6260755D-34)	! Planck's constant [Js]
	PARAMETER	(HBAR =1.05457266D-34)	! Planck's constant/2Pi [Js]
	PARAMETER	(MUZ  =1.2566370614D-6)	! Permeability of vacuum [NA-2]
	PARAMETER	(EPSZ =8.854187817D-12)	! Permittivity of vacuum [Fm-1]

C  Conversion factors:
	REAL*8		C_EVANG,C_ANG_M,C_M2_MM2
	PARAMETER	(C_EVANG=H*C/EC*1.0D10,C_ANG_M=1.0D-10) ! 12398.42
	PARAMETER	(C_M2_MM2=1.0D6)

C  Labeled constants:
	REAL*8		PI,PIHALF,TWOPI
	PARAMETER	(PI    =3.1415 92653 58979 32384 62643D0)
	PARAMETER	(PIHALF=1.5707 96326 79489 66192 31322D0)
	PARAMETER	(TWOPI= 6.2831 85307 17958 64769 25287D0)
	REAL*8		ZERO,ONE,TWO,EPS
	PARAMETER	(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,EPS=1.0D-2)
	REAL*8		CV
	PARAMETER	(CV=ONE/8.0D0/PI/PI*C_ANG_M*C_M2_MM2) ! 1.2665D-6

C  Common blocks:
	INTEGER*4	MODE,METHOD,IHARM,N,NSIGMA
	INTEGER*4	NALPHA
	INTEGER*4	NE,NE1,NE2,NEU
	REAL*8		K,KX,KY,K3,NPI,PE,GAMMA,ER,EW,DEW,LEN
	REAL*8		FAC,C1
	REAL*8		AP2MIN,AP2MAX,AP2CNT,ARGMAX,FU,FV,SIGX2,SIGY2
	REAL*8		CALPHA2
	REAL*8		DE
	INTEGER*4	I1(E_SZ),I2(E_SZ)
	REAL*8		E(E_SZ),EU(E_SZ)
	REAL*8		SPEC0(E_SZ),SPEC1(E_SZ),SPEC2(E_SZ),SPEC3(E_SZ)
	REAL*8		BR0(B_SZ,B_SZ),BR1(B_SZ,B_SZ)
	REAL*8		BR2(B_SZ,B_SZ),BR3(B_SZ,B_SZ)
	REAL*8		RA0,RA1,RA2,RA3

	COMMON		/CALCI/     MODE,METHOD,IHARM,N,NSIGMA
	COMMON		/CALCR/	    K,KX,KY,K3,NPI,PE,GAMMA,ER,EW,DEW,LEN
	COMMON		/FACTOR/    FAC,C1
	COMMON		/BEAM/	    AP2MIN,AP2MAX,AP2CNT,ARGMAX,
     1			    FU,FV,SIGX2,SIGY2
	COMMON		/ANGLE_ALPI/ NALPHA
	COMMON		/ANGLE_ALPR/ CALPHA2
	COMMON		/ENERGY1/    E,EU,SPEC0,SPEC1,SPEC2,SPEC3
	COMMON		/ENERGY2/    I1,I2,NE,NE1,NE2,NEU,DE
	COMMON		/SPECTRA/   BR0,BR1,BR2,BR3,RA0,RA1,RA2,RA3

	IERROR = 0
	LE1 = .FALSE.
	LE2 = .FALSE.

C+ Loop over energies
C  -----------------------------------------------------------------------------
	DO IE=1,NE

	  SPEC0(IE) = ZERO
	  SPEC1(IE) = ZERO
	  SPEC2(IE) = ZERO
	  SPEC3(IE) = ZERO
	  I1(IE)    = 0
	  I2(IE)    = 0
	  CONST     = C1

	  R         = ER/E(IE)

C+ Find range of harmonics that contributes at energy E(IE)
C  -----------------------------------------------------------------------------
	  IF (METHOD .EQ. 1) THEN ! Finite-N
	    DA2  = CALPHA2*R/N
	    IMIN = (AP2MIN+K3-DA2)/R+ONE
	    IMAX = (AP2MAX+K3+DA2)/R
	  ELSE ! Infinite-N
	    IMIN = (AP2MIN+K3)/R+ONE
	    IMAX = (AP2MAX+K3)/R
	  END IF ! METHOD
	  IF (IMAX .LT. IMIN) GO TO 810 ! Next energy
	  LE1 = .TRUE.		
	  IF (IHARM .GT. 0 .AND. (IHARM .LT. IMIN .OR. IHARM .GT. IMAX))
     1     GO TO 810 ! Next energy
	  LE2 = .TRUE.
C-
C  -----------------------------------------------------------------------------
	  RA0 = ZERO
	  RA1 = ZERO
	  RA2 = ZERO
	  RA3 = ZERO

	  IF (IHARM .GT. 0) THEN
	    IH1 = IHARM
	    IH2 = IH1
	  ELSE IF (IHARM .LT. 0) THEN ! Lowest order
	    IH1 = IMIN
	    IH2 = IH1
	  ELSE ! IHARM = 0 ! All contributing harmonics
	    IH1 = IMIN
	    IH2 = IMAX
	  END IF ! IHARM
	  I1(IE)  = IH1
	  I2(IE)  = IH2
	  ICOUNT = 0
C+ Loop over harmonics
C  -----------------------------------------------------------------------------
	  DO I=IH1,IH2
	    ICOUNT  = ICOUNT+1
	    ALPHA2I = R*I-K3

	    IF (METHOD .EQ. 1) THEN ! Finite-N
	      IF (ALPHA2I .GT. ZERO) THEN
	        ALPHAI = SQRT(ALPHA2I)
	      ELSE
	        ALPHAI = ZERO
	      END IF ! ALPHA2I
	      ALPHA2MIN = ALPHA2I-DA2
	      IF (ALPHA2MIN .LT. ZERO) ALPHA2MIN = ZERO
	      ALPHAMIN  = SQRT(ALPHA2MIN)
	      ALPHA2MAX = ALPHA2I+DA2
	      ALPHAMAX  = SQRT(ALPHA2MAX)
	      DALPHA    = (ALPHAMAX-ALPHAMIN)/NALPHA
	      SL = ALPHAMIN+DALPHA/TWO
	      DO IC=1,NALPHA
	        ALPHA(IC) = SL+(IC-1)*DALPHA
	        THETA(IC) = ALPHA(IC)/GAMMA
	      END DO ! IC
	    ELSE ! Infinite-N
	      ALPHAI   = SQRT(ALPHA2I)
	      ALPHA(1) = ALPHAI
	      THETA(1) = ALPHA(1)/GAMMA
	      DALPHA   = R/(TWO*N)
	    END IF ! METHOD

	    IF (ALPHA2I .LT. ZERO) THEN ! may be neg for finite-N only
	      DEL = ZERO
	    ELSE                        ! Walker's approach (infinite-N)
	      DEL = ALPHA2I*N/R
	    END IF ! ALPHA2I

C  Estimate diffraction limited source size
	    IF (DEL .LT. 2.15D0) THEN
	      SIGR2 = (1.29D0+(1.229D0*(DEL-0.8D0)*(DEL-0.8D0)))**2
	    ELSE
	      SIGR2 = 5.81D0*DEL
	    END IF ! DEL

	    SIGR2 = SIGR2*CV*C_EVANG/E(IE)*LEN
	    CONST = C1/(TWOPI*SQRT((SIGX2+SIGR2)*(SIGY2+SIGR2)))

C  Brightness
	    CALL BRIGHTNESS_ARRAY(METHOD,I,R,ALPHAI,ALPHA2I,ALPHA)

C  Two-dimensional convolution of the brightness with the electron distribution
	    CALL CONVOLUTE_DISTRIBUTION(METHOD,CONST,ALPHA,THETA,DALPHA,
     1			        EPS,BR0,BR1,BR2,BR3,
     2			        RA0,RA1,RA2,RA3,ICOUNT)

	    IF (ICOUNT .GT. 1) THEN ! If higher harmonics do not contribute
	      I2(IE) = I
	      GO TO 800
	    END IF ! ICOUNT
	  END DO ! IH
C- Endloop harmonics
C  -----------------------------------------------------------------------------
800	  CONTINUE

C+ Save spectra
C  -----------------------------------------------------------------------------
	  SPEC0(IE) = FAC*RA0
	  SPEC1(IE) = FAC*RA1
	  SPEC2(IE) = ZERO 		! symmetry
	  SPEC3(IE) = FAC*RA3
C-
C  -----------------------------------------------------------------------------

810	  CONTINUE
	END DO ! IE
C- Endloop energy
C  -----------------------------------------------------------------------------
	
	IF (.NOT. LE1 .AND. .NOT. LE2) GO TO 900
	IF (      LE1 .AND. .NOT. LE2) GO TO 910

	IF (METHOD .EQ. 4) THEN ! Infinite-N with convolution
	  CALL CONVOLUTE_ENERGY_VSTEP(SPEC0)
	  CALL CONVOLUTE_ENERGY_VSTEP(SPEC1)
	  CALL CONVOLUTE_ENERGY_VSTEP(SPEC2)
	  CALL CONVOLUTE_ENERGY_VSTEP(SPEC3)
	ELSE IF (METHOD .EQ. 14) THEN ! Infinite-N with convolution (Walker's)
	  CALL CONVOLUTE_ENERGY_ESTEP(SPEC0)
	  CALL CONVOLUTE_ENERGY_ESTEP(SPEC1)
	  CALL CONVOLUTE_ENERGY_ESTEP(SPEC2)
	  CALL CONVOLUTE_ENERGY_ESTEP(SPEC3)
	END IF 

	RETURN

C+ Error returns
C  -----------------------------------------------------------------------------
900	CONTINUE
	PRINT 200,'&SPECTRAL_DISTRIBUTION-E-HARMERR, Harmonic errror'
	PRINT 210,'- no harmonics reachable in the range ',E(NE1),
     1	  ' to ',E(NE2),' eV.'
	GO TO 999

910	CONTINUE
	PRINT 200,'&SPECTRAL_DISTRIBUTION-E-HARMERR, Harmonic errror'
	PRINT 220,'- Harmonic number ',IHARM,' not in the range ',E(NE1),
     1	  ' to ',E(NE2),' eV.'
	GO TO 999

200	FORMAT(' ',8A)
210	FORMAT(' ',A,2(F10.3,A))
220	FORMAT(' ',A,I3,A,2(F10.3,A))

999	CONTINUE
	IERROR = -1
C-
C  -----------------------------------------------------------------------------
	RETURN
	END ! SPECTRAL_DISTRIBUTION

	SUBROUTINE BRIGHTNESS_ARRAY(ICALC,I,R,ALPHAI,ALPHA2I,ALPHA)
 
C[subroutine_header_comments]

C  Size parameters:
	INTEGER*4	E_SZ,A_SZ,B_SZ
	PARAMETER	(E_SZ=40001,A_SZ=400,B_SZ=A_SZ/4+1)

C  Declarations of scalars:
	INTEGER*4	ICALC,I,IC,ID
	REAL*8		ALPHAI,ALPHA2I,R,ARG
	REAL*8		ALPHA2,S0,S1,S2,S3,AXR,AYR,AXI,AYI
	REAL*8		H,SINC

C  Declarations of arrays:
	REAL*8		ALPHA(*)

C  Common blocks:
	INTEGER*4	NPHI4,NPHI,NPHI_BRIGHT
	INTEGER*4	NALPHA
	REAL*8		K,KX,KY,K3,NPI,PE,GAMMA,ER,EW,DEW,LEN
	REAL*8		DPHI
	REAL*8		CALPHA2
	INTEGER*4	INDEX_PHI(A_SZ)
	REAL*8		COSPHI(A_SZ),SINPHI(A_SZ),S2SIGN(A_SZ)
	REAL*8		BR0(B_SZ,B_SZ),BR1(B_SZ,B_SZ)
	REAL*8		BR2(B_SZ,B_SZ),BR3(B_SZ,B_SZ)
	REAL*8		RA0,RA1,RA2,RA3

	COMMON		/CALCR/	    K,KX,KY,K3,NPI,PE,GAMMA,ER,EW,DEW,LEN
	COMMON		/ANGLE_PHII/ INDEX_PHI,NPHI4,NPHI,NPHI_BRIGHT
	COMMON		/ANGLE_PHIR/ COSPHI,SINPHI,S2SIGN,DPHI
	COMMON		/ANGLE_ALPI/ NALPHA
	COMMON		/ANGLE_ALPR/ CALPHA2
	COMMON		/SPECTRA/   BR0,BR1,BR2,BR3,RA0,RA1,RA2,RA3

	IF (ICALC .EQ. 1) THEN ! Finite-N calculation
	  DO IC=1,NALPHA
	    ALPHA2 = ALPHA(IC)*ALPHA(IC)
	    ARG    = NPI*(ALPHA2-ALPHA2I)/R
	    H      = SINC(ARG)
	    DO ID=1,NPHI_BRIGHT
	      CALL BRIGHTE(I,ALPHA(IC),COSPHI(ID),SINPHI(ID),
     1                  S0,S1,S2,S3,AXR,AYR,AXI,AYI)
	      BR0(ID,IC) = H*S0
	      BR1(ID,IC) = H*S1
	      BR2(ID,IC) = H*S2
	      BR3(ID,IC) = H*S3
	    END DO ! ID
	  END DO ! IC
	ELSE ! Infinite-N approximation (H=1.0)
	  DO ID=1,NPHI_BRIGHT
	    CALL BRIGHTE(I,ALPHAI,COSPHI(ID),SINPHI(ID),
     1                S0,S1,S2,S3,AXR,AYR,AXI,AYI)
	    BR0(ID,1) = S0
	    BR1(ID,1) = S1
	    BR2(ID,1) = S2
	    BR3(ID,1) = S3
	  END DO ! ID
	END IF ! ICALC

	RETURN
	END ! BRIGHTNESS_ARRAY

	SUBROUTINE CONVOLUTE_DISTRIBUTION(ICALC,CONST,ALPHA,THETA,DALPHA,
     1				  EPS,BR0,BR1,BR2,BR3,
     2				  RA0,RA1,RA2,RA3,ICOUNT)
 
C[subroutine_header_comments]
 
C  Size parameters:
	INTEGER*4	A_SZ,B_SZ
	PARAMETER	(A_SZ=400,B_SZ=A_SZ/4+1)

C  Declarations of scalars:
	INTEGER*4	ICALC,ICOUNT,IC,ID
	REAL*8		CONST,DALPHA,EPS
	REAL*8		U,V,ARG,SL,P
	REAL*8		SUM0,SUM1,SUM2,SUM3
	REAL*8		DELTA0,DELTA1,DELTA2,DELTA3
 
C  Declarations of arrays:
	REAL*8		ALPHA(*),THETA(*)
	REAL*8		BR0(B_SZ,*),BR1(B_SZ,*)
	REAL*8		BR2(B_SZ,*),BR3(B_SZ,*)
	REAL*8		RA0,RA1,RA2,RA3

C  Labeled constants:
	REAL*8		ZERO
	PARAMETER	(ZERO=0.0D0)
 
C  Common blocks:
	INTEGER*4	NPHI4,NPHI,NPHI_BRIGHT
	INTEGER*4	NALPHA
	REAL*8		AP2MIN,AP2MAX,AP2CNT,ARGMAX,FU,FV,SIGX2,SIGY2
	REAL*8		DPHI
	REAL*8		CALPHA2
	INTEGER*4	INDEX_PHI(A_SZ)
	REAL*8		COSPHI(A_SZ),SINPHI(A_SZ),S2SIGN(A_SZ)

	COMMON		/BEAM/	    AP2MIN,AP2MAX,AP2CNT,ARGMAX,
     1			    FU,FV,SIGX2,SIGY2
	COMMON		/ANGLE_PHII/ INDEX_PHI,NPHI4,NPHI,NPHI_BRIGHT
	COMMON		/ANGLE_PHIR/ COSPHI,SINPHI,S2SIGN,DPHI
	COMMON		/ANGLE_ALPI/ NALPHA
	COMMON		/ANGLE_ALPR/ CALPHA2

	SL = CONST*DALPHA*DPHI
	IF (ICALC .EQ. 1) THEN ! Finite-N calculation
	  SUM0 = ZERO
	  SUM1 = ZERO
	  SUM2 = ZERO
	  SUM3 = ZERO
	  DO IC=1,NALPHA
	    DELTA0 = ZERO
	    DELTA1 = ZERO
	    DELTA2 = ZERO
	    DELTA3 = ZERO
	    DO ID=1,NPHI4
	      U   = THETA(IC)*COSPHI(ID)
	      V   = THETA(IC)*SINPHI(ID)
	      ARG = U*U*FU+V*V*FV
	      IF (ARG .LT. ARGMAX) THEN
		P      = EXP(-ARG)
		DELTA0 = DELTA0+P*BR0(INDEX_PHI(ID),IC)
		DELTA1 = DELTA1+P*BR1(INDEX_PHI(ID),IC)
		DELTA2 = DELTA2+P*BR2(INDEX_PHI(ID),IC)*S2SIGN(ID)
		DELTA3 = DELTA3+P*BR3(INDEX_PHI(ID),IC)
	      END IF ! ARG
	    END DO ! ID
	    SUM0 = SUM0+DELTA0*ALPHA(IC)
	    SUM1 = SUM1+DELTA1*ALPHA(IC)
	    SUM2 = SUM2+DELTA2*ALPHA(IC)
	    SUM3 = SUM3+DELTA3*ALPHA(IC)
	  END DO ! IC
	  SUM0       = SUM0*SL
	  SUM1       = SUM1*SL
	  SUM2       = SUM2*SL
	  SUM3       = SUM3*SL
	  RA0 = RA0 +SUM0
	  RA1 = RA1 +SUM1
	  RA2 = RA2 +SUM2
	  RA3 = RA3 +SUM3
	  IF (SUM0 .GT. EPS*RA0) ICOUNT = 0
	ELSE ! Infinite-N approximation (no alpha integ.)
	  DELTA0 = ZERO
	  DELTA1 = ZERO
	  DELTA2 = ZERO
	  DELTA3 = ZERO
	  DO ID=1,NPHI4
	    U   = THETA(1)*COSPHI(ID)
	    V   = THETA(1)*SINPHI(ID)
	    ARG = U*U*FU+V*V*FV
	    IF (ARG .LT. ARGMAX) THEN
	      P      = EXP(-ARG)
	      DELTA0 = DELTA0+P*BR0(INDEX_PHI(ID),1)
	      DELTA1 = DELTA1+P*BR1(INDEX_PHI(ID),1)
	      DELTA2 = DELTA2+P*BR2(INDEX_PHI(ID),1)*S2SIGN(ID)
	      DELTA3 = DELTA3+P*BR3(INDEX_PHI(ID),1)
	    END IF ! ARG
	  END DO ! ID
	  SUM0       = DELTA0*SL
	  SUM1       = DELTA1*SL
	  SUM2       = DELTA2*SL
	  SUM3       = DELTA3*SL
	  RA0 = RA0 +SUM0
	  RA1 = RA1 +SUM1
	  RA2 = RA2 +SUM2
	  RA3 = RA3 +SUM3
	  IF (SUM0 .GT. EPS*RA0) ICOUNT = 0
	END IF ! ICALC
	
	RETURN
	END ! CONVOLUTE_DISTRIBUTION

	SUBROUTINE CONVOLUTE_ENERGY_ESTEP(SP1)
 
C[subroutine_header_comments]
C  Size parameters:
	INTEGER*4	E_SZ
	PARAMETER	(E_SZ=40001)

C  Declarations of scalars:
	INTEGER*4	IE,IW
 
C  Declarations of arrays:
	REAL*8		SP1(*),SP2(E_SZ)

C  Labeled constants:
	REAL*8		ZERO
	PARAMETER	(ZERO=0.0D0)
 
C  Common blocks:
	INTEGER*4	NE,NE1,NE2,NEU
	INTEGER*4	NW
	REAL*8		DE
	REAL*8		DW
	INTEGER*4	I1(E_SZ),I2(E_SZ)
	REAL*8		E(E_SZ),EU(E_SZ)
	REAL*8		SPEC0(E_SZ),SPEC1(E_SZ),SPEC2(E_SZ),SPEC3(E_SZ)
	REAL*8		HW(E_SZ)

	COMMON		/ENERGY1/    E,EU,SPEC0,SPEC1,SPEC2,SPEC3
	COMMON		/ENERGY2/    I1,I2,NE,NE1,NE2,NEU,DE
	COMMON		/LINE_SHAPEI/ NW
	COMMON		/LINE_SHAPER/ HW,DW

	DO IE=NE1,NE2
	  SP2(IE) = ZERO
	  DO IW=1,NW
	    SP2(IE) = SP2(IE)+SP1(IE+NE1-IW)*HW(IW)
	  END DO ! IW
	END DO ! IE
c  Return in original array
	DO IE=NE1,NE2
	  SP1(IE) = SP2(IE)*DW
	END DO ! IE

	RETURN
	END ! CONVOLUTE_ENERGY_ESTEP

	SUBROUTINE CONVOLUTE_ENERGY_VSTEP(SP1)
 
C[subroutine_header_comments]
C  Size parameters:
	INTEGER*4	E_SZ
	PARAMETER	(E_SZ=40001)

C  Declarations of scalars:
	INTEGER*4	IE,IEU,J1,J2
	REAL*8		EUP,EU1,EU2
	REAL*8		ARG,P,H1,H2,S1,S2,SUM
	REAL*8		SINC
 
C  Declarations of arrays:
	REAL*8		SP1(*),SP2(E_SZ)

C  Labeled constants:
	REAL*8		ZERO,ONE,TWO
	PARAMETER	(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
 
C  Common blocks:
	INTEGER*4	NE,NE1,NE2,NEU
	REAL*8		K,KX,KY,K3,NPI,PE,GAMMA,ER,EW,DEW,LEN
	REAL*8		DE
	INTEGER*4	I1(E_SZ),I2(E_SZ)
	REAL*8		E(E_SZ),EU(E_SZ)
	REAL*8		SPEC0(E_SZ),SPEC1(E_SZ),SPEC2(E_SZ),SPEC3(E_SZ)
	REAL*8		HE(E_SZ),HA(E_SZ)

	COMMON		/CALCR/	    K,KX,KY,K3,NPI,PE,GAMMA,ER,EW,DEW,LEN
	COMMON		/ENERGY1/    E,EU,SPEC0,SPEC1,SPEC2,SPEC3
	COMMON		/ENERGY2/    I1,I2,NE,NE1,NE2,NEU,DE
	COMMON		/STEP/	    HE,HA

	J1 = 1
	J2 = 1
	DO IEU=1,NEU
	  EUP = EU(IEU)
	  EU1 = EUP-EW
	  EU2 = EUP+EW
	  CALL HUNT(E,NE,EU1,J1)
	  CALL HUNT(E,NE,EU2,J2)
	  IF (J2 .EQ. 0 .OR. J1 .EQ. NE) THEN
	    SP2(IEU) = ZERO

	  ELSE IF (J1 .EQ. 0) THEN
	    IF (J2 .LT. NE) THEN
C  Interpolate at upper limit
	      H2  = EU2-E(J2)
	      P   = H2/HE(J2)
	      S2  = (ONE-P)*SP1(J2)+P*SP1(J2+1)
	      ARG = PE*(EUP-EU2)
	      SUM = S2*SINC(ARG)*H2/TWO
	      ARG = PE*(EUP-E(J2))
	      SUM = SUM+SP1(J2)*SINC(ARG)*(H2+HE(J2-1)/TWO)
	      DO IE=1,J2-1
	        ARG = PE*(EUP-E(IE))
	        SUM = SUM+SP1(IE)*SINC(ARG)*HA(IE)
	      END DO ! IE
	      SP2(IEU) = SUM/DEW
	    ELSE ! J2 = NE
	      SUM = ZERO		    
	      DO IE=1,NE
	        ARG = PE*(EUP-E(IE))
	        SUM = SUM+SP1(IE)*SINC(ARG)*HA(IE)
	      END DO ! IE
	      SP2(IEU) = SUM/DEW
	    END IF ! J2

	  ELSE IF (J2 .LT. NE) THEN
	    IF (J1 .EQ. J2) THEN
	      SP2(IEU) = ZERO
	    ELSE
C  Interpolate at lower limit
	      H1  = E(J1+1)-EU1
	      P   = H1/HE(J1)
	      S1  = P*SP1(J1)+(ONE-P)*SP1(J1+1)
	      ARG = PE*(EUP-EU1)
	      SUM = S1*SINC(ARG)*H1/TWO
	      ARG = PE*(EUP-E(J1+1))
	      SUM = SUM+SP1(J1+1)*SINC(ARG)*(H1+HE(J1+1))/TWO
C  Interpolate at upper limit
	      H2  = EU2-E(J2)
	      P   = H2/HE(J2)
	      S2  = (ONE-P)*SP1(J2)+P*SP1(J2+1)
	      ARG = PE*(EUP-EU2)
	      SUM = SUM+S2*SINC(ARG)*H2/TWO
	      ARG = PE*(EUP-E(J2))
	      SUM = SUM+SP1(J2)*SINC(ARG)*(H2+HE(J2-1)/TWO)
	      DO IE=J1+2,J2-1
	        ARG = PE*(EUP-E(IE))
	        SUM = SUM+SP1(IE)*SINC(ARG)*HA(IE)
	      END DO ! IE
	      SP2(IEU) = SUM/DEW
	    END IF

	  ELSE ! J2 = NE
C  Interpolate at lower limit
	    H1  = E(J1+1)-EU1
	    P   = H1/HE(J1)
	    S1  = P*SP1(J1)+(ONE-P)*SP1(J1+1)
	    ARG = PE*(EUP-EU1)
	    SUM = S1*SINC(ARG)*H1/TWO
	    ARG = PE*(EUP-E(J1+1))
	    SUM = SUM+SP1(J1+1)*SINC(ARG)*(H1+HE(J1+1))/TWO
	    DO IE=J1+2,NE
   	      ARG = PE*(EUP-E(IE))
	      SUM = SUM+SP1(IE)*SINC(ARG)*HA(IE)
	    END DO ! IE
	    SP2(IEU) = SUM/DEW
	  END IF
	END DO ! IEU
		
c  Return in original array
	DO IEU=1,NEU
	  SP1(IEU) = SP2(IEU)
	END DO ! IEU

	RETURN
	END ! CONVOLUTE_ENERGY_VSTEP

	REAL*8 FUNCTION SINC(ARG)
C+
C 
C FUNCTIONAL DESCRIPTION:	
C 
C  Routine to calculate the SINC function.
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
C   2-APR-1991
C 
C FUNCTION VALUE:
C 
C  The SINC value
C 
C FORMAL PARAMETERS:
C  
C  Input arguments:
C  ARG		Value of argument
C  
C COMMON BLOCKS:
C  None
C  
C MODIFICATION HISTORY:
C 
C	 Date     | Name  | Description
C ----------------+-------+-----------------------------------------------------
C-
C  Declarations of scalars:
	REAL*8		ARG,SL
 
C  Labeled constants:
	REAL*8		ONE,EPS
	PARAMETER	(ONE=1.0D0,EPS=1.0D-6)
	 
	IF (ABS(ARG) .GT. EPS) THEN
	  SL   = SIN(ARG)/ARG
	  SINC = SL*SL
	ELSE
	  SINC = ONE
	END IF ! ARG

	RETURN
	END ! SINC
