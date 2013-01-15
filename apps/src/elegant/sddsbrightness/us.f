	SUBROUTINE US(ENERGY_I, CUR_I, SIGX_I,
     1    SIGY_I, SIGX1_I, SIGY1_I,
     2    PERIOD_I, N_I, KX_I, KY_I, EMINU_I,
     3    EMAXU_I, NEU_I, DD_I, XPC_I, YPC_I, 
     4    XPS_I, YPS_I, NXP_I, NYP_I,
     5    MODE_I, METHOD_I, IHARM_I,
     6    NPHI_I, NALPHA_I, CALPHA2_I, NOMEGA_I,
     7    COMEGA_I, NSIGMA_I,
     8    E1Z_O, LAMDA1_O, PTOT_O, PD_O, TOTPOWER_O, ISUB_O,
     9    IANG_O, TOTFLUX_O, IMAX_O, IMIN_O,
     1    P1_O, P2_O, P3_O, P4_O, NE_O, XP_O, YP_O, RA0_O,
     2    SPEC0_O, EU_O)
C       +
C       PROGRAM DESCRIPTION:	
C       Program to calculate undulator spectra within the Bessel function 
C       approximation for an ideal planar undulator or an ideal elliptical
C       undulator (including polarization in both cases).  
C       The program may be executed from the xop interface.
C 
C AUTHORS: 
C  Roger J. Dejus Modified by Hairong Shang to be called by c routines.
C  The Advanced Photon Source
C  Experimental Facilities Division
C  Argonne National Laboratory
C 
C CREATION DATE: 
C  25-MAR-1991
C 
C INPUT PARAMETERS:
C  The input parameters are divided into sections related to the storage ring,
C  the undulator device, and the quantity to be calculated. Note: When modifying
C  parameters under the Xus interface, double click the field and make sure to
C  press the RETURN key so that the new parameter is accepted.
C Machine Parameters:
C ENERGY_I     Storage ring energy 			(GeV)
C CUR_I        Storage ring current			(mA)
C SIGX_I      RMS beam size (horizontal)		(mm)
C SIGY_I      RMS beam size (vertical)		(mm)
C SIGX1_I     RMS beam divergence (horizontal)	(mrad)
C SIGY1_I     RMS beam divergence (vertical)	(mrad)
C Undulator Parameters:
C PERIOD_I    Period length			(cm)
C N_I         Number of periods
C KX_I        Deflection parameter (hor.  field) Kx (= 0.0 for a regular planar device)
C KY_I        Deflection parameter (vert. field) Ky
C Scan Parameters:
C EMINU_I     Minimum energy			(eV)
C EMAXU_I     Maximum energy			(eV)
C NEU_I       Number of energy points
C Pinhole Parameters:
C DD_I Distance from the source		(m)
C    (d=0.0 => angular units)
C XPC_I   X-coordinate for center of pinhole	(mm) or (mrad)
C YPC_I   Y-coordinate for center of pinhole	(mm) or (mrad)
C XPS_I   X-size of pinhole (full width)	(mm) or (mrad)
C YPS_I   Y-size of pinhole (full width)	(mm) or (mrad)
C         (for angular units (d=0.0) values are entered in mrad)
C    (X is for horizontal direction)
C    (Y is for the vertical direction)
C NXP_I    Number of subdivisions of pinhole in X (max 500)
C NYP_I    Number of subdivisions of pinhole in Y (max 500)
C    (for plotting 3d results with Xus, the X-size, Y-size, and the number of
C     of subdivisions in the two directions should be equal)
C
C Mode:
C  Depending on the mode selected, some of the pinhole parameters may be
C  set to different values by the program; see the output file us.plt.
C  MODE    1    Angular/spatial flux density distribution
C  MODE    2    Angular/spatial flux density spectrum
C  MODE    3    On-axis brilliance spectrum
C  MODE    4    Flux spectrum through a pinhole
C  MODE    5    Flux spectrum integrated over all angles
C  MODE    6    Power density and integrated power
C
C  Angular/spatial flux density distribution
C    - Flux distribution at the energy chosen as minimum energy.
C  Angular/spatial flux density spectrum
C    - Spectrum at any given point in space as selected by the X and Y
C      coordinate for the center of the pinhole. X is horizontal and Y is
C      vertical.
C  On-axis brilliance spectrum
C  Flux spectrum through a pinhole
C    - Spectrum through a pinhole centered at X-center and Y-center with
C      size X-size and Y-size.  The energy range is from the minimum to the
C      maximum energy.
C  Flux spectrum integrated over all angles.
C    -  The pinhole parameters have no significance here.
C  Power density and integrated power
C    -  Integrated over all energies, thus the energy parameters have no
C       significance here.
C
C Method:
C  METHOD  1    Non-zero emittance; finite-N
C  METHOD  2    Non-zero emittance; infinite-N
C  METHOD  3    Zero emittance;     finite-N
C  METHOD  4    Non-zero emittance; infinite-N + convolution (Dejus' approach)
C  METHOD 14    Non-zero emittance; infinite-N + convolution (Walker's approach)
C
C  Non-zero emittance; finite-N
C    - Use only for "Angular/spatial flux density distribution" and for
C      "Power density and integrated power".
C  Non-zero emittance; infinite-N
C    - For test purposes; do not use (will be removed from menu).
C  Zero emittance; finite-N
C    - Use for zero emittance calculations.
C  Non-zero emittance; infinite-N/convolution
C    - Generally, use for cases where emittance should be included.
C
C Harmonic Number:
C  IHARM   0    All harmonics
C  IHARM  -1    Lowest order harmonic (except MODE=6, include to -IHARM)
C  IHARM   I    I'th harmonic
C
C  All harmonics
C    - Selects all contributing harmonics (generally used).
C  Lowest order harmonic
C    - Selects the lowest order contributing harmonic.
C  Harmonic #
C    - Selects the harmonic number displayed.
C  Edit harmonic number
C    - Modifies the displayed harmonic number.
C   
C Intrinsic Parameters:
C  Several parameters used in the calculations.  Usually not modified by the
C  user.  Please see me (RJD) for further information.
C
C Polarization:
C  The normalized Stokes parameters are calculated including the 
C  unpolarized component.
C  
C DESIGN ISSUES:
C  Program is based on the Bessel function approximation and is valid in the
C  far-field for an ideal sinusoidal magnetic field profile.
C  
C COPYRIGHT:
C  This routine must only be used at The Advanced Photon Source and must not
C  be tranferred or used at any other location without the written consent
C  of the author.
C  
C FILES USED:
C  Input file - us.dat  File in the user's current directory containing the
C                       input parameters.
C  Output file - us.plt File in the user's current directory containing the 
C                       results of the calculation.  The header contains
C                       all input parameters and the calculated zero emittance
C			on-axis first harmonic energy (e1), corresponding
C                       wavelength (l1), total power (ptot), and the on-axis
C                       power density (pd).
C KEYWORDS:
C  Undulator Spectrum, Bessel Function Approximation
C  
C LINK/LIBRARY ISSUES:
C  Calls routines BRIGHTE and HUNT.  BRIGHTE calculates the brilliance and HUNT
C  searches an array of real numbers (from Numerical Recipes).
C  
C PORTABILITY ISSUES:
C  Runs on DEC 3000/400 AXP alpha (Tru64Unix v5.0), SUN (Solaris: SunOS
C  Release v5.6), and Windows 95/98/NT (Pentium and higher).
C  
C TIMING:
C  Execution times vary considerably depending on computer and the 
C  quantity being calculated.  The zero emittance calculations are fast
C  (few seconds), whereas the non-zero emittance calculations may range from
C  seconds (on-axis brilliance) to an hour (flux spectrum through a pinhole).
C
C VERSION:
C  1.91
C  
C MODIFICATION HISTORY:
C 
C	 Date     | Name  | Description
C ----------------+-------+-----------------------------------------------------
C 06-JUL-1994     | RJD   | Modified value for E1MIN for angle-integrated
C                 |       | spectrum (MODE=5) to be non-zero; gamma*theta 
C                 |       | corresponds to sqrt(200) (somewhat arbitrarily
C                 |       | chosen)
C ----------------+-------+-----------------------------------------------------
C 04-OCT-1994     | RJD   | Modified program to include polarization properties.
C                 |       | The four Stokes parameters are now calculated.
C                 |       | Program is for an ideal planar undulator or an ideal
C                 |       | elliptical undulator. Many other changes. The value
C                 |       | of the parameter IHARM has a different meaning.
C                 |       | IHARM=0 now gives 'all harmonics' and IHARM= <0
C                 |       | gives the lowest order harmonic except for the power
C                 |       | option. For the power option, a negative IHARM means
C                 |       | include all harmonics up to and including -IHARM.
C                 |       | This version is 1.6.
C ----------------+-------+-----------------------------------------------------
C 21-JUN-1995     | RJD   | Modified print-out of "Contributing harmonics" in
C		  |       | subroutine PRINT_OUT. Routine incorrectly calculated 
C		  |       | IMIN and IMAX for METHOD 4 (Dejus' method) for
C		  |       | "Spectral distributions". The spectra and integrated
C		  |       | quantities were calculated correctly and are 
C		  |       | unaffected by this modification.
C                 |       | The current version is 1.7.
C ----------------+-------+-----------------------------------------------------
C 04-JAN-1996     | RJD   | Modified the number of decimal places for the sigx1
C                 |       | and sigy1 variables to four in the printout. Added
C                 |       | one more digit for the emax variable to avoid
C                 |       | overflow on rare occasions. Formats 260 and 256 were
C                 |       | changed accordingly.
C                 |       | The current version is 1.8.
C ----------------+-------+-----------------------------------------------------
C 11-NOV-1997     | RJD   | Changed notation: Brightness -> Brilliance.
C                 |       | The current version is 1.9.
C ----------------+-------+-----------------------------------------------------
C 16-JUL-2000     | RJD   | Minor change in the code to compile error-free on
C                 |       | Unix and Windows (no change in results vs. v1.9).
C                 |       | Current version is v1.91.
C ----------------+-------+-----------------------------------------------------
C [change_entry]
C-
C  input parameters
	INTEGER*4       MODE_I,METHOD_I,IHARM_I,N_I,NSIGMA_I
	INTEGER*4	NXP_I,NYP_I,NPHI_I,NALPHA_I, NEU_I, NOMEGA_I
	REAL*8          ENERGY_I, CUR_I, SIGX_I, SIGY_I
	REAL*8          SIGX1_I, SIGY1_I, PERIOD_I, KX_I, KY_I
	REAL*8          EMINU_I, EMAXU_I, XPS_I, YPS_I, DD_I
	REAL*8          CALPAH2_I, COMEGA_I, XPC_I, YPC_I

C  output parameters
	INTEGER*4       ISUB_O, IANG_O, IMAX_O, IMIN_O, NE_O
	REAL*8          E1Z_O, LAMDA1_O, PTOT_O, PD_O, TOTFLUX_O
	REAL*8          TOTPOWER_O
	REAL*8          P1_O(*), P2_O(*), P3_O(*), P4_O(*)
	REAL*8          XP_O(*), YP_O(*), RA0_O(*), SPEC0_O(*)
	REAL*8          EU_O(*)
C       Size parameters:
	INTEGER*4	E_SZ,A_SZ,B_SZ,P_SZ
	PARAMETER	(E_SZ=50001,A_SZ=400,B_SZ=A_SZ/4+1,P_SZ=501)
	INTEGER*4	RECL_SZ,BLOCK_SZ,BUFFER_SZ
	PARAMETER	(BLOCK_SZ=12288,BUFFER_SZ=1) ! Blocksize in bytes
	
C       Declarations of scalars:
	INTEGER*4	IERROR,J,L,N1,N2,NOMEGA,NEI,NALPHAP
	INTEGER*4	NCH,NCH1,NCH2
	INTEGER*4	I,IA,IB,IE,ID,IW,ISIGN,ISUB,IE1,IE2,IDW
	INTEGER*4       IMIN, IMAX, IEU
	REAL*8          ENERGY,CUR,PERIOD,SIGX,SIGY,SIGX1,SIGY1
	REAL*8		EMIN,EMAX,DU,COMEGA
	REAL*8		G2,LAMDAR,LAMDA1,E1Z,ARG,D2,SL
	REAL*8		SIGU2,SIGV2,SIGU,SIGV
	REAL*8		XE,XE1,XE2,YE,YE1,YE2,DEU
	REAL*8		XPMIN,YPMIN,E1MIN,E1MAX,EP,EF,DEF
	REAL*8		PHI1,PHI2,PTOT,PD,GK,K2
	REAL*8          FLUX, POWER, WE,P1,P2,P3,P4

C  Declarations of arrays:
	REAL*4          TD(2)

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
	PARAMETER       (VN='v. 1.91')
	INTEGER*4	IDWMAX
	PARAMETER	(IDWMAX=128)
	REAL*8		PI,PIHALF,TWOPI
	PARAMETER	(PI    =3.1415 92653 58979 32384 62643D0)
	PARAMETER	(PIHALF=1.5707 96326 79489 66192 31322D0)
	PARAMETER	(TWOPI= 6.2831 85307 17958 64769 25287D0)
	REAL*8		FINE_STRUCTURE_CONST
	PARAMETER       (FINE_STRUCTURE_CONST=1.0D19*EC*1.0D19*EC
     1   /(4.0D0*PI*EPSZ*1.0D38*HBAR*C))
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
	DATA            TD /0.0,0.0/
	
C  Common blocks:
	LOGICAL*4	LANG
	INTEGER*4       MODE,METHOD,IHARM,N,NSIGMA
	INTEGER*4	NXP,NYP,NPHI4,NPHI,NPHI_BRIGHT
	INTEGER*4	NALPHA,NE,NE1,NE2,NEU,NW
	REAL*8		D,EMINU,EMAXU,XPC,YPC,XPS,YPS
	REAL*8		K,KX,KY,K3,NPI,PE,GAMMA,ER,EW,DEW,LEN
	REAL*8		FAC,C1,C2,C3,C4,C5
	REAL*8		AP2MIN,AP2MAX,AP2CNT,ARGMAX,FU,FV,SIGX2,SIGY2,GSIGUV
	REAL*8		DXP,DYP,DPHI,CALPHA2,DE,DW

	INTEGER*4	INDEX_PHI(A_SZ)
	INTEGER*4	I1(E_SZ),I2(E_SZ)
	REAL*8		COSPHI(A_SZ),SINPHI(A_SZ),S2SIGN(A_SZ)
	REAL*8		E(E_SZ),EU(E_SZ)
	REAL*8		SPEC0(E_SZ),SPEC1(E_SZ),SPEC2(E_SZ),SPEC3(E_SZ)
	REAL*8		HE(E_SZ),HA(E_SZ),HW(E_SZ)
	REAL*8		BR0(B_SZ,B_SZ),BR1(B_SZ,B_SZ)
	REAL*8		BR2(B_SZ,B_SZ),BR3(B_SZ,B_SZ)
	REAL*8		RA0(P_SZ,P_SZ),RA1(P_SZ,P_SZ)
	REAL*8		RA2(P_SZ,P_SZ),RA3(P_SZ,P_SZ)
	REAL*8		XP(P_SZ),YP(P_SZ),CX(P_SZ),CY(P_SZ)

	COMMON		/PRTI/      LANG
	COMMON		/PRTR/      D,EMINU,EMAXU,XPC,YPC,XPS,YPS
	COMMON		/CALCI/     MODE,METHOD,IHARM,N,NSIGMA
	COMMON		/CALCR/	    K,KX,KY,K3,NPI,PE,GAMMA,ER,EW,DEW,LEN
	COMMON		/FACTOR/    FAC,C1,C2,C3,C4,C5
	COMMON		/BEAM1/	    AP2MIN,AP2MAX,AP2CNT,ARGMAX,
     1			    FU,FV,SIGX2,SIGY2,GSIGUV
	COMMON		/PINHOLE/   XP,YP,CX,CY,NXP,NYP,DXP,DYP
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

C+ Read parameters from data file
	ENERGY=ENERGY_I
	CUR=CUR_I
	SIGX=SIGX_I
	SIGY=SIGY_I
	SIGX1=SIGX1_I
	SIGY1=SIGY1_I
	PERIOD=PERIOD_I
	N=N_I
	KX=KX_I
	KY=KY_I
	DU=DD_I
	XPC=XPC_I
	YPC=YPC_I
	XPS=XPS_I
	YPS=YPS_I
	NXP=NXP_I
	NYP=NYP_I
	MODE=MODE_I
	METHOD=METHOD_I
	IHARM=IHARM_I
	NPHI=NPHI_I
	NALPHA=NALPHA_I
	CALPHA2=CALPHA2_I
	NOMEGA=NOMEGA_I
	COMEGA=COMEGA_I
	NSIGMA=NSIGMA_I
	EMINU=EMINU_I
	EMAXU=EMAXU_I
	NEU=NEU_I
	
C Check for valid input
C (MODE .EQ. 1 .AND. METHOD .EQ. 2) limit use to finite-N
C MODE .EQ. 1 .AND. METHOD .EQ. 4 and zero emittance
C MODE .EQ. 5 .AND. METHOD .EQ. 3 limit use to infinite-N + convolution
C -----------------------------------------------------------------------------
	IF (XPC .LT. ZERO .OR. YPC .LT. ZERO) GO TO 920
	IF ((MODE .EQ. 1 .AND. METHOD .EQ. 2) .OR.
     1    (MODE .EQ. 1 .AND. METHOD .EQ. 4) .OR.
     2    (MODE .EQ. 1 .AND. METHOD .EQ.14)) GO TO 930 
	IF ((MODE .EQ. 5 .AND. METHOD .EQ. 1) .OR.
     1    (MODE .EQ. 5 .AND. METHOD .EQ. 3)) GO TO 940
                                                     
C+ Default values
C  -----------------------------------------------------------------------------
	IF (NPHI    .EQ. 0) NPHI   = 20
	IF (NALPHA  .EQ. 0) NALPHA = 15
	IF (NOMEGA  .EQ. 0) NOMEGA = 16
	IF (NSIGMA  .EQ. 0) NSIGMA = 3
	IF (CALPHA2 .EQ. ZERO) CALPHA2 = 2.0D0
	IF (COMEGA  .EQ. ZERO) COMEGA  = 2.0D0
	
	L = 80

	NCH1 = L
	
C       Determine units (angular or spatial units)
	IF (MODE .EQ. 3 .OR. DU .EQ. ZERO) THEN
	   IANG_O = 1
	   LANG = .TRUE.
	   D    = ONE
	ELSE
	   IANG_O = 0
	   LANG = .FALSE.
	   D    = DU
	END IF			!
	
	IF (MODE .EQ. 2) THEN
	   XPS = ZERO
	   YPS = ZERO
	   NXP = 0
	   NYP = 0
	ELSE IF (MODE .EQ. 3 .OR. MODE .EQ. 5) THEN
	   XPC = ZERO
	   YPC = ZERO
	   XPS = ZERO
	   YPS = ZERO
	   NXP = 0
	   NYP = 0
	END IF			! MODE
	
C+ Definition of constants and calculation of power/power density
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
	GK     = ZERO
	IF (KX .LT. EPSK .OR. KY .LT. EPSK) THEN
	    K  = KX +KY
            GK = K*((K**6)+(24.0D0*(K**4)/7.0D0)+
     1        (4.0D0*K*K)+(16.0D0/7.0D0))/((1.0D0+(K*K))**3.5D0)
	END IF 
	IF (ABS(KX-KY) .LT. EPSK) THEN
	    K  = KX
	    GK = 32.0D0/7.0D0*K/((1.0D0+(K*K))**3.0D0)
	END IF
        PTOT   = PTOT_FAC*N*K2*  (ENERGY**2)*CUR*C_MA_A/(PERIOD*C_CM_M)![W]
        PD     = PD_FAC  *N*K*GK*(ENERGY**4)*CUR*C_MA_A/(PERIOD*C_CM_M)![W/mrad^2]

C+ Beam emittance
C  -----------------------------------------------------------------------------
	FU    = ZERO
	FV    = ZERO
	SIGX2 = SIGX*SIGX
	SIGY2 = SIGY*SIGY
	IF (LANG) THEN ! Angular units
	    SIGU2 = SIGX1*SIGX1*C_MRAD_RAD*C_MRAD_RAD
	    SIGV2 = SIGY1*SIGY1*C_MRAD_RAD*C_MRAD_RAD
	ELSE ! Spatial
	    SIGU2 = (SIGX1*SIGX1+SIGX2/D2)*C_MRAD_RAD*C_MRAD_RAD
	    SIGV2 = (SIGY1*SIGY1+SIGY2/D2)*C_MRAD_RAD*C_MRAD_RAD
	END IF ! LANG
	SIGU  = SQRT(SIGU2) ! [rad]
	SIGV  = SQRT(SIGV2) ! [rad]
	IF (SIGU2 .NE. ZERO) FU = 0.5D0/SIGU2
	IF (SIGV2 .NE. ZERO) FV = 0.5D0/SIGV2
	GSIGUV = GAMMA*MIN(SIGU,SIGV)
C-
C  -----------------------------------------------------------------------------

C+ Determine min and max emission angles and center of pinhole
C  -----------------------------------------------------------------------------
	XE  = (XPC-XPS/TWO)*C_MM_M/D-NSIGMA*SIGU ! Cartesian angle in x-dir.
	YE  = (YPC-YPS/TWO)*C_MM_M/D-NSIGMA*SIGV ! Cartesian angle in y-dir.
	XE1 = XE	    
	YE1 = YE
	IF (XE .LT. ZERO) XE = ZERO
	IF (YE .LT. ZERO) YE = ZERO
	AP2MIN = G2*(XE*XE+YE*YE)

	XE  = (XPC+XPS/TWO)*C_MM_M/D+NSIGMA*SIGU ! Cartesian angle in x-dir.
	YE  = (YPC+YPS/TWO)*C_MM_M/D+NSIGMA*SIGV ! Cartesian angle in y-dir.
	XE2 = XE	    
	YE2 = YE
	AP2MAX = G2*(XE*XE+YE*YE)
	ARGMAX = NSIGMA*NSIGMA/TWO

	XE = XPC*C_MM_M/D			! Cartesian angle in x-dir.
	YE = YPC*C_MM_M/D			! Cartesian angle in y-dir.
	AP2CNT = G2*(XE*XE+YE*YE)
C-
C  -----------------------------------------------------------------------------

C+ Modify NALPHA for power option
C MODE .EQ. 6 .AND. METHOD .NE. 3 !power option and non-zero emitt.
C  -----------------------------------------------------------------------------
	IF (MODE .EQ. 6 .AND. METHOD .NE. 3) THEN
	    NALPHAP  = (SQRT(AP2MAX)-SQRT(AP2MIN))/GSIGUV/0.8D0
	    IF (NALPHAP .LT. MAX(15,NALPHA)) THEN
		NALPHA = MAX(15,NALPHA)
	    ELSE
		NALPHA = NALPHAP
	    END IF ! NALPHAP
	    IF (NALPHA .GT. B_SZ) NALPHA = B_SZ
	END IF
C-
C  -----------------------------------------------------------------------------


C+ Define energy scale and set up array for the line shape function
C  -----------------------------------------------------------------------------
	IF (METHOD .EQ. 4 .AND. MODE .NE. 1 .AND. MODE .NE. 6) THEN
	   IF (MODE .EQ. 5) THEN
	      E1MIN = E1Z*K3/(K3+200.0d0) ! use large value for gamma*theta
	      E1MAX = E1Z
	      DEW   = E1Z/N
	   ELSE
	      E1MIN = E1Z    *K3/(K3+AP2MAX)
	      E1MAX = E1Z    *K3/(K3+AP2MIN)
	      DEW   = (E1Z/N)*K3/(K3+AP2CNT)
	   END IF		! MODE
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
	   END DO		! WHILE
	   IE1  = I
	   EMIN = MAX(IE1*E1MIN,EMIN)
	   
	   IF (MODE .NE. 5) THEN
	      EP = IE1*E1MIN
	      DO WHILE (EP .LT. EMAX)
		 I  = I+1
		 EP = I*E1MIN
	      END DO		! WHILE
	      IE2  = I-1
	      EMAX = MIN(IE2*E1MAX,EMAX)
	   END IF		! MODE
	   
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
		 END IF		! IDW
	      END DO		! WHILE
	      
	      DO WHILE (E(IE) .LT. (EP-EPSE) .AND. 
     1	                E(IE) .LT. (EMAX-EPSE))
		 IF (E(IE) .GT. EP-EF) THEN
		    IDW = IDW*2
		    IF (IDW .GT. IDWMAX) THEN
		       EF = ZERO
		    ELSE
		       DEF = DW /IDW
		       EF  = DEW/IDW
		    END IF	! IDW
		 END IF		! E(IE)
		 
		 IE       = IE+1
		 IF (IE .GT. E_SZ) GO TO 900
		 HE(IE-1) = DEF
		 IF (IE .EQ. 2) THEN
		    HA(IE-1) = HE(IE-1)/TWO
		 ELSE
		    HA(IE-1) = (HE(IE-2)+HE(IE-1))/TWO ! Average
		 END IF		! IE
		 E(IE) = E(IE-1)+HE(IE-1)
	      END DO		! WHILE
	      

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
	      
	      I   = I+1		! Next harmonic
	      DEF = I*E1MIN-E(IE)
	      IF (DEF .GT. ZERO) THEN
		 IE	     = IE+1
		 IF (IE .GT. E_SZ) GO TO 900
		 HE(IE-1) = DEF
		 HA(IE-1) = (HE(IE-2)+HE(IE-1))/TWO ! Average
		 E (IE)   = E(IE-1)+HE(IE-1)
	      END IF		! DEF
	      
	      EP  = I*E1MAX
	   END DO		! WHILE

C  Adjust end points
	   IF (DEF .GT. ZERO) THEN
	      NE = IE-1
	   ELSE
	      NE = IE
	   END IF		! DEF
	   NE1 = 1
	   NE2 = NE
	   HE(NE) = ZERO
	   HA(NE) = HE(NE-1)/TWO
	   
	   IF (NEU .EQ. 0) NEU = INT(NE2/100.0D0+ONE)*100.0D0 ! Default
C	    DEU = (EMAXU-EMINU)/NEU
	   DE  = (EMAXU-EMINU)/NEU
	   NEU = NEU+1
	   DO IE=1,NEU		! Energy scale (user's selection)
C		EU(IE) = EMINU+(IE-1)*DEU
	      EU(IE) = EMINU+(IE-1)*DE
	   END DO		! IE
   	   
	ELSE IF (METHOD .EQ. 14 .AND. MODE .NE. 1 .AND. MODE .NE. 6) THEN 
C   ! Walker's
	   DEW = (E1Z/N)*K3/(K3+AP2CNT)
	   EW  = COMEGA*DEW	! [eV]
	   IF (NEU .GT. 0) THEN ! Redefine NOMEGA
	      SL     = (EMAXU-EMINU)/NEU
	      NOMEGA = TWO*EW/SL+ONE
	      NOMEGA = NOMEGA/2*2 ! Make even
	      IF (NOMEGA .LT. 16) THEN ! Reset to 16
		 NOMEGA = 16
	      END IF		! NOMEGA
	   END IF		! NEU
	   
	   DW   = TWO*COMEGA/NOMEGA
	   DE   = DW*DEW
	   EMIN = EMINU-EW	! Extend lower limit of energy scale
	   IF (EMIN .LT. EPS) EMIN = EPS
	   EMAX = EMAXU+EW	! Extend upper limit of energy scale
	   NE   = (EMAX-EMIN)/DE+ONE ! Redefine the # of intervals
	   NE   = NE+1		! # of points
	   NE1  = NOMEGA/2+1	! Smallest index; odd
	   NE2  = NE+1-NE1	! Largest index
	   NEI  = NE2-NE1	! # of intervals

C  Generate line shape function (used in energy convolution)
	   NW = NOMEGA+1	! # of points
	   DO IW=1,NW
	      ARG = (-COMEGA+(IW-1)*DW)*PI	    
	      IF (ABS(ARG) .GT. EPS) THEN
		 SL = SIN(ARG)/ARG
		 HW(IW) = SL*SL
	      ELSE
		 HW(IW) = ONE
	      END IF		! ARG
	   END DO		! IW
C       Generate energy scale
	   IF (NE .GT. E_SZ) GO TO 905
	   DO IE=1,NE
	      E(IE) = EMIN+(IE-1)*DE ! [eV]
	   END DO		! IE
	ELSE
	   NE  = NEU	    
	   IF (MODE .EQ. 1 .OR. MODE .EQ. 6) NE = 0 ! no energy dependence
	   IF (NE .GT. 0) DE  = (EMAXU-EMINU)/NE
	   NE  = NE+1		! # of points
	   NE1 = 1		! Smallest index
	   NE2 = NE+1-NE1	! Largest index
	   NEI = NE2-NE1	! # of intervals
C  Generate energy scale
	   IF (NE .GT. E_SZ) GO TO 905
	   DO IE=1,NE
	      E(IE) = EMINU+(IE-1)*DE ! [eV]
	   END DO		! IE
	END IF			! METHOD & MODE
C       -
C  -----------------------------------------------------------------------------
	
C       + Pinhole parameters
C  -----------------------------------------------------------------------------
	IF (XPC .EQ. ZERO .AND. YPC .EQ. ZERO) THEN ! Pinhole centered
	   FAC   = 4.0D0
	   XPMIN = ZERO
	   YPMIN = ZERO
	   IF (NXP .GT. 0) DXP = XPS/TWO/NXP
	   IF (NYP .GT. 0) DYP = YPS/TWO/NYP
	ELSE
	   FAC   = ONE
	   XPMIN = XPC-XPS/TWO
	   YPMIN = YPC-YPS/TWO
	   IF (NXP .GT. 0) DXP = XPS/NXP
	   IF (NYP .GT. 0) DYP = YPS/NYP
	END IF
	NXP = NXP+1
	NYP = NYP+1
C-
C  -----------------------------------------------------------------------------

C+ Set up positions within pinhole and associated cartesian angles
C  -----------------------------------------------------------------------------
	DO IA=1,NXP
	   XP(IA) = XPMIN+(IA-1)*DXP ! Position [mm]
	   CX(IA) = XP(IA)*C_MM_M/D ! Angle   [rad]
	END DO			! IA
	DO IB=1,NYP
	   YP(IB) = YPMIN+(IB-1)*DYP ! Position [mm]
	   CY(IB) = YP(IB)*C_MM_M/D ! Angle   [rad]
	END DO ! IB
C-
C  -----------------------------------------------------------------------------

C+ Set up arrays of cos(phi) and sin(phi), and index_phi
C  -----------------------------------------------------------------------------
	IF (METHOD .NE. 3) THEN ! non-zero emittance
	   IF (XE1 .LE. ZERO .AND. YE1 .LE. ZERO) THEN
               IF ((MODE .EQ. 2 .AND. XPC .EQ. 
     1              ZERO .AND. YPC .EQ. ZERO) .OR.
     2              MODE .EQ. 3) THEN
		 NPHI4 = NPHI   ! Use only 0 ... 90 deg.
	      ELSE
		 NPHI4 = 4*NPHI ! Use full range: 0 ... 360 deg.
	      END IF		! MODE
	      NPHI_BRIGHT = NPHI 
C           ! # of indices for the brilliance array (use sym.)
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
		 END IF		! ID
		    IF (ISIGN .NE. 0) THEN
		       S2SIGN(ID) = ISIGN
		    ELSE	! ISIGN = 0
		       S2SIGN(ID) = +S2SIGN(ID-1)
		    END IF
		 END DO		! ID
	    ELSE 
	       IF (XE1 .GT. ZERO) THEN
		  IF (YE1 .GT. ZERO) THEN
		     PHI1 = ATAN(YE1/XE2)
		     IF (XE1 .GT. EPS) THEN
			PHI2 = ATAN(YE2/XE1)
		     ELSE
			PHI2 = PIHALF
		     END IF
		     DPHI  = (PHI2-PHI1)/NPHI
		     NPHI4 = NPHI
		     NPHI_BRIGHT = NPHI
		     ISIGN = +1
		     SL    = PHI1+ISIGN*DPHI/TWO
		  ELSE		! YE1 le 0.0
		     IF (XE1 .GT. EPS) THEN
			PHI1 = ATAN(YE1/XE1)
			PHI2 = ATAN(YE2/XE1)
		     ELSE
			PHI1 = -PIHALF
			PHI2 = +PIHALF
		     END IF
		     DPHI  =  PHI2/NPHI
		     N1    = -PHI1/DPHI+ONE
		     N2    =  NPHI+1
		     NPHI4 = N1+N2
		     NPHI_BRIGHT = N2 ! N2 >= N1
		     ISIGN = +1
		     SL    = ISIGN*DPHI/TWO
		  END IF	! YE1 gt 0
	       ELSE		! XE1 le 0.0
		  PHI1 = ATAN(YE1/XE2)
		  IF (XE1 .LT.- EPS) THEN
		     PHI2 = ATAN(YE1/XE1)+PI
		  ELSE
		     PHI2 = PIHALF
		  END IF
		  DPHI  = (PIHALF-PHI1)/NPHI
		  N1    = NPHI+1
		  N2    = (PHI2-PIHALF)/DPHI+ONE
		  NPHI4 = N1+N2
		  NPHI_BRIGHT = N1 ! N1 >= N2
		  ISIGN = -1
		  SL    = PIHALF+ISIGN*DPHI/TWO
	       END IF		! XE1 gt 0
	       DO ID =1,NPHI4
		  IF (ID .EQ. NPHI_BRIGHT+1) THEN
		     ISIGN = -ISIGN
		     SL    = ARG
		  END IF	! ID
		  ARG        = SL+ISIGN*(ID-1)*DPHI
		  COSPHI(ID) = COS(ARG)
		  SINPHI(ID) = SIN(ARG)
		  IF (ID .LE. NPHI_BRIGHT) THEN
		     INDEX_PHI(ID) = ID
		     S2SIGN(ID)    = +ONE
		  ELSE
		     INDEX_PHI(ID) = ID-NPHI_BRIGHT
		     S2SIGN(ID)    = -ONE
		  END IF	! ID
	       END DO		! ID
	    END IF
	 END IF			! METHOD
C-
C  -----------------------------------------------------------------------------

C+ Scale factors
C  -----------------------------------------------------------------------------
         IF (SIGU .NE. ZERO .AND. SIGV .NE. ZERO)
     1      C1 =N*N*FINE_STRUCTURE_CONST*BW*CUR*C_MA_A/EC
     2          /(TWOPI*SIGU*SIGV*D2*C_M_MM*C_M_MM)
	 C2 = N*N*   FINE_STRUCTURE_CONST*BW*CUR*C_MA_A/EC
	 C3 = N*N*G2*FINE_STRUCTURE_CONST*BW*CUR*C_MA_A/EC
     1        /(D2*C_M_MM*C_M_MM)
	 C4 = PDH_FAC*N*N*(G2**2)*CUR*C_MA_A/(LEN*D2)
	 IF (SIGU .NE. ZERO .AND. SIGV .NE. ZERO)
     1       C5 = PDH_FAC*N*N*G2*CUR*C_MA_A/(TWOPI*SIGU*SIGV*LEN*D2)
C       -
C       -----------------------------------------------------------------------------


C       + Call analysis routine
C       -----------------------------------------------------------------------------
	IF (MODE .EQ. 1) THEN
	   ISUB = 1
	ELSE IF (MODE .GE. 2 .AND. MODE .LE. 4) THEN
	   ISUB = 2
	ELSE IF (MODE .EQ. 5) THEN
	   ISUB = 3
	ELSE IF (MODE .EQ. 6) THEN
	   ISUB = 4
	END IF			! MODE
	IF (METHOD .EQ. 3 .AND. MODE .NE. 6) ISUB = 5
	
	IF (ISUB .EQ. 1) THEN
	   CALL SPACE_DISTRIBUTION(IERROR)
	ELSE IF (ISUB .EQ. 2) THEN
	   CALL SPECTRAL_DISTRIBUTION(IERROR)
	ELSE IF (ISUB .EQ. 3) THEN
	   CALL ANGLE_INTEGRATION(IERROR)
	ELSE IF (ISUB .EQ. 4) THEN
	   CALL POWER_DISTRIBUTION(IERROR)
	ELSE IF (ISUB .EQ. 5) THEN
	   CALL NO_EMITTANCE(IERROR)
	END IF			! ISUB
	PTOT_O = PTOT
	PD_O = PD
	E1Z_O = E1Z
	LAMDA1_O = LAMDA1
C-
C  -----------------------------------------------------------------------------

C+ Print results
C  -----------------------------------------------------------------------------
	IF (IERROR .NE. 0) THEN
	   WRITE(0,*) '&US-F-SUBERR, Subroutine error'
	   WRITE(0,*) '- unsuccessful completion due to error status'
	   STOP
	END IF			! IERROR	    
C  write to output arguments
	DD_I = D
	XPC_I = XPC
	YPC_I = YPC
	XPS_I = XPS
	YPS_I = YPS
	NXP_I = NXP
	NYP_I = NYP
	IF (IANG_O.EQ.1) THEN
		NXP_I=0
		NYP_I=0
	END IF	
	NALPHA_I = NALPHA
	CALPHA2_I = CALPHA2
	NOMEGA_I = NOMEGA
	COMEGA_I = COMEGA
	ISUB_O = ISUB
	FLUX  = ZERO
	POWER = ZERO
	NE_O = 1
C       Space-distribution (non-zero emittance and zero emittance case)
	IF ((ISUB .EQ. 1) .OR. (ISUB .EQ. 5 .AND. MODE .EQ. 1)) THEN
	   CALL TRAPZ2(RA0,FLUX) ! get integrated flux over observation area
	   FLUX  = FAC*FLUX	! ph/s/0.1%bw
	   POWER = FLUX/BW*EC	! W/eV
	   TOTFLUX_O = FLUX
	   TOTPOWER_O = POWER
	   IMIN_O = I1(1)
	   IMAX_O = I2(1)
	   NE_O = 1
	   DO IA=1,NXP
	      DO IB=1,NYP
		 IF (RA0(IA,IB) .GT. ZERO) THEN
		    P1 = RA1(IA,IB)/RA0(IA,IB)	
		    P2 = RA2(IA,IB)/RA0(IA,IB)	
		    P3 = RA3(IA,IB)/RA0(IA,IB)	
		    P4 = ONE -SQRT(P1*P1 +P2*P2 +P3*P3) ! unpolarized
		 ELSE
		    P1 = ZERO
		    P2 = ZERO
		    P3 = ZERO
		    P4 = ZERO
		 END IF 
		 P1_O(NE_O) = P1
		 P2_O(NE_O) = P2
		 P3_O(NE_O) = P3
		 P4_O(NE_O) = P4
		 XP_O(NE_O) = XP(IA)
		 YP_O(NE_O) = YP(IB)
		 RA0_O(NE_O) = RA0(IA,IB)
		 NE_O = NE_O + 1
	      END DO		! IB
	   END DO		! IA   
        ELSE IF ((ISUB .EQ. 2) .OR. (ISUB .EQ. 3) 
     1           .OR. (ISUB .EQ. 5)) THEN
C  Spectral distributions (non-zero emittance and zero emittance case)
C       Get frequency-integrated power (flux) (no meaning for brilliance)
	   IMIN = 10000
	   IMAX = 0
	   DO IE=NE1,NE2
	      IF (I1(IE) .LT. IMIN .AND. I1(IE) .GT. 0) IMIN = I1(IE)
	      IF (I2(IE) .GT. IMAX) IMAX = I2(IE)
	   END DO		! IE
	   IMIN_O = IMIN
	   IMAX_O = IMAX
	   IF (METHOD .EQ. 4) THEN ! Dejus' method
	      DO IEU=1,NEU
		 IF (IEU .EQ. 1 .OR. IEU .EQ. NEU) THEN
		    WE = HALF
		 ELSE
		    WE = ONE
		 END IF		! IEU
		 FLUX  = FLUX +WE*SPEC0(IEU)/EU(IEU) ! -> ps/s/mm^2  /eV or
				! -> ph/s/mrad^2/eV
		 POWER = POWER+WE*SPEC0(IEU) ! -> W/mm^2/eV or W/mrad^2/eV
	      END DO		! IEU
	   ELSE			! METHOD={1,2,3,14}
	      DO IE=NE1,NE2
		 IF (IE .EQ. NE1 .OR. IE .EQ. NE2) THEN
		    WE = HALF
		 ELSE
		    WE = ONE
		 END IF		! IE
		 FLUX  = FLUX +WE*SPEC0(IE)/E(IE) ! -> ps/s/mm^2  /eV or
				! -> ph/s/mrad^2/eV
		 POWER = POWER+WE*SPEC0(IE) ! -> W/mm^2/eV or W/mrad^2/eV
	      END DO		! IE
	   END IF		! METHOD
	   FLUX  = FLUX /BW   *DE ! ph/s/mm^2 or ph/s/mrad^2
	   POWER = POWER/BW*EC*DE ! W or W/mm^2 or W/mrad^2
	   TOTFLUX_O = FLUX
	   TOTPOWER_O = POWER
	   XP_O(1) = XP(1)
	   YP_O(1) = YP(1)
	   IF (METHOD .EQ. 4) THEN ! Dejus'
	      NE_O = NEU
	      DO IEU=1,NEU
		 IF (SPEC0(IEU) .GT. ZERO) THEN
		    P1 = SPEC1(IEU)/SPEC0(IEU)
		    P2 = SPEC2(IEU)/SPEC0(IEU)
		    P3 = SPEC3(IEU)/SPEC0(IEU)
		    P4 = ONE -SQRT(P1*P1 +P2*P2 +P3*P3) ! unpolarized
		 ELSE
		    P1 = ZERO
		    P2 = ZERO
		    P3 = ZERO
		    P4 = ZERO
		 END IF 
		 EU_O(IEU) = EU(IEU)
		 SPEC0_O(IEU) = SPEC0(IEU)
		 P1_O(IEU) = P1
		 P2_O(IEU) = P2
		 P3_O(IEU) = P3
		 P4_O(IEU) = P4
	      END DO		! IEU
	   ELSE			! METHOD={1,2,3,14}
	      DO IE=NE1,NE2
		 IF (SPEC0(IE) .GT. ZERO) THEN
		    P1 = SPEC1(IE)/SPEC0(IE)
		    P2 = SPEC2(IE)/SPEC0(IE)
		    P3 = SPEC3(IE)/SPEC0(IE)
		    P4 = ONE -SQRT(P1*P1 +P2*P2 +P3*P3) ! unpolarized
		 ELSE
		    P1 = ZERO
		    P2 = ZERO
		    P3 = ZERO
		    P4 = ZERO
		 END IF
		 EU_O(NE_O) = E(IE)
		 SPEC0_O(NE_O) = SPEC0(IE)
		 P1_O(NE_O) = P1
		 P2_O(NE_O) = P2
		 P3_O(NE_O) = P3
		 P4_O(NE_O) = P4
		 NE_O = NE_O + 1
	      END DO
	   END IF		! METHOD
	ELSE IF (ISUB .EQ. 4) THEN
C       Power density (non-zero emittance and zero emittance case)
	   CALL TRAPZ2(RA0,POWER)
	   POWER = FAC*POWER	! W
	   TOTPOWER_O = POWER
	   IMIN = I1(1)
	   IMAX = I2(1)
	   IMIN_O = IMIN
	   IMAX_O = IMAX
	   DO IA=1,NXP
	      DO IB=1,NYP
		 IF (RA0(IA,IB) .GT. ZERO) THEN
		    P1 = RA1(IA,IB)/RA0(IA,IB)	
		    P2 = RA2(IA,IB)/RA0(IA,IB)	
		    P3 = RA3(IA,IB)/RA0(IA,IB)	
		    P4 = ONE -SQRT(P1*P1 +P2*P2 +P3*P3) ! unpolarized
		 ELSE
		    P1 = ZERO
		    P2 = ZERO
		    P3 = ZERO
		    P4 = ZERO
		 END IF 
		 XP_O(NE_O) = XP(IA)
		 YP_O(NE_O) = YP(IB)
		 RA0_O(NE_O) = RA0(IA,IB)
		 P1_O(NE_O) = P1
		 P2_O(NE_O) = P2
		 P3_O(NE_O) = P3
		 P4_O(NE_O) = P4
		 NE_O = NE_O + 1
	      END DO		! IB
	   END DO		! IA
	END IF			! ISUB
	NE_O = NE_O - 1
	GO TO 1000
	
C+ Error returns
C  -----------------------------------------------------------------------------
900	WRITE (0,*) 'US-F-BNDERR, Boundary error'
        WRITE (0,*) '- energy array out of bounds; number of points',IE,
     1     ' is greater than ',E_SZ
	GO TO 999

 905	WRITE (0,*) 'US-F-BNDERR, Boundary error'
	WRITE (0,*) '- energy array out of bounds; number of points ',NE,
     1    ' is greater than ',E_SZ
	GO TO 999

 910	WRITE (0,*) '&US-E-HARMERR, Harmonic errror'
	WRITE (0,*) '- no harmonics reachable in the range ',EMINU,
     1    ' to ',EMAXU,' eV.'
	GO TO 999

 920	WRITE(0,*) '&US-E-INVDAT, Invalid data'
        WRITE(0,*) '- check input data file; center of pinhole must ',
     1    'lie in the first quadrant.'
	GO TO 999
	
 930	WRITE(0,*) '&US-E-INVDAT, Invalid data'
	WRITE(0,*) '- check input data file; method ',METHOD,
     1    ' not valid for the flux density distribution.'
	GO TO 999

940	WRITE(0,*) '&US-E-INVDAT, Invalid data'
	WRITE(0,*) '- check input data file; method ',METHOD,
     1    ' not valid for angle-integrated spectrum.'
	GO TO 999

999	WRITE(0,*)
	WRITE(0,*) '&US-F-PRGERR, Program error'
	WRITE(0,*) '- unsuccessful completion due to error status'
	WRITE(0,*)
	STOP
C-
C  -----------------------------------------------------------------------------
1000	END			! US

	SUBROUTINE SPACE_DISTRIBUTION(IERROR)
 
C[subroutine_header_comments]

C  Size parameters:
	INTEGER*4	E_SZ,A_SZ,B_SZ,P_SZ
	PARAMETER	(E_SZ=50001,A_SZ=400,B_SZ=A_SZ/4+1,P_SZ=501)

C  Declarations of scalars:
	INTEGER*4	IERROR,ICOUNT
	INTEGER*4	I,IA,IB,IC,IMIN,IMAX,IH1,IH2

	REAL*8		R,DA2,DALPHA,SL,CONST
	REAL*8		ALPHAI,ALPHA2I,ALPHAMIN,ALPHA2MIN,ALPHAMAX,ALPHA2MAX

C  Declarations of arrays:
	REAL*8		ALPHA(B_SZ),THETA(B_SZ)

C  Labeled constants:
	REAL*8		ZERO,ONE,TWO,EPS
	PARAMETER	(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,EPS=1.0D-2)

C  Common blocks:
	INTEGER*4	MODE,METHOD,IHARM,N,NSIGMA
	INTEGER*4	NXP,NYP
	INTEGER*4	NALPHA
	INTEGER*4	NE,NE1,NE2,NEU
	REAL*8		K,KX,KY,K3,NPI,PE,GAMMA,ER,EW,DEW,LEN
	REAL*8		FAC,C1,C2,C3,C4,C5
	REAL*8		AP2MIN,AP2MAX,AP2CNT,ARGMAX,FU,FV,SIGX2,SIGY2,GSIGUV
	REAL*8		DXP,DYP
	REAL*8		CALPHA2
	REAL*8		DE
	INTEGER*4	I1(E_SZ),I2(E_SZ)
	REAL*8		E(E_SZ),EU(E_SZ)
	REAL*8		SPEC0(E_SZ),SPEC1(E_SZ),SPEC2(E_SZ),SPEC3(E_SZ)
	REAL*8		BR0(B_SZ,B_SZ),BR1(B_SZ,B_SZ)
	REAL*8		BR2(B_SZ,B_SZ),BR3(B_SZ,B_SZ)
	REAL*8		RA0(P_SZ,P_SZ),RA1(P_SZ,P_SZ)
	REAL*8		RA2(P_SZ,P_SZ),RA3(P_SZ,P_SZ)
	REAL*8		XP(P_SZ),YP(P_SZ),CX(P_SZ),CY(P_SZ)

	COMMON		/CALCI/     MODE,METHOD,IHARM,N,NSIGMA
	COMMON		/CALCR/	    K,KX,KY,K3,NPI,PE,GAMMA,ER,EW,DEW,LEN
	COMMON		/FACTOR/    FAC,C1,C2,C3,C4,C5
	COMMON		/BEAM1/	    AP2MIN,AP2MAX,AP2CNT,ARGMAX,
     1			    FU,FV,SIGX2,SIGY2,GSIGUV
	COMMON		/PINHOLE/   XP,YP,CX,CY,NXP,NYP,DXP,DYP
	COMMON		/ANGLE_ALPI/ NALPHA
	COMMON		/ANGLE_ALPR/ CALPHA2
	COMMON		/ENERGY1/    E,EU,SPEC0,SPEC1,SPEC2,SPEC3
	COMMON		/ENERGY2/    I1,I2,NE,NE1,NE2,NEU,DE
	COMMON		/SPECTRA/   BR0,BR1,BR2,BR3,RA0,RA1,RA2,RA3

	IERROR = 0
	CONST  = C1

	I1(1) = 0
	I2(1) = 0

	R   = ER/E(1) ! Set at fixed energy

C+ Find range of harmonics that contributes at this energy
C  -----------------------------------------------------------------------------
	IF (METHOD .EQ. 1) THEN ! Finite-N
	    DA2  = CALPHA2*R/N
	    IMIN = (AP2MIN+K3-DA2)/R+ONE
	    IMAX = (AP2MAX+K3+DA2)/R
	ELSE ! Infinite-N
	    IMIN = (AP2MIN+K3)/R+ONE
	    IMAX = (AP2MAX+K3)/R
	END IF ! METHOD
	IF (IMAX .LT. IMIN) GO TO 900
        IF (IHARM .GT. 0 .AND. (IHARM .LT. IMIN .OR. IHARM .GT. IMAX))
     1   GO TO 910
C-
C  -----------------------------------------------------------------------------
	DO IB=1,NYP
	    DO IA=1,NXP
		RA0(IA,IB) = ZERO
		RA1(IA,IB) = ZERO
		RA2(IA,IB) = ZERO
		RA3(IA,IB) = ZERO
	    END DO ! IA
	END DO ! IB

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
	I1(1)  = IH1
	I2(1)  = IH2
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

C  Brilliance
	    CALL BRIGHTNESS_ARRAY(METHOD,I,R,ALPHAI,ALPHA2I,ALPHA)

C  Two-dimensional convolution of the brilliance with the electron distribution
	    CALL CONVOLUTE_DISTRIBUTION(METHOD,CONST,ALPHA,THETA,DALPHA,
     1				EPS,BR0,BR1,BR2,BR3,
     2				    RA0,RA1,RA2,RA3,ICOUNT)

	    IF (ICOUNT .GT. 1) THEN ! If higher harmonics do not contribute
		I2(1) = I
		GO TO 800
	    END IF ! ICOUNT
		
	END DO ! IH
C- Endloop harmonics
C  -----------------------------------------------------------------------------
800	CONTINUE
	RETURN

C+ Error returns
C  -----------------------------------------------------------------------------
900	CONTINUE
	WRITE(0,200) '&SPACE_DISTRIBUTION-E-HARMERR, Harmonic errror'
	WRITE(0,210) '- no harmonics reachable at ',E(1),' eV.'
	GO TO 999

910	CONTINUE
	WRITE(0,200) '&SPACE_DISTRIBUTION-E-HARMERR, Harmonic errror'
        WRITE(0,220) '- Harmonic number ',IHARM,' not in reachable ',
     1	   IMIN,' range to ',IMAX,' at ',E(1), ' eV.'
	GO TO 999

200	FORMAT(' ',8A)
210	FORMAT(' ',A,F10.3,A)
220	FORMAT(' ',3(A,I3),A,F10.3,A)

999	CONTINUE
	IERROR = -1
C-
C  -----------------------------------------------------------------------------
	RETURN
	END ! SPACE_DISTRIBUTION

	SUBROUTINE SPECTRAL_DISTRIBUTION(IERROR)
 
C[subroutine_header_comments]

C  Size parameters:
	INTEGER*4	E_SZ,A_SZ,B_SZ,P_SZ
	PARAMETER	(E_SZ=50001,A_SZ=400,B_SZ=A_SZ/4+1,P_SZ=501)

C  Declarations of scalars:
	LOGICAL*4	LE1,LE2
	INTEGER*4	IERROR,ICOUNT
	INTEGER*4	I,IA,IB,IC,IE,IMIN,IMAX,IH1,IH2

	REAL*8		R,DA2,DALPHA,SL,CONST
	REAL*8		ALPHAI,ALPHA2I,ALPHAMIN,ALPHA2MIN,ALPHAMAX,ALPHA2MAX
	REAL*8		AREA0,AREA1,AREA2,AREA3
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
	INTEGER*4	NXP,NYP
	INTEGER*4	NALPHA
	INTEGER*4	NE,NE1,NE2,NEU
	REAL*8		K,KX,KY,K3,NPI,PE,GAMMA,ER,EW,DEW,LEN
	REAL*8		FAC,C1,C2,C3,C4,C5
	REAL*8		AP2MIN,AP2MAX,AP2CNT,ARGMAX,FU,FV,SIGX2,SIGY2,GSIGUV
	REAL*8		DXP,DYP
	REAL*8		CALPHA2
	REAL*8		DE
	INTEGER*4	I1(E_SZ),I2(E_SZ)
	REAL*8		E(E_SZ),EU(E_SZ)
	REAL*8		SPEC0(E_SZ),SPEC1(E_SZ),SPEC2(E_SZ),SPEC3(E_SZ)
	REAL*8		BR0(B_SZ,B_SZ),BR1(B_SZ,B_SZ)
	REAL*8		BR2(B_SZ,B_SZ),BR3(B_SZ,B_SZ)
	REAL*8		RA0(P_SZ,P_SZ),RA1(P_SZ,P_SZ)
	REAL*8		RA2(P_SZ,P_SZ),RA3(P_SZ,P_SZ)
	REAL*8		XP(P_SZ),YP(P_SZ),CX(P_SZ),CY(P_SZ)

	COMMON		/CALCI/     MODE,METHOD,IHARM,N,NSIGMA
	COMMON		/CALCR/	    K,KX,KY,K3,NPI,PE,GAMMA,ER,EW,DEW,LEN
	COMMON		/FACTOR/    FAC,C1,C2,C3,C4,C5
	COMMON		/BEAM1/	    AP2MIN,AP2MAX,AP2CNT,ARGMAX,
     1			    FU,FV,SIGX2,SIGY2,GSIGUV
	COMMON		/PINHOLE/   XP,YP,CX,CY,NXP,NYP,DXP,DYP
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
            IF (IHARM .GT. 0 .AND. (IHARM .LT. IMIN .OR. 
     1        IHARM .GT. IMAX))
     2	        GO TO 810 ! Next energy
	    LE2 = .TRUE.
C-
C  -----------------------------------------------------------------------------
	    DO IB=1,NYP
		DO IA=1,NXP
		    RA0(IA,IB) = ZERO
		    RA1(IA,IB) = ZERO
		    RA2(IA,IB) = ZERO
		    RA3(IA,IB) = ZERO
		END DO ! IA
	    END DO ! IB

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

	        IF (MODE .EQ. 3) THEN ! Brilliance
		    IF (ALPHA2I .LT. ZERO) THEN ! may be neg for finite-N only
		        DEL = ZERO
		    ELSE                        ! Walker's approach (infinite-N)
		        DEL = ALPHA2I*N/R
		    END IF ! ALPHA2I

C                   Estimate diffraction limited source size
		    IF (DEL .LT. 2.15D0) THEN
		        SIGR2 = (1.29D0+(1.229D0*(DEL-0.8D0)*(DEL-0.8D0)))**2
		    ELSE
                        SIGR2 = 5.81D0*DEL
		    END IF ! DEL

		    SIGR2 = SIGR2*CV*C_EVANG/E(IE)*LEN
		    CONST = C1/(TWOPI*SQRT((SIGX2+SIGR2)*(SIGY2+SIGR2)))
	        END IF ! MODE

C  Brilliance
		CALL BRIGHTNESS_ARRAY(METHOD,I,R,ALPHAI,ALPHA2I,ALPHA)

C  Two-dimensional convolution of the brilliance with the electron distribution
		CALL CONVOLUTE_DISTRIBUTION(METHOD,CONST,ALPHA,THETA,DALPHA,
     1				    EPS,BR0,BR1,BR2,BR3,
     2				        RA0,RA1,RA2,RA3,ICOUNT)

		IF (ICOUNT .GT. 1) THEN ! If higher harmonics do not contribute
		    I2(IE) = I
		    GO TO 800
		END IF ! ICOUNT
	    END DO ! IH
C- Endloop harmonics
C  -----------------------------------------------------------------------------
800	    CONTINUE

C+ Save spectra
C  -----------------------------------------------------------------------------
	    IF (MODE .EQ. 4) THEN ! Integrate over pinhole
		CALL TRAPZ2(RA0,AREA0)
		CALL TRAPZ2(RA1,AREA1)
		CALL TRAPZ2(RA2,AREA2)
		CALL TRAPZ2(RA3,AREA3)
		SPEC0(IE) = FAC*AREA0
		SPEC1(IE) = FAC*AREA1
		SPEC2(IE) = FAC*AREA2
		SPEC3(IE) = FAC*AREA3
	    ELSE ! Fixed position or brilliance
		SPEC0(IE) = FAC*RA0(1,1)
		SPEC1(IE) = FAC*RA1(1,1)
		SPEC2(IE) = FAC*RA2(1,1)
		SPEC3(IE) = FAC*RA3(1,1)
	    END IF ! MODE
	    IF (FAC .EQ. 4.0D0) SPEC2(IE) = ZERO ! symmetry
C-
C  -----------------------------------------------------------------------------

810	    CONTINUE
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
	WRITE(0,200) '&SPECTRAL_DISTRIBUTION-E-HARMERR, Harmonic errror'
	WRITE(0,210) '- no harmonics reachable in the range ',E(NE1),
     1	  ' to ',E(NE2),' eV.'
	GO TO 999

910	CONTINUE
	WRITE(0,200) '&SPECTRAL_DISTRIBUTION-E-HARMERR, Harmonic errror'
        WRITE(0,220) '- Harmonic number ',IHARM,' not in the range ',
     1	  E(NE1),' to ',E(NE2),' eV.'
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

	SUBROUTINE ANGLE_INTEGRATION(IERROR)
 
C[subroutine_header_comments]

C  Size parameters:
	INTEGER*4	E_SZ,A_SZ,B_SZ,P_SZ
	PARAMETER	(E_SZ=50001,A_SZ=400,B_SZ=A_SZ/4+1,P_SZ=501)

C  Declarations of scalars:
	INTEGER*4	IERROR,ICOUNT
	INTEGER*4	I,ID,IE,IMIN

	REAL*8		R,SL,CONST
	REAL*8		ALPHAI,ALPHA2I
	REAL*8		SUM0,SUM1,SUM3

C  Labeled constants:
	REAL*8		ZERO,ONE,TWO,FOUR,HALF,EPS
	PARAMETER	(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
	PARAMETER	(FOUR=4.0D0,HALF=0.5D0,EPS=1.0D-2)

C  Common blocks:
	INTEGER*4	MODE,METHOD,IHARM,N,NSIGMA
	INTEGER*4	NPHI4,NPHI,NPHI_BRIGHT
	INTEGER*4	NE,NE1,NE2,NEU
	REAL*8		K,KX,KY,K3,NPI,PE,GAMMA,ER,EW,DEW,LEN
	REAL*8		FAC,C1,C2,C3,C4,C5
	REAL*8		DPHI
	REAL*8		DE
	INTEGER*4	INDEX_PHI(A_SZ)
	INTEGER*4	I1(E_SZ),I2(E_SZ)
	REAL*8		COSPHI(A_SZ),SINPHI(A_SZ),S2SIGN(A_SZ)
	REAL*8		E(E_SZ),EU(E_SZ)
	REAL*8		SPEC0(E_SZ),SPEC1(E_SZ),SPEC2(E_SZ),SPEC3(E_SZ)
	REAL*8		BR0(B_SZ,B_SZ),BR1(B_SZ,B_SZ)
	REAL*8		BR2(B_SZ,B_SZ),BR3(B_SZ,B_SZ)
	REAL*8		RA0(P_SZ,P_SZ),RA1(P_SZ,P_SZ)
	REAL*8		RA2(P_SZ,P_SZ),RA3(P_SZ,P_SZ)

	COMMON		/CALCI/     MODE,METHOD,IHARM,N,NSIGMA
	COMMON		/CALCR/	    K,KX,KY,K3,NPI,PE,GAMMA,ER,EW,DEW,LEN
	COMMON		/FACTOR/    FAC,C1,C2,C3,C4,C5
	COMMON		/ANGLE_PHII/ INDEX_PHI,NPHI4,NPHI,NPHI_BRIGHT
	COMMON		/ANGLE_PHIR/ COSPHI,SINPHI,S2SIGN,DPHI
	COMMON		/ENERGY1/    E,EU,SPEC0,SPEC1,SPEC2,SPEC3
	COMMON		/ENERGY2/    I1,I2,NE,NE1,NE2,NEU,DE
	COMMON		/SPECTRA/   BR0,BR1,BR2,BR3,RA0,RA1,RA2,RA3

	IERROR = 0
	SL     = ZERO

C+ Loop over energies
C  -----------------------------------------------------------------------------
	DO IE=1,NE

	    SPEC0(IE) = ZERO
	    SPEC1(IE) = ZERO
	    SPEC2(IE) = ZERO
	    SPEC3(IE) = ZERO
	    I1(IE)    = 0
	    I2(IE)    = 0

	    R     = ER/E(IE)
	    CONST = FOUR*C2*R/(TWO*N)*DPHI

	    IMIN = K3/R+ONE
	    IF (IHARM .GT. 0 .AND. IHARM .LT. IMIN) GO TO 810 ! Next energy

	    IF (IHARM .GT. 0) THEN
		I = IHARM
	    ELSE
		I = IMIN
	    END IF ! IHARM

	    I1(IE) = I
	    ICOUNT = 0

C+ Loop over harmonics
C  -----------------------------------------------------------------------------
	    DO WHILE (.TRUE.)
		ICOUNT  = ICOUNT+1
		ALPHA2I = R*I-K3
		ALPHAI  = SQRT(ALPHA2I)
C  Brilliance
		CALL BRIGHTNESS_ARRAY(METHOD,I,R,ALPHAI,ALPHA2I,SL)

C  Integrate
		SUM0 = ZERO
		SUM1 = ZERO
		SUM3 = ZERO
		DO ID=1,NPHI_BRIGHT
		    SUM0 = SUM0+BR0(ID,1)
		    SUM1 = SUM1+BR1(ID,1)
		    SUM3 = SUM3+BR3(ID,1)
		END DO ! NPHI_BRIGHT

		SUM0 = SUM0*CONST
		SUM1 = SUM1*CONST
		SUM3 = SUM3*CONST
		SPEC0(IE) = SPEC0(IE)+SUM0
		SPEC1(IE) = SPEC1(IE)+SUM1
		SPEC3(IE) = SPEC3(IE)+SUM3

		IF (IHARM .NE. 0) GO TO 800

		IF (SUM0 .GT. EPS*SPEC0(IE)) ICOUNT = 0
		IF (ICOUNT .GT. 1) GO TO 800 ! If higher harmonics do not contribute

		I = I+1
	    END DO ! WHILE
C- Endloop harmonics
C  -----------------------------------------------------------------------------
800	    CONTINUE
	    I2(IE) = I

810	    CONTINUE
	END DO ! IE
C- Endloop energy
C  -----------------------------------------------------------------------------
	
	IF (METHOD .EQ. 4) THEN ! Infinite-N with convolution
	    CALL CONVOLUTE_ENERGY_VSTEP(SPEC0)
	    CALL CONVOLUTE_ENERGY_VSTEP(SPEC1)
	    CALL CONVOLUTE_ENERGY_VSTEP(SPEC3)
	ELSE IF (METHOD .EQ. 14) THEN ! Infinite-N with convolution (Walker's)
	    CALL CONVOLUTE_ENERGY_ESTEP(SPEC0)
	    CALL CONVOLUTE_ENERGY_ESTEP(SPEC1)
	    CALL CONVOLUTE_ENERGY_ESTEP(SPEC3)
	END IF 

	RETURN

C+ Error returns
C  -----------------------------------------------------------------------------
C900	CONTINUE
C	GO TO 999
C
C999	CONTINUE
C	IERROR = -1
C-
C  -----------------------------------------------------------------------------
C	RETURN
	END ! ANGLE_INTEGRATION

	SUBROUTINE POWER_DISTRIBUTION(IERROR)
     
C[subroutine_header_comments]

C  Size parameters:
	INTEGER*4	E_SZ,A_SZ,B_SZ,P_SZ
	PARAMETER	(E_SZ=50001,A_SZ=400,B_SZ=A_SZ/4+1,P_SZ=501)

C  Declarations of scalars:
	INTEGER*4	IERROR,ICOUNT
	INTEGER*4	I,IA,IB,IC

	REAL*8		DALPHA,SL,CONST,DELTAP,PTOT
	REAL*8		ALPHAP,ALPHA2,ALPHAMIN,ALPHAMAX
	REAL*8		XG,YG,COSPHI,SINPHI
	REAL*8		S0,S1,S2,S3,AXR,AYR,AXI,AYI
	REAL*8		AREA0,DELTA0,DELTA1,DELTA2,DELTA3

C  Declarations of arrays:
	REAL*8		ALPHA(B_SZ),THETA(B_SZ)

C  Labeled constants:
	REAL*8		ZERO,ONE,TWO,EPS,EPSP
	PARAMETER	(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,EPS=1.0D-6,EPSP=2.0D-3)

C  Common blocks:
	INTEGER*4	MODE,METHOD,IHARM,N,NSIGMA
	INTEGER*4	NXP,NYP
	INTEGER*4	NALPHA
	REAL*8		K,KX,KY,K3,NPI,PE,GAMMA,ER,EW,DEW,LEN
	REAL*8		FAC,C1,C2,C3,C4,C5
	REAL*8		AP2MIN,AP2MAX,AP2CNT,ARGMAX,FU,FV,SIGX2,SIGY2,GSIGUV
	REAL*8		DXP,DYP
	REAL*8		CALPHA2
	REAL*8		DE
	INTEGER*4	NE,NE1,NE2,NEU
	INTEGER*4	I1(E_SZ),I2(E_SZ)
	REAL*8		BR0(B_SZ,B_SZ),BR1(B_SZ,B_SZ)
	REAL*8		BR2(B_SZ,B_SZ),BR3(B_SZ,B_SZ)
	REAL*8		RA0(P_SZ,P_SZ),RA1(P_SZ,P_SZ)
	REAL*8		RA2(P_SZ,P_SZ),RA3(P_SZ,P_SZ)
	REAL*8		XP(P_SZ),YP(P_SZ),CX(P_SZ),CY(P_SZ)

	COMMON		/CALCI/     MODE,METHOD,IHARM,N,NSIGMA
	COMMON		/CALCR/	    K,KX,KY,K3,NPI,PE,GAMMA,ER,EW,DEW,LEN
	COMMON		/FACTOR/    FAC,C1,C2,C3,C4,C5
	COMMON		/BEAM1/	    AP2MIN,AP2MAX,AP2CNT,ARGMAX,
     1			    FU,FV,SIGX2,SIGY2,GSIGUV
	COMMON		/PINHOLE/   XP,YP,CX,CY,NXP,NYP,DXP,DYP
	COMMON		/ANGLE_ALPI/ NALPHA
	COMMON		/ANGLE_ALPR/ CALPHA2
	COMMON		/ENERGY2/    I1,I2,NE,NE1,NE2,NEU,DE
	COMMON		/SPECTRA/   BR0,BR1,BR2,BR3,RA0,RA1,RA2,RA3

	IERROR = 0
	PTOT   = ZERO
        DO IB=1,NYP
            DO IA=1,NXP
                RA0(IA,IB) = ZERO
	        RA1(IA,IB) = ZERO
	        RA2(IA,IB) = ZERO
	        RA3(IA,IB) = ZERO
	    END DO ! IA
	END DO ! IB

	IF (METHOD .NE. 3) THEN ! non-zero emittance; generate alpha,theta
	    CONST    = C5
	    ALPHAMAX = SQRT(AP2MAX)
	    ALPHAMIN = SQRT(AP2MIN)
	    DALPHA   = (ALPHAMAX-ALPHAMIN)/NALPHA
	    SL       = ALPHAMIN+DALPHA/TWO
	    DO IC=1,NALPHA
	        ALPHA(IC) = SL+(IC-1)*DALPHA
		THETA(IC) = ALPHA(IC)/GAMMA
	    END DO ! IC
	ELSE ! zero emittance
	   CONST = C4
	END IF ! METHOD

	IF (IHARM .GT. 0) THEN
	    I = IHARM
	ELSE
	    I = 1
	END IF ! IHARM
	I1(1)  = I
	ICOUNT = 0

C+ Loop over harmonics
C  -----------------------------------------------------------------------------
	DO WHILE (.TRUE.)
	    ICOUNT = ICOUNT+1
	    IF (METHOD .NE. 3) THEN 
C           ! non-zero emittance; convolve with e-distr.
C  Power density
	        CALL PDF_ARRAY(I,ALPHA)

C  Two-dimensional convolution of the pdf with the electron distribution
		CALL CONVOLUTE_DISTRIBUTION(1,CONST,ALPHA,THETA,DALPHA,
     1				    EPSP,BR0,BR1,BR2,BR3,
     2				         RA0,RA1,RA2,RA3,ICOUNT)

	    ELSE ! zero emittance; direct calculation
	        DO IB=1,NYP
	            DO IA=1,NXP
		        XG     = GAMMA*CX(IA)
		        YG     = GAMMA*CY(IB)
		        ALPHA2 = XG*XG+YG*YG
		        ALPHAP = SQRT(ALPHA2)
		        IF (ALPHAP .LT. EPS) THEN
		            COSPHI = ZERO
			    SINPHI = ONE
		        ELSE
		            COSPHI = XG/ALPHAP
		            SINPHI = YG/ALPHAP
		        END IF ! ALPHAP
		        CALL BRIGHTE(I,ALPHAP,COSPHI,SINPHI,
     1			     S0,S1,S2,S3,AXR,AYR,AXI,AYI)
			DELTA0 = CONST*S0/(K3+ALPHA2)
			DELTA1 = CONST*S1/(K3+ALPHA2)
			DELTA2 = CONST*S2/(K3+ALPHA2)
			DELTA3 = CONST*S3/(K3+ALPHA2)
			RA0(IA,IB) = RA0(IA,IB)+DELTA0
			RA1(IA,IB) = RA1(IA,IB)+DELTA1
			RA2(IA,IB) = RA2(IA,IB)+DELTA2
			RA3(IA,IB) = RA3(IA,IB)+DELTA3
			IF (DELTA0 .GT. EPSP*RA0(IA,IB)) ICOUNT = 0
		    END DO ! IA
		END DO ! IB
	    END IF ! METHOD

	    CALL TRAPZ2(RA0,AREA0) 
C         ! integrated power for harmonics including i (W)
	    DELTAP = FAC*AREA0-PTOT
C         ! incremental change in total power (W)
	    IF (DELTAP .GT. EPSP*PTOT) ICOUNT = 0
	    PTOT = PTOT+DELTAP

	    IF (IHARM .GT. 0) GO TO 800
	    IF (IHARM .LT. 0 .AND. I .GE. -IHARM) GO TO 800
	    IF (IHARM .EQ. 0 .AND. ICOUNT .GT. 1) GO TO 800

	    I = I+1
	END DO ! WHILE
C+ Endloop harmonics
C  -----------------------------------------------------------------------------
800	CONTINUE
        I2(1) = I

	RETURN
	END ! POWER_DISTRIBUTION

	SUBROUTINE NO_EMITTANCE(IERROR)
 
C[subroutine_header_comments]

C  Size parameters:
	INTEGER*4	E_SZ,A_SZ,B_SZ,P_SZ
	PARAMETER	(E_SZ=50001,A_SZ=400,B_SZ=A_SZ/4+1,P_SZ=501)

C  Declarations of scalars:
	LOGICAL*4	LE
	INTEGER*4	IERROR
	INTEGER*4	I,IA,IB,IE

	REAL*8		R,DA2,CONST,ARG
	REAL*8		ALPHA,ALPHA2,ALPHA2I,ALPHA2MIN,ALPHA2MAX,SR2
	REAL*8		XG,YG,COSPHI,SINPHI
	REAL*8		S0,S1,S2,S3,AXR,AYR,AXI,AYI
	REAL*8		AREA0,AREA1,AREA2,AREA3
	REAL*8		SINC_DEJUS

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
	PARAMETER	(C_EVANG=H*C/EC*1.0D10,C_ANG_M=1.0D-10)
	PARAMETER	(C_M2_MM2=1.0D6)

C  Labeled constants:
	REAL*8		PI,PIHALF,TWOPI
	PARAMETER	(PI    =3.1415 92653 58979 32384 62643D0)
	PARAMETER	(PIHALF=1.5707 96326 79489 66192 31322D0)
	PARAMETER	(TWOPI= 6.2831 85307 17958 64769 25287D0)
	REAL*8		ZERO,ONE,TWO,EPS
	PARAMETER	(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,EPS=1.0D-6)
	REAL*8		CV,CW
	PARAMETER	(CV=ONE/8.0D0/PI/PI*C_ANG_M*C_M2_MM2,CW=3.57D0)!1.2665D-6

C  Common blocks:
	INTEGER*4	MODE,METHOD,IHARM,N,NSIGMA
	INTEGER*4	NXP,NYP
	INTEGER*4	NALPHA
	INTEGER*4	NE,NE1,NE2,NEU
	REAL*8		K,KX,KY,K3,NPI,PE,GAMMA,ER,EW,DEW,LEN
	REAL*8		FAC,C1,C2,C3,C4,C5
	REAL*8		DXP,DYP
	REAL*8		CALPHA2
	REAL*8		DE
	INTEGER*4	I1(E_SZ),I2(E_SZ)
	REAL*8		E(E_SZ),EU(E_SZ)
	REAL*8		SPEC0(E_SZ),SPEC1(E_SZ),SPEC2(E_SZ),SPEC3(E_SZ)
	REAL*8		BR0(B_SZ,B_SZ),BR1(B_SZ,B_SZ)
	REAL*8		BR2(B_SZ,B_SZ),BR3(B_SZ,B_SZ)
	REAL*8		RA0(P_SZ,P_SZ),RA1(P_SZ,P_SZ)
	REAL*8		RA2(P_SZ,P_SZ),RA3(P_SZ,P_SZ)
	REAL*8		XP(P_SZ),YP(P_SZ),CX(P_SZ),CY(P_SZ)

	COMMON		/CALCI/     MODE,METHOD,IHARM,N,NSIGMA
	COMMON		/CALCR/	    K,KX,KY,K3,NPI,PE,GAMMA,ER,EW,DEW,LEN
	COMMON		/FACTOR/    FAC,C1,C2,C3,C4,C5
	COMMON		/PINHOLE/   XP,YP,CX,CY,NXP,NYP,DXP,DYP
	COMMON		/ANGLE_ALPI/ NALPHA
	COMMON		/ANGLE_ALPR/ CALPHA2
	COMMON		/ENERGY1/    E,EU,SPEC0,SPEC1,SPEC2,SPEC3
	COMMON		/ENERGY2/    I1,I2,NE,NE1,NE2,NEU,DE
	COMMON		/SPECTRA/   BR0,BR1,BR2,BR3,RA0,RA1,RA2,RA3

	IERROR = 0
	LE     = .FALSE.

C+ Loop over energies
C  -----------------------------------------------------------------------------
	DO IE=1,NE

	    SPEC0(IE) = ZERO
	    SPEC1(IE) = ZERO
	    SPEC2(IE) = ZERO
	    SPEC3(IE) = ZERO
	    I1(IE)    = 0
	    I2(IE)    = 0

	    IF (MODE .EQ. 3) THEN ! Brilliance
		SR2   = CW*CV*C_EVANG/E(IE)*LEN
		CONST = C3/(TWOPI*SR2)
	    ELSE
		CONST = C3
	    END IF ! MODE

	    R   = ER/E(IE)
	    DA2 = CALPHA2*R/N

	    DO IB=1,NYP
		DO IA=1,NXP
		    RA0(IA,IB) = ZERO
		    RA1(IA,IB) = ZERO
		    RA2(IA,IB) = ZERO
		    RA3(IA,IB) = ZERO

		    XG     = GAMMA*CX(IA)
		    YG     = GAMMA*CY(IB)
		    ALPHA2 = XG*XG+YG*YG

		    IF (IHARM .GT. 0) THEN
			I	  = IHARM
			ALPHA2I   = R*I-K3
			ALPHA2MAX = ALPHA2I+DA2
			IF (ALPHA2MAX .LT. ALPHA2) GO TO 800
		    ELSE ! Find harmonic #
			I         = 1
			ALPHA2I   = R*I-K3
			ALPHA2MAX = ALPHA2I+DA2
			DO WHILE (ALPHA2MAX .LT. ALPHA2)
			    I         = I+1
			    ALPHA2I   = R*I-K3
			    ALPHA2MAX = ALPHA2I+DA2
			END DO ! WHILE
		    END IF ! IHARM
			
		    ALPHA2MIN = ALPHA2I-DA2
		    IF (ALPHA2MIN .LE. ALPHA2) THEN
			LE = .TRUE.
			IF (I1(IE) .EQ. 0) THEN
			    I1(IE) = I
			ELSE IF (I .LT. I1(IE)) THEN
			    I1(IE) = I
			END IF ! I1
			I2(IE) = I
			ALPHA  = SQRT(ALPHA2)
			IF (ALPHA .LT. EPS) THEN
			    COSPHI = ZERO
			    SINPHI = ONE
			ELSE
			    COSPHI = XG/ALPHA
			    SINPHI = YG/ALPHA
			END IF ! ALPHA
	        	CALL BRIGHTE(I,ALPHA,COSPHI,SINPHI,
     1	   	             S0,S1,S2,S3,AXR,AYR,AXI,AYI)
			ARG        = NPI*(ALPHA2-ALPHA2I)/R
			RA0(IA,IB) = CONST*SINC_DEJUS(ARG)*S0
			RA1(IA,IB) = CONST*SINC_DEJUS(ARG)*S1
			RA2(IA,IB) = CONST*SINC_DEJUS(ARG)*S2
			RA3(IA,IB) = CONST*SINC_DEJUS(ARG)*S3
		    END IF ! ALPHA2MIN

800		    CONTINUE
		END DO ! IA
	    END DO ! IB

C+ Save spectra
C  -----------------------------------------------------------------------------
	    IF (MODE .EQ. 4) THEN ! Integrate over pinhole
		CALL TRAPZ2(RA0,AREA0)
		CALL TRAPZ2(RA1,AREA1)
		CALL TRAPZ2(RA2,AREA2)
		CALL TRAPZ2(RA3,AREA3)
		SPEC0(IE) = FAC*AREA0
		SPEC1(IE) = FAC*AREA1
		SPEC2(IE) = FAC*AREA2
		SPEC3(IE) = FAC*AREA3
		IF (FAC .EQ. 4.0D0) SPEC2(IE) = ZERO
	    ELSE IF (MODE .EQ. 2 .OR. MODE .EQ. 3) THEN 
C            ! Fixed position or brilliance
		SPEC0(IE) = RA0(1,1)
		SPEC1(IE) = RA1(1,1)
		SPEC2(IE) = RA2(1,1)
		SPEC3(IE) = RA3(1,1)
	    END IF ! MODE (else MODE = 1)
C-
C  -----------------------------------------------------------------------------

	END DO ! IE
C- Endloop energy
C  -----------------------------------------------------------------------------

	IF (.NOT. LE .AND. MODE .EQ. 1) GO TO 900
	IF (.NOT. LE .AND. MODE .NE. 1) GO TO 910

	RETURN

C+ Error returns
C  -----------------------------------------------------------------------------
900	CONTINUE
	WRITE(0,200) '&NO_EMITTANCE-E-HARMERR, Harmonic errror' 
	IF (IHARM .GT. 0) THEN
	    WRITE(0,220) '- Harmonic number ',IHARM,' not reachable at ',
     1             E(1),' eV.'
	ELSE
	    WRITE(0,210) '- no harmonics reachable at ',E(1),' eV.'
	END IF ! IHARM
	GO TO 999

910	CONTINUE
	WRITE(0,200) '&NO_EMITTANCE-E-HARMERR, Harmonic errror' 
	IF (IHARM .GT. 0) THEN
	    WRITE(0,220) '- Harmonic number ',IHARM,' not in the range ',
     1             E(NE1),' to ',E(NE2),' eV.'
	ELSE
	    WRITE(0,210) '- no harmonics reachable in the range ',E(NE1),
     1	      ' to ',E(NE2),' eV.'
	END IF ! IHARM
	GO TO 999

200	FORMAT(' ',8A)
210	FORMAT(' ',A,2(F10.3,A))
220	FORMAT(' ',A,I3,A,2(F10.3,A))

999	CONTINUE
	IERROR = -1
C-
C  -----------------------------------------------------------------------------
	RETURN
	END ! NO_EMITTANCE

	SUBROUTINE BRIGHTNESS_ARRAY(ICALC,I,R,ALPHAI,ALPHA2I,ALPHA)
 
C[subroutine_header_comments]

C  Size parameters:
	INTEGER*4	E_SZ,A_SZ,B_SZ,P_SZ
	PARAMETER	(E_SZ=50001,A_SZ=400,B_SZ=A_SZ/4+1,P_SZ=501)

C  Declarations of scalars:
	INTEGER*4	ICALC,I,IC,ID
	REAL*8		ALPHAI,ALPHA2I,R,ARG
	REAL*8		ALPHA2,S0,S1,S2,S3,AXR,AYR,AXI,AYI
	REAL*8		H,SINC_DEJUS

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
	REAL*8		RA0(P_SZ,P_SZ),RA1(P_SZ,P_SZ)
	REAL*8		RA2(P_SZ,P_SZ),RA3(P_SZ,P_SZ)

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
		H      = SINC_DEJUS(ARG)
		DO ID=1,NPHI_BRIGHT
		    CALL BRIGHTE(I,ALPHA(IC),COSPHI(ID),SINPHI(ID),
     1			 S0,S1,S2,S3,AXR,AYR,AXI,AYI)
		    BR0(ID,IC) = H*S0
		    BR1(ID,IC) = H*S1
		    BR2(ID,IC) = H*S2
		    BR3(ID,IC) = H*S3
		END DO ! ID
	    END DO ! IC
	ELSE ! Infinite-N approximation (H=1.0)
	    DO ID=1,NPHI_BRIGHT
	        CALL BRIGHTE(I,ALPHAI,COSPHI(ID),SINPHI(ID),
     1	   	     S0,S1,S2,S3,AXR,AYR,AXI,AYI)
		BR0(ID,1) = S0
		BR1(ID,1) = S1
		BR2(ID,1) = S2
		BR3(ID,1) = S3
	    END DO ! ID
	END IF ! ICALC

	RETURN
	END ! BRIGHTNESS_ARRAY

	SUBROUTINE PDF_ARRAY(I,ALPHA)
 
C[subroutine_header_comments]

C  Size parameters:
	INTEGER*4	E_SZ,A_SZ,B_SZ,P_SZ
	PARAMETER	(E_SZ=50001,A_SZ=400,B_SZ=A_SZ/4+1,P_SZ=501)

C  Declarations of scalars:
	INTEGER*4	I,IC,ID
	REAL*8		ALPHA2
	REAL*8		S0,S1,S2,S3,AXR,AYR,AXI,AYI

C  Declarations of arrays:
	REAL*8		ALPHA(*)

C  Labeled constants:
	REAL*8		ZERO
	PARAMETER	(ZERO=0.0D0)

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
	REAL*8		RA0(P_SZ,P_SZ),RA1(P_SZ,P_SZ)
	REAL*8		RA2(P_SZ,P_SZ),RA3(P_SZ,P_SZ)

	COMMON		/CALCR/	    K,KX,KY,K3,NPI,PE,GAMMA,ER,EW,DEW,LEN
	COMMON		/ANGLE_PHII/ INDEX_PHI,NPHI4,NPHI,NPHI_BRIGHT
	COMMON		/ANGLE_PHIR/ COSPHI,SINPHI,S2SIGN,DPHI
	COMMON		/ANGLE_ALPI/ NALPHA
	COMMON		/ANGLE_ALPR/ CALPHA2
	COMMON		/SPECTRA/   BR0,BR1,BR2,BR3,RA0,RA1,RA2,RA3

        DO IC=1,NALPHA
	    ALPHA2 = ALPHA(IC)*ALPHA(IC)
	    DO ID=1,NPHI_BRIGHT
	        CALL BRIGHTE(I,ALPHA(IC),COSPHI(ID),SINPHI(ID),
     1	   	     S0,S1,S2,S3,AXR,AYR,AXI,AYI)
		BR0(ID,IC) = S0/(K3+ALPHA2)
		BR1(ID,IC) = S1/(K3+ALPHA2)
		BR2(ID,IC) = S2/(K3+ALPHA2)
		BR3(ID,IC) = S3/(K3+ALPHA2)

c		BR1(ID,IC) = ZERO
c		BR2(ID,IC) = ZERO
c		BR3(ID,IC) = ZERO
	    END DO ! ID
        END DO ! IC

	RETURN
	END ! PDF_ARRAY

	SUBROUTINE CONVOLUTE_DISTRIBUTION(ICALC,CONST,ALPHA,THETA,DALPHA,
     1				  EPS,BR0,BR1,BR2,BR3,
     2				      RA0,RA1,RA2,RA3,ICOUNT)
 
C[subroutine_header_comments]
 
C  Size parameters:
	INTEGER*4	A_SZ,B_SZ,P_SZ
	PARAMETER	(A_SZ=400,B_SZ=A_SZ/4+1,P_SZ=501)

C  Declarations of scalars:
	INTEGER*4	ICALC,ICOUNT,IA,IB,IC,ID
	REAL*8		CONST,DALPHA,EPS
	REAL*8		U,V,ARG,SL,P
	REAL*8		SUM0,SUM1,SUM2,SUM3
	REAL*8		DELTA0,DELTA1,DELTA2,DELTA3
 
C  Declarations of arrays:
	REAL*8		ALPHA(*),THETA(*)
	REAL*8		BR0(B_SZ,*),BR1(B_SZ,*)
	REAL*8		BR2(B_SZ,*),BR3(B_SZ,*)
	REAL*8		RA0(P_SZ,*),RA1(P_SZ,*)
	REAL*8		RA2(P_SZ,*),RA3(P_SZ,*)

C  Labeled constants:
	REAL*8		ZERO
	PARAMETER	(ZERO=0.0D0)
 
C  Common blocks:
	INTEGER*4	NXP,NYP
	INTEGER*4	NPHI4,NPHI,NPHI_BRIGHT
	INTEGER*4	NALPHA
	REAL*8		AP2MIN,AP2MAX,AP2CNT,ARGMAX,FU,FV,SIGX2,SIGY2,GSIGUV
	REAL*8		DXP,DYP
	REAL*8		DPHI
	REAL*8		CALPHA2
	INTEGER*4	INDEX_PHI(A_SZ)
	REAL*8		COSPHI(A_SZ),SINPHI(A_SZ),S2SIGN(A_SZ)
	REAL*8		XP(P_SZ),YP(P_SZ),CX(P_SZ),CY(P_SZ)

	COMMON		/BEAM1/	    AP2MIN,AP2MAX,AP2CNT,ARGMAX,
     1			    FU,FV,SIGX2,SIGY2,GSIGUV
	COMMON		/PINHOLE/   XP,YP,CX,CY,NXP,NYP,DXP,DYP
	COMMON		/ANGLE_PHII/ INDEX_PHI,NPHI4,NPHI,NPHI_BRIGHT
	COMMON		/ANGLE_PHIR/ COSPHI,SINPHI,S2SIGN,DPHI
	COMMON		/ANGLE_ALPI/ NALPHA
	COMMON		/ANGLE_ALPR/ CALPHA2

	SL = CONST*DALPHA*DPHI
	IF (ICALC .EQ. 1) THEN ! Finite-N calculation
	  DO IB=1,NYP
	    DO IA=1,NXP
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
		  U   = CX(IA)-THETA(IC)*COSPHI(ID)
		  V   = CY(IB)-THETA(IC)*SINPHI(ID)
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
	      RA0(IA,IB) = RA0(IA,IB)+SUM0
	      RA1(IA,IB) = RA1(IA,IB)+SUM1
	      RA2(IA,IB) = RA2(IA,IB)+SUM2
	      RA3(IA,IB) = RA3(IA,IB)+SUM3
	      IF (SUM0 .GT. EPS*RA0(IA,IB)) ICOUNT = 0
	    END DO ! IA
	  END DO ! IB
	ELSE ! Infinite-N approximation (no alpha integ.)
	  DO IB=1,NYP
	    DO IA=1,NXP
	      DELTA0 = ZERO
	      DELTA1 = ZERO
	      DELTA2 = ZERO
	      DELTA3 = ZERO
	      DO ID=1,NPHI4
		U   = CX(IA)-THETA(1)*COSPHI(ID)
		V   = CY(IB)-THETA(1)*SINPHI(ID)
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
	      RA0(IA,IB) = RA0(IA,IB)+SUM0
	      RA1(IA,IB) = RA1(IA,IB)+SUM1
	      RA2(IA,IB) = RA2(IA,IB)+SUM2
	      RA3(IA,IB) = RA3(IA,IB)+SUM3
	      IF (SUM0 .GT. EPS*RA0(IA,IB)) ICOUNT = 0
	    END DO ! IA
	  END DO ! IB
	END IF ! ICALC
	
	RETURN
	END ! CONVOLUTE_DISTRIBUTION

	SUBROUTINE CONVOLUTE_ENERGY_ESTEP(SP1)
 
C[subroutine_header_comments]
C  Size parameters:
	INTEGER*4	E_SZ
	PARAMETER	(E_SZ=50001)

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
	PARAMETER	(E_SZ=50001)

C  Declarations of scalars:
	INTEGER*4	IE,IEU,J1,J2
	REAL*8		EUP,EU1,EU2
	REAL*8		ARG,P,H1,H2,S1,S2,SUM
	REAL*8		SINC_DEJUS
 
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
		    SUM = S2*SINC_DEJUS(ARG)*H2/TWO
		    ARG = PE*(EUP-E(J2))
		    SUM = SUM+SP1(J2)*SINC_DEJUS(ARG)*(H2+HE(J2-1)/TWO)
		    DO IE=1,J2-1
			ARG = PE*(EUP-E(IE))
			SUM = SUM+SP1(IE)*SINC_DEJUS(ARG)*HA(IE)
		    END DO ! IE
		    SP2(IEU) = SUM/DEW
		ELSE ! J2 = NE
		    SUM = ZERO		    
		    DO IE=1,NE
			ARG = PE*(EUP-E(IE))
			SUM = SUM+SP1(IE)*SINC_DEJUS(ARG)*HA(IE)
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
		    SUM = S1*SINC_DEJUS(ARG)*H1/TWO
		    ARG = PE*(EUP-E(J1+1))
		    SUM = SUM+SP1(J1+1)*SINC_DEJUS(ARG)*(H1+HE(J1+1))/TWO
C  Interpolate at upper limit
		    H2  = EU2-E(J2)
		    P   = H2/HE(J2)
		    S2  = (ONE-P)*SP1(J2)+P*SP1(J2+1)
		    ARG = PE*(EUP-EU2)
		    SUM = SUM+S2*SINC_DEJUS(ARG)*H2/TWO
		    ARG = PE*(EUP-E(J2))
		    SUM = SUM+SP1(J2)*SINC_DEJUS(ARG)*(H2+HE(J2-1)/TWO)
		    DO IE=J1+2,J2-1
			ARG = PE*(EUP-E(IE))
			SUM = SUM+SP1(IE)*SINC_DEJUS(ARG)*HA(IE)
		    END DO ! IE
		    SP2(IEU) = SUM/DEW
		END IF

	    ELSE ! J2 = NE
C  Interpolate at lower limit
		H1  = E(J1+1)-EU1
		P   = H1/HE(J1)
		S1  = P*SP1(J1)+(ONE-P)*SP1(J1+1)
		ARG = PE*(EUP-EU1)
		SUM = S1*SINC_DEJUS(ARG)*H1/TWO
		ARG = PE*(EUP-E(J1+1))
		SUM = SUM+SP1(J1+1)*SINC_DEJUS(ARG)*(H1+HE(J1+1))/TWO
		DO IE=J1+2,NE
		    ARG = PE*(EUP-E(IE))
		    SUM = SUM+SP1(IE)*SINC_DEJUS(ARG)*HA(IE)
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

	REAL*8 FUNCTION SINC_DEJUS(ARG)
C+
C 
C FUNCTIONAL DESCRIPTION:	
C 
C  Routine to calculate the SINC_DEJUS function.
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
C  The SINC_DEJUS value
C 
C FORMAL PARAMETERS:
C  
C  Input arguments:
C  ALPHA2		Square of gamma*theta; Theta is the polar angle [rad]
C  ALPHA2I		Value of ALPHA2 for peak of the i-th harmonic off-axis
C  R			(1+K**2/2)/(E/E1(0)); K is the deflection parameter;
C			E is the energy, and E1(0) is the energy for the first
C			harmonic on-axis.
C  
C COMMON BLOCKS:
C  
C  /CALCI/     Contains general integer parameters; unchanged on exit
C  /CALCR/     Contains general real    parameters; unchanged on exit
C  
C MODIFICATION HISTORY:
C 
C	 Date     | Name  | Description
C ----------------+-------+-----------------------------------------------------
C [change_entry]
C-
C  Declarations of scalars:
	REAL*8		ARG,SL
 
C  Labeled constants:
	REAL*8		ONE,EPS
	PARAMETER	(ONE=1.0D0,EPS=1.0D-6)
	 
	IF (ABS(ARG) .GT. EPS) THEN
	    SL   = SIN(ARG)/ARG
	    SINC_DEJUS = SL*SL
	ELSE
	    SINC_DEJUS = ONE
	END IF ! ARG

	RETURN
	END ! SINC_DEJUS

	SUBROUTINE TRAPZ2(RA,AREA)
 
C[subroutine_header_comments]
C  Size parameters:
	INTEGER*4	A_SZ,B_SZ,P_SZ
	PARAMETER	(A_SZ=400,B_SZ=A_SZ/4+1,P_SZ=501)

C  Declarations of scalars:
	INTEGER*4	IA,IB
	REAL*8		AREA,SUM,WX,WY
 
C  Declarations of arrays:
	REAL*8		RA(P_SZ,*)

C  Labeled constants:
	REAL*8		ZERO,ONE,HALF
	PARAMETER	(ZERO=0.0D0,ONE=1.0D0,HALF=0.5D0)

C  Common blocks:
	INTEGER*4	NXP,NYP
	REAL*8		DXP,DYP
	REAL*8		XP(P_SZ),YP(P_SZ),CX(P_SZ),CY(P_SZ)

	COMMON		/PINHOLE/   XP,YP,CX,CY,NXP,NYP,DXP,DYP

	SUM = ZERO
	DO IB=1,NYP
	    IF (IB .EQ. 1 .OR. IB .EQ. NYP) THEN
		WY = HALF
	    ELSE
		WY = ONE
	    END IF ! IB
	    DO IA=1,NXP
		IF (IA .EQ. 1 .OR. IA .EQ. NXP) THEN
		    WX = HALF
		ELSE
		    WX = ONE
		END IF ! IA
	    SUM = SUM+WX*WY*RA(IA,IB)
	    END DO ! IA
	END DO ! IB
	AREA = SUM*DXP*DYP

	RETURN
	END ! TRAPZ2
