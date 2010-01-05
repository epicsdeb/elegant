	SUBROUTINE HUNT(X,N,XP,JLO)
C+
C 
C FUNCTIONAL DESCRIPTION:	
C 
C  Subroutine to find the index JLO in an ordered array X with N elements.
C  The index JLO satisfies: X(JLO) < XP < X(JLO+1).
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
C  21-NOV-1990
C 
C FORMAL PARAMETERS:
C  
C  Input arguments:
C  X			Array of ordered values (increasing or decreasing)
C  N			Number of elements in X
C  XP			Value XP is searched in X
C  JLO	  		Initial guess of JLO (also output; see below)
C
C  Output arguments:
C  JLO	  		Final value of JLO; When JLO = 0 or JLO = N the value
C			of XP is outside the range of array X
C
C DESIGN ISSUES:
C  
C  Adapted from "Numerical Recipes" by W.H. Press, B.P. Flannery, S.A. Teukolsky
C  and W.T. Vetterling, Cambridge University Press (1988), p. 91.
C  
C KEYWORDS:
C  
C  Search, ordered table
C  
C [optional_header_tags]...
C 
C MODIFICATION HISTORY:
C 
C	 Date     | Name  | Description
C ----------------+-------+-----------------------------------------------------
C [change_entry]
C-

C  Declarations of scalars:
	LOGICAL*4	ASCEND
	INTEGER*4	N,JLO,JHI,JM,INC
	REAL*8		XP

C  Declarations of arrays:
	REAL*8		X(N)

	ASCEND = X(N) .GT. X(1) ! Indicates whether increasing or decreasing
				! values

	IF (JLO .LT. 1 .OR. JLO .GT. N) THEN
	    JLO = 0
	    JHI = N+1
	    GOTO 10
	END IF ! JLO

	INC = 1 ! Set hunting increment
	IF (XP .GE. X(JLO) .EQV. ASCEND) THEN ! Hunt up
20	    JHI = JLO+INC
	    IF (JHI .GT. N) THEN
		JHI = N+1
	    ELSE IF (XP .GE. X(JHI) .EQV. ASCEND) THEN
		JLO = JHI
		INC = INC+INC
		GOTO 20
	    END IF ! JHI
	ELSE ! Hunt down
	    JHI = JLO
30	    JLO = JHI-INC
	    IF (JLO .LT. 1) THEN
		JLO = 0
	    ELSE IF (XP .LT. X(JLO) .EQV. ASCEND) THEN
		JHI = JLO
		INC = INC+INC
		GOTO 30
	    END IF ! JLO
	END IF ! XP

C  The value XP is now bracketed by X(JLO) < XP < X(JHI)
C  Begin final bisection phase
10	CONTINUE
	DO WHILE (JHI-JLO .GT. 1)
	    JM = (JLO+JHI)/2
	    IF (XP .GT. X(JM) .EQV. ASCEND) THEN
		JLO = JM
	    ELSE
		JHI = JM
	    END IF
	END DO ! WHILE

	RETURN
	END ! HUNT
