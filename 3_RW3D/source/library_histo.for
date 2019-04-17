      module library_histo

	contains

C***********************************************************************
      SUBROUTINE HPLUGIN(X,N,H)                      
C-----------------------------------------------------------------------
C   Algorithm based on Plug-in Methods 
C   Reference: Engel, Herrmann, Gasser (1994), An iterative bandwidth 
C   selector for kernel estimation of densities and their derivatives,
C   Nonparametric Statistics, Vol. 4, pp. 21-34.
C-----------------------------------------------------------------------                                             
C       Purpose:                                                        
C                                                                       
C       Simple  Subroutine for kernel density estimation 
c       with iterative plug-in bandwidth selection
C       
C       This version only uses the gauss kernel and estimates only
c       the density itself and not its derivatives.
C 
C  INPUT    MASS(N)  DOUBLE PRECISION                                                                        
C  INPUT    X(N)   DOUBLE PRECISION   sorted data
C  INPUT    N      INTEGER            length of  X
C  INPUT    Z(M)   DOUBLE PRECISION   output grid (sorted)
C  INPUT    M      INTEGER            length of Z
C  OUTPUT   F(M)   DOUBLE PRECISION   estimated density
C  OUTPUT   H      DOUBLE PRECISION   estimated iterative plugin bandwidth                                                                                                                                          
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N)
      XN=DBLE(N)
      XIQR=X(NINT(.75*XN))-X(NINT(.25*XN)+1)	
      PI=3.141592645D0
      RT2PI = DSQRT(2.D0*PI)
      ITER=5

C-------  Estimate inflation constant C
          H2=(.920*XIQR)/(XN**(1.D0/7.D0))
          H3=(.912*XIQR)/(XN**(1.D0/9.D0))
          S2 = 0.D0
          S3 = 0.D0
          DO 20 I = 1, N-1
           DO 30 J = I+1, N
                     D2 = ((X(I) - X(J))/H2)**2
                     D3 = ((X(I) - X(J))/H3)**2
                     IF(D2.GT.50.AND.D3.GT.60) GOTO 20
                     E2 = DEXP(-0.5D0*D2)
                     E3 = DEXP(-0.5D0*D3)
                     S2 = S2 + (D2**2 - 6.D0*D2 + 3.D0)*E2
               S3 = S3 + (D3**3 - 15.D0*D3**2 + 45.D0*D3 - 15.D0)*E3
30          CONTINUE
20         CONTINUE
          RHAT2 = (2.D0*S2)/((XN**2)*(H2**5)*RT2PI)
          RHAT2 = RHAT2 + 3.D0/(RT2PI*XN*(H2**5))
          RHAT3 = (-2.D0*S3)/((XN**2)*(H3**7)*RT2PI)
          RHAT3 = RHAT3 + 15.D0/(RT2PI*XN*(H3**7))
          CO1 = 1.357D0*(RHAT2/RHAT3)**(1.D0/7.D0)
C-
C-------  Compute constant of asymptotic formula
C-
       CONST=1.D0/(2.D0*DSQRT(PI))
       A=1.132795764/RHAT3**(1.D0/7.D0)*XN**(-1.D0/2.D0)
C-
C------  Loop over iterations
C-
       DO 100 IT=1,ITER
C-                                                                     
C-------  Estimate functional
C
        S=0.D0                                                    
        DO 40 I = 1, N-1
           DO 50 J = I+1, N
                     D2 = ((X(I) - X(J))/A)**2
                     IF(D2.GT.50) GOTO 40
                     E2 = DEXP(-0.5D0*D2)
                     S = S + (D2**2 - 6.D0*D2 + 3.D0)*E2
50          CONTINUE
40       CONTINUE
        R2 = (2.D0*S)/((XN**2)*(A**5)*RT2PI)
        R2 = R2 + 3.D0/(RT2PI*XN*(A**5))
C-                                                                     
C-------  Estimate bandwidth by asymptotic formula
C-
         H=(CONST/(R2*XN))**(.2D0)                                    
         A=CO1*H**(5./7.)
100     CONTINUE
                                                            
      END SUBROUTINE




C***********************************************************************
      SUBROUTINE PLUGIN(X,N,Z,M,F,H,MASS)                      
C-----------------------------------------------------------------------
C   Algorithm based on Plug-in Methods 
C   Reference: Engel, Herrmann, Gasser (1994), An iterative bandwidth 
C   selector for kernel estimation of densities and their derivatives,
C   Nonparametric Statistics, Vol. 4, pp. 21-34.
C-----------------------------------------------------------------------                                             
C       Purpose:                                                        
C                                                                       
C       Simple  Subroutine for kernel density estimation 
c       with iterative plug-in bandwidth selection
C       
C       This version only uses the gauss kernel and estimates only
c       the density itself and not its derivatives.
C 
C  INPUT    MASS(N)  DOUBLE PRECISION                                                                        
C  INPUT    X(N)   DOUBLE PRECISION   sorted data
C  INPUT    N      INTEGER            length of  X
C  INPUT    Z(M)   DOUBLE PRECISION   output grid (sorted)
C  INPUT    M      INTEGER            length of Z
C  OUTPUT   F(M)   DOUBLE PRECISION   estimated density
C  OUTPUT   H      DOUBLE PRECISION   estimated iterative plugin bandwidth                                                                                                                                          
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),Z(M),F(M)
	REAL*8  MASS(N),TOTMASS
      XN=DBLE(N)
      XIQR=X(NINT(.75*XN))-X(NINT(.25*XN)+1)	
      PI=3.141592645D0
      RT2PI = DSQRT(2.D0*PI)
      ITER=5
C-------  Estimate inflation constant C
          H2=(.920*XIQR)/(XN**(1.D0/7.D0))
          H3=(.912*XIQR)/(XN**(1.D0/9.D0))
          S2 = 0.D0
          S3 = 0.D0
          DO 20 I = 1, N-1
           DO 30 J = I+1, N
                     D2 = ((X(I) - X(J))/H2)**2
                     D3 = ((X(I) - X(J))/H3)**2
                     IF(D2.GT.50.AND.D3.GT.60) GOTO 20
                     E2 = DEXP(-0.5D0*D2)
                     E3 = DEXP(-0.5D0*D3)
                     S2 = S2 + (D2**2 - 6.D0*D2 + 3.D0)*E2
               S3 = S3 + (D3**3 - 15.D0*D3**2 + 45.D0*D3 - 15.D0)*E3
30          CONTINUE
20         CONTINUE
          RHAT2 = (2.D0*S2)/((XN**2)*(H2**5)*RT2PI)
          RHAT2 = RHAT2 + 3.D0/(RT2PI*XN*(H2**5))
          RHAT3 = (-2.D0*S3)/((XN**2)*(H3**7)*RT2PI)
          RHAT3 = RHAT3 + 15.D0/(RT2PI*XN*(H3**7))
          CO1 = 1.357D0*(RHAT2/RHAT3)**(1.D0/7.D0)
C-
C-------  Compute constant of asymptotic formula
C-
       CONST=1.D0/(2.D0*DSQRT(PI))
       A=1.132795764/RHAT3**(1.D0/7.D0)*XN**(-1.D0/2.D0)
C-
C------  Loop over iterations
C-
       DO 100 IT=1,ITER
C-                                                                     
C-------  Estimate functional
C
        S=0.D0                                                    
        DO 40 I = 1, N-1
           DO 50 J = I+1, N
                     D2 = ((X(I) - X(J))/A)**2
                     IF(D2.GT.50) GOTO 40
                     E2 = DEXP(-0.5D0*D2)
                     S = S + (D2**2 - 6.D0*D2 + 3.D0)*E2
50          CONTINUE
40       CONTINUE
        R2 = (2.D0*S)/((XN**2)*(A**5)*RT2PI)
        R2 = R2 + 3.D0/(RT2PI*XN*(A**5))
C-                                                                     
C-------  Estimate bandwidth by asymptotic formula
C-
         H=(CONST/(R2*XN))**(.2D0)                                    
         A=CO1*H**(5./7.)
100     CONTINUE
C-
C------- Estimate total mass
C
        TOTMASS = SUM(MASS)
C
C------- Estimate density with plugin bandwidth
C-             
      JBEGIN=1
      JEND=1                                             
      DO 200 I=1,M                                                      
         S=0.D0    
         DO 210 J=JBEGIN,JEND
                T=(Z(I)-X(J))/H
                IF(T.GT.5.0.AND.JBEGIN.LT.N) THEN
                   JBEGIN=JBEGIN+1
                   GOTO 210
                END IF
                S=S+DEXP(-T*T/2.)*MASS(J)
210             CONTINUE         
         DO 220 JEND=J,N
                T=(Z(I)-X(JEND))/H
                IF(T.LT.-5.0) GOTO 230
                S=S+DEXP(-T*T/2.)*MASS(JEND)
220             CONTINUE
C-230         F(I)=S/(TOTMASS*H*RT2PI)
230         F(I)=S/(H*RT2PI)
            JEND=JEND-1
C-                         
200      CONTINUE                                             
C      RETURN                                                            
      END SUBROUTINE
C***********************************************************************	
C     other useful subroutines nicely provided by Engel
C***********************************************************************
      SUBROUTINE HOPTDE(X,N,NKE,NUE,KORD,NBO,NY,ISMO,BO,Z,M,G,   
     .                   ALPHA,W,F,BOP)                      
C-----------------------------------------------------------------------
C       VERSION: MARCH 12, 1992        
C       FORTRAN77 ANSI X3.9-1978/ISO 1539-1980                          
C       COPYRIGHT: STATCOM                                              
C                                                                       
C       PURPOSE:                                                        
C                                                                       
C       GENERAL SUBROUTINE FOR KERNEL ESTIMATION OF A DENSITY OR ITS    
C       DERIVATIVE.                                                     
C       COMPUTATION OF OPTIMAL BANDWIDTH USING THE ASYMPTOTIC FORMULA   
C       WITH ITERATIVE REFINEMENT, FOR OPTIMAL AND NORMAL KERNELS     
C       NUE=0,1,2;
C       IN ACCORDANCE WITH:
C	ENGEL, HERRMANN, GASSER (1992), AN ITERATIVE BANDWIDTH SELECTOR
C	     FOR KERNEL ESTIMATION OF DENSITIES AND THEIR DERIVATIVES
C                                                                       
C  PARAMETERS :                                                         
C                                                                       
C  INPUT    X(N)         DATA (MUST BE SORTED)                          
C  INPUT    N            LENGTH OF X                                    
C  INPUT    BOP          BANDWIDTH FOR ISMO~=0 (HALF-SIDED)             
C  INPUT    NKE          TYPE OF KERNEL (1: MIN. NORMAL
C					 2: OPTIMAL)
C					 0: NORMAL    
C  INPUT    NUE          ORDER OF DERIVATIVE (0-2)                      
C  INPUT    KORD         ORDER OF KERNEL (<=6; NUE+2, NUE+4 OR NUE+6)   
C				     IF NKE=0 KORD=NUE+2
C	                 ONLY IF NKE~=0:
C  INPUT    ISMO         0: OPTIMIZATION OF BANDWIDTH,                  
C                         ELSE NO OPTIMIZATION (SMOOTHING WITH B)       
C  INPUT    NBO          TREATMENT OF BOUNDARY
C                        (0: NO BOUNDARY TREATMENT,                     
C                         1: BIAS AS IN INTERIOR, BANDWIDTH SHRINKING,  
C                         2: BIAS AS IN INTERIOR, BANDWIDTH STATIONARY, 
C                         3: JUST NORMALIZING, BANDWIDTH SHRINKING,     
C                         4: JUST NORMALIZING, BANDWIDTH STATIONARY)    
C  INPUT    NY           0: FIX BANDWIDTH
C                        1: VARIABLE BANDWIDTH (B(I)=B*G(I), I=1,...,M) 
C  INPUT    BO(2)        BO(1) LEFT BOUNDARY, BO(2) RIGHT BOUNDARY      
C                        BO(1)<=X(1), BO(2)>=X(N) (DUMMY FOR NBO=0)     
C
C  INPUT    Z(M)         OUTPUT GRID (SHOULD BE EQUIDISTANT)            
C  INPUT    M            NUMBER OF POINTS TO ESTIMATE                   
C  INPUT    G(M)         BANDW. FACTORS (FOR NY=1 AND ALPHA=0.)         
C  OUTPUT   G(M)         BANDWIDTH FACTORS (FOR NY=1)                   
C                        G(M) IS DUMMY FOR NY=0                         
C  INPUT    ALPHA        0.: G IS INPUT, ~=0,>=-1.,<=1.: G IS PROPORT.  
C                        DENSITY**-ALPHA (COMPUTED IN THE PROGRAM)      
C                        (FOR ALPHA=1. KNN-ESTIMATE)                    
C                        DUMMY FOR NY=0                                 
C  WORK     W(0:N)       WORK ARRAY                                     
C  OUTPUT   F(M)         ESTIMATED DENSITY OR DERIVATIVE OF DENSITY     
C  OUTPUT   BOP          ESTIMATED OPTIMAL BANDWIDTH (ONLY FOR ISMO~=0) 
C                                                                       
C  EXTERNALS :  BANVDE, KERDEN, HHDD, NORDEN
C                                                                       
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),BO(2),Z(M),G(M),F(M),W(0:N)                        
      DIMENSION BIAS(2,2,0:2),VAR(2,2,0:2),FAK2(2:4)                    
      DATA BIAS/.3333,.2,.08571,.04762,.6,.4286,.2381,.1515,1.714,1.333,
     .          .9091,.6293/,                                           
     .  VAR/.5,.6,1.125,1.25,1.5,2.143,9.375,11.93,22.5,35.,275.6,381.6/
      DATA FAK2/4.,36.,576./                                            
      BMAX = 1.d0
	XN=DBLE(N)
      XIQR=X(NINT(.75*XN))-X(NINT(.25*XN)+1)	
C-                                                                      
C-------  NO OPTIMISATION FOR ISMO~=0                                   
      IF(ISMO.NE.0) GOTO 100                                            
C-
	IF(NKE.EQ.0) THEN
        	CALL HHDD(X,N,XIQR,NUE,BOP)
	        GOTO 100
	ELSE
C-------  COMPUTE CONSTANTS FOR LATER USE                               
      KOR2=NUE+2                                                        
      KOR3=NUE+3	
      KK=1                                                              
      IF(KORD.EQ.4.AND.N.GE.100) THEN                                   
         KOR2=4                                                          
         IF(NUE.EQ.0) KK=2                                              
         END IF                                                         
       KK2=KOR2+KOR2                                                       
       EXT=1./(2*KORD+1.)                                                    
       ITER=KORD+1
       NBO1=NBO                                                          
       IF(N.LE.200) NBO1=0                                               
       IF(NBO1.GT.2) NBO1=NBO-2                                          
       NYY=NY                                                            
       IF(KOR2.EQ.4) NYY=0                                                
       IF(NYY.EQ.1) ALPHA2=ALPHA*.2                                      
C-                                                                      
c-      **ESTIMATE INFLATION CONSTANT C
c-         DEPENDING ON ORDER OF DERIVATIVE NUE
c
c	
	IF(NUE.EQ.0) THEN
                A2=(2.57925*XIQR)/(XN**(1.d0/7.d0))
	        A3=(2.901*XIQR)/(XN**(1.d0/9.d0)) 
	   ELSEIF (NUE.EQ.1) THEN
		A2=(2.901*XIQR)/(XN**(1.d0/9.d0)) 
		A3=(3.2323*XIQR)/(XN**(1.d0/11.d0))
	   ELSEIF (NUE.EQ.2) THEN
		A2=(3.2323*XIQR)/(XN**(1.d0/11.d0))
		A3=(3.5569*XIQR)/(XN**(1.d0/13.d0))
	ENDIF	
C	
         CALL KERDEN(X,N,A2,NKE,KOR2,KOR2+2,NBO1,NYY,BO,Z,M,F)   
               R2=.5*(F(1)*F(1)+F(M)*F(M)) 
         DO 16 I=2,M-1
             R2=R2+F(I)*F(I) 	    
16	continue	
	       R2=R2*DBLE(Z(M)-Z(1))/DBLE(M-1)                                
         CALL KERDEN(X,N,A3,NKE,KOR3,KOR3+2,NBO1,NYY,BO,Z,M,F)   
               R3=.5*(F(1)*F(1)+F(M)*F(M)) 
         DO 17 I=2,M-1
17             R3=R3+F(I)*F(I) 	    
	       R3=R3*DBLE(Z(M)-Z(1))/DBLE(M-1)                                
	  IF (NUE.EQ.0) THEN
               CO = 1.459058157d0*(R2**(1.d0/5.d0))/(R3**(1.d0/7.d0))
	  ELSEIF (NUE.EQ.1) THEN
               CO = 1.489839163d0*(R2**(1.d0/7.d0))/(R3**(1.d0/9.d0))
	  ELSEIF(NUE.EQ.2) THEN
               CO = 1.679304212d0*(R2**(1.d0/9.d0))/(R3**(1.d0/11.d0))
	
	  ENDIF
          CO1=CO*XN**(2.d0/(2*KORD+1)/(2*KORD+3.d0))
c
C-------  COMPUTE INTEGRAL(F(T)**ALPHA)DT                               
c         VI=M-1.                                                          
          IF(NY.EQ.1) THEN                                                  
          CALL BANVDE(X,N,BMAX*.1,-ALPHA,Z,G,M,W,F)                      
          VI=.5/G(1)+.5/G(M)                                           
           DO 14 I=2,M-1                                              
14           VI=VI+1./G(I)                                               
c         PRINT *, VI                                                  
          END IF                                                         
C-                                                                      
C-------  COMPUTE CONSTANT OF ASYMPTOTIC FORMULA                        
          CONST=(2*NUE+1)*FAK2(KOR2)*VAR(NKE,KK,NUE)                      
     .      /(2*(KOR2-NUE)*BIAS(NKE,KK,NUE)**2*XN)                   
C-                                                                      
	  BOP=(CONST/R2)**EXT
c          PRINT *, 'Bp=',BOP,'    R2=',R2,'  CHT=',CO
C------  LOOP OVER ITERATIONS                                           
          DO 10 IT=1,ITER                                                   
C-                                                                      
C-------  ESTIMATE DERIVATIVE OF ORDER KOR2                              
          B2=BOP*CO1                                                     
          IF(NYY.EQ.1) CALL BANVDE(X,N,B2,ALPHA2,Z,G,M,W,F)              
          CALL KERDEN(X,N,B2,NKE,KOR2,KOR2+2,NBO1,NYY,BO,Z,M,F)   
C-                                                                      
C-------  ESTIMATE INTEGRAL((F(2)(T)*G(T)**KOR2)**2)DT                   
          IF(NYY.EQ.0) THEN               
               F2=.5*(F(1)*F(1)+F(M)*F(M)) 
               DO 15 I=2,M-1
15             F2=F2+F(I)*F(I) 	    
C-                                                                      
           ELSE                                                     
C-                                                                      
               F2=.5*(F(1)*F(1)*G(1)**KK2+F(M)*F(M)*G(M)**KK2)       
               DO 74 I=2,M-1                                           
74             F2=F2+F(I)*F(I)*G(I)**KK2                                
           END IF     
	   F2=F2*DBLE(Z(M)-Z(1))/DBLE(M-1) 
C-                                                                      
C-------  ESTIMATE BANDWIDTH BY ASYMPTOTIC FORMULA                      
          BOP=(CONST/F2)**EXT                                             
c          PRINT *, IT,BOP,F2                                        
10        CONTINUE                                                       
C-                                                                      
C-------  TRANSFORM BANDWIDTH FOR KORD > KOR2                            
20      CONTINUE 
c	IF(KORD.GT.KOR2) THEN                                              
cc         BOP=BOP*(VAR(NKE,2,NUE)/VAR(NKE,1,NUE))**(1./(2.*NUE+1.))      
c         IF(BOP.GT.BMAX) BOP=BMAX                                       
         END IF                                                         
c	ENDIF
C-                                                                      
C-------  COMPUTE DENSITY OR ITS DERIVATIVE WITH OPTIMAL BANDWIDTH      
100       IF(NKE.EQ.0) CALL NORDEN(X,N,BOP,NUE,Z,M,F)
	  IF(NKE.NE.0) THEN
          IF(NY.EQ.1) CALL BANVDE(X,N,BOP,ALPHA,Z,G,M,W,F)                  
          CALL KERDEN(X,N,BOP,NKE,NUE,KORD,NBO,NY,BO,Z,M,F)                 
          ENDIF
C-                                                                      
          RETURN                                                            
          END SUBROUTINE                                                               
      SUBROUTINE BANVDE(X,N,B,ALPHA,Z,G,M,W,F)                          
C-----------------------------------------------------------------------
C       VERSION: MAY 22, 1988                                           
C       FORTRAN77 ANSI X3.9-1978/ISO 1539-1980                          
C       COPYRIGHT: STATCOM                                              
C                                                                       
C       PURPOSE:                                                        
C                                                                       
C       COMPUTATION OF BANDWIDTH SEQUENCE FOR KERNEL DENSITY ESTIMATION 
C       WITH VARIABLE BANDWIDTH                                         
C                                                                       
C  PARAMETERS :                                                         
C                                                                       
C  INPUT    X(N)         DATA (ORDERED)                                 
C  INPUT    N            LENGTH OF X                                    
C  INPUT    B            ONE-SIDED MEAN BANDWIDTH                       
C  INPUT    ALPHA        0: BANDWIDTH FACTORS ARE IN G                  
C                        ~0: BANDWIDTH FACTORS ARE COMPUTED IN BANVDE   
C  INPUT    Z(M)         OUTPUT GRID                                    
C  INPUT    G(M)         BANDWIDTH FACTORS (ONLY FOR ALPHA=0.)          
C  INPUT    M            LENGTH OF Z, G AND Y                           
C  WORK     W(0:N)       WORK ARRAY                                     
C  OUTPUT   F(M)         BANDWIDTH SEQUENCE                             
C                                                                       
C  EXTERNALS :  KNN, FALPHA                                             
C                                                                       
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION X(N),Z(M),G(M),F(M),W(0:N)                              
      IF(ABS(ALPHA).GE..001) THEN                                       
         IU=1                                                           
         DO 2 I=1,M                                                     
            IF(Z(I).LT.X(1)) IU=I+1                                     
2           IF(Z(I).LE.X(N)) IO=I                                       
         MM=IO-IU+1                                                     
         CALL KNN(X,N,B,W,Z(IU),MM,F(IU))                               
         CALL FALPHA(F(IU),MM,ALPHA,G(IU))                              
         DO 3 I=1,IU                                                    
3           G(I)=G(IU)                                                  
         DO 4 I=IO,M                                                    
4           G(I)=G(IO)                                                  
         END IF                                                         
      DO 1 I=1,M                                                        
1        F(I)=B*G(I)                                                    
      RETURN                                                            
      END SUBROUTINE                                                              
      SUBROUTINE COEFFB(NUE,NKE,KORD,Q,NB,C)                            
C-----------------------------------------------------------------------
C       VERSION: SEPTEMBER 24, 1987                                     
C       FORTRAN77 ANSI X3.9-1978/ISO 1539-1980                          
C       COPYRIGHT: STATCOM                                              
C                                                                       
C       PURPOSE:                                                        
C                                                                       
C       COMPUTES COEFFICIENTS OF POLYNOMIAL BOUNDARY KERNELS,           
C       FOLLOWING GASSER + MUELLER PREPRINT 38 SFB 123 HEIDELBERG       
C       AND UNPUBLISHED RESULTS                                         
C                                                                       
C  PARAMETERS:                                                          
C                                                                       
C  INPUT  NUE        ORDER OF DERIVATIVE (0-4)                          
C  INPUT  NKE        1: MINIMUM VARIANCE KERNEL, 2: OPTIMAL KERNEL      
C  INPUT  KORD       ORDER OF KERNEL (NUE+I, I=2,4,6;  KORD<=6)         
C  INPUT  Q          PERCENTAGE OF WID AT BOUNDARY                      
C  INPUT  NB         < 0 RIGHT BOUNDARY OF DATA                         
C                    > 0 LEFT BOUNDARY OF DATA                          
C  OUTPUT C(7)       POLYNOMIAL KERNEL COEFFICIENTS                     
C                                                                       
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DOUBLE PRECISION C(7)                                             
C                                                                       
      DO 10 I=1,7                                                       
10       C(I)=0.                                                        
      P=-Q                                                              
      P1=1.+Q                                                           
      P3=P1*P1*P1                                                       
C                                                                       
      IF(NKE.EQ.1.AND.NUE.EQ.0.AND.KORD.EQ.2) THEN                      
       D=1./P3                                                          
       C(1)=4.*(1.+P*(1.+P))*D                                          
       C(2)=3.*(1.+P)*D                                                 
       END IF                                                           
C                                                                       
      IF(NKE.EQ.1.AND.NUE.EQ.0.AND.KORD.EQ.4) THEN                      
       D=1./(P3*P3*P1)                                                  
       P11=(1.+P)*D                                                     
       C(1)=16.*(1.+P*(9.+P*(45.+P*(65.+P*(45.+P*(9.+P))))))*D          
       C(2)=60.*(1.+P*(5.+P))*(1.+P*(3.+P))*P11                         
       C(3)=80.*(1.+P*(9.+P*(15.+P*(9.+P))))*D                          
       C(4)=35.*(1.+P*(8.+P))*P11                                       
       END IF                                                           
C                                                                       
      IF(NKE.EQ.1.AND.NUE.EQ.0.AND.KORD.EQ.6) THEN                      
       P6=P3*P3                                                         
       D=P1/(P6*P6)                                                     
       P11=(1.+P)*D                                                     
       C(1)=36.*(1.+P*(25.+P*(325.+P*(1700.+P*(4550.+P*(6202.+P         
     .      *(4550.+P*(1700.+P*(325.+P*(25.+P))))))))))*D               
       C(2)=315.*(1.+P*(14.+P*(36.+P*(14.+P))))*(1.+P*(10.+P*(20.+P     
     .      *(10.+P))))*P11                                             
       C(3)=560.*(2.+P*(50.+P*(335.+P*(992.+P*(1400.+P*(992.+P*(335.+P  
     .      *(50.+P+P))))))))*D                                         
       C(4)=1890.*(1.+P*(24.+P*(112.+P*(188.+P*(112.+P*(24.+P))))))*P11 
       C(5)=1512.*(1.+P*(25.+P*(115.+P*(180.+P*(115.+P*(25.+P))))))*D   
       C(6)=462.*(1.+P*(24.+P*(76.+P*(24.+P))))*P11                     
       END IF                                                           
C                                                                       
      IF(NKE.EQ.2.AND.NUE.EQ.0.AND.KORD.EQ.2) THEN                      
       D=1./(P3*P1)                                                     
       C(1)=(6.+P*(12.+P*18.))*D                                        
       C(2)=9.*(1.+P)*(1.+P)*D                                          
       C(3)=(4.+P*8.)*D                                                 
       END IF                                                           
C                                                                       
      IF(NKE.EQ.2.AND.NUE.EQ.0.AND.KORD.EQ.4) THEN                      
       D=P1/(P3*P3*P3)                                                  
       P12=(1.+P)*(1.+P)*D                                              
       C(1)=20.*(1.+P*(12.+P*(78.+P*(164.+P*(165.+P*(60.+P*10.))))))*D  
       C(2)=100.*(1.+P*(5.+P))**2*P12                                   
       C(3)=200.*(1.+P*(12.+P*(33.+P*(36.+P*(14.+P+P)))))*D             
       C(4)=175.*(1.+P*(10.+P*3.))*P12                                  
       C(5)=56.*(1.+P*(12.+P*(18.+P*4.)))*D                             
       END IF                                                           
C                                                                       
      IF(NKE.EQ.2.AND.NUE.EQ.0.AND.KORD.EQ.6) THEN                      
       P6=P3*P3                                                         
       D=1./(P6*P6)                                                     
       P12=(1.+P)*(1.+P)*D                                              
       C(1)=42.*(1.+P*(30.+P*(465.+P*(3000.+P*(10050.+P*(17772.+P       
     .      *(17430.+P*(9240.+P*(2625.+P*(350.+P*21.))))))))))*D        
       C(2)=441.*(1.+P*(14.+P*(36.+P*(14.+P))))**2*P12                  
       C(3)=1960.*(1.+P*(30.+P*(255.+P*(984.+P*(1902.+P*(1956.+P        
     .      *(1065.+P*(300.+P*(39.+P+P)))))))))*D                       
       C(4)=4410.*(1.+P*(28.+P*(156.+P*(308.+P*(188.+P*(42.+P*3.))))))  
     .      *P12                                                        
       C(5)=5292.*(1.+P*(30.+P*(185.+P*(440.+P*(485.+P*(250.+P*(57.+P*4.
     .      )))))))*D                                                   
       C(6)=3234.*(1.+P*(28.+P*(108.+P*(56.+P*5.))))*P12                
       C(7)=792.*(1.+P*(30.+P*(150.+P*(200.+P*(75.+P*6.)))))*D          
       END IF                                                           
C                                                                       
      IF(NKE.EQ.1.AND.NUE.EQ.1.AND.KORD.EQ.3) THEN                      
       D=-P1/(P3*P3)                                                    
       P11=(1.+P)*D                                                     
       C(1)=36.*(1.+P*(3.+P))*P11                                       
       C(2)=(96.+P*(168.+P*96.))*D                                      
       C(3)=60.*P11                                                     
       END IF                                                           
C                                                                       
      IF(NKE.EQ.1.AND.NUE.EQ.1.AND.KORD.EQ.5) THEN                      
       D=-1./(P3*P3*P3)                                                 
       P11=(1.+P)*D                                                     
       P2P1=(1.+P*(5.+P))*P11                                           
       C(1)=300.*(1.+P*(10.+P*(20.+P*(10.+P))))*P2P1                    
       C(2)=600.*(4.+P*(39.+P*(144.+P*(214.+P*(144.+P*(39.+P*4.))))))*D 
       C(3)=2100.*(3.+P*(8.+P*3.))*(1.+P*(4.+P))*P11                    
       C(4)=840.*(8.+P*(53.+P*(88.+P*(53.+P*8.))))*D                    
       C(5)=2520.*P2P1                                                  
       END IF                                                           
C                                                                       
      IF(NKE.EQ.2.AND.NUE.EQ.1.AND.KORD.EQ.3) THEN                      
       D=-1./(P3*P3)                                                    
       P12=(1.+P)*(1.+P)*D                                              
       C(1)=(60.+P*240.)*P12                                            
       C(2)=120.*(2.+P*(6.+P*(6.+P)))*D                                 
       C(3)=300.*P12                                                    
       C(4)=(120.+P*180.)*D                                             
       END IF                                                           
C                                                                       
      IF(NKE.EQ.2.AND.NUE.EQ.1.AND.KORD.EQ.5) THEN                      
       D=-1./(P3*P3*P3*P1)                                              
       P12=(1.+P)*(1.+P)*D                                              
       C(1)=420.*(1.+P*(18.+P*(98.+P*(176.+P*(75.+P*10.)))))*P12        
       C(2)=2100.*(2.+P*(25.+P*(120.+P*(245.+P*(238.+P*(105.+P*(20.+P)  
     .      ))))))*D                                                    
       C(3)=14700.*(1.+P*(4.+P))**2*P12                                 
       C(4)=5880.*(4.+P*(35.+P*(90.+P*(95.+P*(40.+P*6.)))))*D           
       C(5)=17640.*(1.+P*(6.+P+P))*P12                                  
       C(6)=2520.*(2.+P*(15.+P*(20.+P*5.)))*D                           
       END IF                                                           
C                                                                       
      IF(NKE.EQ.1.AND.NUE.EQ.2.AND.KORD.EQ.4) THEN                      
       D=1./(P3*P3*P1)                                                  
       P11=(1.+P)*D                                                     
       C(1)=480.*(1.+P*(9.+P*(15.+P*(9.+P))))*D                         
       C(2)=900.*(3.+P*(8.+P*3.))*P11                                   
       C(3)=480.*(9.+P*(17.+P*9.))*D                                    
       C(4)=2100.*P11                                                   
       END IF                                                           
C                                                                       
      IF(NKE.EQ.1.AND.NUE.EQ.2.AND.KORD.EQ.6) THEN                      
       P6=P3*P3                                                         
       D=P1/(P6*P6)                                                     
       P11=(1.+P)*D                                                     
       P2P1=(1.+P*(4.+P))*P11                                           
       C(1)=3360.*(2.+P*(50.+P*(335.+P*(992.+P*(1400.+P*(992.+P*(335.+P 
     .      *(50.+P+P))))))))*D                                         
       C(2)=88200.*(1.+P*(8.+P*(15.+P*(8.+P))))*P2P1                    
       C(3)=94080.*(4.+P*(36.+P*(120.+P*(175.+P*(120.+P*(36.+P*4.)))))) 
     .      *D                                                          
       C(4)=176400.*(2.+P*(7.+P+P))*(1.+P+P)*(2.+P)*P11                 
       C(5)=60480.*(10.+P*(58.+P*(95.+P*(58.+P*10.))))*D                
       C(6)=194040.*P2P1                                                
       END IF                                                           
C                                                                       
      IF(NKE.EQ.2.AND.NUE.EQ.2.AND.KORD.EQ.4) THEN                      
       D=P1/(P3*P3*P3)                                                  
       P12=(1.+P)*(1.+P)*D                                              
       C(1)=840.*(1.+P*(12.+P*(28.+P*(24.+P*5.))))*D                    
       C(2)=2100.*(3.+P*(10.+P))*P12                                    
       C(3)=1680.*(9.+P*(28.+P*(27.+P*6.)))*D                           
       C(4)=14700.*P12                                                  
       C(5)=(5040.+P*6720.)*D                                           
       END IF                                                           
C                                                                       
      IF(NKE.EQ.2.AND.NUE.EQ.2.AND.KORD.EQ.6) THEN                      
       P6=P3*P3                                                         
       D=1./(P6*P6)                                                     
       P12=(1.+P)*(1.+P)*D                                              
       C(1)=5040.*(2.+P*(60.+P*(489.+P*(1786.+P*(3195.+P*(2952.+P       
     .      *(1365.+P*(294.+P*21.))))))))*D                             
       C(2)=52920.*(3.+P*(42.+P*(188.+P*(308.+P*(156.+P*(28.+P))))))*P12
       C(3)=141120.*(6.+P*(68.+P*(291.+P*(570.+P*(555.+P*(264.+P        
     .      *(57.+P*4.)))))))*D                                         
       C(4)=529200.*(2.+P*(7.+P+P))**2*P12                              
       C(5)=90720.*(30.+P*(228.+P*(559.+P*(582.+P*(255.+P*40.)))))*D    
       C(6)=582120.*(3.+P*(14.+P*5.))*P12                               
       C(7)=221760.*(2.+P*(12.+P*(15.+P*4.)))*D                         
       END IF                                                           
C                                                                       
      IF(NKE.EQ.1.AND.NUE.EQ.3.AND.KORD.EQ.5) THEN                      
       D=-1./(P3*P3*P3)                                                 
       P11=(1.+P)*D                                                     
       C(1)=8400.*(1.+P*(15.+P*(31.+P*(15.+P))))*P11                    
       C(2)=10080.*(8.+P*(53.+P*(88.+P*(53.+P*8.))))*D                  
       C(3)=117600.*(1.+P+P)*(2.+P)*P11                                 
       C(4)=16800.*(16.+P*(31.+P*16.))*D                                
       C(5)=105840.*P11                                                 
       END IF                                                           
C                                                                       
      IF(NKE.EQ.2.AND.NUE.EQ.3.AND.KORD.EQ.5) THEN                      
       D=-1./(P3*P3*P3*P1)                                              
       P12=(1.+P)*(1.+P)*D                                              
       C(1)=15120.*(1.+P*(18.+P*(38.+P*6.)))*P12                        
       C(2)=45360.*(4.+P*(35.+P*(80.+P*(70.+P*(20.+P)))))*D             
       C(3)=352800.*(2.+P*(6.+P))*P12                                   
       C(4)=151200.*(8.+P*(25.+P*(24.+P*6.)))*D                         
       C(5)=952560.*P12                                                 
       C(6)=70560.*(4.+P*5.)*D                                          
       END IF                                                           
C                                                                       
      IF(NKE.EQ.1.AND.NUE.EQ.4.AND.KORD.EQ.6) THEN                      
       P6=P3*P3                                                         
       D=P1/(P6*P6)                                                     
       P11=(1.+P)*D                                                     
       C(1)=181440.*(1.+P*(25.+P*(115.+P*(180.+P*(115.+P*(25.+P))))))*D 
       C(2)=264600.*(10.+P*(96.+P*(184.+P*(96.+P*10.))))*P11            
       C(3)=1209600.*(10.+P*(58.+P*(95.+P*(58.+P*10.))))*D              
       C(4)=2381400.*(10.+P*(24.+P*10.))*P11                            
       C(5)=211680.*(100.+P*(196.+P*100.))*D                            
       C(6)=6985440.*P11                                                
       END IF                                                           
C                                                                       
      IF(NKE.EQ.2.AND.NUE.EQ.4.AND.KORD.EQ.6) THEN                      
       P6=P3*P3                                                         
       D=1./(P6*P6)                                                     
       P12=(1.+P)*(1.+P)*D                                              
       C(1)=332640.*(1.+P*(30.+P*(171.+P*(340.+P*(285.+P*(90.+P*7.      
     .      ))))))*D                                                    
       C(2)=1164240.*(5.+P*(56.+P*(108.+P*(28.+P))))*P12                
       C(3)=6652800.*(5.+P*(38.+P*(85.+P*(76.+P*(25.+P+P)))))*D         
       C(4)=17463600.*(5.+P*(14.+P*3.))*P12                             
       C(5)=4656960.*(25.+P*(78.+P*(75.+P*20.)))*D                      
       C(6)=76839840.*P12                                               
       C(7)=3991680.*(5.+P*6.)*D                                        
       END IF                                                           
C                                                                       
      IF(NB.GT.0) RETURN                                                
      J=2                                                               
      IF(NUE.EQ.1.OR.NUE.EQ.3) J=1                                      
      DO 2 I=J,KORD,2                                                   
2        C(I)=-C(I)                                                     
      RETURN                                                            
      END SUBROUTINE                                                              
      SUBROUTINE COEFFI(NUE,NKE,KORD,C)                                 
C-----------------------------------------------------------------------
C       VERSION: SEPTEMBER 24, 1987                                     
C       FORTRAN77 ANSI X3.9-1978/ISO 1539-1980                          
C       COPYRIGHT: STATCOM                                              
C                                                                       
C       PURPOSE:                                                        
C                                                                       
C       DEFINES POLYNOMIAL KERNEL COEFFICIENTS FOR INTERIOR.            
C                                                                       
C  PARAMETERS:                                                          
C                                                                       
C  INPUT  NUE        ORDER OF DERIVATIVE (0-4)                          
C  INPUT  NKE        1: MINIMUM VARIANCE KERNEL, 2: OPTIMAL KERNEL      
C  INPUT  KORD       ORDER OF KERNEL (NUE+I, I=2,4,6;  KORD<=6)         
C  OUTPUT C(7)       POLYNOMIAL KERNEL COEFFICIENTS                     
C                                                                       
C-----------------------------------------------------------------------
      DOUBLE PRECISION C(7)                                             
C                                                                       
      DO 10 I=1,7                                                       
10       C(I)=0.                                                        
C                                                                       
      IF(NKE.EQ.1.AND.NUE.EQ.0.AND.KORD.EQ.2) C(1)=.5                   
C                                                                       
      IF(NKE.EQ.1.AND.NUE.EQ.0.AND.KORD.EQ.4) THEN                      
       C(1)=1.125D0                                                     
       C(3)=-0.625D0                                                    
       END IF                                                           
C                                                                       
      IF(NKE.EQ.1.AND.NUE.EQ.0.AND.KORD.EQ.6) THEN                      
       C(1)=1.7578125D0                                                 
       C(3)=-2.734375D0                                                 
       C(5)=1.4765625D0                                                 
       END IF                                                           
C                                                                       
      IF(NKE.EQ.2.AND.NUE.EQ.0.AND.KORD.EQ.2) THEN                      
       C(1)=0.75D0                                                      
       C(3)=-0.25D0                                                     
       END IF                                                           
C                                                                       
      IF(NKE.EQ.2.AND.NUE.EQ.0.AND.KORD.EQ.4) THEN                      
       C(1)=1.40625D0                                                   
       C(3)=-1.5625D0                                                   
       C(5)=0.65625D0                                                   
       END IF                                                           
C                                                                       
      IF(NKE.EQ.2.AND.NUE.EQ.0.AND.KORD.EQ.6) THEN                      
       C(1)=2.05078125D0                                                
       C(3)=-4.78515625D0                                               
       C(5)=5.16796875D0                                                
       C(7)=-1.93359375D0                                               
       END IF                                                           
C                                                                       
      IF(NKE.EQ.1.AND.NUE.EQ.1.AND.KORD.EQ.3) C(2)=-.75D0               
C                                                                       
      IF(NKE.EQ.1.AND.NUE.EQ.1.AND.KORD.EQ.5) THEN                      
       C(2)=-4.6875D0                                                   
       C(4)=3.28125D0                                                   
       END IF                                                           
C                                                                       
      IF(NKE.EQ.2.AND.NUE.EQ.1.AND.KORD.EQ.3) THEN                      
       C(2)=-1.875D0                                                    
       C(4)=0.9375D0                                                    
       END IF                                                           
C                                                                       
      IF(NKE.EQ.2.AND.NUE.EQ.1.AND.KORD.EQ.5) THEN                      
       C(2)=-8.203125D0                                                 
       C(4)=11.484375D0                                                 
       C(6)=-4.921875D0                                                 
       END IF                                                           
C                                                                       
      IF(NKE.EQ.1.AND.NUE.EQ.2.AND.KORD.EQ.4) THEN                      
       C(1)=-3.75D0                                                     
       C(3)=3.75D0                                                      
       END IF                                                           
C                                                                       
      IF(NKE.EQ.1.AND.NUE.EQ.2.AND.KORD.EQ.6) THEN                      
       C(1)=-16.40625D0                                                 
       C(3)=45.9375D0                                                   
       C(5)=-29.53125D0                                                 
       END IF                                                           
C                                                                       
      IF(NKE.EQ.2.AND.NUE.EQ.2.AND.KORD.EQ.4) THEN                      
       C(1)=-6.5625D0                                                   
       C(3)=13.125D0                                                    
       C(5)=-6.5625D0                                                   
       END IF                                                           
C                                                                       
      IF(NKE.EQ.2.AND.NUE.EQ.2.AND.KORD.EQ.6) THEN                      
       C(1)=-24.609375D0                                                
       C(3)=103.359375D0                                                
       C(5)=-132.890625D0                                               
       C(7)=54.140625D0                                                 
       END IF                                                           
C                                                                       
      IF(NKE.EQ.1.AND.NUE.EQ.3.AND.KORD.EQ.5) THEN                      
       C(2)=39.375D0                                                    
       C(4)=-32.8125D0                                                  
       END IF                                                           
C                                                                       
      IF(NKE.EQ.2.AND.NUE.EQ.3.AND.KORD.EQ.5) THEN                      
       C(2)=88.59375D0                                                  
       C(4)=-147.65625D0                                                
       C(6)=68.90625D0                                                  
       END IF                                                           
C                                                                       
      IF(NKE.EQ.1.AND.NUE.EQ.4.AND.KORD.EQ.6) THEN                      
       C(1)=177.1875D0                                                  
       C(3)=-590.625D0                                                  
       C(5)=413.4375D0                                                  
       END IF                                                           
C                                                                       
      IF(NKE.EQ.2.AND.NUE.EQ.4.AND.KORD.EQ.6) THEN                      
       C(1)=324.84375D0                                                 
       C(3)=-1624.21875D0                                               
       C(5)=2273.90625D0                                                
       C(7)=-974.53125D0                                                
       END IF                                                           
C                                                                       
      RETURN                                                            
      END  SUBROUTINE                                                             
      SUBROUTINE DENCON(X,N,B,NKE,NUE,KORD,NBO,NY,BO,Z,M,F)             
C-----------------------------------------------------------------------
C       VERSION: APRIL 12, 1988                                         
C       FORTRAN77 ANSI X3.9-1978/ISO 1539-1980                          
C       COPYRIGHT: STATCOM                                              
C                                                                       
C       PURPOSE:                                                        
C                                                                       
C       SUBROUTINE FOR KERNEL ESTIMATION OF DENSITY OR DERIV. OF DENS.  
C       USES CONVENTIONAL ALGORITHM FOLLOWING ROSENBLATT (1956,AMS) AND 
C       PARZEN (1962,AMS)                                               
C       FIX OR VARIABLE BANDWIDTH                                       
C                                                                       
C  PARAMETERS :                                                         
C                                                                       
C  INPUT    X(N)         DATA (ORDERED)                                 
C  INPUT    N            LENGTH OF X                                    
C  INPUT    B            ONE-SIDED BANDWIDTH (DUMMY FOR NY~=0)          
C  INPUT    NKE          TYPE OF KERNEL (1: MIN. VAR., 2: OPTIMAL)      
C  INPUT    NUE          ORDER OF DERIVATIVE (0-4)                      
C  INPUT    KORD         ORDER OF KERNEL (<=6)                          
C  INPUT    NBO          TREATMENT OF BOUNDARY                          
C                        (0: NO BOUNDARY TREATMENT,                     
C                         1: BIAS AS IN INTERIOR, BANDWIDTH SHRINKING,  
C                         2: BIAS AS IN INTERIOR, BANDWIDTH STATIONARY, 
C                         3: JUST NORMALIZING, BANDWIDTH SHRINKING,     
C                         4: JUST NORMALIZING, BANDWIDTH STATIONARY)    
C  INPUT    NY           0: FIX BANDWIDTH, ELSE VARIABLE BANDWIDTH IN F 
C  INPUT    BO(2)        BO(1) LEFT BOUNDARY, BO(2) RIGHT BOUNDARY      
C                        BO(1)<=X(1), BO(2)>=X(N) (DUMMY FOR NBO=0)     
C  INPUT    Z(M)         OUTPUT GRID (ORDERED) WHERE DENS. IS ESTIMATED 
C  INPUT    M            NUMBER OF POINTS TO ESTIMATE                   
C  INPUT    F(M)         BANDWITH SEQUENCE FOR NY~=0, DUMMY FOR NY=0    
C  OUTPUT   F(M)         ESTIMATED DENSITY OR DERIVATIVE OF DENSITY     
C                                                                       
C  EXTERNALS :  COEFFI, COEFFB                                          
C                                                                       
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION X(N),Z(M),F(M),BO(2)                     
      DOUBLE PRECISION C(7),C1(7)
C-                                                                      
C------  COMPUTE COEFFICIENTS FOR INTERIOR AND CONSTANTS FOR LATER USE  
      CALL COEFFI(NUE,NKE,KORD,C)                                       
      IORDB=KORD+NKE-1                                                  
      IORD=IORDB+NKE-2                                                  
      DO 1 K=2,IORD                                                     
1        IF(C(K).NE.0.) C(K)=K*C(K)                                     
C-                                                                      
      if(nbo.eq.0) then
      bmax=(x(n)-x(1))*.5
      else
      BMAX=(bo(2)-bo(1))*.5                                               
      end if
      NB=0                                                              
      LOW=1                                                             
C-                                                                      
C-------  LOOP OVER OUTPUT GRID                                         
      DO 100 I=1,M                                                      
         BB=B                                                           
         IF(NY.NE.0) BB=F(I)                                            
         IF(BB.GT.BMAX) BB=BMAX                                         
         IF(BB.LE.0.) THEN                                              
            F(I)=0.                                                     
            GOTO 100                                                    
            END IF                                                      
         S=0.                                                           
C-                                                                      
C--------- BOUNDARY TREATMENT IF NBO ~= 0                               
         IF(NBO.NE.0) THEN                                              
            NB=0                                                        
            IF(Z(I).LT.BO(1).OR.Z(I).GT.BO(2)) GOTO 210                 
            IF(Z(I).LT.BO(1)+BB) THEN                                   
               Q=Z(I)-BO(1)                                             
               NB=1                                                     
               END IF                                                   
            IF(Z(I)+BB.GT.BO(2)) THEN                                   
               Q=BO(2)-Z(I)                                             
               NB=-1                                                    
               END IF                                                   
            IF(NB.NE.0) THEN                                            
               IF(NBO.EQ.2.OR.NBO.EQ.4) BB=BB+BB-Q                      
               IF(NBO.LE.2) THEN                                        
C-                                                                      
C---------- COMPUTE BOUNDARY KERNEL COEFFICIENTS IF NBO = 1,2           
                  CALL COEFFB(NUE,NKE,KORD,Q/BB,NB,C1)                  
                  DO 2 K=2,IORDB                                        
2                    C1(K)=K*C1(K)                                      
C-                                                                      
                     ELSE                                               
C-                                                                      
C---------- NORMALISATION OF KERNEL COEFFICIENTS IF NBO = 3,4           
                  Q1=Q/BB                                               
                  Q2=Q1*Q1                                              
                  A=C(1)*(1.+Q1)                                        
                  DO 3 K=3,IORD,2                                       
                     Q1=Q1*Q2                                           
3                    A=A+C(K)*(1.+Q1)/K                                 
                  DO 4 K=1,IORD,2                                       
4                    C1(K)=C(K)/A                                       
                  END IF                                                
               END IF                                                   
            END IF                                                      
C-                                                                      
C--------- COMPUTE FIRST POINT WHICH FALLS INTO THE SMOOTHING INTERVAL  
         IF(NY.NE.0.AND.LOW.GT.1) THEN                                  
10          IF(Z(I).LT.X(LOW-1)+BB) THEN                                
               LOW=LOW-1                                                
               IF(LOW.GT.1) GOTO 10                                     
               END IF                                                   
            END IF                                                      
C-                                                                      
C--------- LOOP OVER DATA POINTS WHICH FALL INTO THE SMOOTHING INTERVAL 
         LOW1=LOW                                                       
         DO 200 J=LOW1,N                                                
            U=(Z(I)-X(J))/BB                                            
            IF(U.GT.1.) THEN                                            
               LOW=LOW+1                                                
               GOTO 200                                                 
               END IF                                                   
            IF(U.LT.-1.) GOTO 210                                       
C-                                                                      
C--------- COMPUTE SUM OF KERNEL WEIGHTS                                
            IF(NB.NE.0) THEN                                            
               U1=U                                                     
               S=S+C1(1)+U*C1(2)                                        
               DO 40 K=3,IORDB                                          
                  U1=U1*U                                               
40                S=S+U1*C1(K)                                          
C-                                                                      
                  ELSE                                                  
C-                                                                      
               U2=U*U                                                   
               IF(IORD.EQ.1) S=S+C(1)                                   
               IF(IORD.EQ.2) S=S+U*C(2)                                 
               IF(IORD.EQ.3) S=S+C(1)+U2*C(3)                           
               IF(IORD.EQ.4) S=S+U*(C(2)+U2*C(4))                       
               IF(IORD.EQ.5) S=S+C(1)+U2*(C(3)+U2*C(5))                 
               IF(IORD.EQ.6) S=S+U*(C(2)+U2*(C(4)+U2*C(6)))             
               IF(IORD.EQ.7) S=S+C(1)+U2*(C(3)+U2*(C(5)+U2*C(7)))       
               END IF                                                   
200         CONTINUE                                                    
C-                                                                      
C--------- COMPUTE KERNEL ESTIMATION AT POINT I                         
210      IF(NUE.EQ.0) THEN                                              
            F(I)=S/(N*BB)                                               
cc            IF(F(I).LT.0.) F(I)=0.                                      
               ELSE                                                     
            F(I)=S/(N*BB**(NUE+1))                                      
            END IF                                                      
100      CONTINUE                                                       
C-                                                                      
      RETURN                                                            
      END SUBROUTINE                                                              
      SUBROUTINE DENLEG(X,N,B,NKE,NUE,KORD,NBO,NY,BO,Z,M,F)             
C-----------------------------------------------------------------------
C       VERSION: MAY 5, 1988                                            
C       FORTRAN77 ANSI X3.9-1978/ISO 1539-1980                          
C       COPYRIGHT: STATCOM                                              
C                                                                       
C       PURPOSE:                                                        
C                                                                       
C       COMPUTATION OF DENSITY ESTIMATE USING O(N) ALGORITHM BASED ON   
C       LEGENDRE POLYNOMIALS FOR FIX OR VARIABLE BANDWIDTH.             
C       (NEW INITIALISATIONS OF THE LEGENDRE SUMS BECAUSE OF NUMERICAL  
C       REASONS)                                                        
C                                                                       
C  PARAMETERS :                                                         
C                                                                       
C  INPUT    X(N)         DATA SORTED IN ASCENDING ORDER                 
C  INPUT    N            LENGTH OF X                                    
C  INPUT    B            ONE-SIDED BANDWIDTH (DUMMY FOR NY~=0)          
C  INPUT    NKE          TYPE OF KERNEL (1: MIN. VAR., 2: OPTIMAL)      
C  INPUT    NUE          ORDER OF DERIVATIVE (0-4)                      
C  INPUT    KORD         ORDER OF KERNEL (<=6)                          
C  INPUT    NBO          TREATMENT OF BOUNDARY                          
C                        (0: NO BOUNDARY TREATMENT,                     
C                         1: BIAS AS IN INTERIOR, BANDWIDTH SHRINKING,  
C                         2: BIAS AS IN INTERIOR, BANDWIDTH STATIONARY, 
C                         3: JUST NORMALIZING, BANDWIDTH SHRINKING,     
C                         4: JUST NORMALIZING, BANDWIDTH STATIONARY)    
C  INPUT    NY           0: FIX BANDWIDTH, ELSE VARIABLE BANDWIDTH IN F 
C  INPUT    BO(2)        BO(1) LEFT BOUNDARY, BO(2) RIGHT BOUNDARY      
C                        BO(1)<=X(1), BO(2)>=X(N) (DUMMY FOR NBO=0)     
C  INPUT    Z(M)         OUTPUT GRID (ORDERED) WHERE DENS. IS ESTIMATED 
C  INPUT    M            NUMBER OF POINTS TO ESTIMATE                   
C  INPUT    F(M)         BANDWITH SEQUENCE FOR NY~=0, DUMMY FOR NY=0    
C  OUTPUT   F(M)         ESTIMATED DENSITY OR DERIVATIVE OF DENSITY     
C                                                                       
C  EXTERNALS :  COEFFB, COEFFI, LEGDEN, FINDEN, DSUDEN                  
C                                                                       
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION X(*),Z(*),F(*),BO(2)                                    
      DOUBLE PRECISION C(7),SW(0:6),A1(6),A2(6)                         
C-                                                                      
C------ COMPUTE CONSTANTS FOR LATER USE                                 
      if(nbo.eq.0) then
      BMAX=(X(N)-X(1))*.5
      else
      bmax=(bo(2)-bo(1))*.5                                               
      end if
      IORD=KORD+NKE-2                                                   
      IF(NBO.GT.2) CALL COEFFI(NUE,NKE,KORD,C)                          
      DO 2 K=3,IORD                                                     
         A1(K)=(K+K-1.D0)/K                                             
 2       A2(K)=(1.D0-K)/K                                               
C-                                                                      
      INIT=0                                                            
      ICALL=0                                                           
C-                                                                      
C                                                                       
C    ************************************                               
C    * S M O O T H I N G   -   L O O P  *                               
C    ************************************                               
C                                                                       
C-                                                                      
      DO 100 I=1,M                                                      
         BB=B                                                           
         IF(NY.NE.0) BB=F(I)                                            
         IF(BB.GT.BMAX) BB=BMAX                                         
         IF(BB.LE.0.) THEN                                              
            F(I)=0.                                                     
            INIT=0                                                      
            GOTO 100                                                    
            END IF                                                      
         IBOUN=0                                                        
         IF(NBO.EQ.0) GOTO 123                                          
C-                                                                      
C----------------------------                                           
C------ BOUNDARY TREATMENT  |                                           
C----------------------------                                           
C-                                                                      
C------ COMPUTE LEFT BOUNDARY KERNEL                                    
         IF(Z(I).LT.BO(1)+BB) THEN                                      
            WWL=BO(1)                                                   
            WWR=BO(1)+BB+BB                                             
            IF(NBO.EQ.1.OR.NBO.EQ.3) WWR=Z(I)+BB                        
            WID=WWR-Z(I)                                                
            IF(NBO.LE.2) CALL COEFFB(NUE,NKE,KORD,(Z(I)-BO(1))/WID,1,C) 
            IBOUN=1                                                     
            END IF                                                      
C-                                                                      
C------ COMPUTE RIGHT BOUNDARY KERNEL                                   
         IF(Z(I)+BB.GT.BO(2)) THEN                                      
            WWL=BO(2)-(BB+BB)                                           
            IF(NBO.EQ.1.OR.NBO.EQ.3) WWL=Z(I)-BB                        
            WWR=BO(2)                                                   
            WID=Z(I)-WWL                                                
            IF(NBO.LE.2) CALL COEFFB(NUE,NKE,KORD,(BO(2)-Z(I))/WID,-1,C)
            IBOUN=-1                                                    
            END IF                                                      
C-                                                                      
 123     IF(IBOUN.EQ.0) THEN                                            
            WID=BB                                                      
            WWL=Z(I)-BB                                                 
            WWR=Z(I)+BB                                                 
            IF(NBO.GT.2) XNOR=1.D0                                      
            END IF                                                      
C-                                                                      
C------- COMPUTE NORMALIZING CONSTANT FOR NBO > 2                       
         IF(NBO.GT.2.AND.IBOUN.NE.0) THEN                               
            IF(IBOUN.EQ.1) Q=(Z(I)-BO(1))/WID                           
            IF(IBOUN.EQ.-1) Q=(BO(2)-Z(I))/WID                          
            QQ=Q*Q                                                      
            XNOR=C(1)*(1.D0+Q)                                          
            DO 4001 K=3,IORD+1,2                                        
               Q=Q*QQ                                                   
 4001          XNOR=XNOR+C(K)*(1.D0+Q)                                  
            IBOUN=0                                                     
            END IF                                                      
C-                                                                      
C------ INITIALISATION FOR INIT=0                                       
         IF(INIT.EQ.0) THEN                                             
C-                                                                      
C------------------------------------------------------                 
C------ INITIALISE LEFT AND RIGHT SUM                                   
C------------------------------------------------------                 
C-                                                                      
            DO 44 K=0,IORD                                              
 44            SW(K)=0.                                                 
            JL=1                                                        
            DO 48 J=1,N                                                 
               IF(X(J).LT.WWL) THEN                                     
                  JL=J+1                                                
                     ELSE                                               
                  IF(X(J).GT.WWR) GOTO 488                              
                  CALL DSUDEN(SW,A1,A2,IORD,X(J),Z(I),WID,1)            
                  END IF                                                
 48            CONTINUE                                                 
 488        JR=J-1                                                      
            WR=WWR                                                      
            INIT=1                                                      
            GOTO 6666                                                   
            END IF                                                      
C-                                                                      
C-------------------------------------------------------------------    
C------ COMPARE OLD SUM WITH NEW SMOOTHING INTERVALL Z(I)-B,Z(I)+B      
C-------------------------------------------------------------------    
C-                                                                      
C------ IS IT POSSIBLE TO USE SOME OF THE OLD TERMS ?                   
         IF(X(JR).GE.WWL) THEN                                          
C-                                                                      
            JNR=JR                                                      
            JNL=JL                                                      
            IF(X(JR).GT.WWR) THEN                                       
               DO 201 J=JR,JL,-1                                        
                  IF(X(J).LE.WWR) GOTO 2011                             
                  CALL DSUDEN(SW,A1,A2,IORD,X(J),Z(I-1),WIDO,-1)        
                  JNR=J-1                                               
 201              CONTINUE                                              
 2011          END IF                                                   
C-                                                                      
            IF(X(JL).LT.WWL) THEN                                       
               DO 301 J=JL,JR                                           
                  IF(X(J).GE.WWL) GOTO 3011                             
                  CALL DSUDEN(SW,A1,A2,IORD,X(J),Z(I-1),WIDO,-1)        
                  JNL=J+1                                               
 301              CONTINUE                                              
 3011          END IF                                                   
C-                                                                      
C------ UPDATING OF SW                                                  
            CALL LEGDEN(SW,IORD,(Z(I)-Z(I-1))/WID,WIDO/WID)             
C-                                                                      
            IF(JNR.EQ.JR) THEN                                          
               DO 401 J=JR+1,N                                          
                  IF(X(J).GT.WWR) GOTO 4011                             
                  CALL DSUDEN(SW,A1,A2,IORD,X(J),Z(I),WID,1)            
                  JNR=J                                                 
 401              CONTINUE                                              
 4011          END IF                                                   
            JR=JNR                                                      
C-                                                                      
            IF(JL.EQ.JNL) THEN                                          
               DO 402 J=JL-1,1,-1                                       
                  IF(X(J).LT.WWL) GOTO 4022                             
                  CALL DSUDEN(SW,A1,A2,IORD,X(J),Z(I),WID,1)            
                  JNL=J                                                 
 402              CONTINUE                                              
 4022          END IF                                                   
            JL=JNL                                                      
C-                                                                      
               ELSE                                                     
C-                                                                      
C------ NEW INITIALISATION OF SW                                        
            DO 22 K=0,IORD                                              
 22            SW(K)=0.                                                 
            DO 202 J=JR,N                                               
               IF(X(J).LT.WWL) THEN                                     
                  JL=J+1                                                
                     ELSE                                               
                  IF(X(J).GT.WWR) GOTO 2022                             
                  CALL DSUDEN(SW,A1,A2,IORD,X(J),Z(I),WID,1)            
                  END IF                                                
 202           CONTINUE                                                 
 2022       JR=J-1                                                      
            WR=WWR                                                      
C-                                                                      
            END IF                                                      
C-                                                                      
 6666    CONTINUE                                                       
C-                                                                      
C------ IF BANDWIDTH IS TOO SMALL NO POINT IN SMOOTHING INTERVALL       
         IF(JL.GT.JR) THEN                                              
            F(I)=0.                                                     
C-                                                                      
               ELSE                                                     
C-                                                                      
C------ NOW THE SUMS ARE BUILT THAT ARE NEEDED TO COMPUTE THE ESTIMATE  
            CALL FINDEN(SW,IORD,NKE,NUE,KORD,IBOUN,F(I),C,ICALL)        
            IF(NUE.EQ.0) THEN                                           
               F(I)=F(I)/(WID*N)                                        
               IF(NBO.GT.2.AND.XNOR.NE.1.D0) F(I)=F(I)/XNOR             
               IF(F(I).LT.0.) F(I)=0.                                   
                  ELSE                                                  
               F(I)=F(I)/(N*WID**(NUE+1))                               
               END IF                                                   
            END IF                                                      
C------ NEW INITIALISATION ?                                            
         IF(JL.GT.JR.OR.WWL.GT.WR) INIT=0                               
         WIDO=WID                                                       
C-                                                                      
 100     CONTINUE                                                       
C-                                                                      
      RETURN                                                            
      END SUBROUTINE                                                              
      SUBROUTINE DSUDEN(SW,A1,A2,IORD,SOBS,T,B,IFLOP)                   
C-----------------------------------------------------------------------
C       VERSION: NOVEMBER 22, 1988                                      
C       FORTRAN77 ANSI X3.9-1978/ISO 1539-1980                          
C       COPYRIGHT: STATCOM                                              
C                                                                       
C       PURPOSE:                                                        
C                                                                       
C       COMPUTES NEW LEGENDRE SUMS (FOR DENSITY ESTIMATION)             
C                                                                       
C       PARAMETERS:                                                     
C                    **************   INPUT   *******************       
C                                                                       
C        SW(0:IORD):   OLD SUM OF DATA WEIGHTS FOR LEGENDRE POLYNOM.    
C        A1(6)     :   CONSTANTS OF RECURSIVE FORMULA FOR LEGENDRE POL. 
C        A2(6)     :                             "                      
C        IORD      :   ORDER OF KERNEL POLYNOMIAL                       
C        SOBS      :   OBSERVATION                                      
C        T         :   POINT WHERE THE DENSITY IS TO BE ESTIMATED       
C        B         :   BANDWIDTH                                        
C        IFLOP     :   1: ADDITION, ELSE SUBTRACTION                    
C                                                                       
C                                                                       
C                    **************   OUTPUT   *******************      
C                                                                       
C        SW(0:IORD):   NEW SUM OF DATA WEIGHTS FOR LEGENDRE POLYNOM.    
C                                                                       
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DOUBLE PRECISION SW(0:IORD),P(6),A1(6),A2(6)                      
C-                                                                      
C------  COMPUTE LEGENDRE POLYNOMIALS                                   
      P(1)=(T-SOBS)/B                                                   
      IF(IORD.GT.1) P(2)=1.5D0*P(1)*P(1)-.5D0                           
      DO 1 K=3,IORD                                                     
 1       P(K)=A1(K)*P(K-1)*P(1)+A2(K)*P(K-2)                            
C-                                                                      
C------  COMPUTE NEW LEGENDRE SUMS                                      
C-                                                                      
      IF(IFLOP.EQ.1) THEN                                               
C-                                                                      
         SW(0)=SW(0)+1.D0                                               
         DO 2 K=1,IORD                                                  
 2          SW(K)=SW(K)+P(K)                                            
C-                                                                      
            ELSE                                                        
C-                                                                      
         SW(0)=SW(0)-1.D0                                               
         DO 3 K=1,IORD                                                  
 3          SW(K)=SW(K)-P(K)                                            
C-                                                                      
         END IF                                                         
C-                                                                      
      RETURN                                                            
      END SUBROUTINE                                                              
      SUBROUTINE FALPHA(F,M,ALPHA,G)                                    
C-----------------------------------------------------------------------
C       VERSION: JULY 17, 1987                                          
C       FORTRAN77 ANSI X3.9-1978/ISO 1539-1980                          
C       COPYRIGHT: STATCOM                                              
C                                                                       
C       PURPOSE:                                                        
C                                                                       
C       COMPUTATION OF BANDWIDTH FACTOR SEQUENCE FOR KERNEL SMOOTHING   
C       WITH BANDWIDTH PROPORTIONAL TO DENSITY**(-ALPHA)                
C       THEN BANDWIDTH SEQUENCE IS B(T)=B*G(T) (SUBROUTINE BANVAR)      
C       KNN MUST BE USED BEFORE TO GET THE F**-1 SEQUENCE               
C       OUTPUT GRID IS ASSUMED TO BE EQUIDISTANT                        
C                                                                       
C  PARAMETERS :                                                         
C                                                                       
C  INPUT    F(M)         KNN-SEQUENCE FROM SUBROUTINE KNN               
C  INPUT    M            NUMBER OF POINTS TO ESTIMATE                   
C  INPUT    ALPHA        EXPONENT                                       
C  OUTPUT   G(M)         BANDWIDTH FACTOR SEQUENCE                      
C                                                                       
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION F(M),G(M)                                               
C-                                                                      
C------  COMPUTE BANDWIDTHS PROP. TO DENSITY**(-ALPHA) AND RIEMANN SUM  
      S=0.                                                              
      DO 2 I=1,M                                                        
         IF(F(I).LE..00001) F(I)=.00001                                 
         G(I)=F(I)**ALPHA                                               
         IF(I.GT.1.AND.I.LT.M) S=S+G(I)                                 
2        CONTINUE                                                       
C-                                                                      
C------  NORMALIZE BANDWIDTHS BY RIEMANN SUM                            
      CO=(M-1.)/(S+.5*(G(1)+G(M)))                                      
      DO 3 I=1,M                                                        
3        G(I)=G(I)*CO                                                   
C-                                                                      
      RETURN                                                            
      END SUBROUTINE                                                              
      SUBROUTINE FINDEN(SW,IORD,NKE,NUE,KORD,IBOUN,F,C,ICALL)           
C------------------------------------------------------------------     
C       VERSION: APRIL 16, 1988                                         
C       FORTRAN77 ANSI X3.9-1978/ISO 1539-1980                          
C       COPYRIGHT: STATCOM                                              
C                                                                       
C       PURPOSE:                                                        
C                                                                       
C-      FINAL COMPUTATION OF A DENSITY VALUE VIA LEGENDRE-POLYNOMIALS   
C-                                                                      
C-      PARAMETERS :                                                    
C-                         **************   INPUT   ******************* 
C-                                                                      
C               SW(0:IORD) : SUM OF DATA WEIGHTS FOR LEGENDRE-POLYN.    
C               IORD       : HIGHEST ORDER OF POLYNOMIAL                
C               NKE        : TYP OF KERNEL                              
C               NUE        : ORDER OF DERIVATIVE                        
C               KORD       : ORDER OF KERNEL                            
C               IBOUN      : ~0 BOUNDARY KERNEL, 0 INTERIOR KERNEL      
C               C(0:IORD)  : SEQUENCE OF POLYN. COEFF. FOR BOUND. KERNEL
C               ICALL      : PARAMETER USED TO INITIALISE COMPUTATION   
C                             OF A MATRIX                               
C-                                                                      
C-                         **************   OUTPUT   ****************** 
C-                                                                      
C               F          :  COMPUTED ETIMATE                          
C-                                                                      
C--------------------------------------------------------------------   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DOUBLE PRECISION SW(0:IORD),C(0:IORD),A(0:6,0:6)                  
      SAVE                                                              
C-                                                                      
C------- DEFINITION OF LEGENDRE COEFFICIENTS FOR BOUNDARY               
      IF(ICALL.EQ.0.AND.IBOUN.NE.0) THEN                                
         A(0,2)=1./3.                                                   
                        A(2,2)=2./3.                                    
               A(1,3)=.6                                                
                              A(3,3)=.4                                 
         A(0,4)=7./35.                                                  
                        A(2,4)=4./7.                                    
                                        A(4,4)=8./35.                   
               A(1,5)=27./63.                                           
                              A(3,5)=28./63.                            
                                              A(5,5)=8./63.             
         A(0,6)=33./231.                                                
                        A(2,6)=110./231.                                
                                        A(4,6)=72./231.                 
                                                       A(6,6)=16./231.  
         ICALL=1                                                        
         END IF                                                         
C-                                                                      
      IF(IBOUN.NE.0) THEN                                               
C-                                                                      
C----->    BOUNDARY! WE HAVE TO EVALUATE THE LEGENDRE REPRESENTATION    
C----->    OF THE BOUNDARY KERNEL                                       
C-                                                                      
C------- COMPUTATION OF THE DENSITY VALUE AT BOUNDARY                   
C-                                                                      
         F=C(0)*SW(0)+2.*C(1)*SW(1)                                     
         DO 1 J=2,IORD                                                  
            WW=A(J,J)*SW(J)                                             
            DO 2 I=J-2,0,-2                                             
2              WW=WW+A(I,J)*SW(I)                                       
            F=F+(J+1.)*C(J)*WW                                          
1           CONTINUE                                                    
C-                                                                      
            ELSE                                                        
C-                                                                      
C------- COMPUTATION OF THE DENSITY VALUE AT INTERIOR                   
C-                                                                      
         IF(NKE.EQ.2) THEN                                              
            IF(NUE.EQ.0) THEN                                           
               IF(KORD.EQ.2) F=.5*(SW(0)-SW(2))                         
               IF(KORD.EQ.4) F=.5*SW(0)-1.25*SW(2)+.75*SW(4)            
               IF(KORD.EQ.6) F=.5*SW(0)-1.25*SW(2)+1.6875*SW(4)         
     .                          -.9375*SW(6)                            
               END IF                                                   
            IF(NUE.EQ.1) THEN                                           
               IF(KORD.EQ.3) F=-1.5*(SW(1)-SW(3))                       
               IF(KORD.EQ.5) F=-1.5*SW(1)+5.25*SW(3)-3.75*SW(5)         
               END IF                                                   
            IF(NUE.EQ.2) THEN                                           
               IF(KORD.EQ.4) F=7.5*(SW(2)-SW(4))                        
               IF(KORD.EQ.6) F=7.5*SW(2)-33.75*SW(4)+26.25*SW(6)        
               END IF                                                   
            IF(NUE.EQ.3) F=-52.5*(SW(3)-SW(5))                          
            IF(NUE.EQ.4) F=472.5*(SW(4)-SW(6))                          
C-                                                                      
               ELSE                                                     
C-                                                                      
            IF(NUE.EQ.0) THEN                                           
               IF(KORD.EQ.2) F=.5*SW(0)                                 
               IF(KORD.EQ.4) F=.5*SW(0)-1.25*SW(2)                      
               IF(KORD.EQ.6) F=.5*SW(0)-1.25*SW(2)+1.6875*SW(4)         
               END IF                                                   
            IF(NUE.EQ.1) THEN                                           
               IF(KORD.EQ.3) F=-1.5*SW(1)                               
               IF(KORD.EQ.5) F=-1.5*SW(1)+5.25*SW(3)                    
               END IF                                                   
            IF(NUE.EQ.2) THEN                                           
               IF(KORD.EQ.4) F=7.5*SW(2)                                
               IF(KORD.EQ.6) F=7.5*SW(2)-33.75*SW(4)                    
               END IF                                                   
            IF(NUE.EQ.3) F=-52.5*SW(3)                                  
            IF(NUE.EQ.4) F=472.5*SW(4)                                  
            END IF                                                      
C-                                                                      
         END IF                                                         
C-                                                                      
      RETURN                                                            
      END SUBROUTINE                                                              
      SUBROUTINE HHDD(X,N,XIQR,K,HHD)
C-----------------------------------------------------------------------
C       VERSION: January 1992
C
C
C                                                                       
C       PURPOSE:                                                        
C                                                                       
C       GENERAL SUBROUTINE FOR AUTOMATIC BANDWIDTH SELECTION WITH 
C       HEIDELBERG METHOD for estimation of k-th derivative
C
C       COMPUTATION OF OPTIMAL BANDWIDTH USING THE ASYMPTOTIC FORMULA   
C       WITH ITERATIVE REFINEMENT, FOR A Normal Kernel
C       ITERATION METHOD,  iterations
C       Bandwidth inflation for estimating functional:
C               A = CO*H*N**2/35 or correspondigly for derivatives
C       (A BANDWIDTH FOR ESTIMATING FUNCTIONAL)                
C                                                       
C  PARAMETERS
C                                                                       
C  INPUT    X(N)         DATA (MUST BE SORTED)                          
C  INPUT    N            LENGTH OF X                                    
C  INPUT    XIQR         INTERQUARTILE RANGE
C  INPUT    K		 ORDER OF DERIVATIVE (0, 1 or 2)
C
C  OUTPUT   HHD          ESTIMATED OPTIMAL BANDWIDTH 
C                                                                       
C                                                                       
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION X(N)
        pi=3.141592645d0
        rt2pi = dsqrt(2.d0*pi)
        xn=dble(n)
	m=k+2
	l=m+2
C-                                                                      
C-                                                                      
C-------  COMPUTE CONSTANTS FOR LATER USE                               
c
        xiqr=x(nint(.75*xn))-x(nint(.25*xn)+1)
c         ITER=3*m-l+1
	iter=3
         EXT=1./(2.*m+1)
         EXINF=2./(2*m+1)/(l+m+1.)       
c
c	if to estimate density itself (k=0)
c
	if (k.eq.0) then
c
c-      **estimate inflation constant c
          h2=(.920*xiqr)/(xn**(1.d0/7.d0))
          h3=(.912*xiqr)/(xn**(1.d0/9.d0))
          s2 = 0.d0
          s3 = 0.d0
          do 20 i = 1, n-1
           do 30 j = i+1, n
                     d2 = ((x(i) - x(j))/h2)**2
                     d3 = ((x(i) - x(j))/h3)**2
                     e2 = dexp(-0.5d0*d2)
                     e3 = dexp(-0.5d0*d3)
                     s2 = s2 + (d2**2 - 6.d0*d2 + 3.d0)*e2
               s3 = s3 + (d3**3 - 15.d0*d3**2 + 45.d0*d3 - 15.d0)*e3
30          continue
20         continue
          rhat2 = (2.d0*s2)/((xn**2)*(h2**5)*rt2pi)
          rhat2 = rhat2 + 3.d0/(rt2pi*xn*(h2**5))
          rhat3 = (-2.d0*s3)/((xn**2)*(h3**7)*rt2pi)
          rhat3 = rhat3 + 15.d0/(rt2pi*xn*(h3**7))
          chat = 1.459058157d0*(rhat2**(1.d0/5.d0))/(rhat3**(1.d0/7.d0))
c	write(*,*) 'rhat2=',rhat2, '  rhat3=',rhat3,'  chat=',chat
c
c
      CO1=chat*(XN**EXINF)
C-
C-------  COMPUTE CONSTANT OF ASYMPTOTIC FORMULA                        
       CONST=1.d0/(2.d0*dsqrt(pi))
	a=1.132795764/rhat3**(1.d0/7.d0)*xn**(-1.d0/7.d0)
C-
C------  LOOP OVER ITERATIONS                                          
      DO 100 IT=1,ITER
C-                                                                     
C-------  ESTIMATE Functional
c
        s=0.d0                                                    
        do 40 i = 1, n-1
           do 50 j = i+1, n
                     d2 = ((x(i) - x(j))/a)**2
                     e2 = dexp(-0.5d0*d2)
                     s = s + (d2**2 - 6.d0*d2 + 3.d0)*e2
50          continue
40       continue
        r2 = (2.d0*s)/((xn**2)*(a**5)*rt2pi)
        r2 = r2 + 3.d0/(rt2pi*xn*(a**5))
C-                                                                     
C-
C-                                                                     
C-------  ESTIMATE BANDWIDTH BY ASYMPTOTIC FORMULA                     
         HHD=(CONST/(R2*XN))**(.2d0)                                    
          A=HHD*CO1
100      continue
c	write(*,*) 'Final rhat2=',r2
	endif
c
c	if to estimate first derivative (k=1)
c
	if (k.eq.1) then
c
c-      **estimate inflation constant c
c
          h3=(.912*xiqr)/(xn**(1.d0/9.d0))
          h4=(.914*xiqr)/(xn**(1.d0/11.d0))
          s3 = 0.d0
          s4 = 0.d0
          do 21 i = 1, n-1
           do 31 j = i+1, n
                     d3 = ((x(i) - x(j))/h3)**2
                     d4 = ((x(i) - x(j))/h4)**2
                     e3 = dexp(-0.5d0*d3)
                     e4 = dexp(-0.5d0*d4)
               s3 = s3 + (d3**3 - 15.d0*d3**2 + 45.d0*d3 - 15.d0)*e3
       s4=s4+ (d4**4-28.d0*d4**3+210.d0*d4**2-420.d0*d4 +105.d0)*e4
31          continue
21         continue
          rhat3 = (-2.d0*s3)/((xn**2)*(h3**7)*rt2pi)
          rhat3 = rhat3 + 15.d0/(rt2pi*xn*(h3**7))
          rhat4 = (2.d0*s4)/((xn**2)*(h4**9)*rt2pi)
          rhat4 = rhat4 + 105.d0/(rt2pi*xn*(h4**9))
        chat = 1.489839163d0*(rhat3**(1.d0/7.d0))/(rhat4**(1.d0/9.d0))
c
c
      CO1=chat*(XN**EXINF)
C-
C-------  COMPUTE CONSTANT OF ASYMPTOTIC FORMULA                        
       CONST=3.d0/4.d0/dsqrt(pi)
	a=1.317592936/rhat4**(1.d0/9.d0)*xn**(-1.d0/9.d0)
C-
C------  LOOP OVER ITERATIONS                                          
      DO 101 IT=1,ITER
C-                                                                     
C-------  ESTIMATE Functional
c
        s=0.d0                                                    
        do 41 i = 1, n-1
           do 51 j = i+1, n
                     d3 = ((x(i) - x(j))/a)**2
                     e3 = dexp(-0.5d0*d3)
                     s = s + (d3**3 -15.d0*d3**2 + 45.d0*d3 -15.d0)*e3
51          continue
41       continue
        r3 = (-2.d0*s)/((xn**2)*(a**7)*rt2pi)
        r3 = r3 + 15.d0/(rt2pi*xn*(a**7))
C-                                                                     
C-
C-                                                                     
C-------  ESTIMATE BANDWIDTH BY ASYMPTOTIC FORMULA                     
         HHD=(CONST/(R3*XN))**(1.d0/7.d0)                               
          A=HHD*CO1
101	continue
c	write(*,*) 'rhat3=',r3
c
	endif
	if (k.eq.2) then
c
c-      **estimate inflation constant c
          h4=(.914*xiqr)/(xn**(1.d0/11.d0))
          h5=(.919*xiqr)/(xn**(1.d0/13.d0))
          s4 = 0.d0
          s5 = 0.d0
          do 22 i = 1, n-1
           do 32 j = i+1, n
                     d4 = ((x(i) - x(j))/h4)**2
                     d5 = ((x(i) - x(j))/h5)**2
                     e4 = dexp(-0.5d0*d4)
                     e5 = dexp(-0.5d0*d5)
	             s4 = s4+ (d4**4 - 28.d0*d4**3 + 210.d0*d4**2 -
     .  			420.d0*d4 +105.d0)*e4
                     s5 = s5 + (d5**5 - 45.d0*d5**4 + 630.d0*d5**3 -
     . 		      		 3150.d0*d5**2 +4725.d0*d5-945.d0)*e5
32          continue
22         continue
          rhat4 = (2.d0*s4)/((xn**2)*(h4**9)*rt2pi)
          rhat4 = rhat4 + 105.d0/(rt2pi*xn*(h4**9))
          rhat5 = (-2.d0*s5)/((xn**2)*(h5**11)*rt2pi)
          rhat5 = rhat5 + 945.d0/(rt2pi*xn*(h5**11))
        chat = 1.679304212d0*(rhat4**(1.d0/9.d0))/(rhat5**(1.d0/11.d0))
c	write(*,*) 'rhat4=',rhat4, '  rhat5=',rhat5,'  chat=',chat
c
      CO1=chat*(XN**EXINF)
C-
C-------  COMPUTE CONSTANT OF ASYMPTOTIC FORMULA                        
       CONST=15.d0/8.d0/dsqrt(pi)
	a=1.495649885/rhat5**(1.d0/11.d0)*xn**(-1.d0/11.d0)
C-
C------  LOOP OVER ITERATIONS                                          
      DO 102 IT=1,ITER
C-                                                                     
C-------  ESTIMATE Functional
c
        s=0.d0                                                    
        do 42 i = 1, n-1
           do 52 j = i+1, n
                     d4 = ((x(i) - x(j))/a)**2
                     e4 = dexp(-0.5d0*d4)
                     s = s + (d4**4 -28.d0*d4**3 + 210.d0*d4**2 -
     . 			          420.d0*d4 + 105.d0)*e4
52          continue
42       continue
        r4 = (2.d0*s)/((xn**2)*(a**9)*rt2pi)
        r4 = r4 + 105.d0/(rt2pi*xn*(a**9))
C-                                                                     
C-
C-                                                                     
C-------  ESTIMATE BANDWIDTH BY ASYMPTOTIC FORMULA                     
         HHD=(CONST/(R4*XN))**(1.d0/9.d0)                               
          A=HHD*CO1
102	continue
c	write(*,*) 'rhat4=' ,r4
c
	endif
C-                    
C-                                                                     
      RETURN                                                           
      END SUBROUTINE                                                             
      SUBROUTINE KERDEN(X,N,B,NKE,NUE,KORD,NBO,NY,BO,Z,M,F)             
C-----------------------------------------------------------------------
C       VERSION: APRIL 29, 1988                                         
C       FORTRAN77 ANSI X3.9-1978/ISO 1539-1980                          
C       COPYRIGHT: STATCOM                                              
C                                                                       
C       PURPOSE:                                                        
C                                                                       
C       SUBROUTINE FOR KERNEL ESTIMATION OF A DENSITY OR ITS DERIVATIVE,
C       CHOOSES BETWEEN CONVENTIONAL AND O(N) ALGORITHM                 
C                                                                       
C  INPUT    X(N)         DATA (MUST BE SORTED)                          
C  INPUT    N            LENGTH OF X                                    
C  INPUT    B            ONE-SIDED BANDWIDTH, FOR NY~=0 MEAN BANDWIDTH  
C  INPUT    NKE          TYPE OF KERNEL (1: MIN. VAR., 2: OPTIMAL)      
C  INPUT    NUE          ORDER OF DERIVATIVE (0-4)                      
C  INPUT    KORD         ORDER OF KERNEL (<=6)                          
C  INPUT    NBO          TREATMENT OF BOUNDARY                          
C                        (0: NO BOUNDARY TREATMENT,                     
C                         1: BIAS AS IN INTERIOR, BANDWIDTH SHRINKING,  
C                         2: BIAS AS IN INTERIOR, BANDWIDTH STATIONARY, 
C                         3: JUST NORMALIZING, BANDWIDTH SHRINKING,     
C                         4: JUST NORMALIZING, BANDWIDTH STATIONARY)    
C  INPUT    NY           0: FIX BANDWIDTH, ELSE VARIABLE BANDWIDTH IN F 
C  INPUT    BO(2)        BO(1) LEFT BOUNDARY, BO(2) RIGHT BOUNDARY      
C                        BO(1)<=X(1), BO(2)>=X(N) (DUMMY FOR NBO=0)     
C  INPUT    Z(M)         OUTPUT GRID (ORDERED) WHERE DENS. IS ESTIMATED 
C  INPUT    M            NUMBER OF POINTS TO ESTIMATE                   
C  INPUT    F(M)         BANDWITH SEQUENCE FOR NY~=0, DUMMY FOR NY=0    
C  OUTPUT   F(M)         ESTIMATED DENSITY OR DERIVATIVE OF DENSITY     
C                                                                       
C  EXTERNALS :  DENCON, DENLEG                                          
C                                                                       
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION X(*),Z(*),F(*),BO(2)                                    
C-                                                                      
C------  COMPUTE CHANGING POINT FOR BANDWIDTH                           
      IF(NY.EQ.0) CHAN=15.*LOG10(1.+SQRT(KORD/FLOAT(N)))/SQRT(FLOAT(M)) 
      IF(NY.NE.0) CHAN=20.*LOG10(1.+SQRT(KORD/FLOAT(N)))/SQRT(FLOAT(M)) 
C-                                                                      
C------  CALL ROUTINE FOR DENSITY ESTIMATION                            
      IF(B.LT.CHAN*(Z(M)-Z(1))) THEN                                    
         CALL DENCON(X,N,B,NKE,NUE,KORD,NBO,NY,BO,Z,M,F)                
            ELSE                                                        
         CALL DENLEG(X,N,B,NKE,NUE,KORD,NBO,NY,BO,Z,M,F)                
         END IF                                                         
C-                                                                      
      RETURN                                                            
      END SUBROUTINE                                                              
      SUBROUTINE KNN(T,N,B,S,TT,M,Y)                                    
C-----------------------------------------------------------------------
C       VERSION: APRIL 4, 1987                                          
C       FORTRAN77 ANSI X3.9-1978/ISO 1539-1980                          
C       COPYRIGHT: STATCOM                                              
C                                                                       
C       PURPOSE:                                                        
C                                                                       
C       COMPUTATION OF BANDWIDTH SEQUENCE FOR K-NEAREST-NEIGHBOUR       
C       KERNEL SMOOTHING                                                
C                                                                       
C  PARAMETERS :                                                         
C                                                                       
C  INPUT    T(N)         INPUT GRID                                     
C  INPUT    N            LENGTH OF T                                    
C  INPUT    B            ONE-SIDED BANDWIDTH (IF GRID WERE EQUIDISTANT) 
C  WORK     S(0:N)       HALF POINT INTERPOLATION SEQUENCE              
C  INPUT    TT(M)        OUTPUT GRID                                    
C  INPUT    M            NUMBER OF POINTS TO ESTIMATE                   
C  OUTPUT   Y(M)         BANDWIDTH SEQUENCE                             
C                                                                       
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION T(N),S(0:N),TT(M),Y(M)                                  
C-                                                                      
C------  COMPUTE S-SEQUENCE                                             
      RA=(T(N)-T(1))*.5/(N-1.)                                          
      S(0)=T(1)-RA                                                      
      DO 1 J=1,N-1                                                      
1        S(J)=(T(J)+T(J+1))*.5                                          
      S(N)=T(N)+RA                                                      
C-                                                                      
C------  COMPUTE K                                                      
      XK=B/RA                                                           
      IF(XK.GT.N) XK=N                                                  
      IF(XK.LT.1.) XK=1.                                                
      K=INT(XK)                                                         
C-                                                                      
      J=0                                                               
C-                                                                      
C------  LOOP OVER OUTPUT GRID                                          
      DO 100 I=1,M                                                      
C-                                                                      
C------  COMPUTE K NEAREST S(J1),...,S(JK) TO TT(I)                     
2        IF(S(J+1).LE.TT(I)) THEN                                       
            J=J+1                                                       
            GOTO 2                                                      
            END IF                                                      
         J1=J                                                           
         IF(TT(I)-S(J).GT.S(J+1)-TT(I)) J1=J+1                          
         JK=J1                                                          
         DO 3 KK=2,K                                                    
            IF(J1.EQ.0) JK=JK+1                                         
            IF(JK.EQ.N) J1=J1-1                                         
            IF(J1.EQ.0.OR.JK.EQ.N) GOTO 3                               
            IF(TT(I)-S(J1-1).LE.S(JK+1)-TT(I)) THEN                     
               J1=J1-1                                                  
                  ELSE                                                  
               JK=JK+1                                                  
               END IF                                                   
3           CONTINUE                                                    
C-                                                                      
C------  COMPUTE BANDWIDTH FOR INTERIOR                                 
         XKK=XK-K+1.                                                    
4        IF(J1.GT.0.AND.JK.LT.N) THEN                                   
            H1=S(J1)-S(J1-1)                                            
            HK=S(JK+1)-S(JK)                                            
            Y(I)=((TT(I)-S(J1))*HK+H1*HK*XKK+(S(JK)-TT(I))*H1)/(HK+H1)  
            IF(TT(I)-Y(I).LT.S(J1-1)) THEN                              
               J1=J1-1                                                  
               XKK=XKK-1.                                               
               GOTO 4                                                   
               END IF                                                   
            IF(TT(I)+Y(I).GT.S(JK+1)) THEN                              
               JK=JK+1                                                  
               XKK=XKK-1.                                               
               GOTO 4                                                   
               END IF                                                   
            END IF                                                      
         IF(XKK.EQ.XK-K.AND.J1.EQ.0) JK=JK-1                            
         IF(XKK.EQ.XK-K.AND.JK.EQ.N) J1=J1+1                            
C-                                                                      
C------  COMPUTE BANDWIDTH FOR BOUNDARY                                 
         IF(J1.EQ.0) THEN                                               
            IF(JK.LT.N) JK=JK+1                                         
            IF(JK.EQ.N) JK=JK-1                                         
            Y(I)=S(JK)-TT(I)+(XK-K)*(S(JK+1)-S(JK))                     
            Y(I)=(Y(I)+TT(I)-S(0))*.5                                   
            END IF                                                      
         IF(JK.EQ.N) THEN                                               
            IF(J1.GT.0) J1=J1-1                                         
            IF(J1.EQ.0) J1=J1+1                                         
            Y(I)=TT(I)-S(J1)+(XK-K)*(S(J1)-S(J1-1))                     
            Y(I)=(Y(I)+S(N)-TT(I))*.5                                   
            END IF                                                      
100      CONTINUE                                                       
      RETURN                                                            
      END SUBROUTINE                                                              
      SUBROUTINE LEGDEN(SW,IORD,D,Q)                                    
C------------------------------------------------------------------     
C       VERSION: APRIL 20, 1988                                         
C       FORTRAN77 ANSI X3.9-1978/ISO 1539-1980                          
C       COPYRIGHT: STATCOM                                              
C                                                                       
C       PURPOSE:                                                        
C                                                                       
C-      UPDATE OF SW-SEQUENCE ACCORDING TO NEW BANDWIDTH AND NEW DATA   
C-      (VERSION FOR DENSITY ESTIMATION)                                
C-                                                                      
C       PARAMETERS :                                                    
C-                         **************   INPUT   ******************* 
C-                                                                      
C               SW(0:IORD) :  SUM OF DATA WEIGHTS FOR LEGENDRE-POLYN.   
C               IORD       :  HIGHEST ORDER OF POLYNOMIAL               
C               D          :  DIST. TO THE NEXT POINT DIVIDED BY BANDW. 
C               Q          :  NEW BANDWIDTH DIVIDED BY OLD BANDWIDTH    
C-                                                                      
C-                                                                      
C-                         **************   OUTPUT   ****************** 
C-                                                                      
C               SW(0:IORD) :  UPDATED VERSION OF SW                     
C-                                                                      
C-  REMARK : SUBROUTINE CHECKS WHETHER D HAS CHANGED SINCE THE PREVIOUS 
C-           CALL, DOLD IS SAVED                                        
C-                                                                      
C---------------------------------------------------------------------  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DOUBLE PRECISION A(0:6,0:6),C(0:6,0:6),SW(0:IORD)                 
      SAVE                                                              
C-                                                                      
C- AENDERUNG 22.2.90: DATA DOLD/0./ WG. PROBLEMEN BEI WATFOR-COMPILER   
      DATA DOLD/0./                                                     
C-                                                                      
C- BUILD UP MATRIX                                                      
      IF(DOLD.NE.D) THEN                                                
         DOLD=D                                                         
         DD=D*D                                                         
C-                                                                      
         IF(IORD.EQ.6) THEN                                             
            C(6,5)=11.D0*D                                              
            C(6,4)=49.5D0*DD                                            
            C(6,3)=(115.5D0*DD+7.D0)*D                                  
            C(6,2)=(144.375D0*DD+45.D0)*DD                              
            C(6,1)=((86.625D0*DD+94.5D0)*DD+3.D0)*D                     
            C(6,0)=((14.4375D0*DD+52.5D0)*DD+10.5D0)*DD                 
            END IF                                                      
C-                                                                      
         IF(IORD.GE.5) THEN                                             
            C(5,4)=9.D0*D                                               
            C(5,3)=31.5D0*DD                                            
            C(5,2)=(52.5D0*DD+5.D0)*D                                   
            C(5,1)=(39.375D0*DD+21.D0)*DD                               
            C(5,0)=((7.875D0*DD+17.5D0)*DD+1.D0)*D                      
            END IF                                                      
C-                                                                      
         IF(IORD.GE.4) THEN                                             
            C(4,3)=7.D0*D                                               
            C(4,2)=17.5D0*DD                                            
            C(4,1)=(17.5D0*DD+3.D0)*D                                   
            C(4,0)=(4.375D0*DD+5.D0)*DD                                 
            END IF                                                      
C-                                                                      
         IF(IORD.GE.3) THEN                                             
            C(3,2)=5.D0*D                                               
            C(3,1)=7.5D0*DD                                             
            C(3,0)=(2.5D0*DD+1.D0)*D                                    
            END IF                                                      
C-                                                                      
         C(2,1)=3.D0*D                                                  
         C(2,0)=1.5D0*DD                                                
C-                                                                      
         C(1,0)=D                                                       
C-                                                                      
         END IF                                                         
C-                                                                      
      IF(Q.LT..9999.OR.Q.GT.1.0001) THEN                                
C-                                                                      
C- BUILT UP MATRIX A=P*Q*P**-1                                          
C-                                                                      
         A(0,0)=1.D0                                                    
         DO 1 K=0,IORD-1                                                
1           A(K+1,K+1)=A(K,K)*Q                                         
C-                                                                      
         WW=Q*Q-1.D0                                                    
         DO 2 K=0,IORD-2                                                
            A(K+2,K)=(K+.5D0)*WW                                        
            WW=WW*Q                                                     
2        CONTINUE                                                       
C-                                                                      
         IF(IORD.GE.4) THEN                                             
            QQ=A(2,2)                                                   
            A(4,0)=.375D0+QQ*(-1.25D0+QQ*.875D0)                        
            END IF                                                      
         IF(IORD.GE.5) A(5,1)=Q*(1.875D0+QQ*(-5.25D0+QQ*3.375D0))       
         IF(IORD.EQ.6) THEN                                             
            A(6,0)=-.3125D0+QQ*(2.1875D0+QQ*(-3.9375D0+QQ*2.0625D0))    
            A(6,2)=QQ*(4.375D0+QQ*(-11.25D0+QQ*6.875D0))                
            END IF                                                      
C-                                                                      
C- COMPUTE A*C AND NEW LEGENDRE SUMS                                    
C-                                                                      
         DO 10 I=IORD,0,-1                                              
            XX=0.                                                       
            DO 20 K=0,I                                                 
               WW=0.                                                    
               DO 30 L=K,I-1,2                                          
 30               WW=WW+A(L,K)*C(I,L)                                   
               IF(MOD(I-K,2).EQ.0) WW=WW+A(I,K)                         
               XX=XX+WW*SW(K)                                           
 20            CONTINUE                                                 
            SW(I)=XX                                                    
 10         CONTINUE                                                    
C-                                                                      
               ELSE                                                     
C-                                                                      
         DO 111 I=IORD,1,-1                                             
            DO 112 K=0,I-1                                              
 112           SW(I)=SW(I)+C(I,K)*SW(K)                                 
 111        CONTINUE                                                    
         END IF                                                         
C-                                                                      
      RETURN                                                            
      END SUBROUTINE                                                              
      SUBROUTINE norden(X,N,H,NUE,Z,M,F)            
C-----------------------------------------------------------------------
C       VERSION: September 1991
C                                                                       
C       PURPOSE:                                                        
C                                                                       
C       SUBROUTINE FOR KERNEL ESTIMATION OF A DENSITY or its derivative
C       based on a normal kernel
C                                                                       
C  INPUT    X(N)         DATA (MUST BE SORTED)                          
C  INPUT    N            LENGTH OF X                                    
C  INPUT    h            ONE-SIDED BANDWIDTH, FOR NY~=0 MEAN BANDWIDTH  
C  INPUT    nue          Derivate (0 - 4)
C  INPUT    Z(M)         OUTPUT GRID (ORDERED) WHERE DENS. IS ESTIMATED 
C  INPUT    M            NUMBER OF POINTS TO ESTIMATE                   
C  OUTPUT   F(M)         ESTIMATED DENSITY 
C                                                                       
C                                                                       
C-----------------------------------------------------------------------
        implicit double precision (a-h,o-z)
        dimension f(*),x(*),z(*)
C-                                                                      
C-      
        c=1./dsqrt(2.d0*3.141592654d0)
        xn=dble(n)
C-                                                                      
C-------  LOOP OVER OUTPUT GRID                                         
      DO 100 I=1,M                                                      
         S=0.d0
         do 200 j=1,n
                t=(z(i)-x(j))/h
                if(nue.eq.0) s=s+dexp(-t*t/2.)
                if(nue.eq.1) s=s-t*dexp(-t*t/2.)
                if(nue.eq.2) s=s+(t*t-1)*dexp(-t*t/2.)
                if(nue.eq.3) s=s+(-t*t*t+3.*t)*dexp(-t*t/2.)
        if(nue.eq.4) s=s+(t*t*t*t-6.*t*t+3.)*dexp(-t*t/2.)
200             continue
                f(i)=c*s/(xn*h**(nue+1))
C-                                                                      
100      CONTINUE                                                       
C-                                                                      
      RETURN                                                            
      END SUBROUTINE                                                      


      SUBROUTINE BSORT(F,N)
C------- SHELL SORT ALGORITHM CACM JULY 1964
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION F(N)
      M=N
20       M=M/2
         K=N-M
         DO 40 J=1,K
            L=J
50             IF(F(L+M).GE.F(L)) GOTO 40
               X=F(L+M)
               F(L+M)=F(L)
               F(L)=X
               IF(L.LE.M) GOTO 40
               L=L-M
               GOTO 50
40          CONTINUE
         IF(M.GT.1) GOTO 20
      RETURN
      END SUBROUTINE
      SUBROUTINE BOPTRU(X,N,tu,to,Z,M,
     .  neu,xtrue,bopt,bmin,bmax,fopt,xerr)
C-----------------------------------------------------------------------
C       VERSION  1995
C       calculates true optmal bandwidth for
c       normal kernel density estimation
C
C-----------------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),Z(M),XTrue(M),Fopt(M),XI(2),V(2)
      DATA MAXIT/60/
      eps=0.0001
C
C-------  COMPUTE BOUNDARY INDICES FOR ERROR COMPUTATION
c
      IU=0
      DO 1 I=1,M
         IF(Z(I).LT.TU.OR.Z(I).GT.TO) GOTO 1
         IF(IU.EQ.0) IU=I
         IO=I
1        CONTINUE
      MM=IO-IU+1
C
C-------  GRID SEARCH
c
      B=BMIN
      grid=1.2
      DO 100 J=1,9999
         CALL norden(X,N,b,NUE,Z,M,Fopt)
         CALL L2RIEM(Z(IU),xTrue(IU),fopt(IU),MM,XI(1))
         IF(J.EQ.1.OR.XI(1).LT.XERR) THEN
            BOPt=B
            XERR=XI(1)
            END If
         B=grid*B
         IF(B.GE.BMAX) GOTO 150
100      CONTINUE
C
C-------  GOLDEN SECTION SEARCH
150   B1=BOPT/grid/1.1
      B2=1.1*grid*BOPT
      IF(B1.LT.BMIN) B1=BMIN
      IF(B2.GT.BMAX) B2=BMAX
      CALL GOLSEC(B1,B2,V,XI,0,IV)
      DO 160 I=1,2
         CALL norden(X,N,V(I),NUE,Z,M,Fopt)
         CALL L2RIEM(Z(IU),xTrue(IU),fopt(IU),MM,XI(I))
160      CONTINUE
      DO 170 IT=1,MAXIT
         CALL GOLSEC(B1,B2,V,XI,IT,I)
         CALL norden(X,N,V(I),NUE,Z,M,Fopt)
         CALL L2RIEM(Z(IU),xTrue(IU),fopt(IU),MM,XI(I))
         IF(((V(2)-V(1).LT.EPS)
     .           .AND.((XI(2)-XI(1))/(XI(2)+XI(1)).LT.EPS))
     .           .OR.IT.EQ.MAXIT) THEN
            XERR=MIN(XI(1),XI(2))
            BOPT=V(1)
            IF(XERR.EQ.XI(2)) BOPT=V(2)
            GOTO 200
            END IF
170      CONTINUE
C
C       SMOOTHING WITH OPTIMAL BANDWIDTH
C
200      CALL norden(X,N,BOPT,NUE,Z,M,Fopt)
C
      RETURN
      END SUBROUTINE
      SUBROUTINE Bcvnor(X,N,Z,M,bopt,bmin,bmax,fopt)
C-----------------------------------------------------------------------
C       VERSION  1995
C       calculates true optmal bandwidth for
c       normal kernel density estimation
C
C-----------------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),Z(M),Fopt(M),XI(2),V(2)
      DATA MAXIT/40/
      eps=0.0001
C
C-------  GRID SEARCH
c
      B=BMIN
      nue=0
      grid=1.3
      DO 100 J=1,9999
         CALL icv(X,n,b,z,m,fopt,xi(1))
         IF(J.EQ.1.OR.XI(1).LT.XERR) THEN
            BOPt=B
            XERR=XI(1)
            END IF
         B=grid*B
         IF(B.GE.BMAX) GOTO 150
100      CONTINUE
C
C-------  GOLDEN SECTION SEARCH
150   B1=BOPT/grid/grid
      B2=grid*grid*BOPT
      IF(B1.LT.BMIN) B1=BMIN
      IF(B2.GT.BMAX) B2=BMAX
      CALL GOLSEC(B1,B2,V,XI,0,IV)
      DO 160 I=1,2
         CALL icv(X,n,V(I),z,m,fopt,xi(I))
160      CONTINUE
      DO 170 IT=1,MAXIT
         CALL GOLSEC(B1,B2,V,XI,IT,I)
         CALL icv(X,n,V(I),z,m,fopt,xi(I))
         IF(((V(2)-V(1).LT.EPS)
     .           .AND.((XI(2)-XI(1))/(XI(2)+XI(1)).LT.EPS))
     .           .OR.IT.EQ.MAXIT) THEN
            XERR=MIN(XI(1),XI(2))
            BOPT=V(1)
            IF(XERR.EQ.XI(2)) BOPT=V(2)
            GOTO 200
            END IF
170      CONTINUE
C
C       SMOOTHING WITH OPTIMAL BANDWIDTH
C
200      CALL norden(X,N,BOPT,NUE,Z,M,Fopt)
C
      RETURN
      END SUBROUTINE
      SUBROUTINE icv(X,N,h,z,m,f,xlscv)
C-----------------------------------------------------------------------
C       VERSION: December 1995
C
C
c      Test of cross-validation function for bandwidth h
c
c      input: x(n) ordered, z(m) equidistant with 
c      z(1) <= x(1) and z(m) >= x(n)
c
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION X(N),z(m),f(m)
        pi=3.141592645d0
        rt2pi = dsqrt(2.d0*pi)
        xn=dble(n)
        nue=0
C-                                                             
	call norden(x,n,h,nue,z,m,f)
        xint=0.0
        xint2=0.0
	sum=0.0
        j0=1            
        w0=2.d0/rt2pi/dble(n-1)/h
c-      
        do 10 i=1,m
           xint=xint+f(i)*f(i)
20         continue
           if(x(j0).le.z(i).and.j0.le.n) then
              j0=j0+1
              sum=sum+f(i)
              if(j0.le.n) goto 20
           end if
10      continue
        xlscv=(xint)*(z(2)-z(1)) - (dble(2)*sum/dble(n-1)-w0)
      RETURN                                                           
      END SUBROUTINE                                                     
      SUBROUTINE L2RIEM(T,X,Y,N,DIST)
C-----------------------------------------------------------------------
C       VERSION  3.4.87  IBM 3090-180
C       FORTRAN77 ANSI X3.9-1978/ISO 1539-1980
C       STATCOM PRODUCTION
C
C       PURPOSE:
C
C       COMPUTATION OF THE MEAN SQUARE DISTANCE BETWEEN TWO FUNCTIONS
C       GRID IS NON-EQUIDISTANT -> RIEMANN SUM
C
C  PARAMETERS :
C
C  INPUT    T(N)         INPUT GRID
C  INPUT    X(N)         FIRST FUNCTION
C  INPUT    Y(N)         SECOND FUNCTION
C  INPUT    N            LENGTH OF T,X,Y
C  OUTPUT   DIST         MEAN SQUARE DISTANCE
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T(N),X(N),Y(N)
      DIST=(X(1)-Y(1))*(X(1)-Y(1))*(T(2)-T(1))
      DO 1 I=2,N-1
1        DIST=DIST+(X(I)-Y(I))*(X(I)-Y(I))*(T(I+1)-T(I-1))
      DIST=.5*(DIST+(X(N)-Y(N))*(X(N)-Y(N))*(T(N)-T(N-1)))/(T(N)-T(1))
      RETURN
      END SUBROUTINE
      SUBROUTINE GOLSEC(A,B,V,F,IT,IV)
C-----------------------------------------------------------------------
C         MINIMUMSUCHE MIT DER METHODE DES GOLDENEN SCHNITTES (EIN STEP)
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(2),F(2)
      C=.382
      IF(IT.EQ.0) THEN
         V(1)=A+C*(B-A)
         V(2)=B-C*(B-A)
         RETURN
         END IF
      IF(F(1).LT.F(2)) THEN
         B=V(2)
         V(2)=V(1)
         F(2)=F(1)
         V(1)=A+C*(B-A)
         IV=1
            ELSE
         A=V(1)
         V(1)=V(2)
         F(1)=F(2)
         V(2)=B-C*(B-A)
         IV=2
         END IF
CC         PRINT *,' V=',V
      RETURN
      END SUBROUTINE




	
	end module                                                   