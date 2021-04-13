C**************************************************************************
C
          PROGRAM GRIDGEN
C 
C   THIS PROGRAM GENERATES BOUNDARY FITTED CURVILINEAR GRID FOR TWO-
C   -DIMENSIONAL ARBITRARY GEOMETRY USING POISSON SOLVER AT A COARSE
C    LEVEL FOLLOWED BY CUBIC SPLINE INTERPOLATION TO OBTAIN THE FINE
C   -GRID FROM THE COARSE GRID
C
C    THE CODE IS DEVELOPED AT THE  CTFD DIVISION , NAL , BANGALORE
C   
C     ENQUIRIES MAY BE DIRECTED TO DR. SEKHAR MAJUMDAR ,
C     SCIENTIST, CTFD DIVISION , NAL , BANGALORE
C                     
C**************************************************************************
C
C	IMPLICIT REAL*8 (A-H,X-Z)
C
         INCLUDE 'comgrid'
      
C
	DIMENSION U(NX,NY),V(NX,NY)     
C

C**************************************************************************
   
C
         if ((inmax.lt.nx).or.(inmax.lt.ny)) stop  
	     PI = 4.D0*DATAN(1.D0)
C
         CALL FILES(1)
C
         CALL INPCOND
C
    	 CALL INDC
	     if((nic.gt.nx).or.(njc.gt.ny))  stop	  
	     if((nif.gt.nx).or.(njf.gt.ny))  stop	  
C
         CALL DEFBND
C
         LEVEL = 1
         CALL BSPEC
C
        IF (BCHECK) THEN
C
             IF( .NOT. HYBRID) THEN
C
           WRITE(30,*)'USER WANTS THE COARSE GRID ONLY'
          WRITE(30,*)'THIS RUN IS MEANT ONLY FOR CHECKING THE BOUNDARY
     >    POINT SPECIFICATIONS'
          WRITE(30,*)'PLEASE CHECK THE FOLLOWING FILES:'
          WRITE(30,*)'FILE  CINFO --> INFORMATION ON GRID SPACINGS AT 
     >    COARSE LEVEL' 
          WRITE(30,*)'FILE XYCOARSE --> X-Y COORDINATES OF 
     >    BOUNDARY POINTS AT COARSE LEVEL'
C
          WRITE(30,*)'     PROGRAM TERMINATED NORMALLY' 
C
	   CALL FILES(0)
C
	   STOP
C     
             ELSE
C
	    CALL INDF
C
            LEVEL = 2
            CALL BSPEC
C
           WRITE(30,*)'USER WANTS BOTH THE COARSE AND FINE
     >     LEVEL GRIDS'
           WRITE(30,*)'THIS RUN IS MEANT ONLY FOR CHECKING THE
     >     BOUNDARY POINT SPECIFICATIONS'
	   WRITE(30,*)'PLEASE CHECK THE FOLLOWING FILES:'
           WRITE(30,*)'FILE  CINFO --> INFORMATION ON GRID SPACINGS
     >     AT COARSE LEVEL'                 
           WRITE(30,*)'FILE XYCOARSE --> X-Y COORDINATES OF
     >     BOUNDARY POINTS AT COARSE LEVEL'
C
           WRITE(30,*)'FILE  FINFO --> INFORMATION ON GRID SPACINGS
     >     AT FINE LEVEL'                   
           WRITE(30,*)'FILE XYFINE --> X-Y COORDINATES OF
     >     BOUNDARY POINTS AT FINE LEVEL'  
C
           WRITE(30,*)' PROGRAM TERMINATED NORMALLY'
C
	   CALL FILES(0)
C
             STOP
C 
             ENDIF
C
         ELSE  
C
         CALL INDAT
C
         CALL PSOLVE

         CALL ORTHCHK

         CALL GRIDCHK  
C
         WRITE(13,*)'TITLE = "GRIDCOARSE"  '
         WRITE(13,*)'VARIABLES = X , Y  '
         WRITE(13,*)'ZONE T = "ZONETITLE" I =  ',NIC,',  J =  ',NJC
C
         DO J=1,NJ
         DO I=1,NI
         U(I,J)=X(I,J)
         V(I,J)=Y(I,J)
         WRITE(13,*)U(I,J),V(I,J)
         END DO
         END DO
C
        CALL PRINT(X,'X-COARSE')
C
        CALL PRINT(Y,'Y-COARSE')
C
         IF(.NOT. HYBRID ) THEN
C
       WRITE(30,*)'GRID GENERATED SUCCESSFULLY AT COARSE LEVEL
     > USING POISSON SOLVER'
       WRITE(30,*)'NOS OF GRID NODES : NI =',NIC, ',   NJ =',NJC
       WRITE(30,*)'PLEASE CHECK THE FILE <RESIDUE> FOR
     > CONVERGENCE HISTORY'
       WRITE(30,*)'GRID COORDINATES ARE STORED IN FILE <TECPC.DAT>
     > IN TECPLOT FORMAT FOR PLOTTING' 
       WRITE(30,*)'FILE <PRINTXY> STORES THE GRID COORDINATES
     > IN A SUITABLE FORMAT FOR OTHER REFERENCES '  
       WRITE(30,*)'FINE LEVEL GRID NOT REQUIRED'
C
	   DO J=1,NJC
		 DO I=1,NIC
		   U(I,J)=X(I,J)
		   V(I,J)=Y(I,J)
         END DO
       END DO

	   WRITE(11) NIC,NJC
	   WRITE(11)((U(I,J),I=1,NIC),J=1,NJC),
     >          ((V(I,J),I=1,NIC),J=1,NJC)
C
C
       REWIND(11)
       WRITE(30,*)'PROGRAM TERMINATED NORMALLY'
C
C
	   CALL FILES(0)
C
		 STOP         
C
		 ELSE
C
		 WRITE(30,*)'GRID GENERATED SUCCESSFULLY AT COARSE LEVEL
     > USING POISSON SOLVER'    
       WRITE(30,*)'NOS OF GRID NODES : NI =',NIC, ',   NJ =',NJC
       WRITE(30,*)'PLEASE CHECK THE FILE <RESIDUE> FOR
     > CONVERGENCE HISTORY'
       WRITE(30,*)'GRID COORDINATES ARE STORED IN FILE <TECPC.DAT>
     > IN TECPLOT FORMAT FOR PLOTTING' 
       WRITE(30,*)'FILE <PRINTXY> STORES THE GRID COORDINATES'
       WRITE(30,*)'CUBIC SPLINE INTERPOLATION USED FOR GENERATION OF
     > FINE LEVEL GRID FROM COARSE LEVEL DATA'
C
	 DO J=1,NJ
	 DO I=1,NI
	 XC(I,J)=X(I,J)
	 YC(I,J)=Y(I,J)
         END DO
	 END DO
C
         CALL INDF
C
	 NI=NIF   
	 NJ=NJF    
C
         DO J = 1,NJF
	 DO I = 1,NIF
	 X(I,J) =0.
	 Y(I,J) =0.
	 END DO
	 END DO
C
         LEVEL = 2
         CALL BSPEC
C
C
         CALL CTOF
C
        WRITE(12,*)'TITLE = "GRIDFINE"  '
		WRITE(12,*)'VARIABLES = X , Y  '
        WRITE(12,*)'ZONE T = "ZONETITLE" I =  ',NIF,',  J =  ',NJF
C
         CALL GRIDCHK  
         
         DO J=1,NJ
         DO I=1,NI
         U(I,J)=X(I,J)
         V(I,J)=Y(I,J)
         END DO
         END DO
         DO J=1,NJ
         DO I=1,NI
         WRITE(12,*)U(I,J),v(i,j)
         END DO
         END DO
C
		 CALL PRINT(X,'X-FINE')
C
		 CALL PRINT(Y,'Y-FINE')
C
       WRITE(30,*)'GRID COORDINATES ARE STORED IN FILE <TECPF.DAT>
     > IN TECPLOT FORMAT FOR PLOTTING' 
       WRITE(30,*)'FILE <PRINTXY> STORES THE GRID COORDINATES'
       WRITE(30,*)'NOS OF FINE GRID NODES : NI =',NIF, ',   NJ =',NJF
C
	   DO J=1,NJF
		 DO I=1,NIF
		   U(I,J)=X(I,J)
		   V(I,J)=Y(I,J)
         END DO
       END DO

	   WRITE(11) NIF,NJF
	   WRITE(11)((U(I,J),I=1,NIF),J=1,NJF),
     >          ((V(I,J),I=1,NIF),J=1,NJF)
C
C
       WRITE(30,*)'PROGRAM TERMINATED NORMALLY'
C
2000     CONTINUE
C
         ENDIF
C
       ENDIF
C
           CALL FILES(0)
C
              STOP
C
	      END
C
C**************************************************************************
C       
        SUBROUTINE BSPEC
C
C**************************************************************************
C
C	IMPLICIT REAL*8 (A-H,O-Z)
C
         INCLUDE 'comgrid'
C
	 DIMENSION XIN(INMAX),YIN(INMAX),AK(INMAX),SARC(INMAX)
         DIMENSION PH(INMAX),AKX(INMAX),AKY(INMAX)
	 DIMENSION XSEG(INMAX,ISEGM,4),YSEG(INMAX,ISEGM,4)
C
C       SET NOS. OF NODES DEPENDING ON LEVEL
C
       IF(LEVEL .EQ. 1) THEN
C
       NI = NIC
       NJ = NJC
C
       ELSE
C
       NI = NIF
       NJ = NJF
C
       ENDIF
C
C    SET INDICES FOR REFERRING TO DIFFERENT SEGMENTS   
C
       DO NBOUND = 1,4

       INDBC(1,NBOUND) = 1
       INDBF(1,NBOUND) = 1
C
       IF(NSEG(NBOUND).EQ.1) GO TO 333
C
       DO NN = 2, NSEG(NBOUND)
       INDBC(NN,NBOUND)=INDBC(NN-1,NBOUND)+NPTC(NN-1,NBOUND)-1
       END DO
C
       IF(LEVEL .EQ. 2) THEN
       DO NN = 2,NSEG(NBOUND)
       INDBF(NN,NBOUND)=INDBF(NN-1,NBOUND)+NPTF(NN-1,NBOUND)-1
       END DO
       ENDIF
C
333    CONTINUE

      END DO
C
        DO NBOUND=1,4
C
C
C  CHECKING COMPATIBILITY OF INPUT DATA FOR EACH BOUNDARY
C
        NTOT=0
C
        DO NN = 1,NSEG(NBOUND)
	IF(LEVEL .EQ. 1) THEN
        NTOT=NTOT+NPTC(NN,NBOUND)
C
	ELSE
C
        NTOT=NTOT+NPTF(NN,NBOUND)
        ENDIF
        END DO
C
        IF(LEVEL .EQ. 1) THEN
        NTOT=NTOT-(NSEG(NBOUND)-1)
        IF(NBOUND.EQ.1.AND.NTOT.NE.NI) GO TO 201
        IF(NBOUND.EQ.3.AND.NTOT.NE.NI) GO TO 202
        IF(NBOUND.EQ.2.AND.NTOT.NE.NJ) GO TO 203
        IF(NBOUND.EQ.4.AND.NTOT.NE.NJ) GO TO 204
	ELSE
        NTOT=NTOT-(NSEG(NBOUND)-1)
        IF(NBOUND.EQ.1.AND.NTOT.NE.NI) GO TO 301
        IF(NBOUND.EQ.3.AND.NTOT.NE.NI) GO TO 302
        IF(NBOUND.EQ.2.AND.NTOT.NE.NJ) GO TO 303
        IF(NBOUND.EQ.4.AND.NTOT.NE.NJ) GO TO 304
	ENDIF
C            
          DO ISEG=1,NSEG(NBOUND)
C
	  IF(LEVEL .EQ. 1) THEN
          NEND=NPTC(ISEG,NBOUND)
	  ELSE
          NEND=NPTF(ISEG,NBOUND)
	  ENDIF
C
C      LAYING OUT POINTS ON LINEAR SEGMENTS
C
             IF (STRAIGHT(ISEG,NBOUND)) THEN
C
             IF(LEVEL .EQ. 1) THEN


	     RSTR=RSTC(ISEG,NBOUND)

	     IF (ISEG .EQ. 1) 
     &       AAFR=AAFBEG(ISEG,NBOUND)

	     ELSE

	     RSTR=RSTF(ISEG,NBOUND)
	     IF (ISEG .EQ. 1) 
     &       AAFR=AAFBEG(ISEG,NBOUND)

	     END IF
C
	     IND=INDX(ISEG,NBOUND)
C	     
             XBEG = XDAT(1,ISEG,NBOUND)
             YBEG = YDAT(1,ISEG,NBOUND)
             XEND = XDAT(2,ISEG,NBOUND)
             YEND = YDAT(2,ISEG,NBOUND)
             XSEG(1,ISEG,NBOUND) = XBEG
             YSEG(1,ISEG,NBOUND) = YBEG
             XSEG(NEND,ISEG,NBOUND) = XEND
             YSEG(NEND,ISEG,NBOUND) = YEND
             DELX = XEND - XBEG
             DELY = YEND - YBEG
             ALTOT = SQRT ( DELX*DELX + DELY* DELY)
             RAL=1./ALTOT
C
             IF( MSTR(ISEG,NBOUND) .EQ. 0 ) THEN
                 CALL STRETCH (NEND,ALTOT,RSTR,PH,0.,IND)
	     ELSEIF( MSTR(ISEG,NBOUND) .EQ. 1 ) THEN
		 IF (ISEG .EQ. 1) THEN
		 WRITE(*,*) 'SEGMENT = ',ISEG,'  AAFR = ',AAFR
                 CALL EXPSTR (NEND,ALTOT,RSTR,PH,0.,AAFR,IND)

	    IF(LEVEL.EQ.1) THEN
		      RAFC(ISEG,NBOUND) = RSTR
	    ELSE
		      RAFF(ISEG,NBOUND) = RSTR
	    ENDIF

	              ABEG = PH(NEND)
		      AAFR = 2.0*PH(NEND) - 3.0*PH(NEND-1)
     &		           + PH(NEND-2)

		 ELSE

		   WRITE(*,*) 'SEGMENT = ',ISEG,'  AAFR = ',AAFR
                      CALL EXPSTR (NEND,ALTOT,RSTR,PH,ABEG,AAFR,IND)

	    IF(LEVEL.EQ.1) THEN
		      RAFC(ISEG,NBOUND) = RSTR
	    ELSE
		      RAFF(ISEG,NBOUND) = RSTR
	    ENDIF

	              ABEG = PH(NEND)
		      AAFR = 2.0*PH(NEND) - 3.0*PH(NEND-1)
     &		           + PH(NEND-2)
		 ENDIF
             ENDIF
C
              IF(LEVEL .EQ. 1) THEN
C
	       SC(1,ISEG,NBOUND)=0.
	       SC(NEND,ISEG,NBOUND)=ALTOT   
	       ELSE
	       SF(1,ISEG,NBOUND)=0.
	       SF(NEND,ISEG,NBOUND)=ALTOT   
	       ENDIF
C
               DO J=2,NEND-1
               XSEG(J,ISEG,NBOUND) = XBEG + DELX*RAL*(PH(J)-PH(1))  
               YSEG(J,ISEG,NBOUND) = YBEG + DELY*RAL*(PH(J)-PH(1)) 
               IF(LEVEL .EQ. 1) THEN
               SC(J,ISEG,NBOUND) = PH(J)        
	       ELSE
               SF(J,ISEG,NBOUND) = PH(J)        
	       ENDIF
               END DO
C
          ELSE
C
C      LAYING OUT POINTS ON NON-LINEAR SEGMENTS
C
             IF(LEVEL .EQ. 1) THEN
	     RSTR=RSTC(ISEG,NBOUND)
	     ELSE
	     RSTR=RSTF(ISEG,NBOUND)
	     ENDIF
	     IND=INDX(ISEG,NBOUND)
C
	     NIN =NDATA(ISEG,NBOUND)
	     DO II = 1,NIN 
	     XIN(II)=XDAT(II,ISEG,NBOUND)
	     YIN(II)=YDAT(II,ISEG,NBOUND)
	     END DO
C
C

C        CALL SPLINE(XIN,YIN,NIN,AK,SARC)  
         CALL HSPLINE(XIN,YIN,NIN,AKX,AKY,SARC)  

C
             ALTOT=SARC(NIN)
C
             XSEG(1,ISEG,NBOUND) = XIN(1)
             YSEG(1,ISEG,NBOUND) = YIN(1)
             XSEG(NEND,ISEG,NBOUND) = XIN(NIN)
             YSEG(NEND,ISEG,NBOUND) = YIN(NIN)
		IF(LEVEL .EQ. 1) THEN
            SC(1,ISEG,NBOUND)=SARC(1)
	    SC(NEND,ISEG,NBOUND)=SARC(NIN)
	    ELSE
	    SF(1,ISEG,NBOUND)=SARC(1)
	    SF(NEND,ISEG,NBOUND)=SARC(NIN)
	    ENDIF
C
             IF( MSTR(ISEG,NBOUND) .EQ. 0 ) THEN
                 CALL STRETCH (NEND,ALTOT,RSTR,PH,0.,IND)
	     ELSEIF( MSTR(ISEG,NBOUND) .EQ. 1 ) THEN
		 IF (ISEG .EQ. 1) THEN
		   WRITE(*,*) 'SEGMENT = ',ISEG,'  AAFR = ',AAFR
                      CALL EXPSTR (NEND,ALTOT,RSTR,PH,0.,AAFR,IND)

	    IF(LEVEL.EQ.1) THEN
		      RAFC(ISEG,NBOUND) = RSTR
	    ELSE
		      RAFF(ISEG,NBOUND) = RSTR
	    ENDIF

	              ABEG = PH(NEND)
		      AAFR = 2.0*PH(NEND) - 3.0*PH(NEND-1)
     &			   + PH(NEND-2)
		 ELSE
		   WRITE(*,*) 'SEGMENT = ',ISEG,'  AAFR = ',AAFR
                      CALL EXPSTR (NEND,ALTOT,RSTR,PH,ABEG,AAFR,IND)

	    IF(LEVEL.EQ.1) THEN
		      RAFC(ISEG,NBOUND) = RSTR
	    ELSE
		      RAFF(ISEG,NBOUND) = RSTR
	    ENDIF
	              ABEG = PH(NEND)
		      AAFR = 2.0*PH(NEND) - 3.0*PH(NEND-1)
     &			   + PH(NEND-2)
		 ENDIF
             ENDIF
C
C
             DO J=2,NEND-1     
	        PHGIV=PH(J)
C
	  CALL SEARCH (NIN,SARC,PHGIV,INTV)
C
	    FACT=(PHGIV-SARC(INTV))/(SARC(INTV+1)-SARC(INTV))
            XSEG(J,ISEG,NBOUND) = FACT*XIN(INTV+1)+(1.-FACT)*XIN(INTV)
            YSEG(J,ISEG,NBOUND) = FACT*YIN(INTV+1)+(1.-FACT)*YIN(INTV)
C	
        	IF(LEVEL .EQ. 1) THEN
                SC(J,ISEG,NBOUND) = PH(J)        
		ELSE
                SF(J,ISEG,NBOUND) = PH(J)        
		ENDIF
             END DO
C
          ENDIF
C
          END DO
C
        END DO
C
C
C      FOR SOUTH BOUNDARY NBOUND = 1 AND J = 1
C      FOR NORTH BOUNDARY NBOUND = 3 AND J = NJ
C
        DO J=1,NJ,NJ-1
        II=0
         IF(J.EQ.1)NBOUND = 1
         IF(J.EQ.NJ)NBOUND = 3
C
            DO I=1,NSEG(NBOUND)
C
 	       IF(LEVEL .EQ. 1) THEN
               KEND=NPTC(I,NBOUND)-1
 	       ELSE
               KEND=NPTF(I,NBOUND)-1
 	       ENDIF
                DO K=1,KEND               
               II=II+1
               X(II,J)=XSEG(K,I,NBOUND)
               Y(II,J)=YSEG(K,I,NBOUND)
               END DO
C
           END DO
C   
          IF(LEVEL .EQ. 1) THEN
         X(NI,J)=XSEG(NPTC(NSEG(NBOUND),NBOUND),NSEG(NBOUND),NBOUND) 
         Y(NI,J)=YSEG(NPTC(NSEG(NBOUND),NBOUND),NSEG(NBOUND),NBOUND)
 	 ELSE
         X(NI,J)=XSEG(NPTF(NSEG(NBOUND),NBOUND),NSEG(NBOUND),NBOUND) 
         Y(NI,J)=YSEG(NPTF(NSEG(NBOUND),NBOUND),NSEG(NBOUND),NBOUND)
 	 ENDIF
C   
C   
 	 END DO
C
C      FOR EAST BOUNDARY NBOUND = 2 AND I = NI
C      FOR WEST BOUNDARY NBOUND = 4 AND I = 1
C
        DO I=1,NI,NI-1
        JJ=0
         IF(I.EQ.1)NBOUND = 4
         IF(I.EQ.NI)NBOUND = 2
C
            DO J=1,NSEG(NBOUND)
 
              IF(LEVEL .EQ. 1) THEN
               KEND = NPTC(J,NBOUND)-1
 	       ELSE
                KEND = NPTF(J,NBOUND)-1
 	       ENDIF
C
               DO K=1,KEND            
               JJ=JJ+1
               X(I,JJ)=XSEG(K,J,NBOUND)
               Y(I,JJ)=YSEG(K,J,NBOUND)
               END DO
C
           END DO
C   
         IF(LEVEL .EQ. 1) THEN
         X(I,NJ)=XSEG(NPTC(NSEG(NBOUND),NBOUND),NSEG(NBOUND),NBOUND) 
         Y(I,NJ)=YSEG(NPTC(NSEG(NBOUND),NBOUND),NSEG(NBOUND),NBOUND)
 	 ELSE
          X(I,NJ)=XSEG(NPTF(NSEG(NBOUND),NBOUND),NSEG(NBOUND),NBOUND) 
          Y(I,NJ)=YSEG(NPTF(NSEG(NBOUND),NBOUND),NSEG(NBOUND),NBOUND)
 	 ENDIF
C   
          END DO
C   
C
C  STORING INFORMATION ABOUT  BOUNDARY COORDINATES ON RELEVANT FILES
C
       	    IF(LEVEL .EQ. 1) THEN
	    WRITE(70,*)'COARSE GRID INFORMATION'
	    WRITE(70,*)'------ ---- -----------'
	    ELSE
	    WRITE(71,*)'FINE GRID INFORMATION'
	    WRITE(71,*)'------ ---- -----------'
	    ENDIF
C
         DO NBOUND = 1,4
C
       	  IF(LEVEL .EQ. 1) THEN
C
C  INFORMATION ON COARSE LEVEL
C
	    IF(NBOUND .EQ. 1) THEN
	    WRITE(70,*)'SOUTH BOUNDARY'
	    WRITE(70,*)'----- --------'
	    ENDIF
C
	    IF(NBOUND .EQ. 2) THEN
	    WRITE(70,*)'EAST BOUNDARY'
	    WRITE(70,*)'---- --------'
	    ENDIF
C
	    IF(NBOUND .EQ. 3) THEN
	    WRITE(70,*)'NORTH BOUNDARY'
	    WRITE(70,*)'----- --------'
	    ENDIF
C
	    IF(NBOUND .EQ. 4) THEN
	    WRITE(70,*)'WEST BOUNDARY'
	    WRITE(70,*)'---- --------'
	    ENDIF
C     
      WRITE(70,*)'SEGMENT   GRID   STRETCH  STRETCH     
     >     GRIDSIZE    GRIDSIZE '
      WRITE(70,*)'  NO     POINTS   RATIO    INDEX     
     >     AT START    AT END '
      WRITE(70,*)'-------  ------  -------  -------    
     >     -------     -------'
C
	 DO ISEG = 1,NSEG(NBOUND)
	 NEND = NPTC(ISEG,NBOUND)
	 AL = SC(2,ISEG,NBOUND)-SC(1,ISEG,NBOUND)
	 AR = SC(NEND,ISEG,NBOUND)-SC(NEND-1,ISEG,NBOUND)
         WRITE(70,10)ISEG,NPTC(ISEG,NBOUND),RSTC(ISEG,NBOUND),
     >              INDX(ISEG,NBOUND),AL,AR
         END DO
C
          ELSE
C
C  INFORMATION ON FINE LEVEL
C
	    IF(NBOUND .EQ. 1) THEN
	    WRITE(71,*)'SOUTH BOUNDARY'
	    WRITE(71,*)'----- --------'
	    ENDIF
C
	    IF(NBOUND .EQ. 2) THEN
	    WRITE(71,*)'EAST BOUNDARY'
	    WRITE(71,*)'---- --------'
	    ENDIF
C
	    IF(NBOUND .EQ. 3) THEN
	    WRITE(71,*)'NORTH BOUNDARY'
	    WRITE(71,*)'----- --------'
	    ENDIF
C
	    IF(NBOUND .EQ. 4) THEN
	    WRITE(71,*)'WEST BOUNDARY'
	    WRITE(71,*)'---- --------'
	    ENDIF
C
      WRITE(71,*)'SEGMENT   GRID   STRETCH  STRETCH      
     >     GRIDSIZE      GRIDSIZE '
      WRITE(71,*)'  NO     POINTS   RATIO    INDEX     
     >     AT START      AT END '
      WRITE(71,*)'-------  ------  -------  -------    
     >     -------       -------'
C
	 DO ISEG = 1,NSEG(NBOUND)
	 NEND = NPTF(ISEG,NBOUND)
	 AL = SF(2,ISEG,NBOUND)-SF(1,ISEG,NBOUND)
	 AR = SF(NEND,ISEG,NBOUND)-SF(NEND-1,ISEG,NBOUND)
         WRITE(71,10)ISEG,NPTF(ISEG,NBOUND),RSTF(ISEG,NBOUND),
     >              INDX(ISEG,NBOUND),AL,AR
         END DO
C
      ENDIF              
C
      END DO
C
C
C
C      STORING X-Y COORDINATES AT THE BOUNDARIES ON RELEVANT FILES
C
C      FOR SOUTH BOUNDARY NBOUND = 1 AND J = 1
C      FOR NORTH BOUNDARY NBOUND = 3 AND J = NJ
C
        DO J=1,NJ,NJ-1
C
		 IF(LEVEL .EQ. 1) THEN
C   
		 IF(J.EQ.1)WRITE(66,*)'SOUTH BOUNDARY'
C   
		 IF(J.EQ.NJ)WRITE(66,*)'NORTH BOUNDARY'
C   
	 WRITE(66,*)'   J-POINTS  I-POINTS      
     > X-COORD             Y-COORD'
C   
	 DO I = 1,NI
	 WRITE(66,30)J,I,X(I,J),Y(I,J)
         END DO 
C   
         ELSE
C   
	 IF(J.EQ.1)WRITE(67,*)'SOUTH BOUNDARY'
	 IF(J.EQ.NJ)WRITE(67,*)'NORTH BOUNDARY'
C   
	 WRITE(67,*)'   J-POINTS  I-POINTS      
     > X-COORD             Y-COORD'
C   
	 DO I = 1,NI
	 WRITE(67,30)J,I,X(I,J),Y(I,J)
         END DO
	 ENDIF
C   
	 END DO
C
C      FOR EAST BOUNDARY NBOUND = 2 AND I = NI
C      FOR WEST BOUNDARY NBOUND = 4 AND I = 1
C
        DO I=1,NI,NI-1
C   
	 IF(LEVEL .EQ. 1) THEN
C   
	 IF(I.EQ.1)WRITE(66,*)'WEST BOUNDARY'
C   
	 IF(I.EQ.NI)WRITE(66,*)'EAST BOUNDARY'
C   
	 WRITE(66,*)'   I-POINTS  J-POINTS      
     > X-COORD             Y-COORD'
C   
	 DO J = 1,NJ
	 WRITE(66,30)I,J,X(I,J),Y(I,J)
         END DO 
C   
         ELSE
C   
	 IF(I.EQ.1)WRITE(67,*)'WEST BOUNDARY'
C   
	 IF(I.EQ.NI)WRITE(67,*)'EAST BOUNDARY'
C   
	 WRITE(67,*)'   I-POINTS  J-POINTS      
     > X-COORD             Y-COORD'
C   
	 DO J = 1,NJ
	 WRITE(67,30)I,J,X(I,J),Y(I,J)
         END DO
C   
	 ENDIF
C   
	 END DO
	 RETURN
C
10        FORMAT(I6,2X,I6,2X,F8.2,2X,I6,2X,E15.5,2X,E15.5)
C
30       FORMAT(I8,2X,I8,2X,E18.5,2X,E18.5)
C
201    WRITE(30,*) 'TOTAL NUMBER OF POINTS INCOMPATIBLE TO SEGMENT
     >POINT DATA ALONG J = 1 LINE FOR COARSE LEVEL'
      WRITE(30,*)'PROGRAM TERMINATED ABNORMALLY'
	 RETURN
202   WRITE(30,*) 'TOTAL NUMBER OF POINTS INCOMPATIBLE TO SEGMENT
     >POINT DATA ALONG J = NJ LINE FOR COARSE LEVEL'
      WRITE(30,*)'PROGRAM TERMINATED ABNORMALLY'
	 RETURN
203   WRITE(30,*) 'TOTAL NUMBER OF POINTS INCOMPATIBLE TO SEGMENT
     >POINT DATA ALONG I = NI LINE FOR COARSE LEVEL'
      WRITE(30,*)'PROGRAM TERMINATED ABNORMALLY'
	 RETURN
204   WRITE(30,*) 'TOTAL NUMBER OF POINTS INCOMPATIBLE TO SEGMENT
     >POINT DATA ALONG I = 1 LINE FOR COARSE LEVEL'
      WRITE(30,*)'PROGRAM TERMINATED ABNORMALLY'
	 RETURN
C
C
301   WRITE(30,*) 'TOTAL NUMBER OF POINTS INCOMPATIBLE TO SEGMENT
     >POINT DATA ALONG J = 1 LINE FOR FINE LEVEL'
      WRITE(30,*)'PROGRAM TERMINATED ABNORMALLY'
	 RETURN
302   WRITE(30,*) 'TOTAL NUMBER OF POINTS INCOMPATIBLE TO SEGMENT
     >POINT DATA ALONG J = NJ LINE FOR FINE LEVEL'
      WRITE(30,*)'PROGRAM TERMINATED ABNORMALLY'
	 RETURN
303   WRITE(30,*) 'TOTAL NUMBER OF POINTS INCOMPATIBLE TO SEGMENT
     >POINT DATA ALONG I = NI LINE FOR FINE LEVEL'
      WRITE(30,*)'PROGRAM TERMINATED ABNORMALLY'
	 RETURN
304   WRITE(30,*) 'TOTAL NUMBER OF POINTS INCOMPATIBLE TO SEGMENT
     >POINT DATA ALONG I = 1 LINE FOR FINE LEVEL'
      WRITE(30,*)'PROGRAM TERMINATED ABNORMALLY'
	 RETURN
C
      END
C  
C*********************************************************************
C
	SUBROUTINE STRETCH(NPTS,ALTOT,RAT,PH,PHBEG,INDEX)
C
C*****************************************************************
C
C	IMPLICIT REAL*8 (A-H,O-Z)
C
         INCLUDE 'comgrid'
	DIMENSION PH(INMAX)
C
	PH(1) = PHBEG
	PH(NPTS)=PHBEG+ALTOT
C
	DELX = ALTOT*(RAT-1)/(RAT**(NPTS-1)-1)
C
	IF (INDEX .EQ. 1) THEN
	  DO I=2,NPTS-1
	    PH(I) = PH(I-1)+DELX
	    DELX = DELX*RAT
	  ENDDO
	  ELSE
C
	IF(INDEX .EQ. -1) THEN
	  DO I=NPTS-1,2,-1
	    PH(I) = PH(I+1)-DELX
	    DELX = DELX*RAT
	  ENDDO
	  ELSE
C
	A = ALTOT/2.
	NP = (NPTS+1)/2.
	PH(NP) = A+PHBEG
        DELX = A*(RAT-1)/(RAT**(NP-1)-1)
	  DO I=2,NP
	    PH(I) = PH(I-1)+DELX
	    DELX = DELX*RAT
	  ENDDO
C
        DELX = A*(RAT-1)/(RAT**(NP-1)-1)
	  DO I=NPTS-1,NP+1,-1
	    PH(I) = PH(I+1)-DELX
	    DELX = DELX*RAT
	  ENDDO
C
	  ENDIF
	  ENDIF
C
	RETURN
	END 
C   
C***********************************************************************
C
	SUBROUTINE EXPSTR(NPTS,ALTOT,RAT,PH,PHBEG,DBEG,INDEX)
C
C*****************************************************************
C
C	IMPLICIT REAL*8 (A-H,O-Z)
C
         INCLUDE 'comgrid'
	DIMENSION PH(INMAX)
C
	PH(1)    = PHBEG
	PH(NPTS) = ALTOT
C
	IF (INDEX .EQ. 1) THEN

	ETA_MY = 1.E0/FLOAT(NPTS-1)
	IF (DBEG .NE. 0) THEN
	CALL STEP(ETA_MY,DBEG,ALTOT,RAT)
	ENDIF

        DO I = 2, NPTS-1
	   PH(I) = ALTOT*(EXP(RAT*ETA_MY*FLOAT(I-1)) - 1.E0)
     &           / (EXP(RAT) - 1.E0)
        ENDDO

	ELSEIF (INDEX .EQ. -1) THEN

	ETA_MY = 1.E0/FLOAT(NPTS-1)
	IF (DBEG .NE. 0) THEN
	CALL STEP(ETA_MY,DBEG,ALTOT,RAT)
	ENDIF

        DO I = NPTS-1,2,-1
	   J     = NPTS + 1 - I
	   PH(I) = ALTOT*(EXP(RAT*ETA_MY*FLOAT(J-1)) - 1.E0)
     &           / (EXP(RAT) - 1.E0)
        ENDDO

	ELSEIF (INDEX .EQ. 2) THEN
C
	ALHALF     = ALTOT/2.
	NPHALF     = (NPTS+1)/2
	PH(NPHALF) = ALHALF

	ETA_MY = 1.E0/FLOAT(NPHALF-1)
	IF (DBEG .NE. 0) THEN
	CALL STEP(ETA_MY,DBEG,ALTOT,RAT)
	ENDIF

	  DO I=2,NPHALF-1
	   PH(I) = ALHALF*(EXP(RAT*ETA_MY*FLOAT(I-1)) - 1.E0)
     &           / (EXP(RAT) - 1.E0)
	  ENDDO
C
	  CONST_T = 2.0*PH(NPHALF)

	  DO I=2,NPHALF-1
	   PH(NPHALF+I) = CONST_T - PH(NPHALF-I)
	  ENDDO

	ELSE

	DPHI = ALTOT/FLOAT(NPTS-1)

	DO I = 2, NPTS-1
	PH(I) = PH(I-1) + DPHI
	ENDDO
C
	ENDIF
C
	RETURN
	END 
C   
C***********************************************************************
	SUBROUTINE STEP(ETA,DYM,H,SK)
C	**********************************************************
C	ETA = Non dimensionalized uniform cell distance
C	DYM = Height of the first cell of the particular begining
C	      grid point
C	H   = Total distance
C	SK  = Stretch coefficient
C	**********************************************************
C	IMPLICIT REAL*8(A-H,O-Z)

	SK_OLD = 1.D0

	DO 10 K=1,1000 
	TERM1 =EXP(SK_OLD)
	TERM2 = EXP(SK_OLD*ETA)
	FUNC  = DYM - H*(TERM2 - 1.E0)/(TERM1 - 1.E0)
	TERM3 = (TERM2 - 1.E0)*TERM1/((TERM1 - 1.E0)*(TERM1 - 1.E0))
	TERM4 = ETA*TERM2/(TERM1 - 1.E0)
	DFUNC =	H*(TERM3 - TERM4) 
	IF (ABS(FUNC) .LT. 1.0E-12) THEN
	SK = SK_OLD 
	WRITE(*,*) 'K = ',K,'     SK = ',SK
	RETURN
	ENDIF
	SK_OLD = SK_OLD - FUNC/DFUNC
 10	CONTINUE
	WRITE(*,*) 'ERROR IN FINDING STRETCH FACTOR K'
	RETURN
	END
C***********************************************************************
C
       SUBROUTINE TDMA(AB,AF,B,IPTS,PHI)
C
C***********************************************************************
C
C	IMPLICIT REAL*8 (A-H,O-Z)
C
	INCLUDE 'comgrid'
C
        DIMENSION AB(INMAX),AF(INMAX),B(INMAX)
     >         ,PP(INMAX),QQ(INMAX),PHI(INMAX)
C
C     TRIANGULATION OF THE MATRIX
C
       PP(1) = 0.
       QQ(1) = 0. 
       DO I=2,IPTS-1
       RDE = 1./(1.+AB(I)*PP(I-1))
       PP(I) = -AF(I)*RDE
       QQ(I) = (B(I)-QQ(I-1)*AB(I))*RDE
       ENDDO	
C
C     BACKWARD SUBSTITUTION
C
       DO I=IPTS-1,2,-1
       PHI(I) = PP(I)*PHI(I+1)+QQ(I)
       ENDDO
       RETURN
       END
C
C
C*********************************************************************
C
	  SUBROUTINE SEARCH(NDAT,PH,PHGIV,INTV)
C
C*********************************************************************
C
C	IMPLICIT REAL*8 (A-H,O-Z)
C
          INCLUDE 'comgrid'
C
	  DIMENSION PH(INMAX)
C
	  DO K=1,NDAT-1
	  PH1=PH(K)-PHGIV
	  PH2=PH(K+1)-PHGIV
	  PROD=PH1*PH2
	  IF(PROD.LE.0.) GO TO 100
	  END DO
  100	  INTV=K
          RETURN
          END
C
C****************************************************
C
         SUBROUTINE PSOLVE
C
C****************************************************
C	
C	THIS PROGRAM GENERATES BODY FITTED CURVILINEAR GRIDS 
C	FOR TWO DIMENSIONAL GEOMETRY SOLVING DIFFERENTIAL EQUATIONS
C	COUPLING THE PHYSICAL (X-Y) AND THE COMPUTATIONAL DOMAIN
C	(ZETA-ETA) 
C
C	IMPLICIT REAL*8 (A-H,O-T,X-Z)
C
         INCLUDE 'comgrid'
C
C       DIMENSION U(NX,NY),V(NX,NY),XOLD(NX,NY),YOLD(NX,NY)
        DIMENSION XOLD(NX,NY),YOLD(NX,NY)
C
C       CALL LINTPJ(X,1,NI,1,NJ)
C       CALL LINTPJ(Y,1,NI,1,NJ)
C
C    GENERATE INITIAL GRID USING TRANSFINITE INTERPOLATION
C
         CALL TRANSNEW(X,1,NI,1,NJ)
         CALL TRANSNEW(Y,1,NI,1,NJ) 
C
C   CHECK THE INITIAL GRID FOR NEGATIVE CONTROL VOLUME
C
         CALL GRIDCHK
C
        DO I=1,NI
        THEI(I)=PI/2.
        THEO(I)=PI/2.
        END DO
        DO J=1,NJ
        THZI(J)=PI/2.
        THZO(J)=PI/2.
        END DO
C
	DO I=1,NI
        DO J=1,NJ
        P(I,J) = 0
        Q(I,J) = 0
	ENDDO
	ENDDO 
C 
	IF (MAXIT .EQ. 0) return
       NITER=0
C
C
	   CALL BFIX
C
C     THE ITERATIVE LOOP STARTS HERE
C
100        CONTINUE
        OPEN(UNIT=51,FILE='moni.dat')
C
	   NITER = NITER +1
C       
	   DO J=1,NJ
	      DO I=1,NI
	         XOLD(I,J)=X(I,J)
	         YOLD(I,J)=Y(I,J)
	      END DO
           END DO
C
	   DO J=2,NJ-1
	      DO I=2,NI-1
	         XETA(I,J) = 0.5*(X(I,J+1)-X(I,J-1))
	         XZETA(I,J) = 0.5*(X(I+1,J)-X(I-1,J))
	         YETA(I,J) = 0.5*(Y(I,J+1)-Y(I,J-1))
	         YZETA(I,J) = 0.5*(Y(I+1,J)-Y(I-1,J))
	      ENDDO
	   ENDDO
C 
           CALL BOUND
C 
	   if (niter .gt. liter) CALL CONTROL
C 
       	   CALL CALCX
C 
	   CALL CALCY
C
	   SUMNX=0.
	   SUMNY=0.
	   SUMDX=0.
	   SUMDY=0.
C
C
          DO J=1,NJ
	      DO I=1,NI
	         XCOR=ABS(XOLD(I,J)-X(I,J))
	         XABS=ABS(X(I,J))
	         YCOR=ABS(YOLD(I,J)-Y(I,J))
	         YABS=ABS(Y(I,J))
	         SUMNX=SUMNX+XCOR
	         SUMNY=SUMNY+YCOR
	         SUMDX=SUMDX+XABS
	         SUMDY=SUMDY+YABS  
               END DO
	 END DO
         RESDU(1)=SUMNX/SUMDX
         RESDU(2)=SUMNY/SUMDY
	WRITE(8,*)NITER,RESDU(1),RESDU(2)
	WRITE(51,*)NITER,RESDU(1),RESDU(2)
	RESM = DMAX1(RESDU(1),RESDU(2))
C
        IF (RESM .GT. 1.) THEN
C     
        WRITE(8,*)'NUMERICAL PROCESS IS DIVERGING'
        WRITE(8,*)'(I) CHECK BOUNDARY POINT COORDINATES'
        WRITE(8,*)'(II) TRY LOWER VALUES OF RELAX AND OMEGA
     >  IF DATA IS CORRECT'
        WRITE(8,*)'PROGRAM TERMINATED ABNORMALLY'
C
	   CALL FILES(0)
C
	STOP
C
	ENDIF
C
	close(51)

        IF (RESM .GT. CC .AND. NITER.LT.MAXIT)  GO TO 100
C
        RETURN
	END
C
C*****************************************************
C
	SUBROUTINE CALCX
C
C	IMPLICIT REAL*8 (A-H,O-Z)
C
C*****************************************************
C
        INCLUDE 'comgrid'
C
C	CALCULATE COEFFICIENTS OF THE DIFFERENCE EQUATION
C
	DO J=2,NJ-1
	  DO I=2,NI-1
	    XETA(I,J) = 0.5*(X(I,J+1)-X(I,J-1))
	    XZETA(I,J) = 0.5*(X(I+1,J)-X(I-1,J))
	    YETA(I,J) = 0.5*(Y(I,J+1)-Y(I,J-1))
	    YZETA(I,J) = 0.5*(Y(I+1,J)-Y(I-1,J))
            ALPHA = XETA(I,J)*XETA(I,J)+YETA(I,J)*YETA(I,J)
	    GAMMA = XZETA(I,J)*XZETA(I,J)+YZETA(I,J)*YZETA(I,J)
	    AE(I,J) = ALPHA 
	    AW(I,J) = ALPHA   
	    AS(I,J) = GAMMA  
	    AN(I,J) = GAMMA   
	  ENDDO
	ENDDO
C
C
C   	CALCULATE SOURCE TERM
C
	DO J=2,NJ-1
	  DO I=2,NI-1
            ALPHA = XETA(I,J)*XETA(I,J)+YETA(I,J)*YETA(I,J)
	    GAMMA = XZETA(I,J)*XZETA(I,J)+YZETA(I,J)*YZETA(I,J)
            BETA=XZETA(I,J)*XETA(I,J)+YZETA(I,J)*YETA(I,J)
	    PHE=X(I+1,J)
	    PHW=X(I-1,J)
	    PHN=X(I,J+1)
	    PHS=X(I,J-1)
	   SU(I,J)=-0.5*BETA*(X(I+1,J+1)-X(I-1,J+1)+X(I-1,J-1)
     >            -X(I+1,J-1))
	   SU(I,J)=SU(I,J)+0.5*ALPHA*P(I,J)*(PHE-PHW)
     >                    +0.5*GAMMA*Q(I,J)*(PHN-PHS)
         ENDDO
	ENDDO
C
C	SOLVE FOR X COORDINATE
C
        NV = 1
C	
	CALL SOLVE(X)
C	
	RETURN
	END
C****************************************************
C
	SUBROUTINE CALCY
C
C	IMPLICIT REAL*8 (A-H,O-Z)
C
C*****************************************************
C
        INCLUDE 'comgrid'
C
C	CALCULATE COEFFICIENTS OF THE DIFFERENCE EQUATION
C
	DO J=2,NJ-1
	  DO I=2,NI-1
	    XETA(I,J) = 0.5*(X(I,J+1)-X(I,J-1))
	    XZETA(I,J) = 0.5*(X(I+1,J)-X(I-1,J))
	    YETA(I,J) = 0.5*(Y(I,J+1)-Y(I,J-1))
	    YZETA(I,J) = 0.5*(Y(I+1,J)-Y(I-1,J))
            ALPHA = XETA(I,J)*XETA(I,J)+YETA(I,J)*YETA(I,J)
	    GAMMA = XZETA(I,J)*XZETA(I,J)+YZETA(I,J)*YZETA(I,J)
	    AE(I,J) = ALPHA 
	    AW(I,J) = ALPHA   
	    AS(I,J) = GAMMA  
	    AN(I,J) = GAMMA    
	  ENDDO
	ENDDO
C
C   	CALCULATE SOURCE TERM
C
C
	DO J=2,NJ-1
	  DO I=2,NI-1
            ALPHA = XETA(I,J)*XETA(I,J)+YETA(I,J)*YETA(I,J)
	    GAMMA = XZETA(I,J)*XZETA(I,J)+YZETA(I,J)*YZETA(I,J)
            BETA=XZETA(I,J)*XETA(I,J)+YZETA(I,J)*YETA(I,J)
	    PHE=Y(I+1,J)
	    PHW=Y(I-1,J)
	    PHN=Y(I,J+1)
	    PHS=Y(I,J-1)
	   SU(I,J)=-0.5*BETA*(Y(I+1,J+1)-Y(I-1,J+1)+Y(I-1,J-1)
     >            -Y(I+1,J-1))
	   SU(I,J)=SU(I,J)+0.5 *ALPHA*P(I,J)*(PHE-PHW)
     >                    +0.5*GAMMA*Q(I,J)*(PHN-PHS)
         ENDDO
	ENDDO
C	
C	SOLVE FOR Y COORDINATE
C	
        NV = 2
C
	CALL SOLVE(Y)
C
	RETURN
	END
C**********************************************************************
C      
         SUBROUTINE SOLVE(PHI)
C
C***********************************************************************
C
C	IMPLICIT REAL*8 (A-H,O-Z)
C
C.....THIS ROUTINE IS A LINEAR EQUATION SOLVER USING
C.....A STRONGLY IMPLICIT PROCEDURE DUE TO STONE
C
        INCLUDE 'comgrid'
C
      DIMENSION ANSX(NX,NY),AEX(NX,NY),SPX(NX,NY)
     *,RESP(NX,NY),PHI(NX,NY)
C
C
C***********************************************************************
C
C ---- PARAMETERS FOR THE SOLVER REQUIRED
C
      ALPHA=.9
      NSWPM=10 
      CRIT1=0.1 
C
      IBEG=2
      IFIN=NI-1
      JBEG=2
      JFIN=NJ-1 
C
C ---- INITIALIZE ARRAYS REQUIRED
C
      RS=0.
      DO 151 J=1,NJ
      DO 151 I=1,NI
      ANSX(I,J)=0.
      AEX(I,J)=0.
      SPX(I,J)=0.
      RESP(I,J)=0.
151   CONTINUE
C
C
C ---- CALCULATES THE RESIDUAL UNBALANCE OF THE CONCERNED
C      VARIABLE OVER THE WHOLE DOMAIN OF COMPUTATION
C
      DO 10 J=JBEG,JFIN
      DO 10 I=IBEG,IFIN
      AEIJ=AE(I,J)
      ANIJ=AN(I,J)
      AWIJ=AW(I,J)
      ASIJ=AS(I,J)
      ASUM = AEIJ + AWIJ + ANIJ + ASIJ
      APIJ=ASUM
      RS=RS+(ANIJ*PHI(I,J+1)+ASIJ*PHI(I,J-1)+
     1        AEIJ*PHI(I+1,J)+AWIJ*PHI(I-1,J)-
     1        APIJ*PHI(I,J)+SU(I,J))**2
1111  CONTINUE
      RELAXM=1.0-RELAX(NV)
      APIJ = APIJ/RELAX(NV)
      SU(I,J) = SU(I,J) + RELAXM*PHI(I,J)*APIJ
  10  CONTINUE
C
      RS=SQRT(RS)
      RESDU(NV) = RS
      RSTART=RS
C
C
C ---- THE SOLUTION ALGORITHM IS A RECURSIVE
C      SUBSTITUTION PROCEDURE THAT FOLLOWS
C
      IBM1 = IBEG - 1
      JBM1=  JBEG - 1
C
C      WRITE(8,*)NV
      ITER = 0
      RELRS=1.
  20  ITER = ITER + 1
C
C ----  SWEEP ALONG SW-NE DIRECTION
C
      IF (ITER .GT. NSWPM.OR.RELRS.LE.CRIT1) RETURN
C
      JSUM=JBEG+JFIN
      ISUM=IBEG+IFIN
201     CONTINUE
   
      DO 21 J=JBEG,JFIN
      DO 21 I=IBEG,IFIN
      AWIJ=AW(I,J)
      ASIJ=AS(I,J)
      ANIJ=AN(I,J)
      AEIJ=AE(I,J)
      APIJ=AWIJ+AEIJ+ASIJ+ANIJ
      APIJ=APIJ/RELAX(NV)
      RESIJ=ANIJ*PHI(I,J+1)+ASIJ*PHI(I,J-1)+
     1        AEIJ*PHI(I+1,J)+AWIJ*PHI(I-1,J)-
     1        APIJ*PHI(I,J)+SU(I,J)
      IF (ABS(RESIJ).LT.1D-10) RESIJ = 0.0
      ASX=-ASIJ/(1.+ALPHA*AEX(I,J-1))
      AWX=-AWIJ/(1.+ALPHA*ANSX(I-1,J))
      APX=APIJ+ASX*(ALPHA*AEX(I,J-1)-ANSX(I,J-1))
     1               +AWX*(ALPHA*ANSX(I-1,J)-AEX(I-1,J))
      RAP=1./APX
      AEX(I,J) =(-AEIJ-ALPHA*ASX*AEX(I,J-1))*RAP
      ANSX(I,J)=(-ANIJ-ALPHA*AWX*ANSX(I-1,J))*RAP
C----'SPX' USED FOR STORING UPPER TRIANGLE MATRIX TIMES SOLUTION VECTOR
      SPX(I,J)=(RESIJ-ASX*SPX(I,J-1)-AWX*SPX(I-1,J))*RAP
      RESP(I,J)=RESIJ*RESIJ
   21 CONTINUE
C
C------ CALCULATION OF SOLUTION VECTOR THROUGH BACK SUBSTITUTION
C
      DO 22 JR=JBEG,JFIN
      J=JSUM-JR
      DO 22 IR=IBEG,IFIN
      I=ISUM-IR
      SPX(I,J)=SPX(I,J)-AEX(I,J)*SPX(I+1,J)-ANSX(I,J)*SPX(I,J+1)
      PHI(I,J)=PHI(I,J)+ SPX(I,J)
  22  CONTINUE
   
C
C ---- STOPPING CRITERION FOR INNER ITERATION PROCESS
C
      RS=0.
      DO 33 J=JBEG,JFIN
      DO 33 I=IBEG,IFIN
      RS=RS+RESP(I,J)
33    CONTINUE
      RS=SQRT(RS)
      IF (RSTART.NE.0.) RELRS=RS/RSTART
      RELDIF=ABS(RELRS-RELAST)
      RELAST=RELRS
C      WRITE(28,7778)ITER,RELRS,RS
C
300   CONTINUE
      GO TO 20
7778  FORMAT(10X,'ITER = ',I4,5X,
     *'RELATIVE DEFECT =',E20.8,5X,'ABS. DEFECT =',E20.8)
      END
C*******************************************************************************
C
       SUBROUTINE BILINT(F,IBEG,IFIN,JBEG,JFIN)
C
C*******************************************************************************
C
C	IMPLICIT REAL*8 (A-H,O-Z)
C
        INCLUDE 'comgrid'
C
       DIMENSION F(NX,NY)
C
       IP1=IBEG+1
       IM1=IFIN-1
       JP1=JBEG+1
       JM1=JFIN-1
       DS=1/FLOAT(IFIN-IBEG)
       DT=1/FLOAT(JFIN-JBEG)
C
       DO 100 I=IP1,IM1
       S=(I-IBEG)*DS
       F1S=(1-S)*F(IBEG,JBEG)+S*F(IFIN,JBEG)
       F1E=(1-S)*F(IBEG,JFIN)+S*F(IFIN,JFIN)
       DO 100 J=JP1,JM1
       F1=(1-S)*F(IBEG,J)+S*F(IFIN,J)    
       T=(J-JBEG)*DT
100    F(I,J)=F1+(1-T)*(F(I,JBEG)-F1S)+T*(F(I,JFIN)-F1E)
       RETURN
       END
C
C*******************************************************************************
C
C
       SUBROUTINE PRINT(PHI,HEAD)
C
C	IMPLICIT REAL*8 (A-H,O-Z)
C******************************************************************************
C     PRINTS THE VALUE OF 2-DIMENSIONAL FUNCTION PHI (I,J)
C******************************************************************************
C
        INCLUDE 'comgrid'
C
       DIMENSION PHI(NX,NY),STORE(INMAX)
       CHARACTER*10  HEAD
C 
C
       WRITE(16,210) HEAD
       ISTR =1
       JSTR=1
       ISTA=ISTR-12
100    ISTA=ISTA+12
       IEND=ISTA+11
       IEND= MIN(NI,IEND)
       WRITE(16,111) (I ,I=ISTA,IEND)
       DO 101 JJ=JSTR,NJ
       J=JSTR+NJ-JJ
       DO 120 I=ISTA,IEND
       A=PHI(I,J)
       IF (ABS(A) .LT.1.D-30) A=0.0
120    STORE(I)=A
101    WRITE(16,113) J,(STORE(I) ,I=ISTA,IEND)
       IF(IEND .LT. NI) GO TO 100
       RETURN
111    FORMAT(/1X,6H         ,I3,11I10)
113    FORMAT(1X,I3,1P12E10.2)
210    FORMAT(/1X,27(2H*-),6X,A10,6X,27(2H-*))
       END
C***********************************************
C
	SUBROUTINE BFIX 
C
C	IMPLICIT REAL*8 (A-H,O-Z)
C
C**********************************************************************
C
C       THIS SUBROUTINE COMPUTES THE METRIC COEFFICIENTS AT TWO
C       OF THE ETA BOUNDARIES (ETA = ETAMIN, ETA = ETAMAX) AND
C       AT TWO OF ZETA BOUNDARIES (ZETA = ZETAMIN,ZETA = ZETAMAX)
C
        INCLUDE 'comgrid'
C 
C       INNER BOUNDARY (ETA = 1 & J=1,ZETA = 1 & I=1)
C       OUTER BOUNDARY (ETA =NJ & J=NJ,ZETA = NI & I=NI)		
C
	DO J=1,NJ,NJ-1
	  DO I=2,NI-1
	    XZETA(I,J) = 0.5*(X(I+1,J)-X(I-1,J))
	    YZETA(I,J) = 0.5*(Y(I+1,J)-Y(I-1,J))
	    XZZETA(I,J) = (X(I+1,J)+X(I-1,J)-2*X(I,J))
	    YZZETA(I,J) = (Y(I+1,J)+Y(I-1,J)-2*Y(I,J))
	  ENDDO
	ENDDO
C
	DO I=1,NI,NI-1
	  DO J=2,NJ-1
	    XETA(I,J) = 0.5*(X(I,J+1)-X(I,J-1))
	    YETA(I,J) = 0.5*(Y(I,J+1)-Y(I,J-1))
	    XEETA(I,J) = (X(I,J+1)+X(I,J-1)-2*X(I,J))
	    YEETA(I,J) = (Y(I,J+1)+Y(I,J-1)-2*Y(I,J))
	  ENDDO
	ENDDO
C
	RETURN
	END
C
C***********************************************************************
C
       SUBROUTINE BOUND
C
C***********************************************************************
C       THIS SUBROUTINE COMPUTES THE METRIC COEFFICIENTS AT TWO
C       OF THE ETA BOUNDARIES (ETA = ETAMIN, ETA = ETAMAX) AND
C       AT TWO OF ZETA BOUNDARIES (ZETA = ZETAMIN,ZETA = ZETAMAX)
C
C	IMPLICIT REAL*8 (A-H,O-T,X-Z)
C
        INCLUDE 'comgrid'
C
       J=1
C
C      FIRST & SECOND DERIVATIVE ALONG ETA=CONSTANT
C
       XZETA(1,J) = X(2,J)-X(1,J)
       YZETA(1,J) = Y(2,J)-Y(1,J)
       XZETA(NI,J)= X(NI,J)-X(NI-1,J)
       YZETA(NI,J)= Y(NI,J)-Y(NI-1,J)
C
C      FIRST DERIVATIVE ACROSS ETA=CONSTANT
C
       DO I=1,NI
       SINTH = SIN(THEI(I))
       COSTH = COS(THEI(I))
       GAMMA = XZETA(I,J)*XZETA(I,J)+YZETA(I,J)*YZETA(I,J)
       SQGM=SQRT(GAMMA)
       DX=X(I,J+1)-X(I,J)
       DY=Y(I,J+1)-Y(I,J)
       DIST=SQRT(DX*DX+DY*DY)
       XETA(I,J)=DIST*(XZETA(I,J)*COSTH-YZETA(I,J)*SINTH)/SQGM         
c      XETA(I,J)=-YZETA(I,J)         
       YETA(I,J)=DIST*(YZETA(I,J)*COSTH+XZETA(I,J)*SINTH)/SQGM          
c      YETA(I,J)=XZETA(I,J)         
       END DO
C
C      SECOND & MIXED DERIVATIVES ACROSS ETA=CONSTANT
C
       DO I=2,NI-1
       XEETA(I,J) = ((-7*X(I,J)-X(I,J+2)+8*X(I,J+1))/2)-3*XETA(I,J)
       YEETA(I,J) = ((-7*Y(I,J)-Y(I,J+2)+8*Y(I,J+1))/2)-3*YETA(I,J)
       XZE(I,J) = 0.5*(XETA(I+1,J) - XETA(I-1,J))
       YZE(I,J) = 0.5*(YETA(I+1,J) - YETA(I-1,J))
       END DO
C     
       J=NJ
C
C      FIRST & SECOND DERIVATIVE ALONG ETA=CONSTANT
C
       XZETA(1,J) = X(2,J)-X(1,J)
       YZETA(1,J) = Y(2,J)-Y(1,J)
       XZETA(NI,J)= X(NI,J)-X(NI-1,J)
       YZETA(NI,J)= Y(NI,J)-Y(NI-1,J)
C
C      FIRST DERIVATIVE ACROSS ETA=CONSTANT
C
       DO I=1,NI
       SINTH = SIN(THEO(I))
       COSTH = COS(THEO(I))
       GAMMA = XZETA(I,J)*XZETA(I,J)+YZETA(I,J)*YZETA(I,J)
       SQGM=SQRT(GAMMA)
       DX=X(I,J)-X(I,J-1)
       DY=Y(I,J)-Y(I,J-1)
       DIST=SQRT(DX*DX+DY*DY)
       XETA(I,J)=DIST*(XZETA(I,J)*COSTH-YZETA(I,J)*SINTH)/SQGM         
c      XETA(I,J)=-YZETA(I,J)         
       YETA(I,J)=DIST*(YZETA(I,J)*COSTH+XZETA(I,J)*SINTH)/SQGM          
c      YETA(I,J)=XZETA(I,J)          
       END DO
C
C      SECOND & MIXED DERIVATIVES ACROSS ETA=CONSTANT
C
       DO I=2,NI-1
       XEETA(I,J) = ((-7*X(I,J)-X(I,J-2)+8*X(I,J-1))/2)+3*XETA(I,J)
       YEETA(I,J) = ((-7*Y(I,J)-Y(I,J-2)+8*Y(I,J-1))/2)+3*YETA(I,J)
       XZE(I,J) = 0.5*(XETA(I+1,J) - XETA(I-1,J))
       YZE(I,J) = 0.5*(YETA(I+1,J) - YETA(I-1,J))
       END DO
C
       I=1
C
C      FIRST & SECOND DERIVATIVES ALONG ZETA=CONSTANT
C
       XETA(I,1) = (X(I,2) - X(I,1))
       YETA(I,1) = (Y(I,2) - Y(I,1))
       XETA(I,NJ) = (X(I,NJ) - X(I,NJ-1))
       YETA(I,NJ) = (Y(I,NJ) - Y(I,NJ-1))
C
C      FIRST DERIVATIVE ACROSS ZETA=CONSTANT
C
       DO J=1,NJ
       SINTH = SIN(THZI(J))
       COSTH = COS(THZI(J))
       ALPHA = XETA(I,J)*XETA(I,J)+YETA(I,J)*YETA(I,J)
       SQAL=SQRT(ALPHA) 
       DX=X(I+1,J)-X(I,J)
       DY=Y(I+1,J)-Y(I,J)
       DIST=SQRT(DX*DX+DY*DY)
       YZETA(I,J) =DIST*(YETA(I,J)*COSTH-XETA(I,J)*SINTH)/SQAL        
c      YZETA(I,J) =-XETA(I,J)        
       XZETA(I,J) =DIST*(XETA(I,J)*COSTH+YETA(I,J)*SINTH)/SQAL        
c      XZETA(I,J) =YETA(I,J)        
       END DO
C
C      SECOND & MIXED DERIVATIVES ACROSS ZETA=CONSTANT
C
       DO J=2,NJ-1
       XZZETA(I,J) = ((-7*X(I,J)-X(I+2,J)+8*X(I+1,J))/2)-3*XZETA(I,J)
       YZZETA(I,J) = ((-7*Y(I,J)-Y(I+2,J)+8*Y(I+1,J))/2)-3*YZETA(I,J)
       XZE(I,J) = 0.5*(XZETA(I,J+1) - XZETA(I,J-1))
       YZE(I,J) = 0.5*(YZETA(I,J+1) - YZETA(I,J-1))
       END DO
C
       I=NI
C 
C      FIRST & SECOND DERIVATIVES ALONG ZETA=CONSTANT
C
       XETA(I,1) = (X(I,2) - X(I,1))
       YETA(I,1) = (Y(I,2) - Y(I,1))
       XETA(I,NJ) = (X(I,NJ) - X(I,NJ-1))
       YETA(I,NJ) = (Y(I,NJ) - Y(I,NJ-1))

C
C      FIRST DERIVATIVE ACROSS ZETA=CONSTANT
C
       DO J=1,NJ
       SINTH =SIN(THZO(J))
       COSTH = COS(THZO(J))
       ALPHA = XETA(I,J)*XETA(I,J)+YETA(I,J)*YETA(I,J)
       SQAL=SQRT(ALPHA)
       DX=X(I,J)-X(I-1,J)
       DY=Y(I,J)-Y(I-1,J)
       DIST=SQRT(DX*DX+DY*DY)
       XZETA(I,J) =DIST*(XETA(I,J)*COSTH+YETA(I,J)*SINTH)/SQAL       
c      XZETA(I,J) =YETA(I,J)       
       YZETA(I,J) =DIST*(YETA(I,J)*COSTH-XETA(I,J)*SINTH)/SQAL         
c      YZETA(I,J) =-XETA(I,J)         
       END DO
C
C      SECOND & MIXED DERIVATIVES ACROSS ZETA=CONSTANT
C
       DO J=2,NJ-1
       XZZETA(I,J) = ((-7*X(I,J)-X(I-2,J)+8*X(I-1,J))/2)+3*XZETA(I,J)
       YZZETA(I,J) = ((-7*Y(I,J)-Y(I-2,J)+8*Y(I-1,J))/2)+3*YZETA(I,J)
       XZE(I,J) = 0.5*(XZETA(I,J+1) - XZETA(I,J-1))
       YZE(I,J) = 0.5*(YZETA(I,J+1) - YZETA(I,J-1))
       END DO
C
       RETURN
       END
C***********************************************************************
C
	SUBROUTINE CONTROL
C
C***********************************************************************
C
C       THIS SUBROUTINE COMPUTES THE CONTROL FUNCTIONS P AND Q
C       ACCORDING TO THE IDEA PROPOSED BY THOMAS AND MIDDELECOFF.
C
C	IMPLICIT REAL*8 (A-H,O-Z)
C
        INCLUDE 'comgrid'
C
        RESM =DMAX1(RESDU(1),RESDU(2))
	DO I=1,NI,NI-1
	  DO J=2,NJ-1
            QOLD=Q(I,J)
	  IF(I.EQ.1)THEN
	    SINTH=SIN(THZI(J))
	    COSTH=COS(THZI(J))
	  ENDIF
	  IF(I.EQ.NI)THEN
	    SINTH=SIN(THZO(J))
	    COSTH=COS(THZO(J))
          ENDIF
	    XE = XETA(I,J)
	    YE = YETA(I,J)
	    XZ = XZETA(I,J)
	    YZ = YZETA(I,J)
	    XEE = XEETA(I,J)
	    YEE = YEETA(I,J)
	    XZZ =XZZETA(I,J)
	    YZZ =YZZETA(I,J)
	    XZET=XZE(I,J)
	    YZET=YZE(I,J)
	    ALPHA=XE*XE+YE*YE
	    GAMMA=XZ*XZ+YZ*YZ
	    BETA=XZ*XE+YZ*YE
	    QTHOM = -((XEE*XE+YEE*YE)/ALPHA)
	    QTHOMC = -((XZZ*XE+YZZ*YE)/GAMMA)
	    QBETA =2*BETA*(XZET*XE+YZET*YE)/(ALPHA*GAMMA)
            QN1 = ((XZZ*YE-YZZ*XE)/GAMMA)+((XEE*YE-YEE*XE)/ALPHA)
	    QN2 = -2*BETA*(XZET*YE-YZET*XE)/(ALPHA*GAMMA)
	    QTHETA= (QN1+QN2)*COSTH/SINTH
	    IF (NITER .EQ. 1) THEN
	       Q(I,J)=QTHOM
	    ELSE
	       QTHOMC_NEW = QTHOMC
	       QTHOMC_OLD=Q(I,J)-QTHOM
	       QTHOMC_DIFF=QTHOMC_NEW-QTHOMC_OLD
	       TERM1=OMEGA*ABS(QTHOMC_DIFF)
	       TERM2=FACT_LIM*MAX(ABS(QTHOMC_OLD),1.D+00)
	       QTHOMC_LIM=MIN(TERM1,TERM2)
	       QTHOMC_COR=SIGN(QTHOMC_LIM,QTHOMC_DIFF)
	       QTHOMC_ADD=QTHOMC_OLD+QTHOMC_COR
  	       Q(I,J)=QTHOM+QTHOMC_ADD
c 	       Q(I,J)=QTHOM+QTHOMC
c	       WRITE(223,*)NITER,J,QTHOMC_OLD,QTHOMC_COR,FACT_LIM
c	       if (niter .eq. MAXIT) then
c	       WRITE(223,*)NITER,I,Q(I,J),QTHOM,QTHOMC,qthomc_add
c	       end if
 	    ENDIF
c 	    QNEW = QTHOM+QTHOMC+QBETA+QTHETA 
c 	    QNEW = QTHOM+QTHOMC 
c	    Q(I,J) = QOLD*(1-OMEGA)+OMEGA*QNEW
	  ENDDO
	ENDDO
C
	DO J=1,NJ,NJ-1
	  DO I=2,NI-1
            POLD=P(I,J)
	    IF(J.EQ.1)THEN
	      SINTH=SIN(THEI(I))
	      COSTH=COS(THEI(I))
            ENDIF
	    IF(J.EQ.NJ)THEN
	      SINTH=SIN(THEO(I))
	      COSTH=COS(THEO(I))
            ENDIF
	    XE = XETA(I,J)
	    YE = YETA(I,J)
	    XZ = XZETA(I,J)
	    YZ = YZETA(I,J)
	    XEE = XEETA(I,J)
	    YEE = YEETA(I,J)
	    XZZ =XZZETA(I,J)
	    YZZ =YZZETA(I,J)
	    XZET=XZE(I,J)
	    YZET=YZE(I,J)
	    ALPHA=XE*XE+YE*YE
	    GAMMA=XZ*XZ+YZ*YZ
	    BETA=XZ*XE+YZ*YE
	    PTHOM= -((XZZ*XZ+YZZ*YZ)/GAMMA)
	    PTHOMC=-((XEE*XZ+YEE*YZ)/ALPHA)
	    PBETA=2*BETA*(XZET*XZ+YZET*YZ)/(ALPHA*GAMMA)
	    PN1 =-((XZZ*YZ-YZZ*XZ)/GAMMA)-((XEE*YZ-YEE*XZ)/ALPHA)
            PN2 = 2*BETA*(XZET*YZ-YZET*XZ)/(ALPHA*GAMMA)
	    PTHETA =(PN1+PN2)*COSTH/SINTH
	    IF (NITER .EQ. 1) THEN
	       P(I,J)=PTHOM
	    ELSE
	       PTHOMC_NEW = PTHOMC
	       PTHOMC_OLD=P(I,J)-PTHOM
	       PTHOMC_DIFF=PTHOMC_NEW-PTHOMC_OLD
	       TERM1=OMEGA*ABS(PTHOMC_DIFF)
	       TERM2=FACT_LIM*MAX(ABS(PTHOMC_OLD),1.D+00)
	       PTHOMC_LIM=MIN(TERM1,TERM2)
	       PTHOMC_COR=SIGN(PTHOMC_LIM,PTHOMC_DIFF)
	       PTHOMC_ADD=PTHOMC_OLD+PTHOMC_COR
  	       P(I,J)=PTHOM+PTHOMC_ADD
c 	       P(I,J)=PTHOM+PTHOMC
c	       if (niter .eq. MAXIT) then
c	       WRITE(224,*)NITER,I,P(I,J),PTHOM,PTHOMC,pthomc_add
c	       WRITE(224,*)NITER,I,PTHOMC_OLD,PTHOMC_NEW,PTHOMC_COR
c	       end if
 	    ENDIF
C           IF(NITER .GT.100 .OR. RESM .LT. 5.D-04)THEN
c 	    PNEW = PTHOM+PTHOMC+PBETA+PTHETA 
c 	    PNEW = PTHOM+PTHOMC 
c           P(I,J) = POLD*(1-OMEGA)+OMEGA*PNEW
	  ENDDO
	ENDDO
C
        CALL LINTPJ(P,1,NI,1,NJ)
        CALL LINTPI(Q,1,NI,1,NJ)
C
C      STRETCHING THE P TERM EXPONENTIALLY
C
C            DO I=2,NI-1
C             DO J=2,NJ-1
C
C	    IF(LEVEL.EQ.1) THEN

C	    IF(J.LT.NPTC(2,4)) ISEG = 1
C	    IF(J.GT.NPTC(2,4)) ISEG = 2
C           ELSE 
C
C	    IF(J.LT.NPTF(2,4)) ISEG = 1
C	    IF(J.GT.NPTF(2,4)) ISEG = 2
C
C	    IF(LEVEL .EQ. 1) THEN
C            NEND=NPTC(ISEG,4)
C            NEND=NPTF(ISEG,4)
C	    ENDIF
C
C            XBEG = XDAT(1,ISEG,4)
C            YBEG = YDAT(1,ISEG,4)
C            XEND = XDAT(2,ISEG,4)
C            YEND = YDAT(2,ISEG,4)
C            DELX = XEND - XBEG
C            DELY = YEND - YBEG
C            ALTOT =SQRT ( DELX*DELX + DELY* DELY)
C
C		 IF (ISEG .EQ. 1) THEN
C		       PBEG = P(I,1)
C	     	       AAFR = P(I,2) - P(I,1)
C                       CALL EXPSTR (NEND,ALTOT,RSTR,P,PBEG,AAFR,IND)
C
C		       AAFR = 2.0*P(I,NEND) - 3.0*P(I,NEND-1)
C     &		           + P(I,NEND-2)
C
C		 ELSE
C
C                       CALL EXPSTR (NEND,ALTOT,RSTR,P,AAFR,AAFR,IND)
C
C
C		 ENDIF
C
C	     ENDDO
C	  ENDDO
C
C	CALL EXPPJ(P,1,NI,1,NJ)
C	CALL EXPPI(Q,1,NI,1,NJ)
C
	RETURN
	END
C*********************************************************************
C
       SUBROUTINE INDAT
C
C	IMPLICIT REAL*8 (A-H,O-Z)
C
        INCLUDE 'comgrid'
C
        READ(80,*)
        READ(80,*)MAXIT,liter,CC
        READ(80,*)
        READ(80,*)(RELAX(I),I=1,2)
        READ(80,*)
        READ(80,*)OMEGA,FACT_LIM
        RETURN
        END
C*******************************************************************************
C
       SUBROUTINE LINTPJ(F,IBEG,IFIN,JBEG,JFIN)
C
C	IMPLICIT REAL*8 (A-H,O-Z)
C
        INCLUDE 'comgrid'
C
       DIMENSION F(NX,NY)
C
       IP1=IBEG+1
       IM1=IFIN-1
       JP1=JBEG+1
       JM1=JFIN-1
       DS=1/FLOAT(JFIN-JBEG)
C
       DO 100 I=IP1,IM1
       DO 100 J=JP1,JM1
       S=(J-JBEG)*DS
       F(I,J)=(1-S)*F(I,JBEG)+S*F(I,JFIN)
100    CONTINUE                                        
       RETURN
       END
C
C*******************************************************************************
C
       SUBROUTINE EXPPJ(F,IBEG,IFIN,JBEG,JFIN)
C*******************************************************************************
C
C	IMPLICIT REAL*8 (A-H,O-Z)
C
        INCLUDE 'comgrid'
C       
       DIMENSION F(NX,NY)
C
       IP1=IBEG+1
       IM1=IFIN-1
       JP1=JBEG+1
       JM1=JFIN-1
C
            DO I=IP1,IM1
            DO J=JP1,JM1
	    IF(LEVEL.EQ.1) THEN
	    IF(J.LE.NPTC(1,4)) MSEG = 1 
	    IF(J.GT.NPTC(2,4)) MSEG = 2
	    RAT = RAFC(MSEG,4) 
	    ETA_MY = 1.D0/(FLOAT(NPTC(MSEG,4)-1))
	    ELSE
	    IF(J.LE.NPTF(1,4)) MSEG = 1 
	    IF(J.GT.NPTF(2,4)) MSEG = 2
	    RAT = RAFF(MSEG,4) 
	    ETA_MY = 1.D0/(FLOAT(NPTC(MSEG,4)-1))
	    ENDIF
	     EVAL= ALTOT*(EXP(RAT*ETA_MY*FLOAT(J-1)) - 1.E0)
     &           / (EXP(RAT) - 1.E0)
             F(I,J) = (1-EVAL)*F(I,JBEG)+EVAL*F(I,JFIN)
            ENDDO
          ENDDO
C
      RETURN
C
      END
C
C
C*******************************************************************************
C
       SUBROUTINE LINTPI(F,IBEG,IFIN,JBEG,JFIN)
C*******************************************************************************
C
C	IMPLICIT REAL*8 (A-H,O-Z)
C
        INCLUDE 'comgrid'
C
       DIMENSION F(NX,NY)
C
       IP1=IBEG+1
       IM1=IFIN-1
       JP1=JBEG+1
       JM1=JFIN-1
       DS=1/FLOAT(IFIN-IBEG)
C
       DO 100 J=JP1,JM1
       DO 100 I=IP1,IM1
       S=(I-IBEG)*DS
       F(I,J)=(1-S)*F(IBEG,J)+S*F(IFIN,J)
100    CONTINUE                                        
       RETURN
       END
C
C**************************************************************************
        SUBROUTINE EXPPI(F,IBEG,IFIN,JBEG,JFIN)
C*******************************************************************************
C
C       IMPLICIT REAL*8 (A-H,O-Z)
 
        INCLUDE 'comgrid'
C
        DIMENSION F(NX,NY)
C

        IP1=IBEG+1
        IM1=IFIN-1
        JP1=JBEG+1
        JM1=JFIN-1
 
          DO J=JP1,JM1
             DO I=IP1,IM1
	    IF(LEVEL.EQ.1) THEN
	    IF(J.LE.NPTC(1,4)) MSEG = 1 
	    IF(J.GT.NPTC(2,4)) MSEG = 2
	    RAT = RAFC(MSEG,4) 
	    ETA_MY = 1.E0/(FLOAT(NPTC(MSEG,4)-1))
	    ELSE
	    IF(J.LE.NPTF(1,4)) MSEG = 1 
	    IF(J.GT.NPTF(2,4)) MSEG = 2
	    RAT = RAFF(MSEG,4) 
	    ETA_MY = 1.E0/(FLOAT(NPTC(MSEG,4)-1))
	    ENDIF
	      EVAL= ALTOT*(EXP(RAT*ETA_MY*FLOAT(I-1)) - 1.E0)
     &          / (EXP(RAT) - 1.E0)
             F(I,J) = (1-EVAL)*F(I,JBEG)+EVAL*F(I,JFIN)
             ENDDO
           ENDDO
C
        RETURN
        END
C
C**************************************************************************
       SUBROUTINE INDC
C**************************************************************************
C	IMPLICIT REAL*8 (A-H,O-T,X-Z)
C
        INCLUDE 'comgrid'
C
        READ(50,*)
        READ(50,*)
        READ(50,*)
        READ(50,*)
        READ(50,*)NIC,NJC                
C
        DO NBOUND= 1,4
        READ(50,*)
        READ(50,*)
		READ(50,*)NSEG(NBOUND)
            DO ISEG = 1, NSEG(NBOUND)
            READ(50,*)
            READ(50,*)NPTC(ISEG,NBOUND)
            READ(50,*)
            READ(50,*)
            READ(50,*)STRAIGHT(ISEG,NBOUND)
            READ(50,*)
            READ(50,*) RSTC(ISEG,NBOUND),INDX(ISEG,NBOUND)
            READ(50,*)
            READ(50,*)
            READ(50,*)
            READ(50,*)MSTR(ISEG,NBOUND),AAFBEG(ISEG,NBOUND)
	     END DO
         END DO
C
         RETURN
	 END
C**************************************************************************
C**************************************************************************
       SUBROUTINE INDF    
C**************************************************************************
C	IMPLICIT REAL*8 (A-H,O-T,X-Z)
C
        INCLUDE 'comgrid'
C
        READ(60,*)
        READ(60,*)
        READ(60,*)
        READ(60,*)
        READ(60,*)NIF,NJF                
C
        DO NBOUND= 1,4
        READ(60,*)
        READ(60,*)
		READ(60,*)NSEG(NBOUND)
C
C
            DO ISEG = 1, NSEG(NBOUND)
            READ(60,*)
            READ(60,*)NPTF(ISEG,NBOUND)
            READ(60,*)
            READ(60,*)
            READ(60,*)STRAIGHT(ISEG,NBOUND)
            READ(60,*)
            READ(60,*) RSTF(ISEG,NBOUND),INDX(ISEG,NBOUND)
            READ(60,*)
            READ(60,*)
            READ(60,*)
            READ(60,*)MSTR(ISEG,NBOUND),AAFBEG(ISEG,NBOUND)
            END DO
         END DO
C
         RETURN
		 END
C**************************************************************************
C
       SUBROUTINE DEFBND
C
C**************************************************************************
C		IMPLICIT REAL*8 (A-H,O-T,X-Z)
C
        INCLUDE 'comgrid'
C
        DO NBOUND= 1,4
            DO ISEG = 1, NSEG(NBOUND)
	     READ(7,*)
	     READ(7,*) NDATA(ISEG,NBOUND)
	     READ(7,*)
		 DO I = 1, NDATA(ISEG,NBOUND)
		 READ(7,*)XDAT(I,ISEG,NBOUND),YDAT(I,ISEG,NBOUND)
		 END DO
	     END DO
        END DO
C
        RETURN
        END
C
C**************************************************************************
C
        SUBROUTINE CTOF     
C
C**************************************************************************
C
C	IMPLICIT REAL*8 (A-H,O-T,X-Z)
C
        INCLUDE 'comgrid'
C
	DIMENSION ZIF(INMAX),ZIC(INMAX)
        DIMENSION ETC(INMAX),ETF(INMAX)
	DIMENSION XIN(INMAX),YIN(INMAX),AK(INMAX),SARC(INMAX)   
	 DIMENSION ZC(INMAX),EC(INMAX),ZF(INMAX),EF(INMAX)
	 DIMENSION XCRS(INMAX,INMAX),YCRS(INMAX,INMAX)
	 DIMENSION XFN(INMAX,INMAX),YFN(INMAX,INMAX)
C
C
C  CALCULATION OF ZETA AND ETA AT COARSE LEVEL
C
         DO I=1,NIC   
	 ZIC(I)=FLOAT(I-1)
	 END DO
C
         DO J=1,NJC   
	 ETC(J)=FLOAT(J-1)
	 END DO
C
C    CALCULATION OF ZETA AND ETA AT FINE LEVEL FROM 
C    THE GIVEN ARCLENGTHS ON J = 1 AND I = 1
C    IE., NBOUND = 1 AND NBOUND = 4 AND ALSO THE COARSE
C    LEVEL COORDINATES ON THE SAME BOUNDARIES
C
	 DO NBOUND=1,4,3
C
         DO NN = 1,NSEG(NBOUND)
C
	 NIN = NPTC(NN,NBOUND)
	 NOUT = NPTF(NN,NBOUND)
C
       	 DO J = 1,NIN
	 JJ = INDBC(NN,NBOUND)+J-1
	 XIN(J) = SC(J,NN,NBOUND)
	 IF(NBOUND.EQ.1)YIN(J) = ZIC(JJ)
	 IF(NBOUND.EQ.4)YIN(J) = ETC(JJ)
         END DO
C
         CALL SPLINE(XIN,YIN,NIN,AK,SARC)  
c        CALL HSPLINE(XIN,YIN,NIN,AKX,AKY,SARC)  
C
	 DO J = 1,NOUT
	 XGIV=SF(J,NN,NBOUND)
	 JJ =INDBF(NN,NBOUND)+J-1
C
	 CALL SPLORD(XIN,YIN,NIN,AK,XGIV,YGIV)
C
	 IF(NBOUND.EQ.1)ZIF(JJ) = YGIV      
	 IF(NBOUND.EQ.4)ETF(JJ) = YGIV      
	 END DO
C
         END DO
C
         END DO
C
C  SEGMENTWISE EVALUATION OF FIELD GRID COORDINATES         
C
           DO NETA = 1,NSEG(4)
             DO NZETA=1,NSEG(1)
	       NILC = NPTC(NZETA,1)
	       NJLC = NPTC(NETA,4)
	       NILF = NPTF(NZETA,1)
	       NJLF = NPTF(NETA,4)
C
             DO I=1,NILC
	       II = INDBC(NZETA,1)+I-1
	       ZC(I) = ZIC(II)
	       END DO
C
            DO I = 1,NILF
	      II = INDBF(NZETA,1)+I-1
	      ZF(I) = ZIF(II)
	      END DO
C
            DO J = 1,NJLC 
	      JJ = INDBC(NETA,4)+J-1
	      EC(J) = ETC(JJ)
	      END DO
C
            DO J = 1,NJLF 
	      JJ = INDBF(NETA,4)+J-1
	      EF(J) = ETF(JJ)
	      END DO
C
           DO J= 1,NJLC
	   DO I = 1,NILC
	   II = INDBC(NZETA,1)+I-1
	   JJ = INDBC(NETA,4)+J-1
	   XCRS(I,J) = XC(II,JJ)
	   YCRS(I,J) = YC(II,JJ)
	   END DO
	   END DO
C

         CALL ALGINT(NILC,NJLC,ZC,EC,ZF,EF,NILF,NJLF,
     >                         XCRS,YCRS,XFN,YFN) 
C
           DO J= 1,NJLF
	   DO I = 1,NILF
	   II = INDBF(NZETA,1)+I-1
	   JJ = INDBF(NETA,4)+J-1
	   X(II,JJ) = XFN(I,J)
	   Y(II,JJ) = YFN(I,J)
	   END DO
	   END DO
C
          END DO
	  END DO
C
	 RETURN
	 END
C  ****************************************************

C       
         SUBROUTINE SPLINE(XIN,YIN,IPTS,AK,SARC)         
C
C  ****************************************************
C
C	     IMPLICIT REAL*8 (A-H,O-Z)
C
        INCLUDE 'comgrid'
C
        DIMENSION XIN(INMAX),YIN(INMAX),AK(INMAX),H(INMAX)
     *   ,AB(INMAX),AF(INMAX),B(INMAX),SARC(INMAX)
       	DIMENSION XX(INMAX),YY(INMAX),TEMP(INMAX)
C       
C       CALCULATION OF H
C
	    H(1)=XIN(2)-XIN(1)
    	DO  I=2,IPTS-1
		 swap=0.
	    H(I) = XIN(I+1)-XIN(I)
	    denm = H(I)+H(I-1)
         if(denm.eq.0.or.h(i).eq.0.or.h(i-1).eq.0.) go to 300
    	ENDDO
	    go to 301
300      continue
        DO I=1,IPTS
		TEMP(I)=XIN(I)
		XIN(I)=YIN(I)
		YIN(I)=TEMP(I)
		END DO

		H(1)=XIN(2)-XIN(1)
		DO I=2,IPTS-1
	    H(i)=XIN(I+1)-XIN(I)
		ENDDO  
		swap =1.
301      continue
C
C     CALCULATION OF END CURVATURE AK(I) FOR THE PIECEWISE CUBIC
C     POLYNOMIAL (SPLINE) BETWEEN XIN(I) & XIN(I+1)
C
         DO I=2,IPTS-1
         AB(I) = H(I-1)/(2.*(H(I)+H(I-1)))
         AF(I) = H(I)/(2.*(H(I)+H(I-1)))
         B1=(YIN(I+1)-YIN(I))/(H(I)*(H(I)+H(I-1)))
         B2=(YIN(I)-YIN(I-1))/(H(I-1)*(H(I)+H(I-1)))
         B(I) = 3.*(B1-B2)
         END DO
C
C     BOUNDARY CONDITIONS
C
        AK(1)=0.
        AK(IPTS)=0.
C
        CALL TDMA(AB,AF,B,IPTS,AK)
        SARC(1)=0.
C
        DO I=1,IPTS-1
    	XDIF=XIN(I+1)-XIN(I)
	    XL=XIN(I)
     	XR=XIN(I+1)
	    XX(1)=XIN(I)
	    YY(1)=YIN(I)
	    XX(NDIV+1)=XIN(I+1)
	    YY(NDIV+1)=YIN(I+1)
C
	    DELX=XDIF/FLOAT(NDIV)
C

	    DO J=2,NDIV
	    XX(J)=XX(J-1)+DELX
	     XDIST=XX(J)
          HI = H(I)
          A1 = 1.-(XDIST-XL)/HI
          A2 = (XDIST-XL)/HI
          A3 = (XDIST-XL)*(XDIST-XR)*(2.-(XDIST-XL)/HI)/6.
          A4 = (XDIST-XL)*(XDIST-XR)*(1.+(XDIST-XL)/HI)/6.
          YDIST  = A1*YIN(I)+A2*YIN(I+1)+A3*AK(I)+A4*AK(I+1)
	    YY(J)=YDIST
        END DO
C
   	    SUM=0.
C
	    DO J=2,NDIV+1
	    DELX=XX(J)-XX(J-1)
	    DELY=YY(J)-YY(J-1)
	    DARC=SQRT(DELX*DELX+DELY*DELY)
	    SUM=SUM+DARC
	    END DO
C
          SARC(I+1)=SARC(I)+SUM
C
	    END DO
		IF(SWAP.NE.1) RETURN
        DO I=1,IPTS
		TEMP(I)=XIN(I)
		XIN(I)=YIN(I)
		YIN(I)=TEMP(I)

C
		END DO
          RETURN

	    END
C
C  ****************************************************

C       
         SUBROUTINE HSPLINE(PX,PY,IPTS,AKX,AKY,SARC)         
C
C  ****************************************************
C
C	     IMPLICIT REAL*8 (A-H,O-Z)
C
        INCLUDE 'comgrid'
C
        DIMENSION AKX(INMAX),AKY(INMAX),PX(INMAX),PY(INMAX),fy(0:inmax)
     * ,fx(0:inmax),AB(INMAX),AF(INMAX),BX(INMAX),BY(INMAX),SARC(INMAX)
C       
C
C     CALCULATION OF END CURVATURE AK(I) FOR THE PIECEWISE CUBIC
C     POLYNOMIAL (SPLINE) BETWEEN XIN(I) & XIN(I+1)
C
         DO I=2,IPTS-1
         AB(I) =.25
         AF(I) =.25
         BX(I) =.75*(PX(I+1)-PX(I-1))
         BY(I) = .75*(PY(I+1)-PY(I-1))
		 enddo
C
C     BOUNDARY CONDITIONS
C
        AKX(1)=PX(2)-PX(1)
        AKX(IPTS)=PX(IPTS)-pX(IPTS-1)
        AKY(1)=PY(2)-PY(1)
        AKY(IPTS)=PY(IPTS)-PY(IPTS-1)
C
        CALL TDMA(AB,AF,BX,IPTS,AKX)
        CALL TDMA(AB,AF,BY,IPTS,AKY)
        SARC(1)=0.
C
C
        DO I=1,IPTS-1
		fx(0)=px(i)
		fy(0)=py(i)
	    fX(NDIV)=PX(I+1) 
	    fY(NDIV)=PY(I+1) 

          DO J=1,NDIV-1
		  U=FLOAT(J)/float(NDIV-1) 
          F1 = 2*U**3-3*U*U+1      
          F2 = -2*U**3+3*U*u
          F3 = U**3-2*U*U+U
          F4 =U**3-U*U
          FX(j)= F1*PX(I)+f2*PX(I+1)+F3*AKX(I)+F4*AKX(I+1)     
          FY(j)= F1*PY(I)+f2*PY(I+1)+F3*AKY(I)+F4*AKY(I+1)     
        END DO
		 
        sum=0.

        DO J=1,NDIV
	    DELX=FX(j)-FX(j-1 )
	    DELY=FY(J)-fY(J-1 )
	    DARC=SQRT(DELX*DELX+DELY*DELY)
	    SUM=SUM+DARC
	    END DO
C
          SARC(I+1)=SARC(I)+SUM

C
	    END DO
          RETURN

	    END
C
C  ****************************************************
C       
         SUBROUTINE SPLORD(XIN,YIN,IPTS,AK,XGIV,YGIV)     
C
C  ****************************************************
C
C	IMPLICIT REAL*8 (A-H,O-Z)
C
        INCLUDE 'comgrid'
C
        DIMENSION XIN(INMAX),YIN(INMAX),AK(INMAX),H(INMAX)
C       
C     CALCULATION OF H
C
	DO  I=1,IPTS-1
	  H(I) = XIN(I+1)-XIN(I)
	ENDDO
C
        CALL SEARCH(IPTS,XIN,XGIV,INTV)
C
	XL = XIN(INTV)
	XR = XIN(INTV+1)
	XDIST = XGIV
          HI = H(INTV)
          A1 = 1.-(XDIST-XL)/HI
          A2 = (XDIST-XL)/HI
          A3 = (XDIST-XL)*(XDIST-XR)*(2.-(XDIST-XL)/HI)/6.
          A4 = (XDIST-XL)*(XDIST-XR)*(1.+(XDIST-XL)/HI)/6.
          YGIV = A1*YIN(INTV)+A2*YIN(INTV+1)+A3*AK(INTV)+A4*AK(INTV+1)
C
          RETURN
	  END
C  
C*********************************************************************
C*********************************************************************
C
       SUBROUTINE INPCOND
C
C	IMPLICIT REAL*8 (A-H,O-Z)
C
        INCLUDE 'comgrid'
C
C
        READ(81,*)
        READ(81,*)
        READ(81,*)
        READ(81,*)BCHECK,HYBRID
C
          RETURN
	  END
C************************************************************************
C
          SUBROUTINE FILES(N)
C
C************************************************************************
C
        IF(N .EQ. 1) THEN
C
C       INPUT FILES
C  
		OPEN(UNIT=81,FILE='task')  
        OPEN(UNIT=80,FILE='itcont')
        OPEN(UNIT=7,FILE='bound.dat')
        OPEN(UNIT=50,FILE='gridc.dat')  
        OPEN(UNIT=60,FILE='gridf.dat') 
C
C       OUTPUT FILES
C      
        OPEN(UNIT=30,FILE='message')
        OPEN(UNIT=12,FILE='tecpf.dat')
        OPEN(UNIT=13,FILE='tecpc.dat')
        OPEN(UNIT=70,FILE='cinfo')    
        OPEN(UNIT=71,FILE='finfo')   
        OPEN(UNIT=66,FILE='xycoarse')
        OPEN(UNIT=67,FILE='xyfine')  
		OPEN(UNIT=8,FILE='residue')
		OPEN(UNIT=16,FILE='printxy') 
		OPEN(UNIT=11,FILE='sbgrid',FORM='UNFORMATTED') 
		OPEN(UNIT=52,FILE='orthchk.dat') 
C
         ELSE
C
         CLOSE(16)
         CLOSE(8)
	     CLOSE(67)
		 CLOSE(66)
		 CLOSE(71)
		 CLOSE(70)
		 CLOSE(13)
		 CLOSE(12)
		 CLOSE(30)
		 CLOSE(60) 
		 CLOSE(50) 
		 CLOSE(7) 
		 CLOSE(80)
         CLOSE(81)
         CLOSE(11)
	 CLOSE(52)
C
         ENDIF
         RETURN
	 END
C
************************************************************************
         SUBROUTINE ALGINT(NILC,NJLC,ZC,EC,ZF,EF,NILF,NJLF,
     >                         XCRS,YCRS,XFN,YFN) 
C
C	IMPLICIT REAL*8 (A-H,O-T,X-Z)
C
         INCLUDE 'comgrid'
C
	 DIMENSION XIN(INMAX),YIN(INMAX),AK(INMAX),SARC(INMAX)
	 DIMENSION ZC(INMAX),EC(INMAX),ZF(INMAX),EF(INMAX)
	 DIMENSION XCRS(INMAX,INMAX),YCRS(INMAX,INMAX)
	 DIMENSION XFN(INMAX,INMAX),YFN(INMAX,INMAX)
C
C  CALCULATION OF X AND Y AT ALL FINE ZETA LEVEL 
C  FOR THE COARSE LEVEL VALUES OF ALL ETA
C  THIS IS THE FIRST STEP OF THE BICUBIC INTERPOLATION
C  ALONG ZETA ONLY
C
       	DO J=1,NJLC
C
       NIN=NILC
C
	   DO I=1,NILC    
	   XIN(I)=ZC(I)
	   YIN(I)=XCRS(I,J)
	   END DO
C
         CALL SPLINE(XIN,YIN,NIN,AK,SARC)                
c        CALL HSPLINE(XIN,YIN,NIN,AKX,AKY,SARC)  
C
	   DO I=1,NILF   
           XGIV=ZF(I) 
C
	 CALL SPLORD(XIN,YIN,NIN,AK,XGIV,YGIV)
C
	   XINT(I,J)=YGIV   
	   END DO
C
	   DO I=1,NILC   
	   XIN(I)=ZC(I) 
	   YIN(I)=YCRS(I,J)
	   END DO
C
       CALL SPLINE(XIN,YIN,NIN,AK,SARC)                
c        CALL HSPLINE(XIN,YIN,NIN,AKX,AKY,SARC)  
C
	   DO I=1,NILF   
           XGIV=ZF(I) 
C
	 CALL SPLORD(XIN,YIN,NIN,AK,XGIV,YGIV)
C
	   YINT(I,J)=YGIV   
	   END DO
C
         END DO
C
C
C  CALCULATION OF X AND Y AT ALL FINE ETA LEVEL 
C  FOR THE FINE LEVEL VALUES OF ALL ZETA
C  THIS IS THE SECOND STEP OF THE BICUBIC INTERPOLATION
C  ALONG ETA ONLY
C
	DO I=1,NILF   
C
	   NIN=NJLC
C
	   DO J=1,NJLC   
	   XIN(J)=EC(J)
	   YIN(J)=XINT(I,J)
	   END DO
C
         CALL SPLINE(XIN,YIN,NIN,AK,SARC)                
c        CALL HSPLINE(XIN,YIN,NIN,AKX,AKY,SARC)  
C
	   DO J=1,NJLF   
	   XGIV=EF(J)
C
	 CALL SPLORD(XIN,YIN,NIN,AK,XGIV,YGIV)
C
	   XFN(I,J)=YGIV
	   END DO
C
	   DO J=1,NJLC   
	   XIN(J)=EC(J)
	   YIN(J)=YINT(I,J)
	   END DO
C
       CALL SPLINE(XIN,YIN,NIN,AK,SARC)                
c        CALL HSPLINE(XIN,YIN,NIN,AKX,AKY,SARC)  
C
	   DO J=1,NJLF   
	   XGIV=EF(J)
C
	 CALL SPLORD(XIN,YIN,NIN,AK,XGIV,YGIV)
C
	   YFN(I,J)=YGIV
	   END DO
	   END DO
C
            RETURN
	    END
C***********************************************************
        SUBROUTINE TRANSNEW(F,IBEG,IFIN,JBEG,JFIN)   
C***********************************************************

         INCLUDE 'comgrid'
C
       DIMENSION F(NX,NY)
   	DIMENSION T1(INMAX),T2(INMAX),S1(INMAX),S2(INMAX),S(INMAX)
C
       IP1=IBEG+1
       IM1=IFIN-1
       JP1=JBEG+1
       JM1=JFIN-1

	I=IBEG
	S(1)=0
        DO J=JP1,JFIN
        DIST = SQRT((X(I,J)-X(I,J-1))**2+(Y(I,J)-Y(I,J-1))**2)
        S(J)=S(J-1)+DIST
        END DO

        DO J=JBEG,JFIN
	     ANUM=S(J)-S(JBEG)
	     ADEN=S(JFIN)-S(JBEG)
	     T1(J)= ANUM/ADEN
         END DO		     

	I=IFIN
	S(1)=0
        DO J=JP1,JFIN
        DIST = SQRT((X(I,J)-X(I,J-1))**2+(Y(I,J)-Y(I,J-1))**2)
        S(J)=S(J-1)+DIST
	END DO

        DO J=JBEG,JFIN
	     ANUM=S(J)-S(JBEG)
	     ADEN=S(JFIN)-S(JBEG)
	     T2(J)= ANUM/ADEN
         END DO		     

	J=JBEG
	S(1)=0
        DO I=IP1,IFIN
        DIST = SQRT((X(I,J)-X(I-1,J))**2+(Y(I,J)-Y(I-1,J))**2)
        S(I)=S(I-1)+DIST
	END DO

        DO I=IBEG,IFIN
           ANUM=S(I)-S(IBEG)
	     ADEN=S(IFIN)-S(IBEG)
	     S1(I)= ANUM/ADEN
         END DO		     

	J=JFIN
	S(1)=0
        DO I=IP1,IFIN
        DIST = SQRT((X(I,J)-X(I-1,J))**2+(Y(I,J)-Y(I-1,J))**2)
        S(I)=S(I-1)+DIST
	END DO
        
	DO I=IBEG,IFIN
	     ANUM=S(I)-S(IBEG)
	     ADEN=S(IFIN)-S(IBEG)
	     S2(I)= ANUM/ADEN
         END DO		     
		

	DO J=JP1,JM1
	 DO I=IP1,IM1
           AN1=(1-T1(J))*S1(I) + T1(J) *S2(I)
	   AD1=1-(S2(I)-S1(I))*(T2(J)-T1(J))
	      U = AN1/AD1

         ALPHA1=1.-U
         ALPHA2= U

           AN2=(1-S1(I))*T1(J) + S1(I) *T2(J)
	   AD2=1-(T2(J)-T1(J))*(S2(I)-S1(I))
	      V = AN2/AD2
	
	BETA1=1.-V
	BETA2=V

	A1V=F(IBEG,J)
	A2V=F(IFIN,J)

	F1UV=ALPHA1*A1V+ALPHA2*A2V

	B1U=F(I,JBEG)
	B2U=F(I,JFIN)

	A1V1=F(IBEG,JBEG)
	A2V1=F(IFIN,JBEG)
	A1V2=F(IBEG,JFIN)
	A2V2=F(IFIN,JFIN)

	FUV=F1UV+BETA1*(B1U-ALPHA1*A1V1-ALPHA2*A2V1)
 	FUV=FUV+BETA2*(B2U-ALPHA1*A1V2-ALPHA2*A2V2)	
	F(I,J)=FUV

	END DO
	END DO

	RETURN
	END
C***********************************************************
	SUBROUTINE GRIDCHK

         INCLUDE 'comgrid'

C                                                                               
C --- VOLUME  OF THE CONTROL VOLUMES
      DO j=1,nj-1
      DO i=1,ni-1
     
      XSW=x(i,j)                                                   
      XSE=x(i+1,j)                                              
      XNW=x(i,j+1)                                              
      XNE=x(i+1,j+1) 
                                        
      YSW=y(i,j)                                                   
      YSE=y(i+1,j)                                              
      YNW=y(i,j+1)                                              
      YNE=y(i+1,j+1)                                         
C                                                                        
C     VOLUME OF THE CONTROL VOLUMES                                                            
        VOLP=0.5*((XNE-XSW)*(YNW-YSE)-(XNW-XSE)*(YNE-YSW))  

      IF (VOLP.LE.0.) then 
	WRITE(30,*)'NEGATIVE VOLUME IN PATCH =',NS
        WRITE(30,*) 'FIRST NEGATIVE VOLUME AT', 'J=',J,'I=',I
 	write(30,*)x(i,j),y(i,j),volp
	WRITE(30,*)'PROGRAM TERMINATE ABNORMALLY DUE TO INTERSECTION'  
	if (maxit .eq. 0) then 
         write(30,*) 'INTERSECTION IN THE ALGEBRAIC GRID'
	else
         write(30,*) 'INTERSECTION IN THE ELLIPTIC GRID'
	end if
	RETURN
      END IF

      END DO
      END DO

	RETURN

      END

C***********************************************************
	SUBROUTINE ORTHCHK

         INCLUDE 'comgrid'
C
C      FIRST DERIVATIVES WRT ZETA AT I=1 AND I=NI
C
	WRITE(52,*)'WEST AND EAST BOUNDARY CHECK'
        DO J=2,NJ-1
          XZETA(1,J) = X(2,J)-X(1,J)
          YZETA(1,J) = Y(2,J)-Y(1,J)
          XZETA(NI,J)= X(NI,J)-X(NI-1,J)
          YZETA(NI,J)= Y(NI,J)-Y(NI-1,J)

	  alfa=(XZETA(1,j)*XZETA(1,J)+YZETA(1,J)*YZETA(1,J))
	  gamma=(XETA(1,j)*XETA(1,J)+YETA(1,J)*YETA(1,J))
	  ORTHW=XZETA(1,J)*XETA(1,J)+YZETA(1,J)*YETA(1,J)
	  ORTHW=ORTHW/SQRT(alfa*gamma)

	  alfa=(XZETA(NI,j)*XZETA(NI,J)+YZETA(NI,J)*YZETA(NI,J))
	  gamma=(XETA(NI,j)*XETA(NI,J)+YETA(NI,J)*YETA(NI,J))
	  ORTHE=XZETA(NI,J)*XETA(NI,J)+YZETA(NI,J)*YETA(NI,J)
	  ORTHE=ORTHE/SQRT(alfa*gamma)

	  WRITE(52,*)J,ORTHW,ORTHE
        END DO
C
C      FIRST DERIVATIVES WRT ETA AT J=1 AND J=NJ
C
	WRITE(52,*)'SOUTH AND NORTH BOUNDARY CHECK'
        DO I=2,NI-1
          XETA(I,1) = X(I,2)-X(I,1)
          YETA(I,1) = Y(I,2)-Y(I,1)
          XETA(I,NJ)= X(I,NJ)-X(I,NJ-1)
          YETA(I,NJ)= Y(I,NJ)-Y(I,NJ-1)

	  alfa=(XZETA(I,1)*XZETA(I,1)+YZETA(I,1)*YZETA(I,1))
	  gamma=(XETA(I,1)*XETA(I,1)+YETA(I,1)*YETA(I,1))
	  ORTHS=XZETA(I,1)*XETA(I,1)+YZETA(I,1)*YETA(I,1)
	  ORTHS=ORTHS/SQRT(alfa*gamma)

	  alfa=(XZETA(I,NJ)*XZETA(I,NJ)+YZETA(I,NJ)*YZETA(I,NJ))
	  gamma=(XETA(I,NJ)*XETA(I,NJ)+YETA(I,NJ)*YETA(I,NJ))
	  ORTHN=XZETA(I,NJ)*XETA(I,NJ)+YZETA(I,NJ)*YETA(I,NJ)
	  ORTHN=ORTHN/SQRT(alfa*gamma)

	  WRITE(52,*)I,ORTHS,ORTHN
        END DO 

	RETURN

        END
C
C***********************************************************************

	SUBROUTINE CONTROL_OLD
C
C***********************************************************************
C
C       THIS SUBROUTINE COMPUTES THE CONTROL FUNCTIONS P AND Q
C       ACCORDING TO THE IDEA PROPOSED BY THOMAS AND MIDDELECOFF.
C
C	IMPLICIT REAL*8 (A-H,O-Z)
C
         INCLUDE 'comgrid'
C
        RESM =DMAX1(RESDU(1),RESDU(2))
	DO I=1,NI,NI-1
	  DO J=2,NJ-1
            QOLD=Q(I,J)
	  IF(I.EQ.1)THEN
	    SINTH=SIN(THZI(J))
	    COSTH=COS(THZI(J))
	  ENDIF
	  IF(I.EQ.NI)THEN
	    SINTH=SIN(THZO(J))
	    COSTH=COS(THZO(J))
          ENDIF
	    XE = XETA(I,J)
	    YE = YETA(I,J)
	    XZ = XZETA(I,J)
	    YZ = YZETA(I,J)
	    XEE = XEETA(I,J)
	    YEE = YEETA(I,J)
	    XZZ =XZZETA(I,J)
	    YZZ =YZZETA(I,J)
	    XZET=XZE(I,J)
	    YZET=YZE(I,J)
	    ALPHA=XE*XE+YE*YE
	    GAMMA=XZ*XZ+YZ*YZ
	    BETA=XZ*XE+YZ*YE
	    QTHOMC = -((XZZ*XE+YZZ*YE)/GAMMA)
	    QTHOM = -((XEE*XE+YEE*YE)/ALPHA)
	    QBETA =2*BETA*(XZET*XE+YZET*YE)/(ALPHA*GAMMA)
            QN1 = ((XZZ*YE-YZZ*XE)/GAMMA)+((XEE*YE-YEE*XE)/ALPHA)
	    QN2 = -2*BETA*(XZET*YE-YZET*XE)/(ALPHA*GAMMA)
	    QTHETA= (QN1+QN2)*COSTH/SINTH
	    QNEW=QTHOM
            IF(NITER .GT.100 .OR. RESM .LT. 5.D-04)THEN
  	    QNEW = QTHOM+QTHOMC+QBETA+QTHETA 
	    ENDIF
	    Q(I,J) = QOLD*(1-OMEGA)+OMEGA*QNEW
	  ENDDO
	ENDDO
C
	DO J=1,NJ,NJ-1
	  DO I=2,NI-1
            POLD=P(I,J)
	    IF(J.EQ.1)THEN
	      SINTH=SIN(THEI(I))
	      COSTH=COS(THEI(I))
            ENDIF
	    IF(J.EQ.NJ)THEN
	      SINTH=SIN(THEO(I))
	      COSTH=COS(THEO(I))
            ENDIF
	    XE = XETA(I,J)
	    YE = YETA(I,J)
	    XZ = XZETA(I,J)
	    YZ = YZETA(I,J)
	    XEE = XEETA(I,J)
	    YEE = YEETA(I,J)
	    XZZ =XZZETA(I,J)
	    YZZ =YZZETA(I,J)
	    XZET=XZE(I,J)
	    YZET=YZE(I,J)
	    ALPHA=XE*XE+YE*YE
	    GAMMA=XZ*XZ+YZ*YZ
	    BETA=XZ*XE+YZ*YE
	    PTHOM= -((XZZ*XZ+YZZ*YZ)/GAMMA)
	    PTHOMC=-((XEE*XZ+YEE*YZ)/ALPHA)
	    PBETA=2*BETA*(XZET*XZ+YZET*YZ)/(ALPHA*GAMMA)
	    PN1 =-((XZZ*YZ-YZZ*XZ)/GAMMA)-((XEE*YZ-YEE*XZ)/ALPHA)
            PN2 = 2*BETA*(XZET*YZ-YZET*XZ)/(ALPHA*GAMMA)
	    PTHETA =(PN1+PN2)*COSTH/SINTH
	    PNEW=PTHOM
            IF(NITER .GT.100 .OR. RESM .LT. 5.D-04)THEN
  	    PNEW = PTHOM+PTHOMC+PBETA+PTHETA 
	    ENDIF
	    P(I,J) = POLD*(1-OMEGA)+OMEGA*PNEW
	  ENDDO
	ENDDO
C
        CALL LINTPJ(P,1,NI,1,NJ)
        CALL LINTPI(Q,1,NI,1,NJ)

C
C      STRETCHING THE P TERM EXPONENTIALLY
C
C            DO I=2,NI-1
C             DO J=2,NJ-1
C
C	    IF(LEVEL.EQ.1) THEN

C	    IF(J.LT.NPTC(2,4)) ISEG = 1
C	    IF(J.GT.NPTC(2,4)) ISEG = 2
C           ELSE 
C
C	    IF(J.LT.NPTF(2,4)) ISEG = 1
C	    IF(J.GT.NPTF(2,4)) ISEG = 2
C
C	    IF(LEVEL .EQ. 1) THEN
C            NEND=NPTC(ISEG,4)
C            NEND=NPTF(ISEG,4)
C	    ENDIF
C
C            XBEG = XDAT(1,ISEG,4)
C            YBEG = YDAT(1,ISEG,4)
C            XEND = XDAT(2,ISEG,4)
C            YEND = YDAT(2,ISEG,4)
C            DELX = XEND - XBEG
C            DELY = YEND - YBEG
C            ALTOT = SQRT ( DELX*DELX + DELY* DELY)
C
C		 IF (ISEG .EQ. 1) THEN
C		       PBEG = P(I,1)
C	     	       AAFR = P(I,2) - P(I,1)
C                       CALL EXPSTR (NEND,ALTOT,RSTR,P,PBEG,AAFR,IND)
C
C		       AAFR = 2.0*P(I,NEND) - 3.0*P(I,NEND-1)
C     &		           + P(I,NEND-2)
C
C		 ELSE
C
C                       CALL EXPSTR (NEND,ALTOT,RSTR,P,AAFR,AAFR,IND)
C
C
C		 ENDIF
C
C	     ENDDO
C	  ENDDO
C
C	CALL EXPPJ(P,1,NI,1,NJ)
C	CALL EXPPI(Q,1,NI,1,NJ)
C
	RETURN
	END
