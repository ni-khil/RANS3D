      SUBROUTINE OUT_NEW(IFLAG) 
C
C     CALCULATE FLUXES OF MASS MOMENTUM
C     AND SCALARS AT INFLOW BOUNDARIES                   
C
C***********************************************************************
C
      INCLUDE 'com2d'                                              
C
	DO I=1,11
	  RNORM(I)=0.
	END DO

      IF(MINLW.NE.0)  THEN                                         
        DO  IM=1,MINLW                                               
              IP=INLEW(IM)                                              
              IB=IP-1                                                   
              CW=(F(ISCW+IP))*mysolv(ip)
        RNORM(1)=RNORM(1)+CW                                            
        RNORM(2)=RNORM(2)+CW*SQRT(F(ISU+IB)**2+F(ISV+IB)**2)            
        RNORM(4)=RNORM(4)+CW*ABS(F(ISWR+IB))                            
        RNORM(5)=RNORM(5)+CW*ABS(F(ISTE+IB))                            
        RNORM(6)=RNORM(6)+CW*ABS(F(ISED+IB))                            
        RNORM(7)=RNORM(7)+CW*ABS(F(ISOM +IB))    
        RNORM(8)=RNORM(8)+CW*ABS(F(ISV2 +IB))    
        RNORM(9)=RNORM(9)+CW*ABS(F(ISF +IB))    
        RNORM(10)=RNORM(10)+CW*ABS(F(ISNUT +IB))    
        RNORM(11)=RNORM(11)+CW*ABS(F(ISS +IB))    
          END DO                        
        RNORM(3)=RNORM(2)
	END IF
C                                               
        IF(MINLE.NE.0) THEN                                          
        DO  IM=1,MINLE                                               
              IP=INLEE(IM)                                            
              IB=IP+1                                            
              CE=-(F(ISCW+IB))*mysolv(ip)
        RNORM(1)=RNORM(1)+CE                                            
        RNORM(2)=RNORM(2)+CE*SQRT(F(ISU+IB)**2+F(ISV+IB)**2)            
        RNORM(4)=RNORM(4)+CE*(F(ISWR+IB))                            
        RNORM(5)=RNORM(5)+CE*(F(ISTE+IB))                            
        RNORM(6)=RNORM(6)+CE*(F(ISED+IB))                            
        RNORM(7)=RNORM(7)+CE*ABS(F(ISOM +IB))    
        RNORM(8)=RNORM(8)+CE*ABS(F(ISV2 +IB))    
        RNORM(9)=RNORM(9)+CE*ABS(F(ISF +IB))    
        RNORM(10)=RNORM(10)+CE*ABS(F(ISNUT +IB))    
        RNORM(11)=RNORM(11)+CE*ABS(F(ISS +IB))    
        END DO
        RNORM(3)=RNORM(2) 
	END IF
C                                              
        IF(MINLS.NE.0) THEN                                          
        DO  IM=1,MINLS                                               
              IP=INLES(IM)                                              
              IB=IP-NI                                                  
              CS=(F(ISCS+IP))*mysolv(ip)
        RNORM(1)=RNORM(1)+CS                                            
        RNORM(2)=RNORM(2)+CS*SQRT(F(ISU+IB)**2+F(ISV+IB)**2)            
        RNORM(4)=RNORM(4)+CS*ABS(F(ISWR+IB))                            
        RNORM(5)=RNORM(5)+CS*ABS(F(ISTE+IB))                            
        RNORM(6)=RNORM(6)+CS*ABS(F(ISED+IB))                            
        RNORM(7)=RNORM(7)+CS*ABS(F(ISOM +IB))    
        RNORM(8)=RNORM(8)+CS*ABS(F(ISV2 +IB))    
        RNORM(9)=RNORM(9)+CS*ABS(F(ISF +IB))    
        RNORM(10)=RNORM(10)+CS*ABS(F(ISNUT +IB))    
        RNORM(11)=RNORM(11)+CS*ABS(F(ISS +IB))    
         END DO
        RNORM(3)=RNORM(2)
	END IF
C                                               
	
        IF(MINLN.NE.0) THEN                                            
        DO  IM=1,MINLN                                               
              IP=INLEN(IM)                                           
              IB=IP+NI                                           
              CN=-(F(ISCS+IB))*mysolv(ip)

              IF (CN .LT. 0) THEN
C	          CALL OUT(IM,IFLAG) ! OUT FLOW CV	
	          CALL OUT(IP,IFLAG) ! OUT FLOW CV	
C		  WRITE(1200+iblk,*)IM,ip,IB,CN
	      ELSE
C		  WRITE(1300+iblk,*)IM,ip,IB,CN
              RNORM(1)=RNORM(1)+CN                                            
              RNORM(2)=RNORM(2)+CN*SQRT(F(ISU+IB)**2+F(ISV+IB)**2)            
              RNORM(3)=RNORM(2)
              RNORM(4)=RNORM(4)+CN*ABS(F(ISWR+IB))                            
              RNORM(5)=RNORM(5)+CN*ABS(F(ISTE+IB))                            
              RNORM(6)=RNORM(6)+CN*ABS(F(ISED+IB))                            
              RNORM(7)=RNORM(7)+CN*ABS(F(ISOM +IB))    
              RNORM(8)=RNORM(8)+CN*ABS(F(ISV2 +IB))    
              RNORM(9)=RNORM(9)+CN*ABS(F(ISF +IB))    
              RNORM(10)=RNORM(10)+CN*ABS(F(ISNUT +IB))    
              RNORM(11)=RNORM(11)+CN*ABS(F(ISS +IB))    
	END IF
        END DO
	END IF

	RETURN
      END                                                               
C
C**********************************************************************
C                            
      SUBROUTINE OUT_NEWPREV(IFLAG)       
C                                                                        
C EXTRAPOLATE FLOW VARIABLES AT OUTFLOW BOUNDARY
C                              
C***********************************************************************
C                                           
      INCLUDE 'com2d'                                               
C                                                                       
C   IFLAG = 1 FOR OUTFLOW AT NORTH                     
C                              
C   IFLAG = 2 FOR OUTFLOW AT SOUTH 
C
C   IFLAG = 3 FOR OUTFLOW AT EAST
C
C   IFLAG = 4 FOR OUTFLOW AT WEST
C                                           
        MFLAG1=0                                                        
        MFLAG2=0                                                        
        MFLAG3=0                                                        
        MFLAG4=0                                                        
      IF(IFLAG.EQ.1) THEN                                               
        MFLAG1=1                                                        
      ELSE IF(IFLAG.EQ.2) THEN                                          
        MFLAG2=1                                                        
      ELSE IF(IFLAG.EQ.3) THEN                                          
        MFLAG3=1                                                        
      ELSE IF(IFLAG.EQ.4) THEN                                          
        MFLAG4=1                                                        
      ENDIF    
                                                                        
      MEXI   =MEXIN*MFLAG1 + MEXIS*MFLAG2 + MEXIE*MFLAG3 + MEXIW*MFLAG4 
      ISANOW = ISAN*MFLAG1 +  ISAS*MFLAG2 +  ISAE*MFLAG3 +  ISAW*MFLAG4 
                                                                        
      IVAR=ISPP+(NN-1)*NIJ
C	write(1200+iblk,*)'SWEEP NO =',IBSWP
C		  WRITE(1200+iblk,*)'IM,ip'
      DO IM=1,MEXI                                                 
         IP  = INEXN(IM)*MFLAG1 + INEXS(IM)*MFLAG2                      
     *       + INEXE(IM)*MFLAG3 + INEXW(IM)*MFLAG4                      

C		  WRITE(1200+iblk,*)IM,ip
	CALL OUT(IP,IFLAG)
      END DO                                                  
      RETURN                                                            
      END                                                               
C
C**********************************************************************
C                            
      SUBROUTINE OUT_OLD(IFLAG)       
C                                                                        
C EXTRAPOLATE FLOW VARIABLES AT OUTFLOW BOUNDARY
C                              
C***********************************************************************
C                                           
      INCLUDE 'com2d'                                               
C                                                                       
C   IFLAG = 1 FOR OUTFLOW AT NORTH                     
C                              
C   IFLAG = 2 FOR OUTFLOW AT SOUTH 
C
C   IFLAG = 3 FOR OUTFLOW AT EAST
C
C   IFLAG = 4 FOR OUTFLOW AT WEST
C                                           
        MFLAG1=0                                                        
        MFLAG2=0                                                        
        MFLAG3=0                                                        
        MFLAG4=0                                                        
      IF(IFLAG.EQ.1) THEN                                               
        MFLAG1=1                                                        
      ELSE IF(IFLAG.EQ.2) THEN                                          
        MFLAG2=1                                                        
      ELSE IF(IFLAG.EQ.3) THEN                                          
        MFLAG3=1                                                        
      ELSE IF(IFLAG.EQ.4) THEN                                          
        MFLAG4=1                                                        
      ENDIF    
                                                                        
      MEXI   =MEXIN*MFLAG1 + MEXIS*MFLAG2 + MEXIE*MFLAG3 + MEXIW*MFLAG4 
      ISANOW = ISAN*MFLAG1 +  ISAS*MFLAG2 +  ISAE*MFLAG3 +  ISAW*MFLAG4 
                                                                        
      IVAR=ISPP+(NN-1)*NIJ
      DO IM=1,MEXI                                                 

         II  = INEXN(IM)*MFLAG1 + INEXS(IM)*MFLAG2                      
     *       + INEXE(IM)*MFLAG3 + INEXW(IM)*MFLAG4                      

         IND_BNDRY =   NI *MFLAG1   -NI *MFLAG2+       MFLAG3    -   MFLAG4  
        
	  INB=IVAR+II
          IB=INB+IND_BNDRY
         F(IB)=F(INB) 

        if (nn .eq. 2) F(IB)=F(INB)
        if (nn .eq. 3) F(IB)=F(INB)

      F(ISANOW+II)=0.  
      END DO                                                  
      RETURN                                                            
      END                                                               
C                              
C***********************************************************************
C                                           
      SUBROUTINE FARFIELD(IFLAG)       
C                                                                        
C EXTRAPOLATE FLOW VARIABLES AT OUTFLOW BOUNDARY
C                              
C***********************************************************************
C                                           
      INCLUDE 'com2d'                                               
C                                                                       
C   IFLAG = 1 FOR OUTFLOW AT NORTH                     
C                              
C   IFLAG = 2 FOR OUTFLOW AT SOUTH 
C
C   IFLAG = 3 FOR OUTFLOW AT EAST
C
C   IFLAG = 4 FOR OUTFLOW AT WEST
C                                           
        MFLAG1=0                                                        
        MFLAG2=0                                                        
        MFLAG3=0                                                        
        MFLAG4=0                                                        
      IF(IFLAG.EQ.1) THEN                                               
        MFLAG1=1                                                        
      ELSE IF(IFLAG.EQ.2) THEN                                          
        MFLAG2=1                                                        
      ELSE IF(IFLAG.EQ.3) THEN                                          
        MFLAG3=1                                                        
      ELSE IF(IFLAG.EQ.4) THEN                                          
        MFLAG4=1                                                        
      ENDIF    

        DO I=1,11
          RNORM(I)=0.
        END DO
                                                                        
      MINL   =MINLN*MFLAG1 + MINLS*MFLAG2 + MINLE*MFLAG3 + MINLW*MFLAG4 
      ISANOW = ISAN*MFLAG1 +  ISAS*MFLAG2 +  ISAE*MFLAG3 +  ISAW*MFLAG4 
                                                                        
      DO IM=1,MINL                                                 
         IP  = INLEN(IM)*MFLAG1 + INLES(IM)*MFLAG2                      
     *       + INLEE(IM)*MFLAG3 + INLEW(IM)*MFLAG4                      

         IND_BNDRY =  NI *MFLAG1 - MFLAG2+  MFLAG3    - MFLAG4
         ICONV_IND_BNDRY =  NI *MFLAG1 -0 *MFLAG2+  MFLAG3    - 0*MFLAG4

	 IB=IP+IND_BNDRY
	 IC_IB=IP+ICONV_IND_BNDRY

         CONV_B=F(ISCS+IC_IB)*(-MFLAG1+MFLAG2)+F(ISCW+IC_IB)*(-MFLAG3+MFLAG4)
	 CONV_B=CONV_B*mysolv(ip)

	 if (CONV_B .LT.0) THEN
            CALL OUT(IP,IFLAG) ! OUT FLOW CV
         ELSE
	    CONV_INF=CONV_B
	    RNORM(1)=RNORM(1)+CONV_INF
            RNORM(2)=RNORM(2)+CONV_INF*SQRT(F(ISU+IB)**2+F(ISV+IB)**2)
            RNORM(3)=RNORM(2)
            RNORM(4)=RNORM(4)+CONV_INF*ABS(F(ISWR+IB))
            RNORM(5)=RNORM(5)+CONV_INF*ABS(F(ISTE+IB))
            RNORM(6)=RNORM(6)+CONV_INF*ABS(F(ISED+IB))
            RNORM(7)=RNORM(7)+CONV_INF*ABS(F(ISOM +IB))
            RNORM(8)=RNORM(8)+CONV_INF*ABS(F(ISV2 +IB))
            RNORM(9)=RNORM(9)+CONV_INF*ABS(F(ISF +IB))
            RNORM(10)=RNORM(10)+CONV_INF*ABS(F(ISNUT +IB))
            RNORM(11)=RNORM(11)+CONV_INF*ABS(F(ISS +IB))

	 end if

      END DO                                                  
      RETURN                                                            
      END                                                               
