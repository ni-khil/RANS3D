C
C*****************************************************************
C
C     COMPUTER PROGRAM FOR PREDICTION OF INCOMPRESSIBLE
C     TWO DIMENSIONAL TURBULENT RECIRCULATIING FLOWS
C
C     ENQUIRIES TO BE DIRECTED TO
C     DR S.MAJUMDAR, CTFD DIVISION, NAL BANGALORE
C
C     PROGRAM STATUS AS ON : FEBRUARY 22, 2005
C  
C     THIS IS THE USER INDEPENDENT PART OF THE PROGRAM    
C                                                                        
C***********************************************************************
C 
      SUBROUTINE BINDEX                                
C                 
C     CALCULATE LOCATION ADDRESS FOR THE NEAR-BOUNDARY NODES
C     AND FILL IN ALL ARRAYS WITH DEFAULT VALUES 
C            
C***********************************************************************
C  
      INCLUDE 'com2d'                                                  
C                                                                        
      DO  I=1,NIJ                                                      
        KW(I)=0                                                         
        KE(I)=0                                                         
        KS(I)=0                                                         
        KN(I)=0                                                         
	mysolv(i)=0
        KBLK(I)=0
        DISTN(I)=1.E20
	APP(I)=1.
      end do	
C                                                         
      IBLOC=0
C                                                           
      DO J=jsolv_b,jsolv_e
      DO I=isolv_b,isolv_e
      II=I+(J-1)*NI
	mysolv(ii)=1
        do indflag=1,4
          tauwal(II,indflag)=0.
          yplus(II,indflag)=0.
        end do

      end do		
      end do		

       DO J=2,NJM                                                      
       DO I=2,NIM                                                      
       KBLK(I+(J-1)*NI)=1 
       end do		
       end do		
C                                               
       DO I=1,100                                                     
       INFDE(I)=0                                                        
       end do		
C
       DO I=1,IDIM
       F(I)=0.    
       END DO
C
       DO I=1,11
       RMOM(I)=0.
       RNORM(I)=0.
       END DO
c
C
       DO I=1,NIJ
       F(ISDEN+I)=DENSIT
       F(ISVIS+I)=VISCOS
       F(ISVISS+I)=VISCOS
       F(ISVISW+I)=VISCOS
       END DO
C                                                                        
         IF(NRBNDA.LT.4.OR.NRBNDA.GT.NRBND) CALL CHECK1(5) 
C     
C     DEFINING KIND OF BOUNDARIES        
C     
      DO N=1,NRBNDA
                                                 
         IF(CBND(N,2)(1:2).EQ.'WA') THEN                            
           N2=11                                                        
           INFDE(6)=INFDE(6)+1                                          
         ELSEIF(CBND(N,2)(1:2).EQ.'SY') THEN                            
           N2=12                                                        
           INFDE(7)=INFDE(7)+1                                          
         ELSEIF(CBND(N,2)(1:2).EQ.'OU') THEN                            
           N2=13                                                        
           INFDE(8)=INFDE(8)+1                                          
         ELSEIF(CBND(N,2)(1:2).EQ.'IN') THEN                            
           N2=14                                                        
           INFDE(9)=INFDE(9)+1                                          
         ELSEIF(CBND(N,2)(1:2).EQ.'CU') THEN                            
           N2=15                                                        
           INFDE(10)=INFDE(10)+1
         ELSEIF(CBND(N,2)(1:2).EQ.'CY') THEN
           N2=16
           INFDE(12)=INFDE(12)+1                                        
         ELSEIF(CBND(N,2)(1:2).EQ.'PC') THEN
           N2=17
           INFDE(13)=INFDE(13)+1                                        
         ELSEIF(CBND(N,2)(1:2).EQ.'BL') THEN                            
           N2=11                                                        
           INFDE(11)=INFDE(11)+1        
         ELSE                                                           
         CALL CHECK1(3)
         ENDIF         
C     
C     DEFINING LOCATION OF BOUNDARIES        
C     
         IF(CBND(N,1)(1:2).EQ.'WE') THEN                            
           I=NBND(N,1)                                                  
           DO  J=NBND(N,3),NBND(N,4)                                 
           KW(I+(J-1)*NI)=N2                                            
	  END DO
           INFDE(1)=INFDE(1)+1                                          
         ELSEIF(CBND(N,1)(1:2).EQ.'EA') THEN                            
           I=NBND(N,2)                                                  
           DO  J=NBND(N,3),NBND(N,4)                                 
           KE(I+(J-1)*NI)=N2                                            
   	   END DO
           INFDE(2)=INFDE(2)+1                                          
         ELSEIF(CBND(N,1)(1:2).EQ.'SO') THEN                            
           J=NBND(N,3)                                                  
           DO  I=NBND(N,1),NBND(N,2)                                 
           KS(I+(J-1)*NI)=N2                                            
	   END DO
           INFDE(3)=INFDE(3)+1                                          
         ELSEIF(CBND(N,1)(1:2).EQ.'NO') THEN                            
           J=NBND(N,4)                                                  
           DO  I=NBND(N,1),NBND(N,2)                                 
           KN(I+(J-1)*NI)=N2                                            
	   END DO
           INFDE(4)=INFDE(4)+1                                          
         ELSEIF(CBND(N,1)(1:2).EQ.'BL') THEN                   
           IBLOC=1                                                      
           DO I=NBND(N,1),NBND(N,2)                                 
           DO J=NBND(N,3),NBND(N,4)                                 
           INFDE(11)=INFDE(11)+1                                        
           KBLK(I+(J-1)*NI)=0        
	   END DO
	   END DO
           INFDE(5)=INFDE(5)+1                                          
         ELSE                                                           
           CALL CHECK1(4)                                               
         ENDIF                                                          
	END DO
C                                                                        
      MWALN=0                                                           
      MWALS=0                                                           
      MWALE=0                                                           
      MWALW=0                                                           
      MSYMN=0                                                           
      MSYMS=0                                                           
      MSYME=0                                                           
      MSYMW=0                                                           
      MEXIN=0                                                           
      MEXIS=0                                                           
      MEXIE=0                                                           
      MEXIW=0                                                           
      MINLN=0                                                           
      MINLS=0                                                           
      MINLE=0                                                           
      MINLW=0                                                           
      MCUTN=0
      MCUTS=0
      MCUTE=0
      MCUTW=0
      MCYCN=0   
      MCYCS=0                                         
      MCYCE=0                               
      MCYCW=0  
      MPCTN=0   
      MPCTS=0                                         
      MPCTE=0                               
      MPCTW=0  
                                                                        
      DO J=1,NJ
      DO I=1,NI	
	II=I+(J-1)*NI
      ISUM=KN(II)+KS(II)+KE(II)+KW(II)                                  
      IF(ISUM.NE.0) THEN 
C
C  WALL BOUNDARY
C  
      IF(KN(II).EQ.11) THEN                        
        MWALN=MWALN+1                                                   
        INDWN(MWALN)=II                                                 
      END IF
      IF(KS(II).EQ.11) THEN                        
        MWALS=MWALS+1                                                   
        INDWS(MWALS)=II                                                 
      END IF
      IF(KE(II).EQ.11) THEN         
        MWALE=MWALE+1                                                   
        INDWE(MWALE)=II                                                 
      END IF
      IF(KW(II).EQ.11) THEN         
        MWALW=MWALW+1                                                   
        INDWW(MWALW)=II                                                 
      END IF
C
C  SYMMETRY BOUNDARY
C  
                                                              
       IF(KN(II).EQ.12) THEN                       
        MSYMN=MSYMN+1                                                   
        INSYN(MSYMN)=II                                                 
      END IF
      IF(KS(II).EQ.12) THEN                       
        MSYMS=MSYMS+1                                                   
        INSYS(MSYMS)=II                                                 
      END IF
      IF(KE(II).EQ.12) THEN                        
        MSYME=MSYME+1                                                   
        INSYE(MSYME)=II                                                 
      END IF
      IF(KW(II).EQ.12) THEN                    
        MSYMW=MSYMW+1                                                   
        INSYW(MSYMW)=II                                                 
      END IF
C
C  EXIT BOUNDARY
C  
                                                                        
      IF(KN(II).EQ.13) THEN                       
        MEXIN=MEXIN+1                                                   
        INEXN(MEXIN)=II                                                 
      END IF
      IF(KS(II).EQ.13) THEN                        
        MEXIS=MEXIS+1                                                   
        INEXS(MEXIS)=II                                                 
      END IF
      IF(KE(II).EQ.13) THEN                       
        MEXIE=MEXIE+1                                                   
        INEXE(MEXIE)=II                                                 
      END IF
      IF(KW(II).EQ.13) THEN                       
        MEXIW=MEXIW+1                                                   
        INEXW(MEXIW)=II                                                 
      END IF
C
C  INFLOW BOUNDARY
C  
      IF(KN(II).EQ.14) THEN                       
        MINLN=MINLN+1                                                   
        INLEN(MINLN)=II                                                 
      END IF
      IF(KS(II).EQ.14) THEN                       
        MINLS=MINLS+1                                                   
        INLES(MINLS)=II                                                 
      END IF
      IF(KE(II).EQ.14) THEN                       
        MINLE=MINLE+1                                                   
        INLEE(MINLE)=II                                                 
      END IF
      IF(KW(II).EQ.14) THEN                     
        MINLW=MINLW+1                                                   
        INLEW(MINLW)=II                                                 
      END IF

C
C  CUT BOUNDARY
C  

      IF(KN(II).EQ.15) THEN                       
        MCUTN=MCUTN+1                                                   
        ICUTN(MCUTN)=II                                                 
      END IF
      IF(KS(II).EQ.15) THEN                       
        MCUTS=MCUTS+1                                                   
        ICUTS(MCUTS)=II                                                 
      END IF
      IF(KE(II).EQ.15) THEN                       
        MCUTE=MCUTE+1                                                   
        ICUTE(MCUTE)=II
      END IF
      IF(KW(II).EQ.15) THEN                       
        MCUTW=MCUTW+1                                                   
        ICUTW(MCUTW)=II
      END IF

C
C  CYCLIC BOUNDARY
C  
          IF(KN(II).EQ.16) THEN
            MCYCN=MCYCN+1
            INCYN(MCYCN)=II
        END IF
          IF(KS(II).EQ.16) THEN
            MCYCS=MCYCS+1
            INCYS(MCYCS)=II
        END IF
          IF(KE(II).EQ.16) THEN
            MCYCE=MCYCE+1
            INCYE(MCYCE)=II
        END IF
          IF(KW(II).EQ.16) THEN
            MCYCW=MCYCW+1
            INCYW(MCYCW)=II
        END IF
C
C  PRESSURE CONSTANT BOUNDARY
C  
          IF(KN(II).EQ.17) THEN
            MPCTN=MPCTN+1
            INPCN(MPCTN)=II
        END IF
          IF(KS(II).EQ.17) THEN
            MPCTS=MPCTS+1
            INPCS(MPCTS)=II
        END IF
          IF(KE(II).EQ.17) THEN
            MPCTE=MPCTE+1
            INPCE(MPCTE)=II
        END IF
          IF(KW(II).EQ.17) THEN
            MPCTW=MPCTW+1
            INPCW(MPCTW)=II
        END IF
        
	END IF
        
	END DO
	END DO
 
      IF(MWALN.GT.IDIM1) CALL CHECK1(2)                                 
      IF(MWALS.GT.IDIM1) CALL CHECK1(2)                                 
      IF(MWALE.GT.IDIM1) CALL CHECK1(2)                                 
      IF(MWALW.GT.IDIM1) CALL CHECK1(2)                                 

      IF(MSYMN.GT.NXF) CALL CHECK1(2)                                 
      IF(MSYMS.GT.NXF) CALL CHECK1(2)                                 
      IF(MSYME.GT.NYF) CALL CHECK1(2)                                 
      IF(MSYMW.GT.NYF) CALL CHECK1(2)                                 

      IF(MEXIN.GT.NXF) CALL CHECK1(2)                                 
      IF(MEXIS.GT.NXF) CALL CHECK1(2)                                 
      IF(MEXIE.GT.NYF) CALL CHECK1(2)                                 
      IF(MEXIW.GT.NYF) CALL CHECK1(2)                                 

      IF(MINLN.GT.NXF) CALL CHECK1(2)                                 
      IF(MINLS.GT.NXF) CALL CHECK1(2)                                 
      IF(MINLE.GT.NYF) CALL CHECK1(2)                                 
      IF(MINLW.GT.NYF) CALL CHECK1(2)                                 

      IF(MCUTN.GT.NXF) CALL CHECK1(2)                                 
      IF(MCUTS.GT.NXF) CALL CHECK1(2)                                 
      IF(MCUTE.GT.NYF) CALL CHECK1(2)                                 
      IF(MCUTW.GT.NYF) CALL CHECK1(2)

      IF(MCYCN.GT.NXF) CALL CHECK1(2)
      IF(MCYCS.GT.NXF) CALL CHECK1(2)                                   
      IF(MCYCE.GT.NYF) CALL CHECK1(2)            
      IF(MCYCW.GT.NYF) CALL CHECK1(2)

      IF(MPCTN.GT.NXF) CALL CHECK1(2)
      IF(MPCTS.GT.NXF) CALL CHECK1(2)                                   
      IF(MPCTE.GT.NYF) CALL CHECK1(2)            
      IF(MPCTW.GT.NYF) CALL CHECK1(2)
                                                                        
      RETURN                                                            
      END                                                               
C                                                                        
C***********************************************************************
C 
             SUBROUTINE BFLUX     
C                                            
C   COMPUTES MASS FLUXES AT BOUNDARIES 
C                                    
C***********************************************************************
C                                   
      INCLUDE 'com2d'
C                                                  
      SQR(VALUE)=VALUE*VALUE                                            
      VALUE=0.
C                                                                       
        IF(MWALW .NE. 0) THEN
           DO IM=1,MWALW                                               
              F(ISCW+INDWW(IM))=0. 
	   END DO
	END IF

        IF(MWALE .NE. 0) THEN
           DO IM=1,MWALE                                               
              F(ISCW+INDWE(IM)+1)=0.                                          
	   END DO
	END IF

       IF(MWALS .NE. 0) THEN

          DO IM=1,MWALS                                               
             IP=INDWS(IM)
             IB=IP-NI
             RS=F(ISFY+IP)*F(ISR+IP)+(1.-F(ISFY+IP))*F(ISR+IB)
             F(ISCS+IP)=F(ISDEN+IB)*RS*(F(ISU+IB)*F(ISB21S+IP)
     >                  +F(ISV+IB)*F(ISB22S+IP))*mysolv(ip)
	  END DO
	END IF

        IF(MWALN .NE. 0) THEN
           DO IM=1,MWALN                                               
              IP=INEXN(IM)
              IB=IP+NI
              RN=F(ISFY+IB)*F(ISR+IB)+(1.-F(ISFY+IB))*F(ISR+IP)
              F(ISCS+IB)=F(ISDEN+IB)*RN*(F(ISU+IB)*F(ISB21S+IB)
     >                   +F(ISV+IB)*F(ISB22S+IB))
	  END DO
	END IF

        IF(MSYMW.NE.0) THEN
           DO IM=1,MSYMW                                               
              F(ISCW+INSYW(IM))=0.
	  END DO
	END IF

        IF(MSYME.NE.0) THEN
           DO IM=1,MSYME                                               
              F(ISCW+INSYE(IM)+1)=0.                                          
	  END DO
	END IF

        IF(MSYMS.NE.0) THEN                                          
           DO IM=1,MSYMS                                               
              F(ISCS+INSYS(IM))=0. 
	  END DO
	END IF

        IF(MSYMN.NE.0) THEN
           DO IM=1,MSYMN                                               
              F(ISCS+INSYN(IM)+NI)=0. 
	  END DO
	END IF
                                                                        
        IF(MEXIW.NE.0) THEN
           DO IM=1,MEXIW                                               
              IP=INEXW(IM)                                             
              IB=IP-1                                                  
              RW=F(ISFX+IP)*F(ISR+IP)+(1.-F(ISFX+IP))*F(ISR+IB)        
              F(ISCW+IP)=F(ISDEN+IB)*RW*(F(ISU+IB)*F(ISB11W+IP)               
     >                   +F(ISV+IB)*F(ISB12W+IP))              
	  END DO
	END IF

        IF(MEXIE.NE.0) THEN
           DO IM=1,MEXIE                                               
              IP=INEXE(IM)                                             
              IB=IP+1                                                  
              RE=F(ISFX+IB)*F(ISR+IB)+(1.-F(ISFX+IB))*F(ISR+IP)        
              F(ISCW+IB)=F(ISDEN+IB)*RE*(F(ISU+IB)*F(ISB11W+IB)               
     >                   +F(ISV+IB)*F(ISB12W+IB))              
	  END DO
	END IF

        IF(MEXIS.NE.0) THEN
           DO IM=1,MEXIS                                               
              IP=INEXS(IM)                                             
              IB=IP-NI                                                 
              RS=F(ISFY+IP)*F(ISR+IP)+(1.-F(ISFY+IP))*F(ISR+IB)        
              F(ISCS+IP)=F(ISDEN+IB)*RS*(F(ISU+IB)*F(ISB21S+IP)               
     >                   +F(ISV+IB)*F(ISB22S+IP))              
	  END DO
	END IF

        IF(MEXIN.NE.0) THEN
           DO IM=1,MEXIN                                               
              IP=INEXN(IM)                                             
              IB=IP+NI                                                 
              RN=F(ISFY+IB)*F(ISR+IB)+(1.-F(ISFY+IB))*F(ISR+IP)        
              F(ISCS+IB)=F(ISDEN+IB)*RN*(F(ISU+IB)*F(ISB21S+IB)               
     >                   +F(ISV+IB)*F(ISB22S+IB))              
	  END DO
	END IF

        IF(MINLW.NE.0) THEN
           DO IM=1,MINLW                                               
              IP=INLEW(IM)                                             
              IB=IP-1                                                  
              RW=F(ISFX+IP)*F(ISR+IP)+(1.-F(ISFX+IP))*F(ISR+IB)        
              F(ISCW+IP)=F(ISDEN+IB)*RW*(F(ISU+IB)*F(ISB11W+IP)               
     >                   +F(ISV+IB)*F(ISB12W+IP))              
	  END DO
	END IF

        IF(MINLE.NE.0) THEN
           DO IM=1,MINLE                                               
              IP=INLEE(IM)                                             
              IB=IP+1                                                  
              RE=F(ISFX+IB)*F(ISR+IB)+(1.-F(ISFX+IB))*F(ISR+IP)        
              F(ISCW+IB)=F(ISDEN+IB)*RE*(F(ISU+IB)*F(ISB11W+IB)               
     >                   +F(ISV+IB)*F(ISB12W+IB))              
	  END DO
	END IF

        IF(MINLS.NE.0) THEN
           DO IM=1,MINLS                                               
              IP=INLES(IM)                                             
              IB=IP-NI                                                 
              RS=F(ISFY+IP)*F(ISR+IP)+(1.-F(ISFY+IP))*F(ISR+IB)        
              F(ISCS+IP)=F(ISDEN+IB)*RS*(F(ISU+IB)*F(ISB21S+IP)               
     >                   +F(ISV+IB)*F(ISB22S+IP))              
	  END DO
	END IF

        IF(MINLN.NE.0) THEN
           DO IM=1,MINLN                                               
              IP=INLEN(IM)                                             
              IB=IP+NI                                                 
              RN=F(ISFY+IB)*F(ISR+IB)+(1.-F(ISFY+IB))*F(ISR+IP)        
              F(ISCS+IB)=F(ISDEN+IB)*RN*(F(ISU+IB)*F(ISB21S+IB)               
     >                   +F(ISV+IB)*F(ISB22S+IB))              
	  END DO
	END IF

        IF(MPCTW.NE.0) THEN 
        DO  IM=1,MPCTW 
               IP=INPCW(IM)
               IB=IP-1                                                  
               RE=F(ISFX+IP)*F(ISR+IP)+(1.-F(ISFX+IP))*F(ISR+IB)        
        F(ISCW+IP)=F(ISDEN+IB)*RE*(F(ISU+IB)*F(ISB11W+IP)               
     >                            +F(ISV+IB)*F(ISB12W+IP))              
	END DO
	END IF

        IF(MPCTE.NE.0) THEN 
        DO  IM=1,MPCTE 
               IP=INPCE(IM)
               IB=IP+1                                                  
               RE=F(ISFX+IB)*F(ISR+IB)+(1.-F(ISFX+IB))*F(ISR+IP)        
        F(ISCW+IB)=F(ISDEN+IB)*RE*(F(ISU+IB)*F(ISB11W+IB)               
     >                            +F(ISV+IB)*F(ISB12W+IB))              
	END DO
	END IF

        IF(MPCTS.NE.0) THEN 
        DO  IM=1,MPCTS 
               IP=INPCS(IM)
               IB=IP-NI                                                  
               RS=F(ISFY+IP)*F(ISR+IP)+(1.-F(ISFY+IP))*F(ISR+IB)
        F(ISCS+IP)=F(ISDEN+IB)*RS*(F(ISU+IB)*F(ISB21S+IP)               
     >                            +F(ISV+IB)*F(ISB22S+IP))              
	END DO
	END IF

        IF(MPCTN.NE.0) THEN 
        DO  IM=1,MPCTN 
               IP=INPCN(IM)
               IB=IP+NI                                                  
               RS=F(ISFY+IB)*F(ISR+IB)+(1.-F(ISFY+IB))*F(ISR+IP)
        F(ISCS+IB)=F(ISDEN+IB)*RS*(F(ISU+IB)*F(ISB21S+IB)               
     >                            +F(ISV+IB)*F(ISB22S+IB))              
	END DO
	END IF

      RETURN                                                            
      END                                                               
C                                                                        
C***********************************************************************
C 
        SUBROUTINE WALL(IFLAG)
C 
C COMPUTES SOURCE TERMS FOR NEAR WALL EFFECTS IN ALL EQUATIONS
C
C***********************************************************************
C
       INCLUDE 'com2d'
C        
C   IFLAG = 1 FOR WALL AT NORTH                     
C                              
C   IFLAG = 2 FOR WALL AT SOUTH 
C
C   IFLAG = 3 FOR WALL AT EAST
C
C   IFLAG = 4 FOR WALL AT WEST
C        
      SQR(VALUE)=VALUE*VALUE                                            
      VALUE=0.
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
        AFLAG1=FLOAT(MFLAG1)                                            
        AFLAG2=FLOAT(MFLAG2)                                            
        AFLAG3=FLOAT(MFLAG3)                                            
        AFLAG4=FLOAT(MFLAG4)                                            
                                                                        
      MWAL   =MWALN*MFLAG1 + MWALS*MFLAG2 + MWALE*MFLAG3 + MWALW*MFLAG4 
      ISANOW = ISAN*MFLAG1 +  ISAS*MFLAG2 +  ISAE*MFLAG3 +  ISAW*MFLAG4 

       IF(IFLAG.EQ.2.OR.IFLAG.EQ.3) THEN    ! WALL AT SOUTH (2) OR AT EAST (3)
       IM_B=MWAL
       IM_E=1
       INTV=-1
       END IF
      
       IF(IFLAG.EQ.1.OR.IFLAG.EQ.4) THEN    ! WALL AT NORTH (1) OR AT WEST (4)
       IM_B=1   
       IM_E=MWAL
       INTV=+1 
       END IF

      DO IM = IM_B,IM_E,INTV

      IP  = INDWN(IM)*MFLAG1 + INDWS(IM)*MFLAG2
     *      + INDWE(IM)*MFLAG3 + INDWW(IM)*MFLAG4
      IB  =IP+( NI *  MFLAG1 -  NI* MFLAG2 +      MFLAG3 -     MFLAG4)
      ID1 =IP+((NI+1)*MFLAG1 +      MFLAG2+(NI+1)*MFLAG3 +  NI*MFLAG4)
      ID2 =IP+( NI *  MFLAG1 +   0* MFLAG2 +      MFLAG3 +  0 *MFLAG4)
C
C     NORMAL VECTOR, AREA AND NORMAL DISTANCE
C
             DX=F(ISCOX+ID1)-F(ISCOX+ID2)
             DY=F(ISCOY+ID1)-F(ISCOY+ID2)
             DL=SQRT(DX*DX+DY*DY)+SMALL
             AN1= DY/DL
             AN2=-DX/DL                           
            
             RNO=F(ISFY+IB)*F(ISR+IB)+(1.-F(ISFY+IB))*F(ISR+IP)
             RS =F(ISFY+IP)*F(ISR+IP)+(1.-F(ISFY+IP))*F(ISR+IB)
             RE =F(ISFX+IB)*F(ISR+IB)+(1.-F(ISFX+IB))*F(ISR+IP)
             RW =F(ISFX+IP)*F(ISR+IP)+(1.-F(ISFX+IP))*F(ISR+IB)
             RN= RNO*AFLAG1 + RS*AFLAG2 + RE*AFLAG3 + RW*AFLAG4
             AREAN=SQRT(SQR(F(ISB21S+IB))+SQR(F(ISB22S+IB)))*RNO
             AREAS=SQRT(SQR(F(ISB21S+IP))+SQR(F(ISB22S+IP)))*RS
             AREAE=SQRT(SQR(F(ISB11W+IB))+SQR(F(ISB12W+IB)))*RE
             AREAW=SQRT(SQR(F(ISB11W+IP))+SQR(F(ISB12W+IP)))*RW
             AREA = AREAN*AFLAG1 + AREAS*AFLAG2 + AREAE*AFLAG3 + AREAW*AFLAG4
         
             DELX= F(ISCOX+ID1)-F(ISCOX+ID2)
             DELY= F(ISCOY+ID1)-F(ISCOY+ID2)
             XU = F(ISCOX+ID2)
             YU = F(ISCOY+ID2)
             XP = F(ISX+IP)
             YP = F(ISY+IP)
             XD = F(ISCOX+ID1)
             YD = F(ISCOY+ID1)

             DELYP=YD-YP
             DELXP=XD-XP
             ADNM =  DELX*DELX+DELY*DELY
             ANUM = (DELYP*DELY+DELXP*DELX)
             ALF = ANUM/ADNM
            
             XN = ALF*XU+(1.-ALF)*XD 
             YN = ALF*YU+(1.-ALF)*YD 
             DEL = (XN-XP)*(XN-XP)+(YN-YP)*(YN-YP)
             DELTA = SQRT(DEL)
        

          
          IF(TURBMOD.EQ.0) CALL LAMINAR(IP,IB,ISANOW,IFLAG,AN1,AN2,DELTA,AREA)

          IF(TURBMOD.EQ.1) CALL KEPSSTD (IP,IB,ISANOW,IFLAG,AN1,AN2,DELTA,AREA)

          IF(TURBMOD.EQ.2) CALL KEPS2LR(IP,IDD2,IB,ISANOW,IFLAG,AN1,AN2,DELTA,AREA)

          IF(TURBMOD.EQ.3) CALL KEPSCHIEN(IP,IB,ISANOW,IFLAG,AN1,AN2,DELTA,AREA)

          IF(TURBMOD.EQ.4) CALL KOMEGA(IP,IB,ISANOW,IFLAG,AN1,AN2,DELTA,AREA)

          IF(TURBMOD.EQ.5) CALL V2F(IP,IB,ISANOW,IFLAG,AN1,AN2,DELTA,AREA)

          IF(TURBMOD.EQ.6) CALL SPL_ALM(IP,IB,ISANOW,IFLAG,AN1,AN2,DELTA,AREA)

          IF(TURBMOD.EQ.7) CALL SST(IP,IB,ISANOW,IFLAG,AN1,AN2,DELTA,AREA)
 
          END DO    ! END OF MWAL LOOP

          RETURN

          END 
C   
C***************************************************************
C   
          SUBROUTINE LAMINAR(IP,IB,ISANOW,IFLAG,AN1,AN2,DELTA,AREA)
C   
C***************************************************************
C   
           INCLUDE 'com2d'
C
          VR=AN2*F(ISU+IP)-AN1*F(ISV+IP)
          UWAL=F(ISU+IB)
          VWAL=F(ISV+IB)
          TWAL=F(ISS+IB)  ! NEEDS TO BE REPLACED BY CONDUCTION FORMULA
          VPARL=AN2*UWAL-AN1*VWAL-VR
          TAUWAL(IP,IFLAG)=VISCOS*VPARL/DELTA
          UTAU=SQRT(ABS(TAUWAL(IP,IFLAG))/F(ISDEN+IP))
          YPLUS(IP,IFLAG)=F(ISDEN+IP)*UTAU*DELTA/VISCOS
          TMULT=VISCOS/DELTA

C         
          A11=1.-AN1*AN1
          A22=1.-AN2*AN2
          A12=AN1*AN2
          SHEAR=TMULT*AREA
           
          IF(NN.EQ.2) THEN
          
C      U-VELOCITY SOURCES
C
       F(ISSP+IP)=F(ISSP+IP)-SHEAR*A11
       F(ISSU+IP)=F(ISSU+IP)+SHEAR*(F(ISV+IP)*A12+UWAL)
C
      ELSEIF(NN.EQ.3) THEN
C
C      V-VELOCITY SOURCES
C
       F(ISSP+IP)=F(ISSP+IP)-SHEAR*A22
       F(ISSU+IP)=F(ISSU+IP)+SHEAR*(F(ISU+IP)*A12+VWAL)
C                                                                
           ELSEIF (NN.EQ.4) THEN 
C
C        WR-VELOCITY SOURCES
C
       RAT=F(ISR+IP)/(RN+SMALL)
       F(ISSP+IP)=F(ISSP+IP)-SHEAR/(RAT+SMALL)
       F(ISSU+IP)=F(ISSU+IP)+SHEAR*RAT*F(ISWR+IB)      
  
          ELSEIF (NN.EQ.11) THEN
C
C        SCALAR  SOURCES
C
C        EVAlUVATE SOLID LIQUID INTERFACE TEMPERATURE
        
        IF(IFLAG .EQ. 4)THEN

        IP1 = 8 + (J-1)*NI

        END IF
        IF(IP .EQ. IP1) THEN
        IP2 = IP1-1
        DXF = F(ISX+IP1)-F(ISCOX+IP1)
        DXS = F(ISCOX+IP1)-F(ISX+IP2)
        TFL = F(ISS+IP)
        TFS = F(ISS+IP1)
        T11 = 1/DXF
        T12 = 2/DXS
        T21 = TFL/DXF
        T22 = 2*TFS/DXS
        TT11 = T11+T12
        TT22 = T21+T22
        TT12 = 1/TT11
        TWAL = TT12*TT22
        
        END IF
         COND=VISCOS*AREA/PRLAM/DELTA
         F(ISSP+IP)=F(ISSP+IP)-COND
         F(ISSU+IP)=F(ISSU+IP)+COND*TWAL
       
         WRITE(987,*) TWAL

           END IF    ! END OF IF LOOP  FOR FLOW VARIABLES
C

C DELINK THE COEFFICIENTS CONNECTING THE WALL NODE
C
            F(ISANOW+IP)=0.
           
            RETURN

          RETURN
          END      
C   
C**********************************************************************
C   
      SUBROUTINE SYMM(IFLAG)     
C                                                                        
C   DELINK COEFFICIENTS AND SET FLUXES TO ZERO
C                              
C***********************************************************************
C                                                                       
      INCLUDE 'com2d'                                               
C                                                                       
C   IFLAG = 1 FOR SYMMETRY AT NORTH                     
C                              
C   IFLAG = 2 FOR SYMMETRY AT SOUTH 
C
C   IFLAG = 3 FOR SYMMETRY AT EAST
C
C   IFLAG = 4 FOR SYMMETRY AT WEST
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
                                                                        
      MSYM   =MSYMN*MFLAG1 + MSYMS*MFLAG2 + MSYME*MFLAG3 + MSYMW*MFLAG4 
      ISANOW = ISAN*MFLAG1 +  ISAS*MFLAG2 +  ISAE*MFLAG3 +  ISAW*MFLAG4 
                                                                        
      DO IM=1,MSYM                                                 
         I   = INSYN(IM)*MFLAG1 + INSYS(IM)*MFLAG2                      
     >       + INSYE(IM)*MFLAG3 + INSYW(IM)*MFLAG4  
    
C      FOR SWIRL ONLY  
              
       IF(NN.EQ.4 .AND.IFLAG.EQ.2) THEN                                 
          F(ISSU+I)=BESUC(I)                                            
     >          -2.*F(ISWR+I)*(F(ISB12P+I)*(F(ISVISW+I+1 )-F(ISVISW+I)) 
     >                        +F(ISB22P+I)*(F(ISVISS+I+NI)-F(ISVISS+I)))
          F(ISSP+I)=BESPC(I)-4.*F(ISVIS+I)*F(ISVOLP+I)/                 
     >             (F(ISR+I)*F(ISR+I)+SMALL)                            
      ELSE                                                              

         IDB =    NI *MFLAG1   -NI *MFLAG2+       MFLAG3    -   MFLAG4  

          INB=ISPP+(NN-1)*NIJ+I
          IB=INB+IDB
         F(IB)=F(INB) 
 
      ENDIF                                                             
      F(ISANOW+I)=0.  
      END DO
                                                  
      RETURN                                                            
      END                                                               
C
C**********************************************************************
C                            
      SUBROUTINE CONVECTIVE_OUT(IFLAG)       
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
         delx=F(ISX+II+IND_BNDRY)-F(ISX+II)
         dely=F(ISy+II+IND_BNDRY)-F(ISy+II)
         deln=sqrt(delx*delx+dely*dely)+small

        ratio=delt/deln
        alfa=(uref*ratio)/(1+ratio*uref)
        if (nn .eq. 2) then
          F(IB)=alfa*f(inb)+(1-alfa)*u(ii+ind_bndry,iblk,nlevel)
        end if
        if (nn .eq. 3) then
          F(IB)=alfa*f(inb)+(1-alfa)*v(ii+ind_bndry,iblk,nlevel)
        end if
        if (nn .eq. 1 .or. nn .gt. 3) f(ib)=f(inb)
      F(ISANOW+II)=0.  
      END DO                                                  
      RETURN                                                            
      END                                                               
C
C**********************************************************************
C                            
      SUBROUTINE PCONST(IFLAG)       
C                                                                        
C EXTRAPOLATE FLOW VARIABLES AT OUTFLOW BOUNDARY
C                              
C***********************************************************************
C                                           
      INCLUDE 'com2d'                                               
C                                                                       
C   IFLAG = 1 FOR PRESSURE CONSTANT AT NORTH                     
C                              
C   IFLAG = 2 FOR PRESSURE CONSTANT AT SOUTH 
C
C   IFLAG = 3 FOR PRESSURE CONSTANT AT EAST
C
C   IFLAG = 4 FOR PRESSURE CONSTANT AT WEST
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
                                                                        
      MPCT   =MPCTN*MFLAG1 + MPCTS*MFLAG2 + MPCTE*MFLAG3 + MPCTW*MFLAG4 
      ISANOW = ISAN*MFLAG1 +  ISAS*MFLAG2 +  ISAE*MFLAG3 +  ISAW*MFLAG4 
       IH=ISPP+(NN-1)*NIJ                                                      
    
       DO IM=1,MPCT                                                 
         IP  = INPCN(IM)*MFLAG1 + INPCS(IM)*MFLAG2                      
     *       + INPCE(IM)*MFLAG3 + INPCW(IM)*MFLAG4                      
         IB  =IP +(NI * MFLAG1 - NI* MFLAG2 +  MFLAG3 -   MFLAG4)  
         F(ISAPU+IB)=2.*F(ISAPU+IP)
         F(ISAPV+IB)=2.*F(ISAPV+IP)

C
C  TREATMENT OF THE ENTRAINMENT BOUNDARY FOR VARIABLES OTHER THAN VELOCITY
C IN CASE THE ENTRAINMENT IS FROM SURROUNDING TO THE DOMAIN THE BOUNDARY
C VALUES ARE LEFT UNALTERED OTHERWISE THE STREAMWISE GRADIENT IS SET TO ZERO
C THROUGH EQUATING THE BOUNDARY AND THE NEAR BOUNDARY VALUES OF THE VARIABLE
  
        IF(NN.GE.5) THEN
          
         CONV_IND_BNDRY =  NI *MFLAG1 -0 *MFLAG2+  MFLAG3    - 0*MFLAG4
	 IC_IB=IP+CONV_IND_BNDRY 
         CONV_B=F(ISCS+IC_IB)*(MFLAG1+MFLAG2)+F(ISCW+IC_B)*(MFLAG3+MFLAG4)
         IF(IFLAG.EQ.1.AND.CONV_B.LT.0.) GO TO 2345
         IF(IFLAG.EQ.2.AND.CONV_B.GT.0.) GO TO 2345
         IF(IFLAG.EQ.3.AND.CONV_B.LT.0.) GO TO 2345
         IF(IFLAG.EQ.4.AND.CONV_B.GT.0.) GO TO 2345

          F(IH+IB)=F(IH+IP)
          F(ISANOW+IP)=0.
    
 2345   CONTINUE
          END IF

      END DO
        
      RETURN                                                            
      END                                                               
C
C**********************************************************************
C                            
      SUBROUTINE BLOCK     
C                            
C     COMPUTES SOURCE TERMS FOR BLOCKING NON-FLUID  
C     ZONES INSIDE THE COMPUTATION DOMAIN 
C                              
C***********************************************************************
C                                                                        
      INCLUDE 'com2d'                                                 
C          
      DO I=IIS,IIE                                                 
      KBLKM=1-KBLK(I)                                              
      F(ISSP+I)=KBLK(I)*F(ISSP+I)-KBLKM*GREAT                           
      F(ISSU+I)=KBLK(I)*F(ISSU+I)+KBLKM*GREAT*F(ISPP+(NN-1)*NIJ+I)
      END DO  
      RETURN                                                            
      END                                                               
C                                                                        
C***********************************************************************
C   
      SUBROUTINE BOUNDS(N)    
C                                          
C     CALL APPROPRIATE BOUNDARY CONDITIONS FOR VARIABLE N
C                              
C***********************************************************************
C
      INCLUDE 'com2d' 
C                                              
      NN=N
     
      IF(NN.EQ.1) GO TO 500
                                          
      IF(MWALW.NE.0) CALL WALL(4)                                       
      IF(MWALN.NE.0) CALL WALL(1)                                       
      IF(MWALE.NE.0) CALL WALL(3)                                       
      IF(MWALS.NE.0) CALL WALL(2)

      IF(MSYME.NE.0.AND.NN.NE.NSYME) CALL SYMM(3)                       
      IF(MSYMW.NE.0.AND.NN.NE.NSYMW) CALL SYMM(4)                       
      IF(MSYMN.NE.0.AND.NN.NE.NSYMN) CALL SYMM(1)                       
      IF(MSYMS.NE.0.AND.NN.NE.NSYMS) CALL SYMM(2)                       
                                                                        
         IF(MEXIE.NE.0) CALL OUT_NEWPREV(3)
         IF(MEXIW.NE.0) CALL OUT_NEWPREV(4)
         IF(MEXIN.NE.0) CALL OUT_NEWPREV(1)
         IF(MEXIS.NE.0) CALL OUT_NEWPREV(2)

C        IF(MINLE.NE.0) CALL OUT_NEW(3)
C        IF(MINLW.NE.0) CALL OUT_NEW(4)         
C        IF(MINLN.NE.0) CALL OUT_NEW(1)                             
C        IF(MINLS.NE.0) CALL OUT_NEW(2)                             

         IF(MINLE.NE.0) CALL FARFIELD(3)
         IF(MINLW.NE.0) CALL FARFIELD(4)         
         IF(MINLN.NE.0) CALL FARFIELD(1)                             
         IF(MINLS.NE.0) CALL FARFIELD(2)                             
                                                                        
  500 IF(IBLOC.NE.0) CALL BLOCK


      IF(MCYCN.NE.0) CALL CYC(1)
      IF(MCYCS.NE.0) CALL CYC(2)
      IF(MCYCE.NE.0) CALL CYC(3)
      IF(MCYCW.NE.0) CALL CYC(4)
                                         
      IF(MPCTN.NE.0) CALL PCONST(1)
      IF(MPCTS.NE.0) CALL PCONST(2)
      IF(MPCTE.NE.0) CALL PCONST(3)
      IF(MPCTW.NE.0) CALL PCONST(4)
                                     
      RETURN                                                            
      END                                                               
C                                                                        
C***********************************************************************
C   
      SUBROUTINE CDFLUX(ISVAR,N)
C                                        
C     COMPUTES CROSS-DERIVATIVE FLUXES FOR VARIABLE N
C                        
C***********************************************************************
C    
          INCLUDE 'com2d'                                                 
C                                               
         PRNL=1.
          IF(N.EQ.11)PRNL=1./PRLAM
         PRNT=1./PR(N)

C                                                
      DO J=jsolv_b,jsolv_e
      DO I=isolv_b,isolv_e
      II=I+(J-1)*NI
C                                                 
C --- INTERPOLATION-FACTORS  
C                                           
      FXP   =F(ISFX+II)                                                 
      FXIP  =F(ISFX+II+1)                                               
      FXJP  =F(ISFX+II+NI)                                              
      FXJM  =F(ISFX+II-NI)                                              
      FXIPJP=F(ISFX+II+1+NI)                                            
      FXIPJM=F(ISFX+II+1-NI)                                            
                                                                        
      FYP   =F(ISFY+II)                                                 
      FYJP  =F(ISFY+II+NI)                                              
      FYIP  =F(ISFY+II+1)                                               
      FYIM  =F(ISFY+II-1)                                               
      FYIPJP=F(ISFY+II+1+NI)                                            
      FYIMJP=F(ISFY+II-1+NI)                                            
C                                                                        
C --- NEIGHBOUR-VALUES OF DEPENDENT VARIABLE
C                            
      UIJK =F(ISVAR+II)                                                 
      UIP  =F(ISVAR+II+1)                                               
      UIPJP=F(ISVAR+II+1+NI)                                            
      UIPJM=F(ISVAR+II+1-NI)                                            
      UIM  =F(ISVAR+II-1)                                               
      UIMJP=F(ISVAR+II-1+NI)                                            
      UIMJM=F(ISVAR+II-1-NI)                                            
      UJP  =F(ISVAR+II+NI)                                              
      UJM  =F(ISVAR+II-NI)                                              
C                                                                        
C --- DIFFERENCES OF DEPENDENT VARIABLE AT CV-SIDES
C                     
      DUEW=FXIP*UIP+(1.-FXIP)*UIJK-(FXP*UIJK+(1.-FXP)*UIM)              
      DUNS=FYJP*UJP+(1.-FYJP)*UIJK-(FYP*UIJK+(1.-FYP)*UJM)              
                                                                        
C     ALONG NORTH-SOUTH ON EAST/WEST 
                                   
      UNSE=(1.-FXIP)*DUNS+    FXIP*((FYIPJP*UIPJP+(1.-FYIPJP)*UIP)      
     >                             -(FYIP  *UIP  +(1.-FYIP  )*UIPJM))   
      UNSW=    FXP  *DUNS+(1.-FXP)*((FYIMJP*UIMJP+(1.-FYIMJP)*UIM)      
     >                             -(FYIM  *UIM  +(1.-FYIM  )*UIMJM))   
                                                                        
C     ALONG EAST-WEST ON NORTH/SOUTH      
                              
      UEWN=(1.-FYJP)*DUEW+    FYJP*((FXIPJP*UIPJP+(1.-FXIPJP)*UJP)      
     >                             -(FXJP  *UJP  +(1.-FXJP  )*UIMJP))   
      UEWS=    FYP  *DUEW+(1.-FYP)*((FXIPJM*UIPJM+(1.-FXIPJM)*UJM)      
     >                             -(FXJM  *UJM  +(1.-FXJM  )*UIMJM))   
                                                                        
C     D'S AT DIFFERENT CELL FACES 
                                      
        RE=FXIP*F(ISR+II+1 )+(1.-FXIP)*F(ISR+II)                        
        RW=FXP *F(ISR+II   )+(1.-FXP )*F(ISR+II-1)                      
        RN=FYJP*F(ISR+II+NI)+(1.-FYJP)*F(ISR+II)                        
        RS=FYP *F(ISR+II   )+(1.-FYP )*F(ISR+II-NI)                     
      D12E=RE*RE*(F(ISB21W+II+1)*F(ISB11W+II+1)+F(ISB22W+II+1)*F(ISB12W+II+1))  
      D12W=RW*RW*(F(ISB21W+II)*F(ISB11W+II)+F(ISB22W+II)*F(ISB12W+II))  
      D21N=RN*RN*(F(ISB21S+II+NI)*F(ISB11S+II+NI)
     >     +F(ISB22S+II+NI)*F(ISB12S+II+NI))  
      D21S=RS*RS*(F(ISB21S+II)*F(ISB11S+II)+F(ISB22S+II)*F(ISB12S+II))  
                                                                        
C    RECIPROCAL OF CELL VOLUMES AROUND DIFFERENT FACES
                 
      RVE=1./(0.5*(F(ISVOLP+II)+F(ISVOLP+II+1 ))+SMALL)                 
      RVW=1./(0.5*(F(ISVOLP+II)+F(ISVOLP+II-1 ))+SMALL)                 
      RVN=1./(0.5*(F(ISVOLP+II)+F(ISVOLP+II+NI))+SMALL)                 
      RVS=1./(0.5*(F(ISVOLP+II)+F(ISVOLP+II-NI))+SMALL)                 
                                                                        
C     CROSS-DERIVATIVE TERMS
C    
C     SPECIAL TREATMENT FOR DIFFUSION COEFFICIENT IN SST MODEL
C
       IF(TURBMOD.EQ.7.AND.N.EQ.(NVTE-NVPP+1)) THEN
C      PR_N=PR1_TE
       PR_N=F1_SST(II)*PR1_TE+(1-F1_SST(II))*PR2_TE
       PRNT=1./PR_N
C      IF(I.EQ.ISOLV_E.AND.IBSWP.EQ.NBSWP_END) THEN 
C      WRITE(2000,*)J,PRNT,F1_SST(II),F2_SST(II)
C      END IF
       END IF

      IF(TURBMOD.EQ.7.AND.N.EQ.(NVOM-NVPP+1)) THEN
C     PR_N=PR1_OM
      PR_N=F1_SST(II)*PR1_OM+(1-F1_SST(II))*PR2_OM
      PRNT=1./PR_N
C     IF(I.EQ.ISOLV_E.AND.IBSWP.EQ.NBSWP_END) THEN 
C     WRITE(2100,*)J,PRNT,F1_SST(II),F2_SST(II)
C     END IF
      END IF

        VISTW=F(ISVISW+II)-VISCOS
        VISTE=F(ISVISW+II+1)-VISCOS
        VISEFW=VISCOS*PRNL+VISTW*PRNT
        VISEFE=VISCOS*PRNL+VISTE*PRNT

        VISTS=F(ISVISS+II)-VISCOS
        VISTN=F(ISVISS+II+NI)-VISCOS
        VISEFS=VISCOS*PRNL+VISTS*PRNT
        VISEFN=VISCOS*PRNL+VISTN*PRNT
C
C  DELINKING THE VISCOSITY AT WALL BOUNDRIES
C
	if (kw(ii) .eq. 11) visefw=0.
	if (ke(ii) .eq. 11) visefe=0.
	if (ks(ii) .eq. 11) visefs=0.
	if (kn(ii) .eq. 11) visefn=0.

       IF(N.EQ.10) THEN

         ANUW=FXP *F(ISNUT+II)+(1.-FXP )*F(ISNUT+II-1)
         ANUE=FXIP*F(ISNUT+II+1)+(1.-FXIP)*F(ISNUT+II)
         VISEFW=(VISCOS+ANUW)*PRNT
         VISEFE=(VISCOS+ANUE)*PRNT

         ANUS=FYP *F(ISNUT+II)+(1.-FYP )*F(ISNUT+II-NI)
         ANUN=FYJP*F(ISNUT+II+NI)+(1.-FYJP)*F(ISNUT+II)
         VISEFS=(VISCOS+ANUS)*PRNT
         VISEFN=(VISCOS+ANUN)*PRNT

        ENDIF

        SUIJ=visefe*RVE*D12E*UNSE-visefw*RVW*D12W*UNSW  
     >      +visefn*RVN*D21N*UEWN-visefs*RVS*D21S*UEWS
        SUIJ=SUIJ*KBLK(II) 
        F(ISSU+II)=F(ISSU+II)+SUIJ
	end do
	end do
C   
      RETURN
C                           
      END                                                               
C                                                                        
C***********************************************************************
C
      SUBROUTINE CHECK1(N)                                              
C
C     CHECKS SYNTAX KIND OF USER-ERRORS IN THE INPUT DATA
C     AT DIFFERENT POINTS OF THE CODE
C
C***********************************************************************
C
      WRITE(6,*) ' '                                                    
      WRITE(6,*) ' USER-ERROR!   JOB STOPS!'                            
      WRITE(6,*) ' ========================'                            
      GOTO (1,2,3,4,5,6) N                                              
    1 WRITE(6,*) ' NX > NXF  OR  NY > NYF  OR  NMAX < MAX(NXF,NYF)'     
      STOP                                                              
    2 WRITE(6,*) ' DIMENSION OF ARRAYS IN COMMON/BND2/ IS NOT ENOUGH'   
      STOP                                                              
    3 WRITE(6,*) ' CHARACTER DEFINED IN CBND(N,2) IS UNRECOGNIZED!'     
      STOP                                                              
    4 WRITE(6,*) ' CHARACTER DEFINED IN CBND(N,1) IS UNRECOGNIZED!'     
      STOP                                                              
    5 WRITE(6,*) ' NRBNDA < 4  OR  NRBNDA > NRBND'                      
      STOP                                                              
    6 WRITE(6,*) ' NUMBER OF GRID POINTS STORED IN THE FILE OF',        
     >           ' CHANNEL 7 IS NOT CORRECT'                            
      STOP                                                              
      END      
C
C***********************************************************************
C     
      SUBROUTINE COFCON                                                 
C
C     CALCULATE COEFFICIENTS FOR PRESSURE-CORRECTION EQUATION AND       
C     CONVECTIVE FLUXES AT CELL FACE USING MOMENTUM INTERPOLATION      
C
C***********************************************************************
C
      INCLUDE 'com2d' 
C  
                                         
         RELAXM=1.-RELAX(NVU-NVPP+1)    
                                         
      DO J=2,NJ  
      DO I=2,NI
      II=I+(J-1)*NI
C                               
C *** X-DIRECTION: EAST AND WEST SIDE                                   
C
       FXP=F(ISFX+II)                                                   
       FXM=1.-FXP                                                       
        RW= FXP*F(ISR  +II)+FXM*F(ISR  +II-1)                           
       DRW=(FXP*F(ISDEN+II)+FXM*F(ISDEN+II-1))*RW                       
      Q11W=(FXP*F(ISAPU+II)+FXM*F(ISAPU+II-1))*F(ISB11W+II)*RW          
      Q12W=(FXP*F(ISAPV+II)+FXM*F(ISAPV+II-1))*F(ISB12W+II)*RW          
      D11W=1.-FXP*F(ISSUMU+II)-FXM*F(ISSUMU+II-1)                       
      D12W=1.-FXP*F(ISSUMV+II)-FXM*F(ISSUMV+II-1)                       
      Q11W=Q11W/(D11W+SMALL)                                            
      Q12W=Q12W/(D12W+SMALL)                                            
C
C     COEFFICIENTS AW AND AE FOR PP-EQUATION                            
C
C
C  FOR CUT BOUNDARY (15), PERIODIC BOUNDARY (16) AND PRESSURE
C CONSTANT BOUNDARY (17), the COEFFICIENTS ARE NOT CUTOFF I.E. IFACT=1
C WHERE AS FOR THE OTHER BOUNDARIES (WALL, INFLOW,OUTFLOW,SYMMETRY,etc.)
C THE BOUNDARY COEFFICIENTS ARE CUT OFF (IFACT=0)
C
	ifact=1-MIN(1,KW(II))
	if (kw(ii) .ge. 15) ifact=1
      F(ISAW+II)=ifact*DRW*(Q11W*F(ISB11W+II)+Q12W*F(ISB12W+II))    
	ifact=1-MIN(1,KE(II-1))
	if (ke(ii-1) .ge. 15) ifact=1
      F(ISAE+II-1)=ifact*F(ISAW+II)                       
C
C     PRESSURE-INTERPOLATION ON EAST/WEST FACE                          
C
c     PW =F(ISFX+II)*F(ISP+II)+(1.-F(ISFX+II))*F(ISP+II-1)  
c     PE  =F(ISFX+II+1)*F(ISP+II+1)+(1.-F(ISFX+II+1))*F(ISP+II) 

      CALL DPRS(II,F(ISP+1),PW,PE,PS,PN)

c     PWW =F(ISFX+II-1)*F(ISP+II-1)+(1.-F(ISFX+II-1))*F(ISP+II-2)

      CALL DPRS(II-1,F(ISP+1),PWW,PPE,PPS,PPN)
C
C     VELOCITY-INTERPOLATIONS ON WEST FACE                              
C
      DPWE=PW-PE                                                        
      DPWWW=PWW-PW                                                      
      DPW=F(ISP+II-1)-F(ISP+II)                                         
      HI_U=F(ISU+II  )-F(ISB11P+II  )*F(ISR+II  )*F(ISAPU+II  )*DPWE  
      HIM1_U=F(ISU+II-1)-F(ISB11P+II-1)*F(ISR+II-1)*F(ISAPU+II-1)*DPWWW 
      UW=FXP*HI_U+FXM*HIM1_U+Q11W*DPW*D11W                                  
      HI_V=F(ISV+II  )-F(ISB12P+II  )*F(ISR+II  )*F(ISAPV+II  )*DPWE  
      HIM1_V=F(ISV+II-1)-F(ISB12P+II-1)*F(ISR+II-1)*F(ISAPV+II-1)*DPWWW 
      VW=FXP*HI_V+FXM*HIM1_V+Q12W*DPW*D12W                                  
C
C UPDATE THE BOUNDARY VELOCITIES AS THE FACE VELOCITIES FOR THAT CV
C FOR PRESSURE CONSTANT BOUNDARY
C
      if (kw(ii) .eq. 17)then 
	  pw=0.
	  pe=f(isp+ii)
	  DPW=PW-PE
          UW=HI_U+Q11W*DPW*D11W                                  
          VW=HI_V+Q12W*DPW*D12W                                  
          f(isu+ii-1)=uw 
          f(isv+ii-1)=vw 
      end if
      if (ke(ii-1) .eq. 17)then
	  pw=f(isp+ii-1)
	  PE=0.
	  DPW=PW-PE
          UE=HIM1_U+Q11W*DPW*D11W                                  
          VE=HIM1_V+Q12W*DPW*D12W 
          f(isu+ii)=ue 
          f(isv+ii)=ve 
      end if
C
C     CONVECTIVE FLUX DUE TO MOMENTUM AND LINEAR INTERPOLATIONS         
C
      CW=DRW*(UW*F(ISB11W+II)+VW*F(ISB12W+II))                          
      UWL=FXP*F(ISU+II)+FXM*F(ISU+II-1)                             
      VWL=FXP*F(ISV+II)+FXM*F(ISV+II-1)                             
      CWL=DRW*(UWL*F(ISB11W+II)+VWL*F(ISB12W+II))                       
C
C     UNDER-RELAXATION                                                  
C
       F(ISCW+II)=CW+RELAXM*(F(ISCW+II)-CWL)                             
C
C     Y-DIRECTION: NORTH AND SOUTH FACES                                
C
       FYP=F(ISFY+II)                                                   
       FYM=1.-FYP                                                       
       RS= FYP*F(ISR  +II)+FYM*F(ISR  +II-NI)                          
       DRS=(FYP*F(ISDEN+II)+FYM*F(ISDEN+II-NI))*RS                      
      Q21S=(FYP*F(ISAPU+II)+FYM*F(ISAPU+II-NI))*F(ISB21S+II)*RS         
      Q22S=(FYP*F(ISAPV+II)+FYM*F(ISAPV+II-NI))*F(ISB22S+II)*RS         
      D21S=1.-FYP*F(ISSUMU+II)-FYM*F(ISSUMU+II-NI)                      
      D22S=1.-FYP*F(ISSUMV+II)-FYM*F(ISSUMV+II-NI)                      
      Q21S=Q21S/(D21S+SMALL)                                            
      Q22S=Q22S/(D22S+SMALL)                                            
C
C     COEFFICIENTS AS AND AN FOR PP-EQUATION                            
C
	jfact=1-MIN(1,KS(II))
	if (ks(ii) .ge. 15) jfact=1
      F(ISAS+II)=jfact*DRS*(Q21S*F(ISB21S+II)+Q22S*F(ISB22S+II))              
	jfact=1-MIN(1,KN(II-NI))
	if (kn(ii-ni) .ge. 15) jfact=1
      F(ISAN+II-NI)=jfact*F(ISAS+II)                     
C
C     PRESSURE-INTERPOLATION ON NORTH/SOUTH FACE                        
C
c     PS =F(ISFY+II)*F(ISP+II)+(1.-F(ISFY+II))*F(ISP+II-NI)    
c     PN =F(ISFY+II+NI)*F(ISP+II+NI)+(1.-F(ISFY+II+NI))*F(ISP+II)

      CALL DPRS(II,F(ISP+1),PW,PE,PS,PN)

c     PSS =F(ISFY+II-NI)*F(ISP+II-NI)+(1.-F(ISFY+II-NI))*F(ISP+II-2*NI) 

      CALL DPRS(II-NI,F(ISP+1),PPW,PPE,PSS,PPN)
C
C      VELOCITY-INTERPOLATIONS ON SOUTH FACE                             
C
           IB=II-NI                                                     
         DPSN=PS-PN                                                     
        DPSSS=PSS-PS                                                    
          DPS=F(ISP+IB)-F(ISP+II)                                       
           HJ_U=F(ISU+II)-F(ISB21P+II)*F(ISR+II)*F(ISAPU+II)*DPSN         
         HJM1_U=F(ISU+IB)-F(ISB21P+IB)*F(ISR+IB)*F(ISAPU+IB)*DPSSS        
      US=FYP*HJ_U+FYM*HJM1_U+Q21S*DPS*D21S                                  
           HJ_V=F(ISV+II)-F(ISB22P+II)*F(ISR+II)*F(ISAPV+II)*DPSN         
         HJM1_V=F(ISV+IB)-F(ISB22P+IB)*F(ISR+IB)*F(ISAPV+IB)*DPSSS        
      VS=FYP*HJ_V+FYM*HJM1_V+Q22S*DPS*D22S                                  
C
C UPDATE THE BOUNDARY VELOCITIES AS THE FACE VELOCITIES FOR THAT CV
C FOR PRESSURE CONSTANT BOUNDARY
C
      if (ks(ii) .eq. 17)then 
	  ps=0.
	  pn=f(isp+ii)
	  DPS=PS-PN
          US=HJ_U+Q21S*DPS*D21S                                  
          VS=HJ_V+Q22S*DPS*D22S                                  
          f(isu+ii-ni)=us 
          f(isv+ii-ni)=vs 
      end if
      if (kn(ii-ni) .eq. 17)then 
	  ps=f(isp+ii-ni)
	  pn=0.
	  DPS=PS-PN
          UN=HJM1_U+Q21S*DPS*D21S                                  
          VN=HJM1_V+Q22S*DPS*D22S                                  
          f(isu+ii)=un 
          f(isv+ii)=vn 
      end if
C
C     CONVECTIVE FLUX DUE TO MOMENTUM AND LINEAR INTERPOLATIONS         
C
      CS=DRS*(US*F(ISB21S+II)+VS*F(ISB22S+II))                          
      USL=FYP*F(ISU+II)+FYM*F(ISU+II-NI)                            
      VSL=FYP*F(ISV+II)+FYM*F(ISV+II-NI)                            
      CSL=DRS*(USL*F(ISB21S+II)+VSL*F(ISB22S+II))                       
C
C     UNDER-RELAXATION                                                  
C
       F(ISCS+II)=CS+RELAXM*(F(ISCS+II)-CSL)

      END DO
      END DO
      
      RETURN                               
      END                                                               
C                                                                        
C***********************************************************************
C
      SUBROUTINE CWCS                                                   
C
C     INITIALISING MASS FLUXES AT CELL FACES                                    
C
C***********************************************************************
C 
      INCLUDE 'com2d'                                              
C
      DO 1040 J= 2, NJ
      DO 1040 I = 2, NI                          
	     II=I+(J-1)*NI
             FXP=F(ISFX+II)                                             
             FXM=1.-FXP                                                 
            DENW=F(ISDEN+II)*FXP+F(ISDEN+II-1)*FXM                      
              RW=F(ISR  +II)*FXP+F(ISR  +II-1)*FXM                      
              UW=F(ISU  +II)*FXP+F(ISU  +II-1)*FXM                      
              VW=F(ISV  +II)*FXP+F(ISV  +II-1)*FXM                      
             FYP=F(ISFY+II)                                             
             FYM=1.-FYP                                                 
            DENS=F(ISDEN+II)*FYP+F(ISDEN+II-NI)*FYM                     
              RS=F(ISR  +II)*FYP+F(ISR  +II-NI)*FYM                     
              US=F(ISU  +II)*FYP+F(ISU  +II-NI)*FYM                     
              VS=F(ISV  +II)*FYP+F(ISV  +II-NI)*FYM                     
      F(ISCW+II)=DENW*(UW*F(ISB11W+II)+VW*F(ISB12W+II))*RW              
 1040 F(ISCS+II)=DENS*(US*F(ISB21S+II)+VS*F(ISB22S+II))*RS              
c     
c
      END                                                               
C
C***********************************************************************
C 
      SUBROUTINE RANSOL                                                 
C
C     SOLVES ALL THE CONSERVATION EQUATIONS NUMBERED "NN"              
C     NN =  1 , 2, 3,  4,  5,  6,  7, 8, 9, 10   11                             
C     FOR   P', U, V, WR, TE, ED, OM, V2, F, NUT,S ( ANY PASSIVE SCALAR)
C                     
C***********************************************************************
C
      INCLUDE 'com2d'                                                 
C
        IIMON=IMON+(JMON-1)*NI                                          
c  
C        
C     SEQUENTIAL SOLUTION OF DIFFERENT EQUATIONS START HERE
C                                          
C *** SOLVING U-MOMENTUM EQUATION (NN=2) 
C                               
      IF(SOLVE(2)) THEN

        IF(NDSCH(2).EQ.1) CALL COEFF1(2)                            
        IF(NDSCH(2).EQ.2) CALL COEFF2(2)                            
        IF(NDSCH(2).EQ.3) CALL COEFF3(2)                            
        IF(NDSCH(2).EQ.4) CALL COEFF4(2)  
        IF(NDSCH(2).EQ.5) CALL COEFF5(2)  

        CALL SORCE2           

        CALL BOUNDS(2)

        IF(SPECL(2)) CALL SPEBC2  

        CALL PRINT(F(ISSU+1),'SU-ARRAY')
        CALL PRINT(F(ISSP+1),'SP-ARRAY')
        CALL SOLVER(ISU,2)


      ENDIF  
C
C *** SOLVING V-MOMENTUM EQUATION (NN=3)                                
C
      IF(SOLVE(3)) THEN                                                  
 
        IF(NDSCH(3).EQ.1) CALL COEFF1(3)                            
        IF(NDSCH(3).EQ.2) CALL COEFF2(3)                            
        IF(NDSCH(3).EQ.3) CALL COEFF3(3)                            
        IF(NDSCH(3).EQ.4) CALL COEFF4(3)
        IF(NDSCH(3).EQ.5) CALL COEFF5(3)

        CALL SORCE3                                                     
                                   
        CALL BOUNDS(3)                                                  

        IF(SPECL(3)) CALL SPEBC3

        CALL PRINT(F(ISSU+1),'SU-ARRAY')
        CALL PRINT(F(ISSP+1),'SP-ARRAY')
        CALL SOLVER(ISV,3)

      ENDIF                                                             
C
C *** SOLVING PRESSURE-CORRECTION-EQUATION (NN=1)                       
C     FOR SATISFACTION OF CONTINUITY
C
      IF(SOLVE(1)) THEN                                                  

        CALL BOUNDS(1)                                                  

        CALL COFCON                                                     

        CALL BFLUX                                                      

   	CALL BLOKCON
	   
        CALL BFLUX                                                      
	
        CALL SORCE1                                                     

        CALL BOUNDS(1)                                                  

        IF(SPECL(1)) CALL SPEBC1                                        

        CALL SOLVER(ISPP,1)

        CALL PVCORR
                                                
        CALL BOUNDS(1)                                                  

        CALL BFLUX    

      ENDIF
C
C *** SOLVING WR-MOMENTUM EQUATION (NN=4)                               
C
      IF(SOLVE(4)) THEN                                                  
 
        IF(NDSCH(4).EQ.1) CALL COEFF1(4)                           
        IF(NDSCH(4).EQ.2) CALL COEFF2(4)                           
        IF(NDSCH(4).EQ.3) CALL COEFF3(4)                           
        IF(NDSCH(4).EQ.4) CALL COEFF4(4)                           
        IF(NDSCH(4).EQ.5) CALL COEFF5(4)                           
                           
        CALL SORCE4                                                     
                                    
        CALL BOUNDS(4)                                                  
 
        IF(SPECL(4)) CALL SPEBC4                                        

        CALL SOLVER(ISWR,4)

        DO 150 II=1,NIJ                                                 
  150   F(ISW+II)=F(ISWR+II)/(F(ISR+II)+SMALL)                          

      ENDIF                                                             

C
C SOLVING FOR K-EQUATION
C
      IF(SOLVE(5)) THEN
 
        IF(NDSCH(5).EQ.1) CALL COEFF1(5)                           
        IF(NDSCH(5).EQ.2) CALL COEFF2(5)                           
        IF(NDSCH(5).EQ.3) CALL COEFF3(5)                           
        IF(NDSCH(5).EQ.4) CALL COEFF4(5)                           
        IF(NDSCH(5).EQ.5) CALL COEFF5(5)                           


        CALL SORCE5 
 
        CALL BOUNDS(5) 

        IF(SPECL(5)) CALL SPEBC5   
    
        CALL SOLVER(ISTE,5)


C        CHECK FOR   ZERO OR NEGATIVE K VALUES

        DO J=jsolv_b,jsolv_e
        DO I=isolv_b,isolv_e
           II = I+(J-1)*NI
         TEIJ = F (ISTE + II)
         IF(TEIJ.LT.0)   THEN
         WRITE(77,*)ntstep,iblk,ibswp,I,J,TEIJ
	  go to 777
         END IF 
	end do
	end do
777	continue
        END IF       ! IF CHECK FOR SOLVING K-EQUATION
C       
C
C SOLVING FOR EPSILON-EQUATION
C
      IF(SOLVE(6)) THEN 
C       
                     
        IF(NDSCH(6).EQ.1) CALL COEFF1(6)                           
        IF(NDSCH(6).EQ.2) CALL COEFF2(6)                           
        IF(NDSCH(6).EQ.3) CALL COEFF3(6)                           
        IF(NDSCH(6).EQ.4) CALL COEFF4(6)                           
        IF(NDSCH(6).EQ.5) CALL COEFF5(6)                           
                           
        CALL SORCE6                                                     
                                    
        CALL BOUNDS(6)                                                  
	
        IF(SPECL(6)) CALL SPEBC6    

        CALL SOLVER(ISED,6)

c       CHECK FOR   ZERO OR NEGATIVE EPSILON VALUES

        DO J=jsolv_b,jsolv_e
        DO I=isolv_b,isolv_e
           II = I+(J-1)*NI
         EDIJ = F (ISED + II)
         IF(EDIJ.LT.0)   THEN
         WRITE(78,*)ntstep,iblk,ibswp,I,J,EDIJ
	 go to 778
         END IF
         ENDDO
         ENDDO
778	 continue

        END IF       ! IF CHECK FOR SOLVING EPSILON-EQUATION
          
C       
C       
C SOLVING OMEGA-EQUATION 
C 
      IF(SOLVE(7)) THEN 

        IF(NDSCH(7).EQ.1) CALL COEFF1(7)                           
        IF(NDSCH(7).EQ.2) CALL COEFF2(7)                           
        IF(NDSCH(7).EQ.3) CALL COEFF3(7)                           
        IF(NDSCH(7).EQ.4) CALL COEFF4(7)                           
        IF(NDSCH(7).EQ.5) CALL COEFF5(7)                           
                           
C       CALL BNDPH(F(ISOM+1))

        CALL SORCE7                                                     
                                    
        CALL BOUNDS(7)                                                  
	
        IF(SPECL(7)) CALL SPEBC7                                        

        CALL SOLVER(ISOM,7)

c       CHECK FOR   ZERO OR NEGATIVE OMEGA VALUES

        DO J=jsolv_b,jsolv_e
        DO I=isolv_b,isolv_e
           II = I+(J-1)*NI
         AOMIJ = F (ISOM + II)
         IF(AOMIJ.LT.0)   THEN
         WRITE(88,*)iblk,ibswp,I,J,AOMIJ
          goto 779
         write(88,*)f(isaw+ii),f(isae+ii),f(isas+ii),f(isan+ii)
         write(88,*)f(issu+ii),f(issp+ii)
         write(6,*)'PROGRAM TERMINATED DUE TO NEGATIVE OMEGA-VALUE'
         stop
         END IF
         ENDDO
         ENDDO
C
779	 continue
        END IF       ! IF CHECK FOR SOLVING OMEGA-EQUATION
C
C       
C SOLVING  V-SQUARE EQUATION 
C 
      IF(SOLVE(8)) THEN 

        IF(NDSCH(8).EQ.1) CALL COEFF1(8)                           
        IF(NDSCH(8).EQ.2) CALL COEFF2(8)                           
        IF(NDSCH(8).EQ.3) CALL COEFF3(8)                           
        IF(NDSCH(8).EQ.4) CALL COEFF4(8)                           
        IF(NDSCH(8).EQ.5) CALL COEFF5(8)                           
                           
        CALL SORCE8                                                     
                                    
        CALL BOUNDS(8)                                                  
	
        IF(SPECL(8)) CALL SPEBC8                                        

        CALL SOLVER(ISV2,8)

        END IF       ! IF CHECK FOR SOLVING V-SQUARE-EQUATION
C       
C SOLVING F EQUATION 
C 
      IF(SOLVE(9)) THEN 

        CALL COEFF_F
                           
        CALL SORCE9                                                     
                                    
        CALL BOUNDS(9)                                                  
	
        IF(SPECL(9)) CALL SPEBC9                                        

        CALL SOLVER(ISF,9)
C
        END IF       ! IF CHECK FOR SOLVING F-EQUATION
C
C *** SOLVING NUTILDE EQUATION FOR SPALART ALMARAS MODEL(NN=10)

C
      IF(SOLVE(10)) THEN
C
       IF(NDSCH(10).EQ.1) CALL COEFF1(10)
       IF(NDSCH(10).EQ.2) CALL COEFF2(10)
       IF(NDSCH(10).EQ.3) CALL COEFF3(10)
       IF(NDSCH(10).EQ.4) CALL COEFF4(10)
       IF(NDSCH(10).EQ.5) CALL COEFF5(10)
C
        CALL SORCE10
C
        CALL BOUNDS(10)

C
        IF(SPECL(10)) CALL SPEBC10

        CALL SOLVER(ISNUT,10)

        END IF       ! IF CHECK FOR SOLVING NUTILDE EQUATION
C
C   COMPUTE EDDY VISCOSITY
C   WHEN F(ISNUT+II) IS USED FOR NUTILDE IN SA-TURBULENCE MODEL
C

         CALL EDYVIS
C
C
C****************************************************************
C
C *** SOLVING SCALAR EQUATION (NN=11)                                    
C
      IF(SOLVE(11)) THEN                                                  
C 
C     RELEASE THE BLOCK FOR TEMPERATURE TO ALLOW DIFFUSION OF HEAT 
C
        KBLK = 1
       IF(NDSCH(11).EQ.1) CALL COEFF1(11)                            
       IF(NDSCH(11).EQ.2) CALL COEFF2(11)                            
       IF(NDSCH(11).EQ.3) CALL COEFF3(11)                            
       IF(NDSCH(11).EQ.4) CALL COEFF4(11)                            
       IF(NDSCH(11).EQ.5) CALL COEFF5(11)                            
C                             
        CALL SORCE11                                                     
                                    
        CALL BOUNDS(11)                                                  
 
        IF(SPECL(11)) CALL SPEBC11                                        
 
        CALL SOLVER(ISS,11)

        END IF       ! IF CHECK FOR SOLVING SCALAR-EQUATION
C
	   CALL BLKRES
C
C
          WRITE(17,*) iblk,IBSWP,(RMOM(I),I=1,11),F(ISP+IIMON)
           
C
        CALL PRINT(F(ISAW+1),'AW FOR T')                            
        CALL PRINT(F(ISAS+1),'AS FOR T')                            
                                                                 
                                                                 
      RETURN                                                            

    3  FORMAT(/1X,'IBSWP' 
     >,'---------FIELD VALUES AT MONITORING '        
     >,'LOCATION(',I3,',',I3,')----',/9X,'P',8X,'U',8X,'V',6X,'W',  
     >6X,'SCAL'/)
                
    4  FORMAT(/1X,'IBSWP' 
     >,'---------FIELD VALUES AT MONITORING '        
     >,'LOCATION(',I3,',',I3,')----',/9X,'K',6X,'EPS',8X,'OM',6X,'V2',  
     >6X,'F',6X,'NUT'/)                

    2  FORMAT(1X,I4,1X,1P6E8.1)                               
                             
      END
C                                                               
C***********************************************************************
C
      SUBROUTINE GEOINP
C                                  
C***********************************************************************
C
c       READING X and Y COORDINATES FROM
C       THE GRID GENERATION OUTPUT FILE
C

      INCLUDE 'com2d'                                                  
C
      READ(7,rec=iblk) nim,njm,((cox(I,J),I=2,ni),j=2,nj),
     >                         ((coy(I,J),I=2,NI),j=2,nj)


      DO 20 J=1,NJ                                                      
      COX(1,J)=COX(2,J)                                                 
   20 COY(1,J)=COY(2,J)                                                 
C 
      DO 30 I=1,NI                                                      
      COX(I,1)=COX(I,2)                                                 
   30 COY(I,1)=COY(I,2)                                                 

      DO J=1,NJ
      DO I=1,NI
        II=I+(J-1)*NI
      ISKP=1
      JSKP=1
      IF(I.EQ.NI) ISKP=0
      IF(J.EQ.NJ) JSKP=0
C
C     CALCULATE CENTER COORDINATE OF CONTROL-VOLUMES
C
      XCENT(II)=0.25*(cox(i,j)+cox(i,j+jskp)+cox(i+iskp,j+jskp)+cox(i+iskp,j))
      yCENT(II)=0.25*(coy(i,j)+coy(i,j+jskp)+coy(i+iskp,j+jskp)+coy(i+iskp,j))


        END DO
        END DO

     
        do j=1,nj
         do i=1,ni
           ii=i+(j-1)*ni
           f(iscox+ii)= cox(i,j)
           f(iscoy+ii)= coy(i,j)
         end do
        end do
   

	return

      END 
	
C                                                               
C                                                                        
C***********************************************************************
C
      SUBROUTINE GEOMET                                                 
C
C     CALCULATE ALL GEOMETRICAL QUANTITIES                                  
C     REQUIRED FOR THE CONTROL VOLUME
C
C***********************************************************************
C
      INCLUDE 'com2d'                                    
C
      SQR(VALUE)=VALUE*VALUE                                            
      VALUE=0.

      DO 1000 J=1,NJ                                                    
      DO 1000 I=1,NI                                                    
        II=I+(J-1)*NI                                                   
      ISKP=1                                                            
      JSKP=NI                                                           
      IF(I.EQ.NI) ISKP=0                                                
      IF(J.EQ.NJ) JSKP=0                                                
                                                                        
      XSW=F(ISCOX+II)                                                   
      XSE=F(ISCOX+II+ISKP)                                              
      XNW=F(ISCOX+II+JSKP)                                              
      XNE=F(ISCOX+II+JSKP+ISKP)                                         
                                                                        
      YSW=F(ISCOY+II)                                                   
      YSE=F(ISCOY+II+ISKP)                                              
      YNW=F(ISCOY+II+JSKP)                                              
      YNE=F(ISCOY+II+JSKP+ISKP)                                         
                                                                        
      XW=0.5*(XNW+XSW)                                                  
      XE=0.5*(XNE+XSE)                                                  
      XS=0.5*(XSW+XSE)                                                  
      XN=0.5*(XNW+XNE)                                                  
                                                                        
      YW=0.5*(YNW+YSW)                                                  
      YE=0.5*(YNE+YSE)                                                  
      YS=0.5*(YSW+YSE)                                                  
      YN=0.5*(YNW+YNE)                                                  
C                                                                        
C     CALCULATE CENTER COORDINATE OF CONTROL-VOLUMES                        
C
      XP=CENTER(ISCOX,II)                                               
      YP=CENTER(ISCOY,II)                                               
      F(ISX+II)=XP                                                      
      F(ISY+II)=YP                                                      
      F(ISR+II)=1.                                                      
      IF(AXISYM_X) F(ISR+II)=YP
      IF(AXISYM_Y) F(ISR+II)=XP
                                                                        
      IF(I.GT.1) THEN                                                   
            ISKP=1                                                      
         XIM=CENTER(ISCOX,II-1)                                         
         YIM=CENTER(ISCOY,II-1)                                         
            IF(I.EQ.NI) ISKP=0                                          
      ELSE                                                              
         XIM=XP                                                         
         YIM=YP                                                         
      ENDIF                                                             
                                                                        
      IF(J.GT.1) THEN                                                   
            JSKP=NI                                                     
        XJM=CENTER(ISCOX,II-NI)                                         
        YJM=CENTER(ISCOY,II-NI)                                         
            IF(J.EQ.NJ) JSKP=0                                          
      ELSE                                                              
        XJM=XP                                                          
        YJM=YP                                                          
      ENDIF                                                             
C                                                                        
C     COMPUTE INTERPOLATION FACTORS                                             
C
        SMH=SQRT(SQR(XW-XIM)+SQR(YW-YIM))                               
      F(ISFX+II)=SMH/(SQRT(SQR(XP-XW)+SQR(YP-YW))+SMH+SMALL)            
        SMH=SQRT(SQR(XS-XJM)+SQR(YS-YJM))                               
      F(ISFY+II)=SMH/(SQRT(SQR(XP-XS)+SQR(YP-YS))+SMH+SMALL)            
C                                                                        
C     PROJECTION AREA BIJ'S                                                               
C
      F(ISB11W+II)=YNW-YSW                                              
      F(ISB12W+II)=XSW-XNW                                              
      F(ISB21W+II)=YIM-YP                                               
      F(ISB22W+II)=XP -XIM                                              
                                                                        
      F(ISB11S+II)=YP -YJM                                              
      F(ISB12S+II)=XJM-XP                                               
      F(ISB21S+II)=YSW-YSE                                              
      F(ISB22S+II)=XSE-XSW                                              
                                                                        
      F(ISB11P+II)=YN -YS                                               
      F(ISB12P+II)=XS -XN                                               
      F(ISB21P+II)=YW -YE                                               
      F(ISB22P+II)=XE -XW                                               
C                                                                        
C     VOLUME OF THE CONTROL VOLUMES                                                            
C
       F(ISVOLP+II)=0.5*F(ISR+II)*ABS((XNE-XSW)*(YNW-YSE)-(XNW-XSE)*(YNE-YSW))

C      write(196,*)J,I,F(ISX+II),F(ISY+II)
        
1000   CONTINUE                                                          

       RETURN
       END  
C
C***********************************************************************
C                                                             
      REAL FUNCTION CENTER(IST,II)                                      
C
C
      INCLUDE 'com2d'                                           
C
      CENTER=0.25*(F(IST+II)+F(IST+II+JSKP)                             
     >            +F(IST+II+JSKP+ISKP)+F(IST+II+ISKP))                  
       RETURN
       END                                                               
C                      
C***********************************************************************
C
      SUBROUTINE INFLUX 
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

      IF(MINLW.EQ.0) GO TO 419                                          
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
C                                               
  419 IF(MINLE.EQ.0) GO TO 429                                          
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
C                                              
  429 IF(MINLS.EQ.0) GO TO 439                                          
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
C                                               
  439 IF(MINLN.EQ.0) GOTO 460                                           
        DO  IM=1,MINLN                                               
              IP=INLEN(IM)                                           
              IB=IP+NI                                           
              CN=-(F(ISCS+IB))*mysolv(ip)
        RNORM(1)=RNORM(1)+CN                                            
        RNORM(2)=RNORM(2)+CN*SQRT(F(ISU+IB)**2+F(ISV+IB)**2)            
        RNORM(4)=RNORM(4)+CN*ABS(F(ISWR+IB))                            
        RNORM(5)=RNORM(5)+CN*ABS(F(ISTE+IB))                            
        RNORM(6)=RNORM(6)+CN*ABS(F(ISED+IB))                            
        RNORM(7)=RNORM(7)+CN*ABS(F(ISOM +IB))    
        RNORM(8)=RNORM(8)+CN*ABS(F(ISV2 +IB))    
        RNORM(9)=RNORM(9)+CN*ABS(F(ISF +IB))    
        RNORM(10)=RNORM(10)+CN*ABS(F(ISNUT +IB))    
        RNORM(11)=RNORM(11)+CN*ABS(F(ISS +IB))    
        END DO
        RNORM(3)=RNORM(2)
460	continue

	RETURN
      END                                                               
C                                                                        
C***********************************************************************
C
      SUBROUTINE INFORM                                                 
C
C     PRINT IMPORTANT INFORMATION ON THE CALCULATION                    
C
C***********************************************************************
C
      INCLUDE 'com2d'                                           
C
      WRITE(6,*) ' '                                                      
      WRITE(6,*) ' '                                                    
      WRITE(6,*) ' TOTAL NUMBER OF BLOCKS = ',NBLOCK          
      WRITE(6,*) ' '                                                    
      WRITE(6,*) ' '                                                    
      WRITE(6,*) ' TIME STEP =  ',NTSTEP,'AND T=',tinst
      WRITE(6,*) ' '                                                    
      WRITE(6,*) ' '                                                    
      WRITE(6,*) ' INFORMATION FOR BLOCK NUMBER =  ',IBLK,'AND SWEEP=',IBSWP
      WRITE(6,*) ' '                                                    

      IF(OVACOR) WRITE(6,*)                                             
     >  ' OVERALL CONTINUITY CORRECTION IS APPLIED AT EXIT'             
      IF(.NOT.OVACOR) WRITE(6,*)                                        
     >  ' NO OVERALL CONTINUITY CORRECTION AT EXIT'                     

      IF(.NOT.AXISYM_Y) WRITE(6,20)                            
   20 FORMAT(/2X,'FLOW IS STEADY TWO-DIMENSIONAL PLANE')                
      IF(AXISYM_Y) WRITE(6,24)                                 
   24 FORMAT(/2X,'FLOW IS STEADY AXISYMMETRIC ABOUT X-AXIS')
       IF(LINKPV.EQ.1) WRITE(6,22)                            
   22 FORMAT(/2X,'SIMPLE STRATEGY USED FOR PRESSURE VELOCITY LINK')                
       IF(LINKPV.EQ.2) WRITE(6,23)                                 
   23 FORMAT(/2X,'SIMPLEC STRATEGY USED FOR PRESSURE VELOCITY LINK')

      WRITE(6,30) (RELAX(I),I=1,10)   
   30 FORMAT(/2X,'RELAXATION FACTORS: ',F5.2,'(PP),',F5.2,'(U),',F5.2,
     >       '(V),',F5.2,'(WR),',F5.2,'(K),',F5.2,'(ED),',F5.2,'(OM),',
     >       F5.2,'(V2),',F5.2,'(F),',F5.2,'(S)')
      WRITE(6,301) (RELAX(I),I=11,13)   
  301 FORMAT(/2X,'RELAXATION FACTORS: ',F5.2,'(P),',F5.2,'(DEN),',F5.2,
     >      '(VIS)')
      WRITE(6,32) (RNORM(I),I=1,6)                                      
   32 FORMAT(/2X,'INLET FLUXES: ',1PE10.3,'(MASS),',1PE10.3,'(U),',     
     >       1PE10.3,'(V),',1PE10.3,'(WR),',1PE10.3,'(K),',1PE10.3,     
     >       '(ED)')
      WRITE(6,321) (RNORM(I),I=7,10)                                      
  321 FORMAT(/2X,'INLET FLUXES: ',1PE10.3,'(OM),',1PE10.3,'(V2),',     
     >       1PE10.3,'(F),',1PE10.3,'(S)')     
                           
      WRITE(6,60) (NDSCH(I),I=2,10)                                      
   60 FORMAT(/2X,'CONVECTION SCHEMES ',                                 
     > '(1-HYBRID, 2-QUICK, 3-HLPA, 4-CENTRAL DIFFERENCE): ',              
     >     I2,'(U),',I2,'(V),',I2,'(WR),',I2,'(K),',I2,'(ED),',I2,'(OM),',
     >     I2,'(V2),',I2,'(F),',I2,'(S)')
       WRITE(6,61) (NSOLV(I),I=1,10)                                      
   61 FORMAT(/2X,'LINEAR SOLVER ',                                 
     > '(1- ALTERNATE DIRECTION TDMA, 2-SIP OF STONE): ',I2,'(PP),',   
     >     I2,'(U),',I2,'(V),',I2,'(WR),',I2,'(K),',I2,'(ED),',I2,'(OM),',
     >     I2,'(V2),',I2'(F),',I2,'(S)')

      IP=0                                                              
      DO 70 I=2,10                                                       
   70 IP=IP+NDSCH(I)                                                    
      IF(IP.GE.10) WRITE(6,72) (FDEFER(I),I=2,10)                         
   72 FORMAT(/2X,'BLENDING FACTORS: ',F5.2,'(U)',F5.2,'(V)',F5.2,'(WR)' 
     >       , F5.2,'(K)',F5.2,'(ED)',F5.2,'(OM)',F5.2,'(V2)',F5.2,'(F)',
     >          F5.2,'(S)')                       

	IF(TURBMOD.GT.0) THEN
	 write(6,*)
	 write(6,*)'THE FLOW IS TURBULENT'
	write(6,*)
	go to (1,2,3,4,5,6,7) turbmod
1	continue
	   write(6,*)'        STANDARD K-EPSILON MODEL USED'
	go to 99
C
2	continue
	   write(6,*)'        RODI TWO LAYER MODEL USED'
	go to 99
C
3	continue
	   write(6,*)'         CHIEN LOW RE VERSION OF K-EPSILON MODEL USED'
	go to 99
C
4	continue
	   write(6,*)'        K- OMEGA  MODEL USED'
	go to 99
C
5	continue
	   write(6,*)'         K-EPSILON V2-F  MODEL USED'
	go to 99
C
6	continue
	   write(6,*)'         SPALART-ALAMARAS MODEL USED'
	go to 99
C
7	continue
	   write(6,*)'         SHEAR STRESS TRANSPORT MODEL  USED'
	go to 99
C
99      continue
	else
	 write(6,*)
	 write(6,*)'THE FLOW IS LAMINAR' 
	end if
	write(6,*)
C
C     INFORMATION FOR DEBUGGING
C
        IF (CONVECTIVE_BC)THEN
           write(6,*)
           write(6,*)'CONVECTIVE BOUNDARY CONDITION APPLIED AT THE EXIT'
        ELSE
           write(6,*)
           write(6,*)'CONVECTIVE BOUNDARY CONDITION NOT APPLIED AT THE EXIT'
        END IF
           write(6,*)
           write(6,*)'FLOW ANALYSIS USING FIXED GRID'

C
C     INFORMATION FOR DEBUGGING                                         
C
      IF(IDEBUG.NE.1) RETURN                                            
      WRITE(6,*) ' '                                                    
      WRITE(6,*) ' INFORMATION FOR DEBUGGING WHEN SETTING IDEBUG=1'     
     >          ,' IN SUBROUTINE TITLE'                                      
      WRITE(6,*) ' ***********************************************'     
     >          ,'***************'                                      
      WRITE(6,*) ' '                                                    
      WRITE(6,*) ' NUMBER OF WEST BOUNDARIES : ',INFDE( 1)              
      WRITE(6,*) ' NUMBER OF EAST BOUNDARIES : ',INFDE( 2)              
      WRITE(6,*) ' NUMBER OF SOUTH BOUNDARIES: ',INFDE( 3)              
      WRITE(6,*) ' NUMBER OF NORTH BOUNDARIES: ',INFDE( 4)              
      WRITE(6,*) ' NUMBER OF BLOCKED REGIONS : ',INFDE( 5)              
      WRITE(6,*) ' '                                                    
      WRITE(6,*) ' NUMBER OF WALL BOUNDARIES          : ',INFDE( 6)     
      WRITE(6,*) ' NUMBER OF SYMMETRY BOUNDARIES      : ',INFDE( 7)     
      WRITE(6,*) ' NUMBER OF OUTFLOW BOUNDARIES       : ',INFDE( 8)     
      WRITE(6,*) ' NUMBER OF INFLOW BOUNDARIES        : ',INFDE( 9)     
      WRITE(6,*) ' NUMBER OF CUT BOUNDARIES           : ',INFDE(10)
      WRITE(6,*) ' NUMBER OF PERIODIC BOUNDARIES      : ',INFDE(12)
      WRITE(6,*) ' NUMBER OF NODES IN BLOCKED REGIONS : ',INFDE(11)     
      WRITE(6,*) ' '                                                    
      WRITE(6,*) ' TOTAL NUMBER OF BOUNDARY-NODES USED FOR'             
      WRITE(6,*) '   WEST  WALL     (MWALW): ',MWALW                    
      WRITE(6,*) '   EAST  WALL     (MWALE): ',MWALE                    
      WRITE(6,*) '   SOUTH WALL     (MWALS): ',MWALS                    
      WRITE(6,*) '   NORTH WALL     (MWALN): ',MWALN                    
      WRITE(6,*) '   WEST  SYMMETRY (MSYMW): ',MSYMW                    
      WRITE(6,*) '   EAST  SYMMETRY (MSYME): ',MSYME                    
      WRITE(6,*) '   SOUTH SYMMETRY (MSYMS): ',MSYMS                    
      WRITE(6,*) '   NORTH SYMMETRY (MSYMN): ',MSYMN                    
      WRITE(6,*) '   WEST  OUTFLOW  (MEXIW): ',MEXIW                    
      WRITE(6,*) '   EAST  OUTFLOW  (MEXIE): ',MEXIE                    
      WRITE(6,*) '   SOUTH OUTFLOW  (MEXIS): ',MEXIS                    
      WRITE(6,*) '   NORTH OUTFLOW  (MEXIN): ',MEXIN                    
      WRITE(6,*) '   WEST  INFLOW   (MINLW): ',MINLW                    
      WRITE(6,*) '   EAST  INFLOW   (MINLE): ',MINLE                    
      WRITE(6,*) '   SOUTH INFLOW   (MINLS): ',MINLS                    
      WRITE(6,*) '   NORTH INFLOW   (MINLN): ',MINLN                    
      WRITE(6,*) '   WEST  CUT      (MCUTW): ',MCUTW                    
      WRITE(6,*) '   EAST  CUT      (MCUTE): ',MCUTE                    
      WRITE(6,*) '   SOUTH CUT      (MCUTS): ',MCUTS                    
      WRITE(6,*) '   NORTH CUT      (MCUTN): ',MCUTN
      WRITE(6,*) '   WEST  PERIODIC (MCYCW): ',MCYCW               
      WRITE(6,*) '   EAST  PERIODIC (MCYCE): ',MCYCE               
      WRITE(6,*) '   SOUTH PERIODIC (MCYCS): ',MCYCS
      WRITE(6,*) '   NORTH PERIODIC (MCYCN): ',MCYCN                     
      WRITE(6,*) '   WEST  PRESSURE BNDRY (MPCTW): ',MPCTW
      WRITE(6,*) '   EAST  PRESSURE BNDRY (MPCTE): ',MPCTE
      WRITE(6,*) '   SOUTH PRESSURE BNDRY (MPCTS): ',MPCTS
      WRITE(6,*) '   NORTH PRESSURE BNDRY (MPCTN): ',MPCTN
      WRITE(6,*) ' ***********************************************'     
     >          ,'***************'                                      
 
      RETURN                                                            
      END      
C
C***********************************************************************
C
      SUBROUTINE ISTADR                                                 
C
C     CALCULATE STARTING ADDRESSES FOR DIFFERENT
C     VARIABLES STORED IN THE ARRAY F                                      
C
C***********************************************************************
C
      INCLUDE 'com2d'    
C                                             
	  DATA NVCOX,NVCOY,NVX,NVY,NVR,NVFX,NVFY,NVVOLP                
     >,NVPP,NVU,NVV,NVWR,NVTE,NVED,NVOM,NVV2,NVF,NVNUT,NVS,NVP,NVDEN,NVVIS,
     >NVVISW,NVVISS,NVW,NVCW,NVCS,
     >NVB11W,NVB12W,NVB21W,NVB22W,
     >NVB11S,NVB12S,NVB21S,NVB22S,NVB11P,NVB12P,NVB21P,NVB22P,
     >NVAW,NVAE,NVAS,NVAN,NVSU,NVSP,NVAP,NVAPU,NVAPV,NVSUMU,NVSUMV 
     >/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     >23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,
     >43,44,45,46,47,48,49,50/
C

      IF((Ni.GT.NXF.OR.Nj.GT.NYF).OR.NMAX.LT.MAX(NXF,NYF))
     >        CALL CHECK1(1)
      NIM=NI-1                                                          
      NJM=NJ-1                                                          
      NIJ=NI*NJ                                                         
      IIS=2+NI                                                          
      IIE=NIJ-NI-1                                                      

      NV(1)=1
C                                                           
      DO 25 IV=2,NVTOT                                                  
   25 NV(IV)=NV(IV-1)+NIJ                                               
                                                                        
      IAUW=0                                                            
      IAUE=IAUW+6*NIJ                                                   
      IAUS=IAUE+6*NIJ                                                   
      IAUN=IAUS+6*NIJ                                                   

      ISCOX =NV(  1)-1                                                  
      ISCOY =NV(  2)-1                                                  
      ISX   =NV(  3)-1                                                  
      ISY   =NV(  4)-1                                                  
      ISR   =NV(  5)-1                                                  
      ISFX  =NV(  6)-1                                                  
      ISFY  =NV(  7)-1                                                  
      ISVOLP=NV(  8)-1                                                  
      ISPP  =NV(  9)-1                                                  
      ISU   =NV( 10)-1                                                  
      ISV   =NV( 11)-1                                                  
      ISWR  =NV( 12)-1                                                  
      ISTE  =NV( 13)-1                                                  
      ISED  =NV( 14)-1    
      ISOM  =NV( 15)-1                                              
      ISV2  =NV( 16)-1                                              
      ISF   =NV( 17)-1                                              
      ISNUT =NV( 18)-1 
      ISS   =NV( 19)-1                                                  
      ISP   =NV( 20)-1                                                  
      ISDEN =NV( 21)-1                                                  
      ISVIS =NV( 22)-1                                                  
      ISVISW=NV( 23)-1                                                  
      ISVISS=NV( 24)-1                                                  
      ISW   =NV( 25)-1                                                  
      ISCW  =NV( 26)-1                                                  
      ISCS  =NV( 27)-1                                                  
      ISB11W=NV( 28)-1                                                  
      ISB12W=NV( 29)-1                                                  
      ISB21W=NV( 30)-1                                                  
      ISB22W=NV( 31)-1                                                  
      ISB11S=NV( 32)-1                                                  
      ISB12S=NV( 33)-1                                                  
      ISB21S=NV( 34)-1                                                  
      ISB22S=NV( 35)-1                                                  
      ISB11P=NV( 36)-1                                                  
      ISB12P=NV( 37)-1                                                  
      ISB21P=NV( 38)-1                                                  
      ISB22P=NV( 39)-1                                                  
      ISAW  =NV( 40)-1                                                  
      ISAE  =NV( 41)-1                                                  
      ISAS  =NV( 42)-1                                                  
      ISAN  =NV( 43)-1                                                  
      ISSU  =NV( 44)-1                                                  
      ISSP  =NV( 45)-1 
      ISAP  =NV( 46)-1                                                 
      ISAPU =NV( 47)-1                                                  
      ISAPV =NV( 48)-1
      ISSUMU=NV( 49)-1
      ISSUMV=NV( 50)-1

      END
C                                                               
C***********************************************************************
C
      SUBROUTINE ADSOR(N)                                               
C
C     CALCULATE ADDITIONAL SOURCE TERMS DUE TO PRESSURE
C     DIFFERENCE AND OTHER DIFFUSIVE TYPE STRESS TERMS
C     IN THE MOMENTUM EQUATIONS                
C
C***********************************************************************
C
      INCLUDE 'com2d'  
C                                            
      DO J=jsolv_b,jsolv_e
      DO I=isolv_b,isolv_e
      II=I+(J-1)*NI
C                                
C     INTERPOLATION-FACTORS
C                                             
      FXP   =F(ISFX+II)                                                 
      FXIP  =F(ISFX+II+1)                                               
      FXJP  =F(ISFX+II+NI)                                              
      FXJM  =F(ISFX+II-NI)                                              
      FXIPJP=F(ISFX+II+1+NI)                                            
      FXIPJM=F(ISFX+II+1-NI)                                            
                                                                        
      FYP   =F(ISFY+II)                                                 
      FYIP  =F(ISFY+II+1)                                               
      FYIM  =F(ISFY+II-1)                                               
      FYJP  =F(ISFY+II+NI)                                              
      FYIPJP=F(ISFY+II+1+NI)                                            
      FYIMJP=F(ISFY+II-1+NI)                                            
C                                                                        
C     NEIGHBOUR-UALUES OF U-VELOCITY  
C                                  
      UIJK =F(ISU+II)                                                   
      UIP  =F(ISU+II+1)                                                 
      UIPJP=F(ISU+II+1+NI)                                              
      UIPJM=F(ISU+II+1-NI)                                              
      UIM  =F(ISU+II-1)                                                 
      UIMJP=F(ISU+II-1+NI)                                              
      UIMJM=F(ISU+II-1-NI)                                              
      UJP  =F(ISU+II+NI)                                                
      UJM  =F(ISU+II-NI)                                                
C                                                                        
C     U-VELOCITY DIFFERENCES AROUND THE CONTROL VOLUME
C                  
      DUDXE=UIP-UIJK                                                    
      DUDXW=UIJK-UIM                                                    
      DUDYN=UJP-UIJK                                                    
      DUDYS=UIJK-UJM                                                    
                                                                        
      DUNS=FYJP*UJP+(1.-FYJP)*UIJK-(FYP*UIJK+(1.-FYP)*UJM)              
      DUEW=FXIP*UIP+(1.-FXIP)*UIJK-(FXP*UIJK+(1.-FXP)*UIM)              
C                                                                        
C     ALONG NORTH-SOUTH FOR EAST/WEST
C                                   
      DUDYE=(1.-FXIP)*DUNS+    FXIP *((FYIPJP*UIPJP+(1.-FYIPJP)*UIP)    
     >                               -(FYIP  *UIP  +(1.-FYIP  )*UIPJM)) 
      DUDYW=    FXP  *DUNS+(1.-FXP )*((FYIMJP*UIMJP+(1.-FYIMJP)*UIM)    
     >                               -(FYIM  *UIM  +(1.-FYIM  )*UIMJM)) 
C                                                                        
C     ALONG EAST-WEST FOR NORTH/SOUTH
C                                   
      DUDXN=(1.-FYJP)*DUEW+    FYJP *((FXIPJP*UIPJP+(1.-FXIPJP)*UJP)    
     >                               -(FXJP  *UJP  +(1.-FXJP  )*UIMJP)) 
      DUDXS=    FYP  *DUEW+(1.-FYP )*((FXIPJM*UIPJM+(1.-FXIPJM)*UJM)    
     >                               -(FXJM  *UJM  +(1.-FXJM  )*UIMJM)) 
C                                                                        
C     NEIGHBOUR-VALUES OF V-VELOCITY 
C                                   
      VIJK =F(ISV+II)                                                   
      VIP  =F(ISV+II+1)                                                 
      VIPJP=F(ISV+II+1+NI)                                              
      VIPJM=F(ISV+II+1-NI)                                              
      VIM  =F(ISV+II-1)                                                 
      VIMJP=F(ISV+II-1+NI)                                              
      VIMJM=F(ISV+II-1-NI)                                              
      VJP  =F(ISV+II+NI)                                                
      VJM  =F(ISV+II-NI)                                                
C                                                                        
C     V-VELOCITY DIFFERENCES AROUND THE CONTROL VOLUME
C                  
      DVDXE=VIP-VIJK                                                    
      DVDXW=VIJK-VIM                                                    
      DVDYN=VJP-VIJK                                                    
      DVDYS=VIJK-VJM                                                    
                                                                        
      DVEW=FXIP*VIP+(1.-FXIP)*VIJK-(FXP*VIJK+(1.-FXP)*VIM)              
      DVNS=FYJP*VJP+(1.-FYJP)*VIJK-(FYP*VIJK+(1.-FYP)*VJM)              
C                                                                        
C     ALONG NORTH-SOUTH FOR EAST/WEST
C                                   
      DVDYE=(1.-FXIP)*DVNS+    FXIP *((FYIPJP*VIPJP+(1.-FYIPJP)*VIP)    
     >                               -(FYIP  *VIP  +(1.-FYIP  )*VIPJM)) 
      DVDYW=    FXP  *DVNS+(1.-FXP )*((FYIMJP*VIMJP+(1.-FYIMJP)*VIM)    
     >                               -(FYIM  *VIM  +(1.-FYIM  )*VIMJM)) 
C                                                                        
C     ALONG EAST-WEST FOR NORTH/SOUTH
C                                   
      DVDXN=(1.-FYJP)*DVEW+    FYJP *((FXIPJP*VIPJP+(1.-FXIPJP)*VJP)    
     >                               -(FXJP  *VJP  +(1.-FXJP  )*VIMJP)) 
      DVDXS=    FYP  *DVEW+(1.-FYP )*((FXIPJM*VIPJM+(1.-FXIPJM)*VJM)    
     >                               -(FXJM  *VJM  +(1.-FXJM  )*VIMJM)) 
C                                                                        
C     RECIPROCAL OF CELL VOLUMES AROUND DIFFERENT FACES
C                 
      RVE=1./(0.5*(F(ISVOLP+II)+F(ISVOLP+II+1 ))+SMALL)                 
      RVW=1./(0.5*(F(ISVOLP+II)+F(ISVOLP+II-1 ))+SMALL)                 
      RVN=1./(0.5*(F(ISVOLP+II)+F(ISVOLP+II+NI))+SMALL)                 
      RVS=1./(0.5*(F(ISVOLP+II)+F(ISVOLP+II-NI))+SMALL)                 
C                                                                        
C     RADII AT CELL FACES 
C                                              
      RE=FXIP*F(ISR+II+1 )+(1.-FXIP)*F(ISR+II)                          
      RW=FXP *F(ISR+II   )+(1.-FXP )*F(ISR+II-1)                        
      RN=FYJP*F(ISR+II+NI)+(1.-FYJP)*F(ISR+II)                          
      RS=FYP *F(ISR+II   )+(1.-FYP )*F(ISR+II-NI)                       
C                                                                        
C     DIFFUSIVE TYPE STRESS TERMS FOR MOMENTUM EQUATIONS
C                
       B1NW=F(ISB11W+II+(N-2)*NIJ)                                      
       B2NW=F(ISB21W+II+(N-2)*NIJ)                                      
      TERMW=(DUDXW*B1NW+DUDYW*B2NW)*F(ISB11W+II)                        
     >     +(DVDXW*B1NW+DVDYW*B2NW)*F(ISB12W+II)                        
       B1NE=F(ISB11W+II+1+(N-2)*NIJ)                                    
       B2NE=F(ISB21W+II+1+(N-2)*NIJ)                                    
      TERME=(DUDXE*B1NE+DUDYE*B2NE)*F(ISB11W+II+1)                      
     >     +(DVDXE*B1NE+DVDYE*B2NE)*F(ISB12W+II+1)                      
       B1NS=F(ISB11S+II+(N-2)*NIJ)                                      
       B2NS=F(ISB21S+II+(N-2)*NIJ)                                      
      TERMS=(DUDXS*B1NS+DUDYS*B2NS)*F(ISB21S+II)                        
     >     +(DVDXS*B1NS+DVDYS*B2NS)*F(ISB22S+II)                        
       B1NN=F(ISB11S+II+NI+(N-2)*NIJ)                                   
       B2NN=F(ISB21S+II+NI+(N-2)*NIJ)                                   
      TERMN=(DUDXN*B1NN+DUDYN*B2NN)*F(ISB21S+II+NI)                     
     >     +(DVDXN*B1NN+DVDYN*B2NN)*F(ISB22S+II+NI) 

      TERM1=TERME*F(ISVISW+II+1 )*RE*RE*RVE-TERMW*F(ISVISW+II)*RW*RW*RVW
     >     +TERMN*F(ISVISS+II+NI)*RN*RN*RVN-TERMS*F(ISVISS+II)*RS*RS*RVS
C                                                                        
C     PRESSURE GRADIENT FORCES AS SOURCE TERM OF THE MOMENTUM EQUATIONS
C  
      CALL DPRS(II,F(ISP+1),PW,PE,PS,PN)
C               
      TERM2=F(ISR+II)*(F(ISB11P+II+(N-2)*NIJ)*(PW-PE)+                  
     >                 F(ISB21P+II+(N-2)*NIJ)*(PS-PN))                  
C
C     ADD SOURCE TERMS DUE TO DIFFUSION AND PRESSURE GRADIENT
C
       SUIJ=(TERM1+TERM2)*KBLK(II)
C
      F(ISSU+II)=F(ISSU+II)+SUIJ 
	end do
	end do
C
      RETURN                              
      END
C                                                               
C***********************************************************************
C
      SUBROUTINE PRINT(PHI,HEAD)                                        
C
C     PRINT THE VALUES OF 2-DIMENSIONAL FIELD VARIABLE PHI(I,J)               
C
C***********************************************************************
C
      INCLUDE 'com2d'
C                                               
      DIMENSION PHI(Ni,nj),STORE(NMAX)                                  
      CHARACTER*10 HEAD
C                                                 
      JS=1                                                              
      WRITE(6,210) HEAD                                                 
      ISTA=-11                                                          
  100 ISTA=ISTA+12                                                      
      IEND=ISTA+11                                                      
      IEND=MIN0(NI,IEND)                                                
      WRITE(6,111) (I,I=ISTA,IEND)                                      
      DO 101 JJ=JS,NJ                                                   
        J=JS+NJ-JJ                                                      
        DO 120 I=ISTA,IEND                                              
          A=PHI(I,J)                                                    
          IF(ABS(A).LT.1.E-30) A=0.0                                    
  120   STORE(I)=A                                                      
  101 WRITE(6,113) J,(STORE(I),I=ISTA,IEND)                             
      IF(IEND.LT.NI) GO TO 100                                          
      RETURN                                                            
  111 FORMAT(/1X,6H      ,I3,11I10)                                     
  113 FORMAT(1X,I3,1P12E10.2)                                           
  210 FORMAT(/1X,27(2H*-),6X,A10,6X,27(2H-*))                           
      END      
C                                                         
C***********************************************************************
C
      SUBROUTINE  PVCORR                                                 
C
C     CORRECT VELOCITY, MASS FLUX AND PRESSURE FIELD
C                         
C***********************************************************************
C
      INCLUDE 'com2d'                                                
C
C     CORRECT VELOCITY AT CELL-CENTER AND FLUXES AT CELL-FACES          
C
	do j=2,njm
	do i=2,nim
C
      II=I+(J-1)*NI        

      CALL DPRS(II,F(ISPP+1),PPW,PPE,PPS,PPN)

      DPPWE=PPW-PPE                                                     
      DPPSN=PPS-PPN                                                     
      F(ISU+II)=F(ISU+II)+(F(ISB11P+II)*DPPWE+F(ISB21P+II)*DPPSN)       
     >         *F(ISR+II)*F(ISAPU+II)*KBLK(II)/(1.-F(ISSUMU+II))        
      F(ISV+II)=F(ISV+II)+(F(ISB12P+II)*DPPWE+F(ISB22P+II)*DPPSN)       
     >         *F(ISR+II)*F(ISAPV+II)*KBLK(II)/(1.-F(ISSUMV+II))        
       F(ISCW+II)=F(ISCW+II)+(F(ISPP+II-1 )-F(ISPP+II))*F(ISAW+II)*APP(II)
       F(ISCS+II)=F(ISCS+II)+(F(ISPP+II-NI)-F(ISPP+II))*F(ISAS+II)*APP(II)

	end do
	end do
C
C CORRECTION OF U & V AT PRESSURE BOUNDARY PLANE
C
     	IF(MPCTW.NE.0) then
           DO  IM=1,MPCTW
               IP=INPCW(IM)
               IB=IP-1

	       PPW=0.
	       PPE=f(ispp+ip)
               DPPWE=PPW-PPE 
            F(ISU+IB)=F(ISU+IB)+(F(ISB11P+IB)*DPPWE)       
     >                *F(ISR+IB)*F(ISAPU+IB)/(1.-F(ISSUMU+IB))        
            F(ISV+IB)=F(ISV+IB)+(F(ISB12P+IB)*DPPWE)       
     >               *F(ISR+IB)*F(ISAPV+IB)/(1.-F(ISSUMV+IB))        
           end do
      	end if

        IF(MPCTE.NE.0) then
           DO  IM=1,MPCTE
               IP=INPCE(IM)
               IB=IP+1

	       PPW=f(ispp+ip)
	       PPE=0
               DPPWE=PPW-PPE
            F(ISU+IB)=F(ISU+IB)+(F(ISB11P+IB)*DPPWE)
     >                *F(ISR+IB)*F(ISAPU+IB)/(1.-F(ISSUMU+IB))
            F(ISV+IB)=F(ISV+IB)+(F(ISB12P+IB)*DPPWE)
     >               *F(ISR+IB)*F(ISAPV+IB)/(1.-F(ISSUMV+IB))
           end do
        end if

        IF(MPCTS.NE.0) then
           DO  IM=1,MPCTS
               IP=INPCS(IM)
               IB=IP-NI

	       PPS=0.
	       PPN=f(ispp+ip)

               DPPSN=PPS-PPN
            F(ISU+IB)=F(ISU+IB)+(F(ISB21P+IB)*DPPSN)
     >                *F(ISR+IB)*F(ISAPU+IB)/(1.-F(ISSUMU+IB))
            F(ISV+IB)=F(ISV+IB)+(F(ISB22P+IB)*DPPSN)
     >               *F(ISR+IB)*F(ISAPV+IB)/(1.-F(ISSUMV+IB))
           end do
        end if

        IF(MPCTN.NE.0) then
           DO  IM=1,MPCTN
               IP=INPCN(IM)
               IB=IP+NI

	       PPS=f(ispp+ip)
	       PPN=0.

               DPPSN=PPS-PPN
            F(ISU+IB)=F(ISU+IB)+(F(ISB21P+IB)*DPPSN)
     >                *F(ISR+IB)*F(ISAPU+IB)/(1.-F(ISSUMU+IB))
            F(ISV+IB)=F(ISV+IB)+(F(ISB22P+IB)*DPPSN)
     >               *F(ISR+IB)*F(ISAPV+IB)/(1.-F(ISSUMV+IB))
           end do
        end if
C                                                                        
C     CORRECT PRESSURE AT CELL-CENTER                                   
C
      DO II=1,NIJ                                                   
      F(ISP+II)=F(ISP+II)+RELAX(NVP-NVPP+1)*F(ISPP+II)           
      F(ISPP+II)=0.                                                     
      END DO	

      RETURN                             
      END
C                                                               
C***********************************************************************
C
      SUBROUTINE RDBLOK
C
C     READ COMPUTED RESULTS FOR EACH BLOCK FROM A DIRECT ACCESS FILE
C
C***********************************************************************
C
      INCLUDE 'com2d'                                            
C
      IPb=9*NIJ+1
      IPe=27*NIJ

      READ(8,rec=iblk) IBSWP_NO,NTSTEP_NO,rnorm,delt_o,tinst_no,
     >                  (f(iphi),iphi=ipb,ipe) 
      return
      END
C
C                                                               
      SUBROUTINE WTBLOK
C                                                 
C     SAVE COMPUTED RESULTS FOR EACH BLOCK ON A DIRECT ACCESS FILE
C              
C***********************************************************************
C
      INCLUDE 'com2d'                                                
C
      IPb=9*nij+1
      IPe=27*NIJ
      WRITE(8,rec=iblk) IBSWP,NTSTEP,rnorm,delt,tinst+delt,
     >                  (f(iphi),iphi=ipb,ipe)
C
C
      return
      END
C                                                               
C***********************************************************************
C
      SUBROUTINE WTTIME
C                                                 
C     SAVE COMPUTED RESULTS FOR EACH BLOCK ON A DIRECT ACCESS FILE
C              
C***********************************************************************
C
      INCLUDE 'com2d'                                                
C                                                               
	ipb=9*nij+1
	ipe=20*nij
	
        ntdiff=ntstep-ntsb
        if (mod(ntdiff,ntsint).eq.0) then
	     ifact=ntdiff/ntsint
	     irnum=iblk+nblock*ifact
             WRITE(10,rec=irnum) IBSWP,NTSTEP,RNORM,delt,tinst+delt,
     >                           (f(ip),ip=ipb,ipe)
        END IF


	return
C
      END
C              
C***********************************************************************
C
      SUBROUTINE SORCE1
C
C     CALCULATION OF CELLWISE MASS IMBALANCE 
C     AS SOURCE FOR CONTINUITY EQUATION
C              
C***********************************************************************
C
      INCLUDE 'com2d'                                                
C
      DO J=jsolv_b,jsolv_e
      DO I=isolv_b,isolv_e
      II=I+(J-1)*NI
      SMP=(F(ISCW+II+1)-F(ISCW+II)+F(ISCS+II+NI)-F(ISCS+II))*KBLK(II)   
      SMP_PREV=(CW_OLD(II+1,iblk)-CW_OLD(II,iblk)
     >            +CS_OLD(II+NI,iblk)-CS_OLD(II,iblk))*KBLK(II)   
 	SMP_PREV=0

         IF( DELT .GT. 1) SMP_PREV=0. ! LARGE DELT (STEADY STATE SOLN)

       F(ISSU+II)=-(SMP-SMP_PREV)                                                   
      F(ISSP+II)=0.                                                     
      end do
      end do
      RETURN 
      END                                                           
C              
C***********************************************************************
C
      SUBROUTINE SORCE2                                                      
C
C     CALCULATION OF SOURCE TERMS FOR U-MOMENTUM EQUATION                
C
C************************************************************************
C
      INCLUDE 'com2d'                                                
C
C  CROSS DERIVATIVE DIFFUSIVE FLUXES
C
      CALL CDFLUX(ISU,2)  
C
C  ADDITIONAL DIFFUSION AND PRESSURE GRADIENT TERMS AS SOURCE
C                                              
      CALL ADSOR(2)                                                     
C
C  INCLUSION OF UNSTEADY TERMS
C
      DO J=jsolv_b,jsolv_e
      DO I=isolv_b,isolv_e
        ii=i+(j-1)*ni
        VOLP=F(isvolp+ii)
        term=volp*kblk(ii)/delt

        if (TEMPSCH .EQ. 1 .OR. NTSTEP .EQ.1) THEN
C
C FIRST ORDER TIME DISCRETISATION
C
         su_unst=term*u(ii,iblk,1)
         sp_unst=term
        else
C
C SECOND ORDER TIME DISCRETISATION
C
         su_unst=term*(2*u(ii,iblk,2)-0.5*u(ii,iblk,1))
         sp_unst=1.5*term
        end if

        F(ISSU+ii)=f(issu+ii)+su_unst
        F(ISSp+ii)=f(issp+ii)-sp_unst
        end do
        end do
  
       IF(AXISYM_Y) THEN  
      
C      EXTRA TERMS FOR AXISYMMETRIC SITUATION ( Y-AXIS IS AXIS OF SYMMETRY)
C
      FAKVW=0.1
C
      DO J=jsolv_b,jsolv_e
      DO I=isolv_b,isolv_e
	II=I+(J-1)*NI
      F(ISSP+II)=F(ISSP+II)-2.*F(ISVIS+II)*F(ISVOLP+II)
     >                    *KBLK(II)/ (F(ISR+II)*F(ISR+II)+SMALL)

      IF(SWIRL) THEN
        SUIJ=F(ISDEN+II)*F(ISW+II)*F(ISW+II)*F(ISVOLP+II)/
     >                                    (F(ISR+II)+SMALL)
        SUIJ=SUIJ*KBLK(II)
      F(ISSU+II)=F(ISSU+II)+SUIJ*(1.+FAKVW*F(ISU+II)/(F(ISW+II)+SMALL))
      F(ISSP+II)=F(ISSP+II)-FAKVW*SUIJ/(F(ISW+II)+SMALL)
      ENDIF
        end do
        end do

       ENDIF

      END 
C              
C***********************************************************************
C
      SUBROUTINE SORCE3                                                      
C
C     CALCULATION OF SOURCE TERMS FOR V-MOMENTUM EQUATION                 
C
C************************************************************************
C                                                
      INCLUDE 'com2d'                                                
C
C  CROSS DERIVATIVE DIFFUSIVE FLUXES
C
      CALL CDFLUX(ISV,3)                                                
C 
C  ADDITIONAL DIFFUSION AND PRESSURE GRADIENT TERMS AS SOURCE
C
      CALL ADSOR(3)                                                     
C
       IF(AXISYM_X) THEN
C
C      EXTRA VISCOUS TERMS FOR AXISYMM. SITUATION ( X-AXIS IS AXIS OF SYMMETRY)
C
      FAKVW=0.1
C
      DO J=jsolv_b,jsolv_e
      DO I=isolv_b,isolv_e
        II=I+(J-1)*NI
      F(ISSP+II)=F(ISSP+II)-2.*F(ISVIS+II)*F(ISVOLP+II)
     >                    *KBLK(II)/ (F(ISR+II)*F(ISR+II)+SMALL)

      IF(SWIRL) THEN
        SUIJ=F(ISDEN+II)*F(ISW+II)*F(ISW+II)*F(ISVOLP+II)/
     >                                    (F(ISR+II)+SMALL)
        SUIJ=SUIJ*KBLK(II)
      F(ISSU+II)=F(ISSU+II)+SUIJ*(1.+FAKVW*F(ISV+II)/(F(ISW+II)+SMALL))
      F(ISSP+II)=F(ISSP+II)-FAKVW*SUIJ/(F(ISW+II)+SMALL)
      ENDIF
        end do
        end do
       
       ENDIF

       IF(SOLVE(11).AND.NATCONV) THEN 

C*    V-MOMENTUM SOURCE TERM DUE TO THERMAL BUOYANCY                      
C     GRAVITATIONAL FORCE ASSUMED TO ACT ALONG -V (DOWNWARD)  DIRECTION
C
        RAFACT=RAYLEIGH*PRLAM                       
        DO J=2,NJ-1
        DO I=2,NI-1 
          II=I+ (J-1)*NI         
          VOLP = F(ISVOLP +II)
          F(ISSU+II)=F(ISSU+II)+RAFACT*F(ISS+II)*VOLP
        END DO
        END DO

       ENDIF     ! END OF NATCONV LOOP
C                                                          
C  INCLUSION OF UNSTEADY TERMS
C
      DO J=jsolv_b,jsolv_e
      DO I=isolv_b,isolv_e
        ii=i+(j-1)*ni
        VOLP=F(isvolp+ii)
        term=volp*kblk(ii)/delt
        if (TEMPSCH .EQ. 1 .OR. NTSTEP .EQ.1) THEN
C
C FIRST ORDER TIME DISCRETISATION
C
         su_unst=term*v(ii,iblk,1)
         sp_unst=term
        else
C
C SECOND ORDER TIME DISCRETISATION
C
         su_unst=term*(2*v(ii,iblk,2)-0.5*v(ii,iblk,1))
         sp_unst=1.5*term
        end if

        F(ISSU+ii)=f(issu+ii)+su_unst
        F(ISSp+ii)=f(issp+ii)-sp_unst
        end do
        end do

C                                                          
      RETURN                                                            
      END 
C
C************************************************************************
C                                                
      SUBROUTINE SORCE4
C
C     CALCULATION OF SOURCE TERMS FOR WR-MOMENTUM EQUATION                    
C
C************************************************************************
C                                                      
      INCLUDE 'com2d'                                                
C
C    CROSS DERIVATIVE DIFFUSIVE FLUXES
C 
      CALL CDFLUX(ISWR,4) 
C                                              
      DO J=jsolv_b,jsolv_e
      DO I=isolv_b,isolv_e
        II=I+(J-1)*NI
        WE=F(ISFX+II+1 )*F(ISWR+II+1 )+(1.-F(ISFX+II+1 ))*F(ISWR+II)        
        WW=F(ISFX+II   )*F(ISWR+II   )+(1.-F(ISFX+II   ))*F(ISWR+II-1)      
        WN=F(ISFY+II+NI)*F(ISWR+II+NI)+(1.-F(ISFY+II+NI))*F(ISWR+II)        
        WS=F(ISFY+II   )*F(ISWR+II   )+(1.-F(ISFY+II   ))*F(ISWR+II-NI)     
       BESUC(II)=F(ISSU+II)                                               
       BESPC(II)=F(ISSP+II)                                               
      F(ISSU+II)=F(ISSU+II)                                               
     >-2.*(F(ISB12W+II+1 )*F(ISVISW+II+1 )*WE-F(ISB12W+II)*F(ISVISW+II)*WW  
     >    +F(ISB22S+II+NI)*F(ISVISS+II+NI)*WN-F(ISB22S+II)*F(ISVISS+II)*WS) 
	end do
	end do
      RETURN                                                            
      END
C
C************************************************************************
C                                                      
      SUBROUTINE SORCE11
C
C     CALCULATION OF SOURCE TERMS FOR SCALAR EQUATION                           
C
C************************************************************************
C                                                      
      INCLUDE 'com2d'
C                                                
C    CROSS DERIVATIVE DIFFUSIVE FLUXES
C
      CALL CDFLUX(ISS,11)
C
C  INCLUSION OF UNSTEADY TERMS
C
      DO J=jsolv_b,jsolv_e
      DO I=isolv_b,isolv_e
        ii=i+(j-1)*ni
        VOLP=F(isvolp+ii)
        term=volp*kblk(ii)/delt
        if (TEMPSCH .EQ. 1 .OR. NTSTEP .EQ.1) THEN
C
C FIRST ORDER TIME DISCRETISATION
C
         su_unst=term*scal(ii,iblk,1)
         sp_unst=term
        else
C
C SECOND ORDER TIME DISCRETISATION
C
         su_unst=term*(2*scal(ii,iblk,2)-0.5*scal(ii,iblk,1))
         sp_unst=1.5*term
        end if

        F(ISSU+ii)=f(issu+ii)+su_unst
        F(ISSP+ii)=f(issp+ii)-sp_unst
      
        end do
        end do

      RETURN                                                            
      END     
C
C************************************************************************
C                                                      
C  
      SUBROUTINE ADI(ISVAR,N)
C
C     THIS IS LINE ITERATIVE SOLVER BASED ON 
C     THOMA'S ALGORITHM COUPLED TO SWEEPING
C     ALONG ALTERNATE DIRECTION
C                       
C     VARIABLE INDEX N: 1(PP), 2(U), 3(V), 4(WR), 5(TE), 6(ED), 7(S)
C    
C***********************************************************************
C
      INCLUDE 'com2d'                                               
C
C
      DIMENSION RES(IDIM1),P(NMAX),Q(NMAX),B(NMAX)                    
C
      DATA FRAC/0.3/    

	do ii=1,nij
	  res(ii)=0.
	  resd(n,ii)=0.
	end do
                   
C
C     COMPUTATION OF STARTING RESIDUALS
C
        RS=0.
C                                                    
        DO J=jsolv_b,jsolv_e
        DO I=isolv_b,isolv_e
	  ii=i+(j-1)*ni
         TERM=F(ISAN+II)*F(ISVAR+II+NI)+F(ISAS+II)*F(ISVAR+II-NI)       
     >       +F(ISAE+II)*F(ISVAR+II+1 )+F(ISAW+II)*F(ISVAR+II-1 )       
     >       +F(ISSU+II)-F(ISAP+II)*F(ISVAR+II)
         RES(II)=TERM*KBLK(II)
           RS=RS+RES(II)*RES(II)  
	   resd(n,ii)=res(ii)
         END DO                          
         END DO                          
C
        RSTART=SQRT(RS)
C
        IF(N.NE.1) RMOM(N)=RSTART
        RELRS=1.
C                                                                  
C     INNER ITERATION SWEEP STARTS HERE                                        
C
      ITER = 0                                            
C
20    CONTINUE
        IF (RELRS .GT. 1.) write(55,*)ntstep,ibswp,iblk,iter,n,relrs
C
C
C       WRITE(6,2)'ST2D IT. ',ITER,' REL.RES.: ',RELRS,'ABS.RES.: ',RS
C                             
         IF(RELRS.LE.FRAC.OR.ITER.GE.NSWP(N)) GO TO 4000
C                             
       ITER=ITER+1
C      
C     COMMENCE W-E SWEEP                                             
C
C     TDMA ALONG I DIRECTION ONLY
C
        DO J=jsolv_b,jsolv_e
C
         P(isolv_b-1)=0.
         Q(isolv_b-1)=F(ISVAR+isolv_b-1+(J-1)*NI)
C
        DO I=isolv_b,isolv_e
          II=I+(J-1)*NI
          B(I)=F(ISAS+II)*F(ISVAR+II-NI)+F(ISAN+II)*F(ISVAR+II+NI)
     >         +F(ISSU+II)
          RDNM=1./(F(ISAP+II)-F(ISAW+II)*P(I-1))
          P(I)=F(ISAE+II)*RDNM
          Q(I)=(B(I)+Q(I-1)*F(ISAW+II))*RDNM
          END DO
C
        DO I=isolv_e,isolv_b,-1
          II=I+(J-1)*NI
          F(ISVAR+II)=F(ISVAR+II+1)*P(I)+Q(I)
          END DO
C
       END DO
C                                                                   
C     COMMENCE S-N SWEEP                                             
C
C     TDMA ALONG J DIRECTION ONLY
C
        DO I=isolv_b,isolv_e
C
         P(jsolv_b-1)=0.
         Q(jsolv_b-1)=F(ISVAR+i+(Jsolv_b-1-1)*NI)
C
        DO J=jsolv_b,jsolv_e
          II=I+(J-1)*NI
          B(J)=F(ISAW+II)*F(ISVAR+II-1)+F(ISAE+II)*F(ISVAR+II+1)
     >         +F(ISSU+II)
          RDNM=1./(F(ISAP+II)-F(ISAS+II)*P(J-1))
          P(J)=F(ISAN+II)*RDNM
          Q(J)=(B(J)+Q(J-1)*F(ISAS+II))*RDNM
          END DO
C
        DO J=jsolv_e,jsolv_b,-1
          II=I+(J-1)*NI
          F(ISVAR+II)=F(ISVAR+II+NI)*P(J)+Q(J)
          END DO
C
       END DO
C             
C        SPECIAL STOPPING CRITERION FOR THE PP-EQUATION ONLY            
C
          RS=0.      
C                                                    
        DO J=jsolv_b,jsolv_e
        DO I=isolv_b,isolv_e
         II=I+(J-1)*NI
         TERM=F(ISAN+II)*F(ISVAR+II+NI)+F(ISAS+II)*F(ISVAR+II-NI)       
     >       +F(ISAE+II)*F(ISVAR+II+1 )+F(ISAW+II)*F(ISVAR+II-1 )       
     >       +F(ISSU+II)-F(ISAP+II)*F(ISVAR+II)  
         RES(II)=TERM*KBLK(II)                               
         RS=RS+RES(II)*RES(II)                            
         END DO       
         END DO 
C         
         RS=SQRT(RS)                                                    
         RELRS=RS/RSTART                                                
C
C        
         GO TO 20                                                         
              
 4000 CONTINUE
C                                
 2    FORMAT(1X,A,I3,1P,1X,2(A,E12.4,1X))
C
      RETURN                                                            
      END                                                               
C    
C***********************************************************************
C
      SUBROUTINE SIP(ISVAR,N)
C
C  THIS SOLVER USES THE STRONGLY IMPLICIT METHOD BY STONE         
C 
C*******************************************************************
C
      INCLUDE 'com2d'                                                     
C
      DIMENSION BN(IDIM1),BS(IDIM1),BE(IDIM1)
     >         ,BW(IDIM1),BP(IDIM1),RES(IDIM1)
C                                
      DATA FRAC/0.1/
      DATA ALFA /0.9/
C
      DO 5567 II=1,NIJ
      BN(II)=0.
      BS(II)=0.
      BE(II)=0.
      BW(II)=0.
      BP(II)=0.
      RES(II)=0.
	resd(n,ii)=0.
5567  CONTINUE
C
C   COMPUTATION OF STARTING RESIDUALS
C
        RS=0.
C                                                    
        DO J=jsolv_b,jsolv_e
        DO I=isolv_b,isolv_e
         II = I+(J-1)*NI
         TERM=F(ISAN+II)*F(ISVAR+II+NI)+F(ISAS+II)*F(ISVAR+II-NI)       
     >       +F(ISAE+II)*F(ISVAR+II+1 )+F(ISAW+II)*F(ISVAR+II-1 )       
     >       +F(ISSU+II)-F(ISAP+II)*F(ISVAR+II)
         RES(II)=TERM*KBLK(II)

	 resd(n,ii)=res(ii)	
         RS=RS+RES(II)*RES(II)   
        END DO                        
        END DO
C
        RS=SQRT(RS)+SMALL
	RSTART=RS
C
        RELRS=1.
C
C      CALCULATION OF COEFICIENTS OF THE MATRIX AFTER
C      DECOMPOSITION INTO UPPER AND LOWER TRIANGLES
C
      DO J=jsolv_b,jsolv_e
      DO I=isolv_b,isolv_e
C
       II=I+(J-1)*NI
C
        ASII=F(ISAS+II)
        AWII=F(ISAW+II)
        ANII=F(ISAN+II)
        AEII=F(ISAE+II)
        APII=F(ISAP+II)  
    
      BW(II)=-AWII/(1.+ALFA*BN(II-1))
      BS(II)=-ASII/(1.+ALFA*BE(II-NI))
      P1=ALFA*BW(II)*BN(II-1)
      P2=ALFA*BS(II)*BE(II-NI)
      BP(II)=1./(APII+P1+P2-BW(II)*BE(II-1)
     *                   -BS(II)*BN(II-NI))
      BN(II)=(-ANII-P1)*BP(II)
      BE(II)=(-AEII-P2)*BP(II)
	end do
	end do
C
C     INNER ITERATION SWEEP STARTS HERE
C
        ITER = 0
C
  20    CONTINUE
	
 	IF (RELRS .GT. 1.) write(55,*)ntstep,ibswp,iblk,iter,n,relrs 
c	IF (n .eq. 1) write(56,*)ntstep,ibswp,iblk,iter,n,relrs 
C
c      WRITE(6,2)'ST2D IT. ',ITER,' REL.RES.: ',RELRS,'ABS.RES.: ',RS
C
        IF(RELRS.LE.FRAC.OR.ITER.GE.NSWP(N)) GO TO 4000
C
        ITER = ITER+1
C
        DO 27 J=jsolv_b,jsolv_e
        DO 27 I=isolv_b,isolv_e
         II=I+(J-1)*NI
        RES(II)=(RES(II)-BS(II)*RES(II-NI)-
     >            BW(II)*RES(II-1))*BP(II)
  27    CONTINUE
C
C       CALCULATION OF SOLUTION VECTOR AND UPDATING OF THE FIELD
C
        DO 28 J=jsolv_b,jsolv_e
        DO 28 I=isolv_b,isolv_e
        JB=jsolv_b+jsolv_e-J
        IB=isolv_b+isolv_e-I
        II=IB+(JB-1)*NI
        RES(II)=RES(II)-BE(II)*RES(II+1)-BN(II)*RES(II+NI)
  28    CONTINUE

	do j=jsolv_b,jsolv_e
	do i=isolv_b,isolv_e
	 ii=i+(j-1)*ni
         FNEW=F(ISVAR+II)+RES(II)*KBLK(II)                       
         F(ISVAR+II)=FNEW 
	end do
	end do
C                                                  
C  SPECIAL ARRANGEMNET FOR INTRA INNER SWEEP TRANSFER
C  OF FIELD VALUES FOR THE CYCLIC BOUNDARY           
C                                                  
        CALL BOUNDS(1)                                          
                          
C                                                                        
C        SPECIAL STOPPING CRITERION FOR THE PP-EQUATION ONLY            
C
   
         RS=0.      
C                                                    
         DO  J=jsolv_b,jsolv_e
         DO  I=isolv_b,isolv_e
         II=I+(J-1)*NI
         TERM=F(ISAN+II)*F(ISVAR+II+NI)+F(ISAS+II)*F(ISVAR+II-NI)       
     >       +F(ISAE+II)*F(ISVAR+II+1 )+F(ISAW+II)*F(ISVAR+II-1 )       
     >       +F(ISSU+II)-F(ISAP+II)*F(ISVAR+II)
         RES(II)=TERM*KBLK(II)
         RS=RS+RES(II)*RES(II)
         END DO                            
         END DO
C        
         RS=SQRT(RS)                                                    
         RELRS=RS/RSTART                                                
C         
         GO TO 20 
C              
C                                                   
4000     CONTINUE

 2    FORMAT(1X,A,I3,1P,1X,2(A,E12.4,1X))
C 
      RETURN
      END
C
C***********************************************************************
C      
      SUBROUTINE SOLVER(ISVAR,N)                                                
C
C     CALLS THE APPROPRIATE SOLVER FOR THE VARIABLE N
C     AND STORES THE REQUIRED SUM AND RECIPROCAL OF 
C     SUM OF COEFFICIENTS FOR MOMENTUM EQUATIONS TO 
C     BE USED IN THE MOMENTUM INTERPOLATION
C                                                    
C***********************************************************************        
C
      INCLUDE 'com2d'
C                    
C FOR IMPLICIT UNDER RELAXATION
C                                                                 
      RINV=1./RELAX(N)
      RELAXM=1.-RELAX(N)                               
C                                                   
C                                                        
C
C     SOURCE TERM LINEARISATION TO ENSURE DIAGONAL DOMINANCE
C
      IF(N.NE.1)  CALL DIAG_DOM(ISVAR,N)

C
      IF(N.GT.1.AND.N.LT.5) THEN 
C                                               
      DO J=jsolv_b,jsolv_e
      DO I=isolv_b,isolv_e
C
            II=I+(J-1)*NI
            ASUM=F(ISAW+II)+F(ISAE+II)+F(ISAN+II)+F(ISAS+II)                  
            APSUM=(ASUM-F(ISSP+II))*RINV                                       
            F(ISAP+II)=APSUM
C
            F(ISAPU+II+(N-2)*NIJ)=1./(APSUM+SMALL)
C
            IF(LINKPV.EQ.2) THEN
	    F(ISSUMU+II+(N-2)*NIJ)=ASUM/(APSUM+SMALL)
C
            ELSE
            F(ISSUMU+II+(N-2)*NIJ)=0.
C
            ENDIF
C
            F(ISSU+II)=F(ISSU+II)+RELAXM*F(ISVAR+II)*APSUM                     
	end do
	end do
C               
       ELSE
C                                                                     
      DO J=jsolv_b,jsolv_e
      DO I=isolv_b,isolv_e
         II=I+(J-1)*NI                                                

        APSUM=(F(ISAW+II)+F(ISAE+II)+F(ISAN+II)+F(ISAS+II)
     >            -F(ISSP+II))*RINV
            F(ISAP+II)=APSUM
C
c           IF(N.EQ.1)  APP(II)=F(ISAP+II) ! NORMCOEFF CALLED
            IF(N.EQ.1)  APP(II)=1.        
C
         F(ISSU+II)=F(ISSU+II)+RELAXM*F(ISVAR+II)*APSUM

	end do
	end do

      ENDIF
C
      CALL COEFF_CHECK(N)
C
       IF (N.NE. 1) CALL NORMCOEFF
C
      IF(NSOLV(N).EQ.1)  CALL ADI(ISVAR,N)

      IF(NSOLV(N).EQ.2)  CALL SIP(ISVAR,N)                        
                                                                               
      RETURN
      END                                                                       
C
C************************************************************************
C
      SUBROUTINE BNDPH(FF)               
C                         
C     COMPUTES PRESSURE AT BOUNDARY NODES USING
C     LINEAR EXTRAPOLATION FROM INTERIOR POINTS
C                         
C***********************************************************************
C
      INCLUDE 'com2d'
C                                                 
      DIMENSION FF(NIJ)

      DO 100 J=2,NJM
      DO 100 I=2,NIM
      II=I+(J-1)*NI
          L1=MIN(1,KW(II)/10)                                           
	  if(kw(ii) .eq. 15) l1=0
      FF(II-1 )=(1-L1)*FF(II-1)                                   
     >   +L1*((1.+F(ISFX+II+1))*FF(II)-F(ISFX+II+1)*FF(II+1))     
          L2=MIN(1,KE(II)/10)                                           
	  if(ke(ii) .eq. 15) l2=0
      FF(II+1 )=(1-L2)*FF(II+1)                                   
     >   +L2*((2.-F(ISFX+II))*FF(II)-(1.-F(ISFX+II))*FF(II-1))    
          L3=MIN(1,KS(II)/10)                                           
	  if(ks(ii) .eq. 15 ) l3=0
      FF(II-NI)=(1-L3)*FF(II-NI)                                  
     >   +L3*((1.+F(ISFY+II+NI))*FF(II)-F(ISFY+II+NI)*FF(II+NI))  
          L4=MIN(1,KN(II)/10)                                           
	  if(kn(ii) .eq. 15 ) l4=0
      FF(II+NI)=(1-L4)*FF(II+NI)                                  
     >   +L4*((2.-F(ISFY+II))*FF(II)-(1.-F(ISFY+II))*FF(II-NI))   
         
100   CONTINUE                                                          
          II=2+(1-1)*NI                                                 
      FF(II-1)=(1.+F(ISFX+II+1))*FF(II)-F(ISFX+II+1)*FF(II+1)  
          II=2+(NJ-1)*NI                                                
      FF(II-1)=(1.+F(ISFX+II+1))*FF(II)-F(ISFX+II+1)*FF(II+1)  
          II=NIM+(1-1)*NI                                               
      FF(II+1)=(2.-F(ISFX+II))*FF(II)-(1.-F(ISFX+II))*FF(II-1) 
          II=NIM+(NJ-1)*NI                                              
      FF(II+1)=(2.-F(ISFX+II))*FF(II)-(1.-F(ISFX+II))*FF(II-1) 
      RETURN
      END                                                               
C
C***********************************************************************
C 
      SUBROUTINE BNDVIS(FF)               
C                         
C     COMPUTES VISCOSITY AT BOUNDARY NODES USING
C     ZERO ORDER EXTRAPOLATION FROM INTERIOR POINTS
C                         
C***********************************************************************
C
      INCLUDE 'com2d'
C                                                 
      DIMENSION FF(NIJ)

      DO  J=2,NJM
      DO  I=2,NIM
      II=I+(J-1)*NI
          L1=MIN(1,KW(II)/10)                                           
	  if(kw(ii) .ge. 14.or.kw(ii).eq.11) l1=0
      FF(II-1 )=(1-L1)*FF(II-1)+L1*ff(II)                         
          L2=MIN(1,KE(II)/10)                                           
	  if(ke(ii) .ge. 14.or.ke(ii).eq.11) l2=0
      FF(II+1 )=(1-L2)*FF(II+1)+L2*FF(II)                                   
          L3=MIN(1,KS(II)/10)                                           
	  if(ks(ii) .ge. 14.or.ks(ii).eq.11) l3=0
      FF(II-NI)=(1-L3)*FF(II-NI) +L3*FF(II)                                 
          L4=MIN(1,KN(II)/10)                                           
	  if(kn(ii) .ge. 14.or.kn(ii).eq.11) l4=0
      FF(II+NI)=(1-L4)*FF(II+NI)+L4*FF(II)                        
      END DO
      END DO
      RETURN
      END              
                                                 
      SUBROUTINE  DIAG_DOM(ISVAR,N)
C
C
C      MODIFY SU AND SP ARRAY TO ENSURE DIAGONAL DOMINANCE       
C
C***********************************************************************
C
      INCLUDE 'com2d'                                              
C
C***********************************************************************
C
C       write(1234,*) 'NTSTEP=',NTSTEP,'IBSWP=',IBSWP,'IBLK=',IBLK,'N=', N
C       write(1235,*) 'NTSTEP=',NTSTEP,'IBSWP=',IBSWP,'IBLK=',IBLK,'N=', N
     
      DO J=JSOLV_B,JSOLV_E
      DO I=ISOLV_B,ISOLV_E
      II=I+(J-1)*NI
      SUII=F(ISSU+II)
      SPII=F(ISSP+II)
      PHII=F(ISVAR+II)
C      write(1234,*) J,I,SUII,SPII,PHII

      IF(SPII.GT.0.) THEN
       F(ISSU+II)=F(ISSU+II)+SPII*PHII
       F(ISSP+II)=0.
       SUII=F(ISSU+II)
       SPII=F(ISSP+II)
      ENDIF
      
C
C     IF(PHII.NE.0.) THEN

C     SU_BYPH=SUII/(PHII+small)
C     IF(SU_BYPH.LT.0) THEN
C      F(ISSU+II)=0.
C      F(ISSP+II)=F(ISSP+II)+SU_BYPH 
C     END IF

C     END IF 
C      write(1235,*) J,I,F(ISSU+II),F(ISSP+II)

      END DO
      END DO
    
      RETURN 
      END
      
C
C***********************************************************************
C
      SUBROUTINE BINTER                                                         
C                                                                               
C     THIS SUBROUTINE SETS THE RIGHT INTERPOLATION FACTORS                      
C     FOR ALL NEAR BOUNDARY CELLS                                               
C                                                                               
C***********************************************************************
C
      INCLUDE 'com2d'                                                           
C                                                                               
	do j=2,njm
	do i=2,nim
	ii=i+(j-1)*ni
      L1=MIN(1,KW(II)/10)                                                       
        if (kw(ii) .ge. 15 .and. kw(ii).lt.17) L1=0
      L2=MIN(1,KE(II)/10)                                                       
        if (ke(ii) .ge. 15 .and. ke(ii).lt.17) L2=0
      L3=MIN(1,KS(II)/10)                                                       
        if (ks(ii) .ge. 15 .and. ks(ii).lt.17) L3=0
      L4=MIN(1,KN(II)/10)                                                       
        if (kn(ii) .ge. 15 .and. kn(ii).lt.17) L4=0
      F(ISFX+II)=(1-L1)*F(ISFX+II)                                              
      F(ISFX+II+1)=(1-L2)*F(ISFX+II+1)+L2                                       
      F(ISFY+II)=(1-L3)*F(ISFY+II)                                              
      F(ISFY+II+NI)=(1-L4)*F(ISFY+II+NI)+L4                                     
	end do
	end do
      RETURN                                                                    
      END                                                                       
C                         
C***********************************************************************
C
      SUBROUTINE CUTRD
C                                                     
C     READ FIELD VALUES ON CUT BOUNDARIES FOR THE NEIGHBOURING BLOCKS
C                                         
C************************************************************************
C
C                                                     
      INCLUDE 'com2d'
C
       dimension temp(nmax*28)
C
	is_beg=isu
	ivar_b=nvu
	ivar_e=nvvis
	ivar_t=ivar_e-ivar_b+1

	do nbound=1,nrbnda
	  if (cbnd(nbound,2)(1:2).eq. 'CU') then
            nrrec=nbcutn(nbound)+(nbblk(nbound)-1)*maxcut

		read(9,rec=nrrec)temp
	    if(cbnd(nbound,1)(1:2).eq.'EA'.or.
     >           cbnd(nbound,1)(1:2).eq.'WE') then
	     if (cbnd(nbound,1)(1:2).eq. 'EA') IFLAG=-1 
	     if (cbnd(nbound,1)(1:2).eq. 'WE') IFLAG=1 
           i=nbnd(nbound,1)-iflag 
	     js=nbnd(nbound,3)
	     je=nbnd(nbound,4)

	     iskip=0
	     lensten=je-js+1
           do ire=1,2
	        do ivar=ivar_b,ivar_e
	          isphi=is_beg+(ivar-ivar_b)*nij
                DO jlim=js,je           
	           j=jlim
	           if (ndir(nbound) .eq. -1) j=js+je-jlim
                 II=I+(j-1)*NI
  	           ii1=ii-iflag
	           ij=jlim-js+1+(ivar-ivar_b)*lensten+iskip
	           if (ire .eq. 1) then
	              f(isphi+ii)=temp(ij)
	           else
	              f(isphi+ii1)=temp(ij)
	           end if
	          end do
	        end do
	        iskip=iskip+ivar_t*lensten
	     end do

	     do ivar=nvapu,nvapv
	        isphi=is_beg+(ivar-ivar_b)*nij
               DO jlim=js,je           
	             j=jlim
	             if (ndir(nbound) .eq. -1) j=js+je-jlim
                II=I+(j-1)*NI
	          ij=jlim-js+1+(ivar-nvapu)*lensten+iskip 
	          f(isphi+ii)=temp(ij)
	        end do
	        iskip=iskip+lensten
	     end do

	    else

	       if (cbnd(nbound,1)(1:2).eq. 'NO') JFLAG=-1 
	       if (cbnd(nbound,1)(1:2).eq. 'SO') JFLAG=1 
	       j=nbnd(nbound,3)-jflag  
	       is=nbnd(nbound,1)
	       ie=nbnd(nbound,2)


	       jskip=0
	       lensten=ie-is+1
	       do jre=1,2
	         do ivar=ivar_b,ivar_e
	           isphi=is_beg+(ivar-ivar_b)*nij
                 DO ilim=is,ie           
	             i=ilim
	             if (ndir(nbound) .eq. -1) i=is+ie-ilim
                   II=I+(j-1)*NI
	             ii1=ii-ni*jflag
	             ij=ilim-is+1+(ivar-ivar_b)*lensten+jskip
	             if (jre .eq. 1) then
	               f(isphi+ii)=temp(ij)
	             else
	               f(isphi+ii1)=temp(ij)
	             end if
	            end do
	         end do
	         jskip=jskip+ivar_t*lensten
	       end do

	       do ivar=nvapu,nvapv
	         isphi=is_beg+(ivar-ivar_b)*nij
               DO ilim=is,ie           
	           i=ilim
	           if (ndir(nbound) .eq. -1) i=is+ie-ilim
                 II=I+(j-1)*NI
	           ij=ilim-is+1+(ivar-nvapu)*lensten+jskip 
	           f(isphi+ii)=temp(ij)
	         end do
	       end do

          end if
      end if
	end do
	
	return
	end
C
C***********************************************************************
C 
      SUBROUTINE CUTWR  
C                                                     
C     WRITE FIELD VALUES ON CUT BOUNDARIES FOR THE NEIGHBOURING BLOCKS
C                                         
C************************************************************************
C                                                     
      INCLUDE 'com2d'
C
	dimension temp(nmax*28)
C                                                     

	is_beg=isu
	ivar_b=nvu
	ivar_e=nvvis
	ivar_t=ivar_e-ivar_b+1
	icut=0

	do nbound=1,nrbnda
	  if (cbnd(nbound,2)(1:2).eq. 'CU') then
	    icut=icut+1
            nwrec=icut+(iblk-1)*maxcut

	    if(cbnd(nbound,1)(1:2).eq.'EA'.or.
     >           cbnd(nbound,1)(1:2).eq.'WE') then
	     if (cbnd(nbound,1)(1:2).eq. 'EA') IFLAG=-1 
	     if (cbnd(nbound,1)(1:2).eq. 'WE') IFLAG=1 
           i=nbnd(nbound,1) 
	     js=nbnd(nbound,3)
	     je=nbnd(nbound,4)

	     iskip=0
	     lensten=je-js+1
           do iwr=1,2
 	         do ivar=ivar_b,ivar_e
	           isphi=is_beg+(ivar-ivar_b)*nij
              DO jlim=js,je           
	             j=jlim
	             if (ndir(nbound) .eq. -1) j=js+je-jlim
                II1=I+(j-1)*NI
	             ii2=ii1+iflag
	             ij=jlim-js+1+(ivar-ivar_b)*lensten+iskip 
	             if (iwr .eq.1 ) then
	               temp(ij)=f(isphi+ii1)
	             else  
	               temp(ij)=0.5*(f(isphi+ii1)+f(isphi+ii2))
	             end if
	           end do
	       end do
	       iskip=iskip+ivar_t*lensten
	     end do

	     do ivar=nvapu,nvapv
	       isphi=is_beg+(ivar-ivar_b)*nij
              DO jlim=js,je           
	             j=jlim
	             if (ndir(nbound) .eq. -1) j=js+je-jlim
               II1=I+(j-1)*NI
	         ij=jlim-js+1+(ivar-nvapu)*lensten+iskip 
	         temp(ij)=f(isphi+ii1)
	       end do
	       iskip=iskip+lensten
	     end do

	  else
	     if (cbnd(nbound,1)(1:2).eq. 'NO') JFLAG=-1 
	     if (cbnd(nbound,1)(1:2).eq. 'SO') JFLAG=1 
	       j=nbnd(nbound,3) 
	       is=nbnd(nbound,1)
	       ie=nbnd(nbound,2)

	       jskip=0
	       lensten=ie-is+1
	       do jwr=1,2
	         do ivar=ivar_b,ivar_e
	           isphi=is_beg+(ivar-ivar_b)*nij
                 DO ilim=is,ie           
	             i=ilim
	             if (ndir(nbound) .eq. -1) i=is+ie-ilim
                   II1=I+(j-1)*NI
	             II2=II1+jflag*ni
	             ij=ilim-is+1+(ivar-ivar_b)*lensten+jskip 
	             if (jwr .eq.1 ) then
	               temp(ij)=f(isphi+ii1)
	             else  
	               temp(ij)=0.5*(f(isphi+ii1)+f(isphi+ii2))
	             end if
	           end do
	         end do
	       jskip=jskip+ivar_t*lensten
	       end do

	       do ivar=nvapu,nvapv
	         isphi=is_beg+(ivar-ivar_b)*nij
               DO ilim=is,ie           
		       i=ilim
	           if (ndir(nbound) .eq. -1) i=is+ie-ilim
                 II1=I+(j-1)*NI
	           ij=ilim-is+1+(ivar-nvapu)*lensten+jskip 
	           temp(ij)=f(isphi+ii1)
	         end do
	       end do
        end if

 	   write(9,rec=nwrec)temp
      end if
	end do
	
	return
	end
C                                         
C************************************************************************
C                                                     
      SUBROUTINE QINTP(DELUU,DELU,DELD,PHUU,PHU,PHD,PHUD)
C                                          
C     THIS SUBROUTINE IS USED TO COMPUTE CELL FACE VALUE (BETWEEN U AND D) OF
C     A VARIABLE USING QUADRATIC INTERPOLATION FROM VALUES AT THREE CONSECUTIVE
C     NODES UU,U,D                                   
C
C*************************************************************************
C
      SMALL = 1.0e-20
                                                  
      P = DELU + DELD     
      Q = DELU + DELUU                          
      R = PHD - PHU                              
      S = PHUU - PHU    
                       

      RDNM = 1./((P*P*Q+P*Q*Q)+SMALL)                   
      A = 4.*(Q*R+P*S)*RDNM                     
      B = 2.*(Q*Q*R-P*P*S)*RDNM                        
      X = DELU/2.                                     
      PHUD = A*X*X+B*X+PHU                          

      RETURN                                      
                                                                        
      END                                                               
C
C*************************************************************************
C                                                        
      FUNCTION DELEW(II)                   
C
C*************************************************************************
C                                                        
C                                                                            
      INCLUDE 'com2d'    
C                                                      
      XSW=F(ISCOX+II)                              
      XSE=F(ISCOX+II+1)                                                 
      XNW=F(ISCOX+II+NI)                          
      XNE=F(ISCOX+II+1+NI)
      YSW=F(ISCOY+II)                           
      YSE=F(ISCOY+II+1)                          
      YNW=F(ISCOY+II+NI)                         
      YNE=F(ISCOY+II+1+NI)                    
      DXEW=0.5*(XSE+XNE-XSW-XNW)                
      DYEW=0.5*(YSE+YNE-YSW-YNW)                 
      DELEW=SQRT(DXEW*DXEW+DYEW*DYEW)                 
      RETURN                                          
      END                                           
C
C*************************************************************************
C
      FUNCTION DELNS(II)                          
C
C*************************************************************************
C
      INCLUDE 'com2d'                                    
C
      XSW=F(ISCOX+II)                                                        
      XSE=F(ISCOX+II+1)                                                        
      XNW=F(ISCOX+II+NI)                             
      XNE=F(ISCOX+II+1+NI)                         
      YSW=F(ISCOY+II)                                                   
      YSE=F(ISCOY+II+1)                           
      YNW=F(ISCOY+II+NI)  
      YNE=F(ISCOY+II+1+NI)                      
      DXNS=0.5*(XNE+XNW-XSE-XSW)                 
      DYNS=0.5*(YNE+YNW-YSE-YSW)                 
      DELNS=SQRT(DXNS*DXNS+DYNS*DYNS)   
C      
      RETURN                                    
      END                                        
C
C*************************************************************************
C
	SUBROUTINE BLKRES
C
C*************************************************************************
C
	include 'com2d'
C
	do ivar=1,11
	 rmom(ivar)=0.
	 do j=jsolv_b,jsolv_e
	   do i=isolv_b,isolv_e
	      ii=i+(j-1)*ni 
	      rmom(ivar)=rmom(ivar)+resd(ivar,ii)*resd(ivar,ii)
	   end do
	 end do
	 rmom(ivar)=sqrt(rmom(ivar))
	end do
	return
	end   
C
C*************************************************************************
C
        SUBROUTINE DIFF_OLD(IGIV,JGIV,ISPHI,DPHDX,DPHDY)
C      
C      EVALUATION OF SPATIAL GRADIENT OF A VARIABLE
C      IN A FIELD PROBLEM USING GREEN'S THEOREM
C
C*************************************************************************
C      
	include 'com2d'
C
      COMMON/deriv/PHC(5)
      DIMENSION PH(NXF,NYF)

C       
c       CALL CORNPHI(IGIV,JGIV,PH)
C       
        
        CALL DIFFPHI(IGIV,JGIV,ISPHI,DPHDX,DPHDY)

        RETURN
        END
C
C*************************************************************************
C      
           SUBROUTINE CORNPHI(IGIV,JGIV,PH)
C      
C THIS SUBROUTINE CALCULATES ALL THE FOUR NEIGHBOURING CORNER PHI VALUES
C FOR A GIVEN CENTRAL POINT IN A CONTROL VOLUME SPECIFIED WITH I,J VALUES
C      
C***********************************************************************
C      
	include 'com2d'
C      
      COMMON/deriv/PHC(5)
      DIMENSION PH(NXF,NYF)
      
        II = IGIV+(JGIV-1)*NI
                                                                          
         IS=0
         JS=0
         
         DO KOUNT=1,4
         
         IF(KOUNT.EQ.2) IS=1
         IF(KOUNT.EQ.3) JS=1
         IF(KOUNT.EQ.4) IS=0
         
         XP=F(ISCOX+IGIV+IS+(JGIV+JS-1)*NI)
         YP=F(ISCOY+IGIV+IS+(JGIV+JS-1)*NI)
         
         X1=F(ISX+IGIV-1+IS+(JGIV+JS-2)*NI)
         Y1=F(ISY+IGIV-1+IS+(JGIV+JS-2)*NI)
         X2=F(ISX+IGIV+IS+(JGIV+JS-2)*NI)
         Y2=F(ISY+IGIV+IS+(JGIV+JS-2)*NI)
         X3=F(ISX+IGIV+IS+(JGIV+JS-1)*NI)
         Y3=F(ISY+IGIV+IS+(JGIV+JS-1)*NI)
         X4=F(ISX+IGIV-1+IS+(JGIV+JS-1)*NI)
         Y4=F(ISY+IGIV-1+IS+(JGIV+JS-1)*NI)
         
         RSW=SQRT((XP-X1)*(XP-X1)+(YP-Y1)*(YP-Y1))+SMALL
         RSE=SQRT((XP-X2)*(XP-X2)+(YP-Y2)*(YP-Y2))+SMALL
                                                                          
         RNE=SQRT((XP-X3)*(XP-X3)+(YP-Y3)*(YP-Y3))+SMALL
         RNW=SQRT((XP-X4)*(XP-X4)+(YP-Y4)*(YP-Y4))+SMALL
         
        ADNM=1./RSW + 1./RSE + 1./RNE + 1./RNW
        
        PH1 = PH(IGIV-1+IS,JGIV-1+JS)
        PH2 = PH(IGIV+IS,JGIV-1+JS)
        PH3 = PH(IGIV+IS,JGIV+JS)
        PH4 = PH(IGIV-1+IS,JGIV+JS)
        
        ANUM=PH1/RSW + PH2/RSE + PH3/RNE + PH4/RNW
        
        PHC(KOUNT)=ANUM/ADNM
C	write(17,*)phc(kount)

        END DO
        
        IF(KW(II) .EQ. 1) THEN
          PHC(1)=0.5*(PH(IGIV-1,JGIV)+PH(IGIV-1,JGIV-1))
          PHC(4)=0.5*(PH(IGIV-1,JGIV)+PH(IGIV-1,JGIV+1))
        END IF
        
        IF(KE(II) .EQ. 1) THEN
          PHC(2)=0.5*(PH(IGIV+1,JGIV)+PH(IGIV+1,JGIV-1))
          PHC(3)=0.5*(PH(IGIV+1,JGIV)+PH(IGIV+1,JGIV+1))
                                                                          
        END IF
        
        IF(KS(II) .EQ.1) THEN
          PHC(1)=0.5*(PH(IGIV,JGIV-1)+PH(IGIV-1,JGIV-1))
          PHC(2)=0.5*(PH(IGIV,JGIV-1)+PH(IGIV+1,JGIV-1))
        END IF
        
        IF(KN(II) .EQ.1) THEN
          PHC(3)=0.5*(PH(IGIV,JGIV+1)+PH(IGIV+1,JGIV+1))
          PHC(4)=0.5*(PH(IGIV,JGIV+1)+PH(IGIV-1,JGIV+1))
        END IF
        
        PHC(5)=PHC(1)
        
        RETURN
        END
C      
C***********************************************************************
C      
        SUBROUTINE DIFFPHI(IGIV,JGIV,ISPHI,DPHDX,DPHDY)
C      
C       SUBROUTINE DIFFPHI CALCULATES DIFFERENTIAL OF PHI W.R.T X OR Y
C      
C************************************************************************
C      
	include 'com2d'
C      
      COMMON/deriv/PHC(5)

      DIMENSION XX(5),YY(5),iiskp(4),ijskp(4),f_int(4)
C     
        iiskp(1)=0
        ijskp(1)=-1
        iiskp(2)=1
        ijskp(2)=0
        iiskp(3)=0
        ijskp(3)=1
        iiskp(4)=-1
        ijskp(4)=0

	iigiv=igiv+(jgiv-1)*ni
	f_int(1)=f(isfy+iigiv)
	f_int(2)=1-f(isfx+iigiv+1)
	f_int(3)=1-f(isfy+iigiv+ni)
	f_int(4)=f(isfx+iigiv)

        XX(1)=F(ISCOX+IGIV+(JGIV-1)*NI)
        XX(2)=F(ISCOX+IGIV+1+(JGIV-1)*NI)
        XX(3)=F(ISCOX+IGIV+1+(JGIV+1-1)*NI)
        XX(4)=F(ISCOX+IGIV+(JGIV+1-1)*NI)
        XX(5)=XX(1)
        
        YY(1)=F(ISCOY+IGIV+(JGIV-1)*NI)
        YY(2)=F(ISCOY+IGIV+1+(JGIV-1)*NI)
        YY(3)=F(ISCOY+IGIV+1+(JGIV+1-1)*NI)
        YY(4)=F(ISCOY+IGIV+(JGIV+1-1)*NI)
        YY(5)=YY(1)
                                                                          
        DPHDX=0.
        DPHDY=0.
        DXAB=XX(1)-XX(3)
        DYAB=YY(1)-YY(3)
        DXCD=XX(2)-XX(4)
        DYCD=YY(2)-YY(4)
        
         AREA=ABS(0.5*(DXAB*DYCD-DXCD*DYAB))
         
        DO ISIDE = 1,4
c       PHAVG=0.5*(PHC(ISIDE)+PHC(ISIDE+1))
        iigiv=igiv+(jgiv-1)*ni
        iinb=igiv+iiskp(iside)+(jgiv+ijskp(iside)-1)*ni
        phavg=0.5*(f(isphi+iigiv)+f(isphi+iinb))
        phavg=f_int(iside)*f(isphi+iigiv)+(1-f_int(iside))*f(isphi+iinb)
        AN1=-(YY(ISIDE)-YY(ISIDE+1))
        AN2=(XX(ISIDE)-XX(ISIDE+1))
        DPHDX=DPHDX+PHAVG*AN1
        DPHDY=DPHDY+PHAVG*AN2
C	write(16,*)area,an1,an2,phavg	
        END DO
        DPHDX=DPHDX/AREA
        DPHDY=DPHDY/AREA
C	write(16,*)dpHdx,dpHdy	
        
        RETURN
        END
C
C                                                                               
      SUBROUTINE DIFF(II,ISPHI,DPHDX,DPHDY)

C     
C     DIFFERENTIATION OF PHI WITH RESPECT TO TWO CARTESIAN DIRECTION
C
C**********************************************************************        
      INCLUDE 'com2d' 
C                    
C --- DIFFERENCES ALONG X1 AND X2 DIRECTION
C     FOR GRADIENT EVALUATION AT CENTRAL NODE              
C
       FXP=F(ISFX+II)                                                         
       FYP=F(ISFY+II)                                                         
       FXM=1.-FXP                                                             
       FYM=1.-FYP                                                             
       FXIP=F(ISFX+II+1)                                                       
       FYJP=F(ISFY+II+NI)                                                      
       FXIM=1.-FXIP                                                            
       FYJM=1.-FYJP                                                            
        
       PHIP=F(ISPHI+II)
 
       DPHEW=FXIP*F(ISPHI+II+1)+FXIM*PHIP
     >                           -(FXP*PHIP+FXM*F(ISPHI+II-1))                
       DPHNS=FYJP*F(ISPHI+II+NI)+FYJM*PHIP
     >                           -(FYP*PHIP+FYM*F(ISPHI+II-NI ))                
      
       B11P=F(ISB11P+II)
       B21P=F(ISB21P+II)
       B12P=F(ISB12P+II)
       B22P=F(ISB22P+II)
       
       RVOLP=1./F(ISVOLP+II)

       DPHDX=(DPHEW*B11P+DPHNS*B21P)*RVOLP               
       DPHDY=(DPHEW*B12P+DPHNS*B22P)*RVOLP               

       RETURN
       END              

         
C***********************************************************************
C 
      SUBROUTINE DPRS(II,FF,PHIW,PHIE,PHIS,PHIN)
C                         
C     COMPUTES PRESSURE OR PRESSURE CORRECTION  AT BOUNDARY NODES USING
C     LINEAR EXTRAPOLATION FROM INTERIOR POINTS
C                         
C***********************************************************************
C
      INCLUDE 'com2d'
C                                                 
      DIMENSION FF(NIJ)

      FXP   =F(ISFX+II)                                                 
      FXIP  =F(ISFX+II+1)                                               
      FYP   =F(ISFY+II)                                                 
      FYJP  =F(ISFY+II+NI)                                              

C     CALCULATING CELL FACE VALUES FOR CELLS WITHOUT ANY BOUNDARY FACE
C      
      PHIW=       FXP*FF(II)+(1.-FXP)*FF(II-1)  
      PHIE=       FXIP*FF(II+1)+(1.-FXIP)*FF(II)                  
      PHIS=       FYP*FF(II)+(1.-FYP)*FF(II-NI)                   
      PHIN=       FYJP*FF(II+NI)+(1.-FYJP)*FF(II)  
C               
      IF(KW(II).GT.0.AND.KW(II).LT.15)                            
     >  PHIW=(1.+FXIP)*FF(II)-FXIP*FF(II+1)     
      IF(KE(II).GT.0.AND.KE(II).LT.15)                            
     >  PHIE=(2.-FXP)*FF(II)-(1.-FXP)*FF(II-1)     
      IF(KS(II).GT.0.AND.KS(II).LT.15)                            
     >  PHIS=(1.+FYJP)*FF(II)-FYJP*FF(II+NI)    
      IF(KN(II).GT.0.AND.KN(II).LT.15)                            
     >  PHIN=(2.-FYP)*FF(II)-(1.-FYP)*FF(II-NI)   
C
C FOR PRESSURE BOUNDARY THE CONDITION THE PRESSURE/PRESSURE CORRECTION
C IS SET TO ZERO
C
      IF(KW(II).EQ.17) PHIW=0.
      IF(KE(II).EQ.17) PHIE=0.
      IF(KS(II).EQ.17) PHIS=0.
      IF(KN(II).EQ.17) PHIN=0.

 
      RETURN
      END               
C      
C***********************************************************************
C 
      SUBROUTINE CYC(IFLAG)       
C                                                                        
C TRANSFER OF FLOW VARIABLES AT PERIODIC BOUNDARY
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
                                                                        
      MCYC   =MCYCN*MFLAG1 + MCYCS*MFLAG2 + MCYCE*MFLAG3 + MCYCW*MFLAG4 


      IH =ISPP+(NN-1)*NIJ                                               
      DO IM=1,MCYC                                                 
         IP   = INCYN(IM)*MFLAG1 + INCYS(IM)*MFLAG2                      
     *       + INCYE(IM)*MFLAG3 + INCYW(IM)*MFLAG4                      

         IB  =IP +(NI * MFLAG1 - NI* MFLAG2 +  MFLAG3 -   MFLAG4)  

         IBB = IB + (NI*MFLAG1 -  NI*MFLAG2 +  MFLAG3  -   MFLAG4 )
         ICYC=(INCYN(IM)-(NJ-2-3)*NI)*MFLAG1
     >       +(INCYS(IM)+(NJ-2-3)*NI)*MFLAG2
     >       +(INCYE(IM)-(NI-2-3))*MFLAG3
     >       +(INCYW(IM)+(NI-2-3))*MFLAG4
         ICYCB = ICYC + (NI*MFLAG1 -  NI*MFLAG2 + MFLAG3 -   MFLAG4)

      F(IH+IB) = F(IH+ICYC)
        F(IH+IBB)=0.5*(F(IH+ICYC)+F(IH+ICYCB))

      IF(NN.EQ.1) THEN
       F(ISU+IB)=F(ISU+ICYC)
       F(ISU+IBB)=0.5*(F(ISU+ICYC)+F(ISU+ICYCB))
       F(ISV+IB)=F(ISV+ICYC)
       F(ISV+IBB)=0.5*(F(ISV+ICYC)+F(ISV+ICYCB))
       F(ISAPU+IB)=F(ISAPU+ICYC)
       F(ISAPV+IB)=F(ISAPV+ICYC)
       F(ISCW+IB)=F(ISCW+ICYC)
       F(ISCS+IB)=F(ISCS+ICYC)
      END IF
      END DO
C
      RETURN                                                            
      END                                                               
C
C**********************************************************************
C                            
         SUBROUTINE NORMCOEFF
C
C**********************************************************************
C
         include 'com2d'
C
         DO J=jsolv_b,jsolv_e
         DO I=isolv_b,isolv_e
         II=I+(J-1)*NI
        F(ISAE+II)=F(ISAE+II)/F(ISAP+II)
        F(ISAW+II)=F(ISAW+II)/F(ISAP+II)
        F(ISAN+II)=F(ISAN+II)/F(ISAP+II)
        F(ISAS+II)=F(ISAS+II)/F(ISAP+II)
        F(ISSU+II)=F(ISSU+II)/F(ISAP+II)
        F(ISSP+II)=F(ISSP+II)/F(ISAP+II)
        F(ISAP+II)=1.

        END DO
        END DO

        RETURN
        END
