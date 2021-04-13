       SUBROUTINE COEFF1(N) 
C                                       
C     CALCULATE CONVECTIVE-DIFFUSIVE COEFFICIENTS 
C     BY USING CENTRAL/UPWIND HYBRID DIFFERENCE SCHEME
C     DEFERRED CORRECTION PROCEDURE USED COUPLED TO PURE UPWIND 
C              
C***********************************************************************
C
      INCLUDE 'com2d'                                                 
C
      SQR(VALUE)=VALUE*VALUE                                            
                                                   
	 VALUE=0.
         PRNL=1.
         IF(N.EQ.11)PRNL=1./PRLAM
         PRNT=1./PR(N)
         
         ISVAR=ISPP+(N-1)*NI*NJ 

      DO J=jsolv_b,jsolv_e
      DO I=isolv_b,isolv_e
	 
      ii=i+(j-1)*ni
      IP=II+1                                                   
      JP=II+NI     
      IJP=ISVAR+II    
    
      SUH=0.
C
C     SPECIAL TREATMENT FOR DIFFUSION COEFFICIENT FOR SST MODEL
C      
      IF(TURBMOD.EQ.7.AND.N.EQ.(NVTE-NVPP+1)) THEN  
      PR_N=F1_SST(II)*PR1_TE+(1-F1_SST(II))*PR2_TE
      PRNT=1./PR_N
      END IF
 
      IF(TURBMOD.EQ.7.AND.N.EQ.(NVOM-NVPP+1)) THEN  
      PR_N=F1_SST(II)*PR1_OM+(1-F1_SST(II))*PR2_OM
      PRNT=1./PR_N
      END IF

C
C     DIFFERENCING ALONG X-DIRECTION
C
      FXP=F(ISFX+II)
      FXIP=F(ISFX+IP)
      D11W=SQR(F(ISB11W+II))+SQR(F(ISB12W+II))
      D11E=SQR(F(ISB11W+IP))+SQR(F(ISB12W+IP))
      VOLW=0.5*(F(ISVOLP+II)+F(ISVOLP+II-1))
      VOLE=0.5*(F(ISVOLP+IP)+F(ISVOLP+II))
      RW=FXP *F(ISR+II)+(1.-FXP )*F(ISR+II-1)
      RE=FXIP*F(ISR+IP)+(1.-FXIP)*F(ISR+II)
      VISTW=F(ISVISW+II)-VISCOS
      VISTE=F(ISVISW+IP)-VISCOS
      VISEFW=VISCOS*PRNL+VISTW*PRNT
      VISEFE=VISCOS*PRNL+VISTE*PRNT

      IF(N.EQ.10) THEN
        ANUW=FXP *F(ISNUT+II)+(1.-FXP )*F(ISNUT+II-1)
        ANUE=FXIP*F(ISNUT+IP)+(1.-FXIP)*F(ISNUT+II)
        VISEFW=(VISCOS+ANUW)*PRNT
        VISEFE=(VISCOS+ANUE)*PRNT
      ENDIF

      DW=VISEFW*RW*RW*D11W/(VOLW+SMALL)
      DE=VISEFE*RE*RE*D11E/(VOLE+SMALL)
      CW=F(ISCW+II)
      CE=F(ISCW+IP) 
C
C     UPWIND DIFFERENCE COEFFICIENTS
C
      AWU=DW+MAX(0.,CW)
      AEU=DE+MAX(0.,-CE)
   
      F(ISAW+II)=AWU
      F(ISAE+II)=AEU
C
C     CENTRAL DIFFERENCE COEFFICIENTS
C
      AWC=DW+CW*(1.-FXP)
      AEC=DE-CE*FXIP
C
C     WEST FACE
C
      PEW=ABS(CW/(DW+SMALL))
      IF(PEW.GT.2) FACT=1.
      IF(PEW.LE.2) FACT=0.

      F(ISAW+II)=FACT*AWU+(1.-FACT)*AWC
C
C     EAST FACE
C
      PEE=ABS(CE/(DE+SMALL))
      IF(PEE.GT.2) FACT=1.
      IF(PEE.LE.2) FACT=0.
      
      F(ISAE+II)=FACT*AEU+(1.-FACT)*AEC
C
C        DIFFERENCING SCHEMES IN Y-DIRECTION
C
        FYP=F(ISFY+II)
        FYJP=F(ISFY+JP)
        D22S=SQR(F(ISB21S+II))+SQR(F(ISB22S+II))
        D22N=SQR(F(ISB21S+JP))+SQR(F(ISB22S+JP))
        VOLS=0.5*(F(ISVOLP+II)+F(ISVOLP+II-NI))
        VOLN=0.5*(F(ISVOLP+JP)+F(ISVOLP+II))
        RS= FYP*F(ISR+II)+(1.-FYP )*F(ISR+II-NI)
        RN=FYJP*F(ISR+JP)+(1.-FYJP)*F(ISR+II)
        VISTS=F(ISVISS+II)-VISCOS
        VISTN=F(ISVISS+JP)-VISCOS
        VISEFS=VISCOS*PRNL+VISTS*PRNT
        VISEFN=VISCOS*PRNL+VISTN*PRNT

          IF(N.EQ.10) THEN

         ANUS=FYP *F(ISNUT+II)+(1.-FYP )*F(ISNUT+II-NI)
         ANUN=FYJP*F(ISNUT+JP)+(1.-FYJP)*F(ISNUT+II)
         VISEFS=(VISCOS+ANUS)*PRNT
         VISEFN=(VISCOS+ANUN)*PRNT
        ENDIF
      
        DS=VISEFS*RS*RS*D22S/(VOLS+SMALL)
        DN=VISEFN*RN*RN*D22N/(VOLN+SMALL)
        CS=F(ISCS+II)
        CN=F(ISCS+JP)    
C
C     UPWIND DIFFERENCE COEFFICIENTS
C
      ASU=DS+MAX(0.,CS)
      ANU=DN+MAX(0.,-CN)
C
C     CENTRAL DIFFERENCE COEFFICIENTS
C
      ASC=DS+CS*(1.-FYP)
      ANC=DN-CN*FYJP
C
C     SOUTH FACE
C
      PES=ABS(CS/(DS+SMALL))
      IF(PES.GT.2) FACT=1.
      IF(PES.LE.2) FACT=0.

      F(ISAS+II)=FACT*ASU+(1.-FACT)*ASC
C
C     NORTH FACE
C
      PEN=ABS(CN/(DN+SMALL))
      IF(PEN.GT.2) FACT=1.
      IF(PEN.LE.2) FACT=0.
      
      F(ISAN+II)=FACT*ANU+(1.-FACT)*ANC 
C
           DELQ=(CE-CW+CN-CS)*KBLK(II)
           SP_CONT  = - MAX(0., DELQ)
           SU_CONT  = - MIN(0.,DELQ)*F(ISVAR+II)
           F(ISSP+II)  = SP_CONT
           F(ISSU+II) = SU_CONT

      END DO       !     DO I=isolv_b,isolv_e
      END DO       !     DO J=jsolv_b,jsolv_e

      RETURN                          
      END                                                               
C
C**********************************************************************
C
      SUBROUTINE COEFF2(N)                 
C                                             
C     CALCULATE CONVECTIVE-DIFFUSIVE COEFFICIENTS
C     BY USING THE QUICK SCHEME 
C     DEFERRED CORRECTION PROCEDURE USED COUPLED TO PURE UPWIND     
C                                                     
C***********************************************************************
C                                                   
       INCLUDE 'com2d'                            
C                                           
      SQR(VALUE)=VALUE*VALUE  
             
	 VALUE=0.
         PRNL=1.
         IF(N.EQ.11)PRNL=1./PRLAM
         PRNT=1./PR(N)
                                            
C                                                  
         ISVAR=ISPP+(N-1)*NI*NJ 
 
                                    
                                                  
      DO J=JSOLV_B,JSOLV_E
      DO I=ISOLV_B,ISOLV_E                
                                                                       
      II=I+(J-1)*NI                              
      IP=II+1                                 
      JP=II+NI                                   
      IJP=ISVAR+II                                              

      IF(TURBMOD.EQ.7.AND.N.EQ.(NVTE-NVPP+1)) THEN  
C     PR_N=PR1_TE
      PR_N=F1_SST(II)*PR1_TE+(1-F1_SST(II))*PR2_TE
      PRNT=1./PR_N
      END IF
 
      IF(TURBMOD.EQ.7.AND.N.EQ.(NVOM-NVPP+1)) THEN  
C     PR_N=PR1_OM
      PR_N=F1_SST(II)*PR1_OM+(1-F1_SST(II))*PR2_OM
      PRNT=1./PR_N
      END IF


      SUH=0.                                          
C                                                     
C     DIFFERENCING SCHEMES IN X-DIRECTION                               
C                                                   
      FXP=F(ISFX+II)                              
      FXIP=F(ISFX+IP)                       
                                  
      D11W=SQR(F(ISB11W+II))+SQR(F(ISB12W+II))
      D11E=SQR(F(ISB11W+IP))+SQR(F(ISB12W+IP))
      VOLW=0.5*(F(ISVOLP+II)+F(ISVOLP+II-1))         
      VOLE=0.5*(F(ISVOLP+IP)+F(ISVOLP+II))           
      RW=FXP *F(ISR+II)+(1.-FXP )*F(ISR+II-1)      
      RE=FXIP*F(ISR+IP)+(1.-FXIP)*F(ISR+II)        
      VISTW=F(ISVISW+II)-VISCOS
      VISTE=F(ISVISW+IP)-VISCOS
      VISEFW=VISCOS*PRNL+VISTW*PRNT
      VISEFE=VISCOS*PRNL+VISTE*PRNT

      IF(N.EQ.10) THEN
        ANUW=FXP *F(ISNUT+II)+(1.-FXP )*F(ISNUT+II-1)
        ANUE=FXIP*F(ISNUT+IP)+(1.-FXIP)*F(ISNUT+II)
        VISEFW=(VISCOS+ANUW)*PRNT
        VISEFE=(VISCOS+ANUE)*PRNT
      ENDIF

      DW=VISEFW*RW*RW*D11W/(VOLW+SMALL)
      DE=VISEFE*RE*RE*D11E/(VOLE+SMALL)
      CW=F(ISCW+II)                             
      CE=F(ISCW+IP)                             

C     
C     QUICK SCHEME WITH DEFERRED CORRECTION
C     
C      PURE UPWIND COEFFICIENTS
C      
      AWU=DW+MAX(0.,CW)
      AEU=DE+MAX(0.,-CE)
      F(ISAW+II)=AWU
      F(ISAE+II)=AEU
C     
C      QUICK COEFFICIENTS
C      
       FP=F(IJP)
       FW=F(IJP-1)
       FE=F(IJP+1)
       FWW=F(IJP-2)
       FEE=F(IJP+2)
       
        IF(KW(II).EQ.0 .OR.KW(II) .EQ.15) THEN
        
        IF(CW.GT.0.) THEN
C
C  IF BLOCK FOR CW POSITIVE (OR) VELOCITY FROM W TO P
C
         FWU=FW
         DELU=DELEW(II-1)*KBLK(II-1)
         DELD=DELEW(II)*KBLK(II)
         DELUU=DELEW(II-2)*KBLK(II-2)
         
         CALL QUICK(DELUU,DELU,DELD,FWW,FW,FP,FWQ)
         
         ELSE
C
C  ELSE BLOCK FOR CW NEGATIVE (OR) VELOCITY FROM P TO W
C
         FWU=FP
         DELU=DELEW(II)*KBLK(II)
         DELD=DELEW(II-1)*KBLK(II-1)
         DELUU=DELEW(II+1)*KBLK(II+1)
         
         CALL QUICK(DELUU,DELU,DELD,FE,FP,FW,FWQ)
         
         END IF
C
C ADD SOURCE TERMS FOR FLUX DIFFERENCE
C
         SUH=SUH+CW*(FWQ-FWU)
        ELSE
          F(ISAW+II)=AWU
          END IF
          
        IF(KE(II).EQ.0 .OR. KE(II).EQ.15) THEN
                                                                        
         
        IF(CE.GT.0.) THEN
C
C  IF BLOCK FOR CE POSITIVE (OR) VELOCITY FROM P TO E
C
         FEU=FP
         DELU=DELEW(II)*KBLK(II)
         DELD=DELEW(II+1)*KBLK(II+1)
         DELUU=DELEW(II-1)*KBLK(II-1)
         
         CALL QUICK(DELUU,DELU,DELD,FW,FP,FE,FEQ)
         
         ELSE
C
C  ELSE BLOCK FOR CE NEGATIVE (OR) VELOCITY FROM E TO P
C
         FEU=FE
         DELU=DELEW(II+1)*KBLK(II+1)
         DELD=DELEW(II)*KBLK(II)
         DELUU=DELEW(II+2)*KBLK(II+2)
         
         CALL QUICK(DELUU,DELU,DELD,FEE,FE,FP,FEQ)
        
         END IF
C
C ADD SOURCE TERMS FOR FLUX DIFFERENCE
C
           SUH=SUH-CE*(FEQ-FEU)
         ELSE
         F(ISAE+II)=AEU                     
                                                                        
         END IF
C                        
C        DIFFERENCING SCHEMES IN Y-DIRECTION
C                      
        FYP=F(ISFY+II)   
        FYJP=F(ISFY+JP)   
        D22S=SQR(F(ISB21S+II))+SQR(F(ISB22S+II))
        D22N=SQR(F(ISB21S+JP))+SQR(F(ISB22S+JP)) 
        VOLS=0.5*(F(ISVOLP+II)+F(ISVOLP+II-NI))
        VOLN=0.5*(F(ISVOLP+JP)+F(ISVOLP+II))
        RS= FYP*F(ISR+II)+(1.-FYP )*F(ISR+II-NI)
        RN=FYJP*F(ISR+JP)+(1.-FYJP)*F(ISR+II)
        VISTS=F(ISVISS+II)-VISCOS
        VISTN=F(ISVISS+JP)-VISCOS
        VISEFS=VISCOS*PRNL+VISTS*PRNT
        VISEFN=VISCOS*PRNL+VISTN*PRNT

          IF(N.EQ.10) THEN

         ANUS=FYP *F(ISNUT+II)+(1.-FYP )*F(ISNUT+II-NI)
         ANUN=FYJP*F(ISNUT+JP)+(1.-FYJP)*F(ISNUT+II)
         VISEFS=(VISCOS+ANUS)*PRNT
         VISEFN=(VISCOS+ANUN)*PRNT
        ENDIF
      
        DS=VISEFS*RS*RS*D22S/(VOLS+SMALL)
        DN=VISEFN*RN*RN*D22N/(VOLN+SMALL)
        CS=F(ISCS+II)
        CN=F(ISCS+JP)                             

C       
C     QUICK SCHEME WITH DEFERRED CORRECTION
C                              
C      PURE UPWIND COEFFICIENTS
C                                           
          ASU=DS+MAX(0.,CS)                 
          ANU=DN+MAX(0.,-CN)
                         
        F(ISAS+II)=ASU    
        F(ISAN+II)=ANU   
C                                                
C      QUICK COEFFICIENTS                      
C                                           
       FP=F(IJP)                                
       FS=F(IJP-NI)                          
       FN=F(IJP+NI)                                   
       FSS=F(IJP-2*NI)                                
       FNN=F(IJP+2*NI)
                                                  
        IF(KS(II).EQ.0 .OR. KS(II).EQ.15) THEN
                                          
        IF(CS.GT.0.) THEN      
C
C  IF BLOCK FOR CS POSITIVE (OR) VELOCITY FROM S TO P
C
         FSU=FS        
         DELU=DELNS(II-NI)*KBLK(II-NI)                  
         DELD=DELNS(II)*KBLK(II)                              
         DELUU=DELNS(II-2*NI)*KBLK(II-2*NI)                        
                                            
         CALL QUICK(DELUU,DELU,DELD,FSS,FS,FP,FSQ)
                         
         ELSE             
C
C  ELSE BLOCK FOR CS NEGATIVE (OR) VELOCITY FROM P TO S
C
         FSU=FP                                  
         DELU=DELNS(II)*KBLK(II)                          
         DELD=DELNS(II-NI)*KBLK(II-NI)                     
         DELUU=DELNS(II+NI)*KBLK(II+NI)                 
                                                
         CALL QUICK(DELUU,DELU,DELD,FN,FP,FS,FSQ)
                                                      
         END IF                                       
C
C  ADD SOURCE TERMS FOR FLUX DIFFERENCE 
C
        SUH=SUH+CS*(FSQ-FSU)                
        ELSE                              
        F(ISAS+II)=ASU         
        END IF                           
                                           
        IF(KN(II).EQ.0 .OR. KN(II).EQ.15)THEN         
                                                     
        IF(CN.GT.0.) THEN                   
C
C  IF BLOCK FOR CN POSITIVE (OR) VELOCITY FROM P TO N
C
         FNU=FP                                    
         DELU=DELNS(II)*KBLK(II)  
         DELD=DELNS(II+NI)*KBLK(II+NI)
         DELUU=DELNS(II-NI)*KBLK(II-NI)                     
                                                 
         CALL QUICK(DELUU,DELU,DELD,FS,FP,FN,FNQ)
                                            
         ELSE                                   
C
C  ELSE BLOCK FOR CN NEGATIVE (OR) VELOCITY FROM N TO P
C
         FNU=FN                                   
         DELU=DELNS(II+NI)*KBLK(II+NI)                            
         DELD=DELNS(II)*KBLK(II)                               
         DELUU=DELNS(II+2*NI)*KBLK(II+2*NI)
                                                  
C
         CALL QUICK(DELUU,DELU,DELD,FNN,FN,FP,FNQ)
                                          
         END IF                
C
C  ADD SOURCE TERMS FOR FLUX DIFFERENCE 
C
        SUH=SUH-CN*(FNQ-FNU)                         
        ELSE                                         
        F(ISAN+II)=ANU                      
        END IF    
                          
C
           DELQ=(CE-CW+CN-CS)*KBLK(II)
           SP_CONT  = - MAX(0., DELQ)
           SU_CONT  = - MIN(0.,DELQ)*F(ISVAR+II)
           F(ISSP+II)  = SP_CONT
           F(ISSU+II)=  SU_CONT+FDEFER(N)*SUH 
c
        END DO       !     DO I=isolv_b,isolv_e
        END DO       !     DO J=jsolv_b,jsolv_e
                                                      
        RETURN                                        
      END                    
C                                                     
C***********************************************************************
C                                                   
      SUBROUTINE COEFF3(N)                 
C                                             
C     CALCULATE CONVECTIVE-DIFFUSIVE COEFFICIENTS
C     BY USING THE HLPA SCHEME 
C     DEFERRED CORRECTION PROCEDURE USED COUPLED TO PURE UPWIND     
C                                                     
C***********************************************************************
C                                                   
       INCLUDE 'com2d'                            
C                                           
      SQR(VALUE)=VALUE*VALUE               
                                            
	 VALUE=0.
         PRNL=1.
          IF(N.EQ.11)PRNL=1./PRLAM
         PRNT=1./PR(N)
C                                                  
         ISVAR=ISPP+(N-1)*NI*NJ 

         PRNINV=1./PR(N)                           
                                                  
      DO J=JSOLV_B,JSOLV_E
      DO I=ISOLV_B,ISOLV_E                
                                                                       
      II=I+(J-1)*NI                              
      IP=II+1                                 
      JP=II+NI                                   
      IJP=ISVAR+II                                              
      SUH=0.                                          
C                                                     
C
C     SPECIAL TREATMENT FOR DIFFUSION COEFFICIENT FOR SST MODEL
C      

      IF(TURBMOD.EQ.7.AND.N.EQ.(NVTE-NVPP+1)) THEN  
      PR_N=F1_SST(II)*PR1_TE+(1-F1_SST(II))*PR2_TE
      PRNT=1./PR_N
      END IF
 
      IF(TURBMOD.EQ.7.AND.N.EQ.(NVOM-NVPP+1)) THEN  
      PR_N=F1_SST(II)*PR1_OM+(1-F1_SST(II))*PR2_OM
      PRNT=1./PR_N
      END IF

C     DIFFERENCING SCHEMES IN X-DIRECTION                               
C                                                   
      FXP=F(ISFX+II)                              
      FXIP=F(ISFX+IP)                       
                                  
      D11W=SQR(F(ISB11W+II))+SQR(F(ISB12W+II))
      D11E=SQR(F(ISB11W+IP))+SQR(F(ISB12W+IP))
      VOLW=0.5*(F(ISVOLP+II)+F(ISVOLP+II-1))         
      VOLE=0.5*(F(ISVOLP+IP)+F(ISVOLP+II))           
      RW=FXP *F(ISR+II)+(1.-FXP )*F(ISR+II-1)      
      RE=FXIP*F(ISR+IP)+(1.-FXIP)*F(ISR+II)        
      VISTW=F(ISVISW+II)-VISCOS
      VISTE=F(ISVISW+IP)-VISCOS
      VISEFW=VISCOS*PRNL+VISTW*PRNT
      VISEFE=VISCOS*PRNL+VISTE*PRNT

      IF(N.EQ.10) THEN
        ANUW=FXP *F(ISNUT+II)+(1.-FXP )*F(ISNUT+II-1)
        ANUE=FXIP*F(ISNUT+IP)+(1.-FXIP)*F(ISNUT+II)
        VISEFW=(VISCOS+ANUW)*PRNT
        VISEFE=(VISCOS+ANUE)*PRNT
      ENDIF

      DW=VISEFW*RW*RW*D11W/(VOLW+SMALL)
      DE=VISEFE*RE*RE*D11E/(VOLE+SMALL)
      CW=F(ISCW+II)                             
      CE=F(ISCW+IP)                             

C     
C     HLPA SCHEME WITH DEFERRED CORRECTION
C     
C      PURE UPWIND COEFFICIENTS
C      
      AWU=DW+MAX(0.,CW)
      AEU=DE+MAX(0.,-CE)
      F(ISAW+II)=AWU
      F(ISAE+II)=AEU
C     
C      HLPA COEFFICIENTS
C      
       FP=F(IJP)
       FW=F(IJP-1)
       FE=F(IJP+1)
       FWW=F(IJP-2)
       FEE=F(IJP+2)
       
        IF(KW(II).EQ.0 .OR.KW(II) .EQ.15) THEN
        
        IF(CW.GT.0.) THEN
C
C  IF BLOCK FOR CW POSITIVE (OR) VELOCITY FROM W TO P
C
         FWU=FW
         DELU=DELEW(II-1)*KBLK(II-1)
         DELD=DELEW(II)*KBLK(II)
         DELUU=DELEW(II-2)*KBLK(II-2)
         
         CALL HLPA(DELUU,DELU,DELD,FWW,FW,FP,FWH)
         
         ELSE
C
C  ELSE BLOCK FOR CW NEGATIVE (OR) VELOCITY FROM P TO W
C
         FWU=FP
         DELU=DELEW(II)*KBLK(II)
         DELD=DELEW(II-1)*KBLK(II-1)
         DELUU=DELEW(II+1)*KBLK(II+1)
         
         CALL HLPA(DELUU,DELU,DELD,FE,FP,FW,FWH)
         
         END IF
C
C ADD SOURCE TERMS FOR FLUX DIFFERENCE
C
         SUH=SUH+CW*(FWH-FWU)
        
          END IF
          
        IF(KE(II).EQ.0 .OR. KE(II).EQ.15) THEN
                                                        
         
        IF(CE.GT.0.) THEN
C
C  IF BLOCK FOR CE POSITIVE (OR) VELOCITY FROM P TO E
C
         FEU=FP
         DELU=DELEW(II)*KBLK(II)
         DELD=DELEW(II+1)*KBLK(II+1)
         DELUU=DELEW(II-1)*KBLK(II-1)
         
         CALL HLPA(DELUU,DELU,DELD,FW,FP,FE,FEH)
         
         ELSE
C
C  ELSE BLOCK FOR CE NEGATIVE (OR) VELOCITY FROM E TO P
C
         FEU=FE
         DELU=DELEW(II+1)*KBLK(II+1)
         DELD=DELEW(II)*KBLK(II)
         DELUU=DELEW(II+2)*KBLK(II+2)
         
         CALL HLPA(DELUU,DELU,DELD,FEE,FE,FP,FEH)
        
         END IF
C
C ADD SOURCE TERMS FOR FLUX DIFFERENCE
C
           SUH=SUH-CE*(FEH-FEU)
                                                             
         END IF
C                        
C        DIFFERENCING SCHEMES IN Y-DIRECTION
C                      
        FYP=F(ISFY+II)   
        FYJP=F(ISFY+JP)   
        D22S=SQR(F(ISB21S+II))+SQR(F(ISB22S+II))
        D22N=SQR(F(ISB21S+JP))+SQR(F(ISB22S+JP)) 
        VOLS=0.5*(F(ISVOLP+II)+F(ISVOLP+II-NI))
        VOLN=0.5*(F(ISVOLP+JP)+F(ISVOLP+II))
        RS= FYP*F(ISR+II)+(1.-FYP )*F(ISR+II-NI)
        RN=FYJP*F(ISR+JP)+(1.-FYJP)*F(ISR+II)
        VISTS=F(ISVISS+II)-VISCOS
        VISTN=F(ISVISS+JP)-VISCOS
        VISEFS=VISCOS*PRNL+VISTS*PRNT
        VISEFN=VISCOS*PRNL+VISTN*PRNT

          IF(N.EQ.10) THEN

         ANUS=FYP *F(ISNUT+II)+(1.-FYP )*F(ISNUT+II-NI)
         ANUN=FYJP*F(ISNUT+JP)+(1.-FYJP)*F(ISNUT+II)
         VISEFS=(VISCOS+ANUS)*PRNT
         VISEFN=(VISCOS+ANUN)*PRNT
        ENDIF
      
        DS=VISEFS*RS*RS*D22S/(VOLS+SMALL)
        DN=VISEFN*RN*RN*D22N/(VOLN+SMALL)
        CS=F(ISCS+II)
        CN=F(ISCS+JP)                             

C       
C     HLPA SCHEME WITH DEFERRED CORRECTION
C                              
C      PURE UPWIND COEFFICIENTS
C                                           
          ASU=DS+MAX(0.,CS)                 
          ANU=DN+MAX(0.,-CN)
                         
        F(ISAS+II)=ASU    
        F(ISAN+II)=ANU   
C                                                
C      HLPA COEFFICIENTS                      
C                                           
       FP=F(IJP)                                
       FS=F(IJP-NI)                          
       FN=F(IJP+NI)                                   
       FSS=F(IJP-2*NI)                                
       FNN=F(IJP+2*NI)
                                                  
        IF(KS(II).EQ.0 .OR. KS(II).EQ.15) THEN
                                          
        IF(CS.GT.0.) THEN    
C
C  IF BLOCK FOR CS POSITIVE (OR) VELOCITY FROM S TO P
C
         FSU=FS        
         DELU=DELNS(II-NI)*KBLK(II-NI)                  
         DELD=DELNS(II)*KBLK(II)                              
         DELUU=DELNS(II-2*NI)*KBLK(II-2*NI)                        
                                            
         CALL HLPA(DELUU,DELU,DELD,FSS,FS,FP,FSH)
                         
         ELSE
C
C  ELSE BLOCK FOR CS NEGATIVE (OR) VELOCITY FROM P TO S
C
         FSU=FP                                  
         DELU=DELNS(II)*KBLK(II)                          
         DELD=DELNS(II-NI)*KBLK(II-NI)                     
         DELUU=DELNS(II+NI)*KBLK(II+NI)                 
                                                
         CALL HLPA(DELUU,DELU,DELD,FN,FP,FS,FSH)
                                                      
         END IF
C
C  ADD SOURCE TERMS FOR FLUX DIFFERENCE
C
        SUH=SUH+CS*(FSH-FSU)                
                
        END IF                           
                                           
        IF(KN(II).EQ.0 .OR. KN(II).EQ.15)THEN         
                                                     
        IF(CN.GT.0.) THEN   
C
C  IF BLOCK FOR CN POSITIVE (OR) VELOCITY FROM P TO N
C
         FNU=FP                                    
         DELU=DELNS(II)*KBLK(II)  
         DELD=DELNS(II+NI)*KBLK(II+NI)
         DELUU=DELNS(II-NI)*KBLK(II-NI)                     
                                                 
         CALL HLPA(DELUU,DELU,DELD,FS,FP,FN,FNH)
                                            
         ELSE 

C  ELSE BLOCK FOR CN NEGATIVE (OR) VELOCITY FROM N TO P
C
         FNU=FN                                   
         DELU=DELNS(II+NI)*KBLK(II+NI)                            
         DELD=DELNS(II)*KBLK(II)                               
         DELUU=DELNS(II+2*NI)*KBLK(II+2*NI)
                                                  
         CALL HLPA(DELUU,DELU,DELD,FNN,FN,FP,FNH)
                                          
         END IF  
C
C  ADD SOURCE TERMS FOR FLUX DIFFERENCE
C
        SUH=SUH-CN*(FNH-FNU)                         
                      
        END IF    
                   
C                            
       DELQ=(CE-CW+CN-CS)*KBLK(II)                         
       SP_CONT=-MAX(0.,DELQ)
       SU_CONT=-MIN(0.,DELQ)*F(ISVAR+II)*KBLK(II)
       F(ISSP+II)=SP_CONT
       F(ISSU+II)=SU_CONT+SUH*FDEFER(N)                        

               
c
        END DO       !     DO I=isolv_b,isolv_e
        END DO       !     DO J=jsolv_b,jsolv_e
                                                      
        RETURN                                        
      END                    
C                                                     
C***********************************************************************
C                                                   
       SUBROUTINE COEFF4(N) 
C                                       
C     CALCULATE CONVECTIVE-DIFFUSIVE COEFFICIENTS 
C     BY USING CENTRAL DIFFERENCE SCHEME
C     DEFERRED CORRECTION PROCEDURE USED COUPLED TO PURE UPWIND 
C              
C***********************************************************************
C
      INCLUDE 'com2d'                                                 
C
      SQR(VALUE)=VALUE*VALUE                                            
                                                   
	 VALUE=0.
         PRNL=1.
          IF(N.EQ.11)PRNL=1./PRLAM
         PRNT=1./PR(N)
         
         ISVAR=ISPP+(N-1)*NI*NJ   
                                   
      DO J=jsolv_b,jsolv_e
      DO I=isolv_b,isolv_e
	 
      ii=i+(j-1)*ni
      IP=II+1                                                   
      JP=II+NI     
      IJP=ISVAR+II    
    
      SUH=0.
C
C     SPECIAL TREATMENT FOR DIFFUSION COEFFICIENT FOR SST MODEL
C      

      IF(TURBMOD.EQ.7.AND.N.EQ.(NVTE-NVPP+1)) THEN  
      PR_N=F1_SST(II)*PR1_TE+(1-F1_SST(II))*PR2_TE
      PRNT=1./PR_N
      END IF
 
      IF(TURBMOD.EQ.7.AND.N.EQ.(NVOM-NVPP+1)) THEN  
      PR_N=F1_SST(II)*PR1_OM+(1-F1_SST(II))*PR2_OM
      PRNT=1./PR_N
      END IF
C
C     DIFFERENCING ALONG X-DIRECTION
C
      FXP=F(ISFX+II)
      FXIP=F(ISFX+IP)
      D11W=SQR(F(ISB11W+II))+SQR(F(ISB12W+II))
      D11E=SQR(F(ISB11W+IP))+SQR(F(ISB12W+IP))
      VOLW=0.5*(F(ISVOLP+II)+F(ISVOLP+II-1))
      VOLE=0.5*(F(ISVOLP+IP)+F(ISVOLP+II))
      RW=FXP *F(ISR+II)+(1.-FXP )*F(ISR+II-1)
      RE=FXIP*F(ISR+IP)+(1.-FXIP)*F(ISR+II)
      VISTW=F(ISVISW+II)-VISCOS
      VISTE=F(ISVISW+IP)-VISCOS
      VISEFW=VISCOS*PRNL+VISTW*PRNT
      VISEFE=VISCOS*PRNL+VISTE*PRNT

      IF(N.EQ.10) THEN
        ANUW=FXP *F(ISNUT+II)+(1.-FXP )*F(ISNUT+II-1)
        ANUE=FXIP*F(ISNUT+IP)+(1.-FXIP)*F(ISNUT+II)
        VISEFW=(VISCOS+ANUW)*PRNT
        VISEFE=(VISCOS+ANUE)*PRNT
      ENDIF

      DW=VISEFW*RW*RW*D11W/(VOLW+SMALL)
      DE=VISEFE*RE*RE*D11E/(VOLE+SMALL)
      CW=F(ISCW+II)
      CE=F(ISCW+IP) 
C
C     UPWIND DIFFERENCE COEFFICIENTS
C
      AWU=DW+MAX(0.,CW)
      AEU=DE+MAX(0.,-CE)
      F(ISAW+II)=AWU
      F(ISAE+II)=AEU
C
C     WEST FACE
C
      FP=F(IJP)
      FW=F(IJP-1)

      IF(CW.GT.0) FWU=FW
      IF(CW.LE.0) FWU=FP
      
      DELU=DELEW(II-1)*KBLK(II-1)
      DELD=DELEW(II)*KBLK(II)    
C           
C      CENTRAL DIFFERENCE COEFFICIENTS
C
      CALL CENTRAL(DELU,DELD,FW,FP,FWC)
C
C ADD SOURCE TERMS FOR FLUX DIFFERENCE
C
      IF (KBLK(II) .NE. 0)SUH=SUH+CW*(FWC-FWU)
      
C
C     EAST FACE
C
      
      FP=F(IJP)
      FE=F(IJP+1)

      IF(CE.GT.0) FEU=FP
      IF(CE.LE.0) FEU=FE
      
      DELU=DELEW(II)*KBLK(II)
      DELD=DELEW(II+1)*KBLK(II+1)    
C           
C      CENTRAL DIFFERENCE COEFFICIENTS
C
      CALL CENTRAL(DELU,DELD,FP,FE,FEC)
C
C ADD SOURCE TERMS FOR FLUX DIFFERENCE
C
      IF (KBLK(II) .NE. 0)SUH=SUH-CE*(FEC-FEU)
C
C        DIFFERENCING SCHEMES IN Y-DIRECTION
C
        FYP=F(ISFY+II)
        FYJP=F(ISFY+JP)
        D22S=SQR(F(ISB21S+II))+SQR(F(ISB22S+II))
        D22N=SQR(F(ISB21S+JP))+SQR(F(ISB22S+JP))
        VOLS=0.5*(F(ISVOLP+II)+F(ISVOLP+II-NI))
        VOLN=0.5*(F(ISVOLP+JP)+F(ISVOLP+II))
        RS= FYP*F(ISR+II)+(1.-FYP )*F(ISR+II-NI)
        RN=FYJP*F(ISR+JP)+(1.-FYJP)*F(ISR+II)
        VISTS=F(ISVISS+II)-VISCOS
        VISTN=F(ISVISS+JP)-VISCOS
        VISEFS=VISCOS*PRNL+VISTS*PRNT
        VISEFN=VISCOS*PRNL+VISTN*PRNT

          IF(N.EQ.10) THEN

         ANUS=FYP *F(ISNUT+II)+(1.-FYP )*F(ISNUT+II-NI)
         ANUN=FYJP*F(ISNUT+JP)+(1.-FYJP)*F(ISNUT+II)
         VISEFS=(VISCOS+ANUS)*PRNT
         VISEFN=(VISCOS+ANUN)*PRNT
        ENDIF
      
        DS=VISEFS*RS*RS*D22S/(VOLS+SMALL)
        DN=VISEFN*RN*RN*D22N/(VOLN+SMALL)
        CS=F(ISCS+II)
        CN=F(ISCS+JP)    
C
C     UPWIND DIFFERENCE COEFFICIENTS
C
      ASU=DS+MAX(0.,CS)
      ANU=DN+MAX(0.,-CN)

      F(ISAS+II)=ASU
      F(ISAN+II)=ANU
C
C     SOUTH FACE
C

      FP=F(IJP)
      FS=F(IJP-NI)

      IF(CS.GT.0) FSU=FS
      IF(CS.LE.0) FSU=FP

      DELU=DELNS(II-NI)*KBLK(II-NI)
      DELD=DELNS(II)*KBLK(II)
C
C      CENTRAL DIFFERENCE COEFFICIENTS
C
      CALL CENTRAL(DELU,DELD,FS,FP,FSC)
C
C ADD SOURCE TERMS FOR FLUX DIFFERENCE
C
      IF (KBLK(II) .NE. 0)SUH=SUH+CS*(FSC-FSU)

C  
C     NORTH FACE
C
      FP=F(IJP)
      FN=F(IJP+NI)

      IF(CN.GT.0) FNU=FP
      IF(CN.LE.0) FNU=FN

      DELU=DELNS(II)*KBLK(II)
      DELD=DELNS(II+NI)*KBLK(II+NI)
C
C      CENTRAL DIFFERENCE COEFFICIENTS
C
      CALL CENTRAL(DELU,DELD,FP,FN,FNC)
C
C ADD SOURCE TERMS FOR FLUX DIFFERENCE
C
      IF (KBLK(II) .NE. 0)SUH=SUH-CN*(FNC-FNU)

C                            
      DELQ=(CE-CW+CN-CS)*KBLK(II)                         
      SP_CONT=-MAX(0.,DELQ)
      SU_CONT=-MIN(0.,DELQ)*F(ISVAR+II)*KBLK(II)
      F(ISSP+II)=SP_CONT
      F(ISSU+II)=SU_CONT+SUH*FDEFER(N)                        
C
        END DO       !     DO I=isolv_b,isolv_e
        END DO       !     DO J=jsolv_b,jsolv_e
                                
        RETURN                                        
      END   
C
C**********************************************************************
C
      SUBROUTINE COEFF5(N)                 
C                                             
C     CALCULATE CONVECTIVE-DIFFUSIVE COEFFICIENTS
C     BY USING THE SECOND ORDER UPWIND SCHEME 
C     DEFERRED CORRECTION PROCEDURE USED COUPLED TO PURE UPWIND     
C                                                     
C***********************************************************************
C                                                   
       INCLUDE 'com2d'                            
C                                           
      SQR(VALUE)=VALUE*VALUE  
             
	 VALUE=0.
         PRNL=1.
          IF(N.EQ.11)PRNL=1./PRLAM
         PRNT=1./PR(N)
                                            
C                                                  
         ISVAR=ISPP+(N-1)*NI*NJ 
                                     
                                                  
      DO J=JSOLV_B,JSOLV_E
      DO I=ISOLV_B,ISOLV_E                
                                                                       
      II=I+(J-1)*NI                              
      IP=II+1                                 
      JP=II+NI                                   
      IJP=ISVAR+II                                              
      SUH=0.                                          
C
C     SPECIAL TREATMENT FOR DIFFUSION COEFFICIENT FOR SST MODEL
C      

      IF(TURBMOD.EQ.7.AND.N.EQ.(NVTE-NVPP+1)) THEN  
      PR_N=F1_SST(II)*PR1_TE+(1-F1_SST(II))*PR2_TE
      PRNT=1./PR_N
      END IF
 
      IF(TURBMOD.EQ.7.AND.N.EQ.(NVOM-NVPP+1)) THEN  
      PR_N=F1_SST(II)*PR1_OM+(1-F1_SST(II))*PR2_OM
      PRNT=1./PR_N
      END IF
C                                                     
C     DIFFERENCING SCHEMES IN X-DIRECTION                               
C                                                   
      FXP=F(ISFX+II)                              
      FXIP=F(ISFX+IP)                       
                                  
      D11W=SQR(F(ISB11W+II))+SQR(F(ISB12W+II))
      D11E=SQR(F(ISB11W+IP))+SQR(F(ISB12W+IP))
      VOLW=0.5*(F(ISVOLP+II)+F(ISVOLP+II-1))         
      VOLE=0.5*(F(ISVOLP+IP)+F(ISVOLP+II))           
      RW=FXP *F(ISR+II)+(1.-FXP )*F(ISR+II-1)      
      RE=FXIP*F(ISR+IP)+(1.-FXIP)*F(ISR+II)        
      VISTW=F(ISVISW+II)-VISCOS
      VISTE=F(ISVISW+IP)-VISCOS
      VISEFW=VISCOS*PRNL+VISTW*PRNT
      VISEFE=VISCOS*PRNL+VISTE*PRNT

      IF(N.EQ.10) THEN
        ANUW=FXP *F(ISNUT+II)+(1.-FXP )*F(ISNUT+II-1)
        ANUE=FXIP*F(ISNUT+IP)+(1.-FXIP)*F(ISNUT+II)
        VISEFW=(VISCOS+ANUW)*PRNT
        VISEFE=(VISCOS+ANUE)*PRNT
      ENDIF

      DW=VISEFW*RW*RW*D11W/(VOLW+SMALL)
      DE=VISEFE*RE*RE*D11E/(VOLE+SMALL)
      CW=F(ISCW+II)                             
      CE=F(ISCW+IP)                             
C     
C     SECOND ORDER UPWIND  SCHEME WITH DEFERRED CORRECTION
C     
C      PURE UPWIND COEFFICIENTS
C      
      AWU=DW+MAX(0.,CW)
      AEU=DE+MAX(0.,-CE)
      F(ISAW+II)=AWU
      F(ISAE+II)=AEU
C     
C      LINEAR UPWINDIG  COEFFICIENTS
C      
       FP=F(IJP)
       FW=F(IJP-1)
       FE=F(IJP+1)
       FWW=F(IJP-2)
       FEE=F(IJP+2)
       
        IF(KW(II).EQ.0 .OR.KW(II) .EQ.15) THEN
        
        IF(CW.GT.0.) THEN
C
C  IF BLOCK FOR CW POSITIVE (OR) VELOCITY FROM W TO P
C
         FWU=FW
         DELU=DELEW(II-1)*KBLK(II-1)
         DELUU=DELEW(II-2)*KBLK(II-2)
         
         CALL LUINTP(DELUU,DELU,FWW,FW,FWLU)
         
         ELSE
C
C  ELSE BLOCK FOR CW NEGATIVE (OR) VELOCITY FROM P TO W
C
         FWU=FP
         DELU=DELEW(II)*KBLK(II)
         DELUU=DELEW(II+1)*KBLK(II+1)
         
         CALL LUINTP(DELUU,DELU,FE,FP,FWLU)
         
         END IF
C
C ADD SOURCE TERMS FOR FLUX DIFFERENCE
C
         SUH=SUH+CW*(FWLU-FWU)
        ELSE
          F(ISAW+II)=AWU
          END IF
          
        IF(KE(II).EQ.0 .OR. KE(II).EQ.15) THEN
                                                                        
         
        IF(CE.GT.0.) THEN
C
C  IF BLOCK FOR CE POSITIVE (OR) VELOCITY FROM P TO E
C
         FEU=FP
         DELU=DELEW(II)*KBLK(II)
         DELUU=DELEW(II-1)*KBLK(II-1)
         
         CALL LUINTP(DELUU,DELU,FW,FP,FELU)
         
         ELSE
C
C  ELSE BLOCK FOR CE NEGATIVE (OR) VELOCITY FROM E TO P
C
         FEU=FE
         DELU=DELEW(II+1)*KBLK(II+1)
         DELUU=DELEW(II+2)*KBLK(II+2)
         
         CALL LUINTP(DELUU,DELU,FEE,FE,FELU)
        
         END IF
C
C ADD SOURCE TERMS FOR FLUX DIFFERENCE
C
           SUH=SUH-CE*(FELU-FEU)
         ELSE
         F(ISAE+II)=AEU                     
                                                                        
         END IF
C                        
C        DIFFERENCING SCHEMES IN Y-DIRECTION
C                      
        FYP=F(ISFY+II)   
        FYJP=F(ISFY+JP)   
        D22S=SQR(F(ISB21S+II))+SQR(F(ISB22S+II))
        D22N=SQR(F(ISB21S+JP))+SQR(F(ISB22S+JP)) 
        VOLS=0.5*(F(ISVOLP+II)+F(ISVOLP+II-NI))
        VOLN=0.5*(F(ISVOLP+JP)+F(ISVOLP+II))
        RS= FYP*F(ISR+II)+(1.-FYP )*F(ISR+II-NI)
        RN=FYJP*F(ISR+JP)+(1.-FYJP)*F(ISR+II)
        VISTS=F(ISVISS+II)-VISCOS
        VISTN=F(ISVISS+JP)-VISCOS
        VISEFS=VISCOS*PRNL+VISTS*PRNT
        VISEFN=VISCOS*PRNL+VISTN*PRNT

          IF(N.EQ.10) THEN

         ANUS=FYP *F(ISNUT+II)+(1.-FYP )*F(ISNUT+II-NI)
         ANUN=FYJP*F(ISNUT+JP)+(1.-FYJP)*F(ISNUT+II)
         VISEFS=(VISCOS+ANUS)*PRNT
         VISEFN=(VISCOS+ANUN)*PRNT
        ENDIF
      
        DS=VISEFS*RS*RS*D22S/(VOLS+SMALL)
        DN=VISEFN*RN*RN*D22N/(VOLN+SMALL)
        CS=F(ISCS+II)
        CN=F(ISCS+JP)  
C        
C     QUIK SCHEME WITH DEFERRED CORRECTION
C                              
C      PURE UPWIND COEFFICIENTS
C                                           
          ASU=DS+MAX(0.,CS)                 
          ANU=DN+MAX(0.,-CN)
                         
        F(ISAS+II)=ASU    
        F(ISAN+II)=ANU   
C                                                
C      QUICK COEFFICIENTS                      
C                                           
       FP=F(IJP)                                
       FS=F(IJP-NI)                          
       FN=F(IJP+NI)                                   
       FSS=F(IJP-2*NI)                                
       FNN=F(IJP+2*NI)
                                                  
        IF(KS(II).EQ.0 .OR. KS(II).EQ.15) THEN
                                          
        IF(CS.GT.0.) THEN  
C
C  IF BLOCK FOR CS POSITIVE (OR) VELOCITY FROM S TO P
C
         FSU=FS        
         DELU=DELNS(II-NI)*KBLK(II-NI)                  
         DELUU=DELNS(II-2*NI)*KBLK(II-2*NI)                        
                                            
         CALL LUINTP(DELUU,DELU,FSS,FS,FSLU)
                         
         ELSE
C
C  ELSE BLOCK FOR CS NEGATIVE (OR) VELOCITY FROM P TO S
C
         FSU=FP                                  
         DELU=DELNS(II)*KBLK(II)                          
         DELUU=DELNS(II+NI)*KBLK(II+NI)                 
                                                
         CALL LUINTP(DELUU,DELU,FN,FP,FSLU)
                                                      
         END IF
C
C  ADD SOURCE TERMS FOR FLUX DIFFERENCE
C
        SUH=SUH+CS*(FSLU-FSU)                
        ELSE                              
        F(ISAS+II)=ASU         
        END IF                           
                                           
        IF(KN(II).EQ.0 .OR. KN(II).EQ.15)THEN         
                                                     
        IF(CN.GT.0.) THEN    
C
C  IF BLOCK FOR CN POSITIVE (OR) VELOCITY FROM P TO N
C
         FNU=FP                                    
         DELU=DELNS(II)*KBLK(II)  
         DELUU=DELNS(II-NI)*KBLK(II-NI)                     
                                                 
         CALL LUINTP(DELUU,DELU,FS,FP,FNLU)
                                            
         ELSE
C
C  ELSE BLOCK FOR CN NEGATIVE (OR) VELOCITY FROM N TO P
C
         FNU=FN                                   
         DELU=DELNS(II+NI)*KBLK(II+NI)                            
         DELUU=DELNS(II+2*NI)*KBLK(II+2*NI)
                                                  
         CALL LUINTP(DELUU,DELU,FNN,FN,FNLU)
                                          
         END IF
C
C  ADD SOURCE TERMS FOR FLUX DIFFERENCE
C
        SUH=SUH-CN*(FNLU-FNU)                         
        ELSE                                         
        F(ISAN+II)=ANU                      
        END IF    
                          
C                            
      DELQ=(CE-CW+CN-CS)*KBLK(II)                         
      SP_CONT=-MAX(0.,DELQ)
      SU_CONT=-MIN(0.,DELQ)*F(ISVAR+II)*KBLK(II)
      F(ISSP+II)=SP_CONT
      F(ISSU+II)=SU_CONT+SUH*FDEFER(N)                        
c
        END DO       !     DO I=isolv_b,isolv_e
        END DO       !     DO J=jsolv_b,jsolv_e
                                                      
        RETURN                                        
      END                    
C                                                     
C***********************************************************************
C                                                                       
      SUBROUTINE QUICK(DELUU,DELU,DELD,PHUU,PHU,PHD,PHUD)
C                                          
C     THIS SUBROUTINE IS USED TO COMPUTE CELL FACE VALUE (BETWEEN U AND D) OF
C     A VARIABLE USING QUADRATIC INTERPOLATION FROM VALUES AT THREE CONSECUTIVE
C     NODES UU,U,D                                   
C
C****************************************************************************  
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
C****************************************************************************  
C
      SUBROUTINE CENTRAL(DELU,DELD,PHU,PHD,PHUD)
C                                          
C     THIS SUBROUTINE IS USED TO COMPUTE CELL FACE VALUE (BETWEEN U AND D) OF
C     A VARIABLE USING LINEAR INTERPOLATION FROM VALUES AT TWO CONSECUTIVE
C     NODES U AND D                                   
C
C******************************************************************************
C
      SMALL = 1.0e-20
                                                  
      FACT=DELD/(DELU+DELD+SMALL)
      PHUD=FACT*PHU+(1.-FACT)*PHD      

      RETURN                                      
      END
C
C*****************************************************************************
C                                            
      SUBROUTINE HLPA(DELUU,DELU,DELD,PHUU,PHU,PHD,PHUD)
C                                          
C   THIS SUBROUTINE IS USED TO COMPUTE CELL FACE VALUE (BETWEEN U AND D) OF
C   A VARIABLE USING HYBRID LINEAR PARABOLIC INTERPOLATION FROM VALUES AT THREE 
C   CONSECUTIVE NODES UU,U,D                                   
C
C******************************************************************************
C
      SMALL=1.E-20        
      PHUHAT=(PHU-PHUU)/(PHD-PHUU+SMALL)
      GAMMA=0.
      IF(ABS(PHUHAT-0.5).lt.0.5) GAMMA=1.
      PHUD=PHU+GAMMA*(PHD-PHU)*PHUHAT
     
      RETURN                                      
      END
C
C****************************************************************************
C                                                                       
      SUBROUTINE LUINTP(DELUU,DELU,PHUU,PHU,PHUD)
C                                          
C   THIS SUBROUTINE IS USED TO COMPUTE CELL FACE VALUE (BETWEEN U AND D) OF
C   A VARIABLE USING LINEAR INTERPOLATION FROM VALUES AT TWO CONSECUTIVE  
C   NODES UU( UP-UPSTREAM) AND U(UPSTREAM)                                  
C
C******************************************************************************
C
        Q=DELU/(DELU+DELUU)
C
	PHUD=(1.+Q)*PHU-Q*PHUU
C
        RETURN
        END
C                              
C***********************************************************************
C
      SUBROUTINE COEFF_CHECK(N)     
C                              
C CHECKING THE POSITIVITY OF COEFFICIENTS 
C                              
C***********************************************************************
C
      INCLUDE 'com2d'
C
      SQR(VALUE)=VALUE*VALUE                                            
                                                   
	 VALUE=0.
         PRNINV=1./PR(N) 

C
C Checking the Positivity 
C
         DO J=jsolv_b,jsolv_e
         DO I=isolv_b,isolv_e
         II=I+(J-1)*NI

           IP=II+1                                                   
           JP=II+NI     
C
C     DIFFERENCING ALONG X-DIRECTION
C
      FXP=F(ISFX+II)
      FXIP=F(ISFX+IP)
      D11W=SQR(F(ISB11W+II))+SQR(F(ISB12W+II))
      D11E=SQR(F(ISB11W+IP))+SQR(F(ISB12W+IP))
      VOLW=0.5*(F(ISVOLP+II)+F(ISVOLP+II-1))
      VOLE=0.5*(F(ISVOLP+IP)+F(ISVOLP+II))
      RW=FXP *F(ISR+II)+(1.-FXP )*F(ISR+II-1)
      RE=FXIP*F(ISR+IP)+(1.-FXIP)*F(ISR+II)
      VISTW=F(ISVISW+II)-VISCOS
      VISTE=F(ISVISW+IP)-VISCOS
      VISEFW=VISCOS+VISTW*PRNINV
      VISEFE=VISCOS+VISTE*PRNINV


      DW=VISEFW*RW*RW*D11W/(VOLW+SMALL)
      DE=VISEFE*RE*RE*D11E/(VOLE+SMALL)
      CW=F(ISCW+II)
      CE=F(ISCW+IP) 
      PECW=ABS(CW/(DW+SMALL))
      PECE=ABS(CE/(DE+SMALL))
C
C  FOR WEST 
C
         IF(F(ISAW+II).LT.0.AND.PECW.GT.2) then

C     UPWIND DIFFERENCE COEFFICIENTS
C
      AWU=DW+MAX(0.,CW)
      AEU=DE+MAX(0.,-CE)
   
C     CENTRAL DIFFERENCE COEFFICIENTS
C
      AWC=DW+CW*(1.-FXP)
      AEC=DE-CE*FXIP
C
C     WEST FACE
C
      PEW=ABS(CW/(DW+SMALL))
      IF(PEW.GT.2) FACT=1.
      IF(PEW.LE.2) FACT=0.

       AWII=FACT*AWU+(1.-FACT)*AWC
c        WRITE(37,*)'WEST',N,J,I,PECW,F(ISAW+II)
         END IF  
C
C  FOR EAST 
C
         IF(F(ISAE+II).LT.0.AND.PECE.GT.2) then
C        WRITE(37,*)'EAST',N,J,I,PECE,F(ISAE+II)
         END IF  
C
C        DIFFERENCING SCHEMES IN Y-DIRECTION
C
        FYP=F(ISFY+II)
        FYJP=F(ISFY+JP)
        D22S=SQR(F(ISB21S+II))+SQR(F(ISB22S+II))
        D22N=SQR(F(ISB21S+JP))+SQR(F(ISB22S+JP))
        VOLS=0.5*(F(ISVOLP+II)+F(ISVOLP+II-NI))
        VOLN=0.5*(F(ISVOLP+JP)+F(ISVOLP+II))
        RS= FYP*F(ISR+II)+(1.-FYP )*F(ISR+II-NI)
        RN=FYJP*F(ISR+JP)+(1.-FYJP)*F(ISR+II)
        VISTS=F(ISVISS+II)-VISCOS
        VISTN=F(ISVISS+JP)-VISCOS
        VISEFS=VISCOS+VISTS*PRNINV
        VISEFN=VISCOS+VISTN*PRNINV  

      
        DS=VISEFS*RS*RS*D22S/(VOLS+SMALL)
        DN=VISEFN*RN*RN*D22N/(VOLN+SMALL)
        CS=F(ISCS+II)
        CN=F(ISCS+JP)    
        PECS=ABS(CS/(DS+SMALL))
        PECN=ABS(CN/(DN+SMALL))
C
C  FOR SOUTH 
C
         IF(F(ISAS+II).LT.0.AND.PECS.GT.2) then
C        WRITE(37,*)'SOUTH',N,J,I,PECS,F(ISAS+II)
         END IF  
C
C  FOR NORTH 
C
         IF(F(ISAN+II).LT.0.AND.PECN.GT.2) then
C        WRITE(37,*)'NORTH',N,J,I,PECN,F(ISAN+II)
         END IF  

         END DO
         END DO
C
         RETURN
         END
C
C *****************************************************************************          
