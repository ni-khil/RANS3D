      SUBROUTINE SORCE5
C
C     CALCULATION OF SOURCE TERMS FOR K EQUATION                           
C
C************************************************************************
C                                                      
      INCLUDE 'com2d'                                                
C
C    CROSS DERIVATIVE DIFFUSIVE FLUXES
C

      CALL CDFLUX(ISTE,5) 

C
C   NEAR WALL DAMPING FUNCTIONS FOR LOW RE CHIEN MODEL ( TURBMOD=3)

	if (turbmod .eq. 3)  CALL CHIEN_MOD
         
      DO J=jsolv_b,jsolv_e
      DO I=isolv_b,isolv_e
       
           II=I+(J-1)*NI
C
C     VELOCITY DIFFERENCES FOR GRADIENT EVALUATION AT CENTRAL NODE      
C
         FXP=F(ISFX+II)                                                 
         FYP=F(ISFY+II)                                                 
         FXM=1.-FXP                                                     
         FYM=1.-FYP                                                     
        FXIP=F(ISFX+II+1)                                               
        FYJP=F(ISFY+II+NI)                                              
        FXIM=1.-FXIP                                                    
        FYJM=1.-FYJP    
        RP2=F(ISR+II)*F(ISR+II)+SMALL                             
        VISTUR=F(ISVIS+II)-VISCOS+SMALL                              
        VOLP=F(ISVOLP+II)+SMALL                                    
C
C EVALUATE THE VELOCITY GRADIENTS FOR COMPUTATION OF PRODUCTION TERM
C
        CALL DIFF(II,ISU,DUDX,DUDY)
        CALL DIFF(II,ISV,DVDX,DVDY)

        VORTIC(II)=ABS(DUDY-DVDX)

        DWDX=0.
        DWDY=0.

        IF(SWIRL) THEN 

          CALL DIFF(II,ISW,DWDX,DWDY)

          IF(AXISYM_X) DWDX=DWDX-F(ISW+II)/RP2
          IF(AXISYM_Y) DWDY=DWDY-F(ISW+II)/RP2

         END IF 
C                                                                        
C    PRODUCTION OF TURBULENT KINETIC ENERGY                            
C
	STRAIN=(2*DUDX*DUDX+2*DVDY*DVDY+(DUDY+DVDX)**2+DWDX*DWDX+DWDY*DWDY)
      
       IF(AXISYM_X) STRAIN=STRAIN+2.*F(ISV+II)*F(ISV+II)/RP2   
       IF(AXISYM_Y) STRAIN=STRAIN+2.*F(ISU+II)*F(ISU+II)/RP2   

          STRN(II)=SQRT(STRAIN)
          
 	  GENR=VISTUR*STRAIN 
C
C COMPUTATION OF VORTICITY MAGNITUDE FOR KATO-LAUNDER MODIFICATION 
C
	IF (KATO) GENR=VISTUR*SQRT(STRAIN)*VORTIC(II)
 	
C
C  ADDITIUONAL SOURCE TERMS DUE TO BOUYANCY GENERATED TURBULENCE
C
        GK(II)=0.0

        IF(BUOY_SORCE) THEN
C
C  EVALUATE THE LOCAL VERTICAL GRADIENT OF TEMPERATURE
C
        CALL DIFF(II,ISS,DTDX,DTDY) 
C
      GK(II)=-RAYLEIGH*PRLAM*DTDY/PR(11)

        END IF


        PROD(II)=GENR  

C                    
C      RETAINING THE SOURCE TERMS BEFORE UPDATING THE SOURCE DUE TO GENERATION
C      TO BE USED FOR THE NEAR WALL CELLS IN THE STANDARD K_EPSILON MODEL
C
       BESUC(II)=F(ISSU+II)                                             
       BESPC(II)=F(ISSP+II) 
       VTERM=VOLP*KBLK(II)                                          
 
       F(ISSU+II)=F(ISSU+II)+(PROD(II)+GK(II))*VTERM                                           
C
      IF(TURBMOD.LE.3) THEN
      TLFACT=ALMU(II)/ALEPS(II)
      EPSBYK=CMU*FMU(II)*F(ISDEN+II)*F(ISTE+II)*TLFACT/VISTUR
      GAMMA_K=F(ISDEN+II)*EPSBYK*VTERM
      END IF
C
       IF (TURBMOD.EQ.4) THEN
       GAMMA_K=F(ISDEN+II)*F(ISOM+II)*BETA_STR*VTERM
       END IF
C
       IF(TURBMOD.EQ.5) THEN
C      EPSBYK=1./TFAC(I,J)
       EPSBYK=F(ISED+II)/(F(ISTE+II)+small)
       GAMMA_K=F(ISDEN+II)*EPSBYK*VTERM
       ENDIF
C
       IF(TURBMOD.EQ.7) THEN
       GAMMA_K=F(ISDEN+II)*F(ISOM+II)*BETA_STR*VTERM
       ENDIF
C
C Changes made by BNR on 5 July, 2002
C
C  A LIMITER FOR THE DISSIPATION TERM (DENSIT*EPSILON = Gamma_k*k) IS APPLIED
C  BASED ON THE SUGGESTION BY A.J.Lew et al, IJCFD,2001, VOl.14, PP. 201-200
C  GAMMA_K IS LIMITED TO ENSURE POSITIVE VALUE SINCE A NEGATIVE VALUE CAUSES
C  THE EXPONTENTIAL GROWTH OF THE SOLUTION
C
       F(ISSP+II)=F(ISSP+II) - max(GAMMA_K,0.)
C
C  THE ADDITIONAL TERM D IN THE K EQUATION IS NON ZERO ONLY FOR
C  THE CHIEN VERSION OF K-EPSILON MODEL (TURBMOD=3); OTHERWISE D=0

       F(ISSP+II)=F(ISSP+II)+D(II)*VTERM

	end do
	end do
C
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
C FIRST ORDER TIME  DISCRETISATION
C
         su_unst=term*te(ii,iblk,1)
         sp_unst=term
        else
C
C SECOND ORDER TIME DISCRETISATION
C
         su_unst=term*(2*te(ii,iblk,2)-0.5*te(ii,iblk,1))
         sp_unst=1.5*term
        end if

        F(ISSU+ii)=f(issu+ii)+su_unst
        F(ISSp+ii)=f(issp+ii)-sp_unst
C                    
C      RETAINING THE SOURCE TERMS BEFORE UPDATING THE SOURCE DUE TO GENERATION
C      TO BE USED FOR THE NEAR WALL CELLS IN HE STANDARD K_EPSILON MODEL
C
        BESUC(II)=BESUC(II)+su_unst
        BESPC(II)=BESPC(II)-sp_unst
    
        end do
        end do

      RETURN                                                            
      END
C
C************************************************************************
C                                                      
      SUBROUTINE SORCE6                                                      
C
C     CALCULATION OF SOURCE TERMS FOR EPSILON EQUATION                         
C
C************************************************************************
C
      INCLUDE 'com2d'                                                
C
C    CROSS DERIVATIVE DIFFUSIVE FLUXES
C
      CALL CDFLUX(ISED,6)                                               
C
C   NEAR WALL DAMPING FUNCTIONS FOR LOW RE CHIEN MODEL ( TURBMOD=3)
C
	IF (TURBMOD.EQ.3)  CALL CHIEN_MOD
C                                        
      DO J=jsolv_b,jsolv_e
      DO I=isolv_b,isolv_e
     
      II=I+(J-1)*NI
      VTERM=F(ISVOLP+II)*KBLK(II)  
      TLFACT=ALMU(II)/ALEPS(II)
      VISTUR=F(ISVIS+II)-VISCOS+SMALL     
      
      IF(TURBMOD.LE.3) THEN                          
      EPSBYK=CMU*FMU(II)*F(ISDEN+II)*(F(ISTE+II))*TLFACT/VISTUR  
c      F_EPSILON=C1*F1(II)*PROD(II)*EPSBYK*VTERM
C
C  EFFECT OF BUOYANCY ON DISSIPATION IS INCLUDED
C
      IF(NATCONV) CEPS_3=TANH(ABS(F(ISV+II)/F(ISU+II)))
      F_EPSILON=C1*F1(II)*(PROD(II)+CEPS_3*GK(II))*EPSBYK*VTERM
C 
      GAMMA_EPSILON=C2*F2(II)*F(ISDEN+II)*EPSBYK*VTERM
      ENDIF

      IF(TURBMOD.EQ.5) THEN
       AEPS=GREAT
       IF(ENN.EQ.6) AEPS=0.00285
       TKEII=ABS(F(ISTE+II))
       RY=SMIN(II)*SQRT(TKEII)/VISCOS
       IF(C_EPS .eq. 1)  
     >  F1(II) = 1.4*(1.+0.045*SQRT(abs(F(ISTE+II))/(F(ISV2+II)+SMALL)))
       IF(C_EPS .eq. 2)
     >   F1(II) = (1.3+0.25/(1+(0.15*SMIN(II)/(ALFAC(I,J)+small))**2)**4)
       IF(C_EPS .eq. 3)  F1(II) = 1.55+exp(-AEPS*RY*RY)
C      EPSBYK=1./TFAC(I,J)
       EPSBYK=F(ISED+II)/(F(ISTE+II)+small)
       F_EPSILON=F1(II)*PROD(II)*EPSBYK*VTERM
       GAMMA_EPSILON=C2*F(ISDEN+II)*EPSBYK*VTERM
       ENDIF
C
C Changes made by BNR on 5 July, 2002
C
C  A LIMITER FOR THE PRODUCTION OF EPSILON IS APPLIED
C  BASED ON THE SUGGESTION BY A.J.Lew et al, IJCFD,2001, VOl.14, PP. 201-200
C  F_EPSILON IS LIMITED TO ENSURE POSITIVE VALUE 
C
       F(ISSU+II)=F(ISSU+II)+max(F_EPSILON,0.)
C
C  A LIMITER FOR THE PRODUCTION OF EPSILON IS APPLIED
C  BASED ON THE SUGGESTION BY A.J.Lew et al, IJCFD,2001, VOl.14, PP. 201-200
C  GAMMA_EPSILON IS LIMITED TO ENSURE POSITIVE VALUE 
C
C
C  THE ADDITIONAL TERM E IN THE EPSILON EQUATION IS NON ZERO ONLY FOR
C  THE CHIEN VERSION OF K-EPSILON MODEL (TURBMOD=3);OTHERWISE E=0
C
       F(ISSP+II)=F(ISSP+II)-max(GAMMA_EPSILON,0.)+E(II)*VTERM

	end do
	end do
C
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
C FIRST ORDER DISCRETISATIOn
C
         su_unst=term*ed(ii,iblk,1)
         sp_unst=term
        else
C
C SECOND ORDER DISCRETISATIOn
C
         su_unst=term*(2*ed(ii,iblk,2)-0.5*ed(ii,iblk,1))
         sp_unst=1.5*term
        end if

        F(ISSU+ii)=f(issu+ii)+su_unst
        F(ISSp+ii)=f(issp+ii)-sp_unst
        end do
        end do

      RETURN                                                            
      END
C
C************************************************************************
C
      SUBROUTINE SORCE7
C
C     CALCULATION OF SOURCE TERMS FOR OMEGA EQUATION                           
C
C************************************************************************
C                                                      
      INCLUDE 'com2d'
C              
C    CROSS DERIVATIVE DIFFUSIVE FLUXES
C
      CALL CDFLUX(ISOM,7)
C
      IF(TURBMOD.EQ.4) THEN 

C
C     NON ZERO LIMITER USED FOR PRODUCTION AND DESTRUCTION TERMS 
C
      DO J=jsolv_b,jsolv_e
      DO I=isolv_b,isolv_e
      II=I+(J-1)*NI
      VTERM=F(ISVOLP+II)*KBLK(II)  
      VISTUR=F(ISVIS+II)-VISCOS+SMALL  
C     SUII=GAMMA*PROD(II)*F(ISDEN+II)/VISTUR 
C     SUII=GAMMA*STRN(II)*STRN(II)
      SUII=GAMMA*PROD(II)*F(ISOM+II)/F(ISTE+II)
    
      SUII=SUII*VTERM 
      F(ISSU+II)=F(ISSU+II)+SUII
C
      SPII=BETA*F(ISDEN+II)*F(ISOM+II)*VTERM
C     SPII=BETA*F(ISDEN+II)*F(ISDEN+II)*F(ISTE+II)*VTERM/VISTUR

      F(ISSP+II)=F(ISSP+II)-SPII
C
      END DO
      END DO
     
       END IF          !  IF (TURBMOD.EQ.4 )  LOOP ENDS HERE
C
C   FOR SST MODEL OF MENTER 
C
C     NON ZERO LIMITER USED FOR PRODUCTION AND DESTRUCTION TERMS 
C
       IF (TURBMOD.EQ.7 )  THEN

      DO J=jsolv_b,jsolv_e
      DO I=isolv_b,isolv_e
      II=I+(J-1)*NI
      TERM1=F(ISVOLP+II)*KBLK(II)
      VISTUR=F(ISVIS+II)-VISCOS+SMALL
      CALL DIFF(II,ISTE,DKDX,DKDY)
      CALL DIFF(II,ISOM,DOMDX,DOMDY)
      CDKW1=(DKDX*DOMDX+DKDY*DOMDY)
      SIGMA_OMG2=(1./PR2_OM)
      CDKW1=CDKW1*2.*F(ISDEN+II)*SIGMA_OMG2/F(ISOM+II)
      CDKW=MAX(CDKW1,SMALL)
C
      GAMA=GAMA1*F1_SST(II)+(1.-F1_SST(II))*GAMA2
C
C      SUII=GAMA*PROD(II)*F(ISOM+II)/F(ISTE+II)+2.*(1.-F1_SST(II))*F(ISDEN+II)
C     >     *SIGMA_OMG2*(DKDX*DOMDX+DKDY*DOMDY)/(F(ISOM+II)+SMALL)
C     SUII=GAMA*PROD(II)*F(ISOM+II)/F(ISTE+II)+2.*(1.-F1_SST(II))*F(ISDEN+II)
C    >     *SIGMA_OMG2*(DKDX*DOMDX+DKDY*DOMDY)/F(ISOM+II)
      SUII=GAMA*PROD(II)*F(ISOM+II)/F(ISTE+II)+(1.-F1_SST(II))*CDKW
C     SUII=GAMA*STRN(II)*STRN(II)+(1.-F1_SST(II))*CDKW

      SUII=SUII*TERM1
C
C     SUII=MAX(0.,SUII)
      F(ISSU+II)=F(ISSU+II)+SUII
C
      BETA=BETA1*F1_SST(II)+(1.-F1_SST(II))*BETA2
C
      SPII=BETA*F(ISDEN+II)*F(ISOM+II)*TERM1
      SPII=MAX(0.,SPII)
      F(ISSP+II)=F(ISSP+II)-SPII
C
      END DO
      END DO
C
       END IF          !  IF (TURBMOD.EQ.7 )  LOOP ENDS HERE
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
         su_unst=term*aom(ii,iblk,1)
         sp_unst=term
        else
C
C SECOND ORDER TIME DISCRETISATION
C
         su_unst=term*(2*aom(ii,iblk,2)-0.5*aom(ii,iblk,1))
         sp_unst=1.5*term
        end if

        F(ISSU+ii)=f(issu+ii)+su_unst
        F(ISSp+ii)=f(issp+ii)-sp_unst

        end do
        end do

      RETURN                                                            
      END
C
C***********************************************************************
C
      REAL FUNCTION F1_SST(II)
C
C***********************************************************************
C
      INCLUDE 'com2d'
C
      ARG_F11=SQRT(abs(F(ISTE+II)))/0.09/(F(ISOM+II)+SMALL)/(SMIN(II)+SMALL)
      ARG_F12=500.*VISCOS/(SMIN(II)+SMALL)/(SMIN(II)+SMALL)/(F(ISOM+II)+SMALL)
      ARG_F1=MAX(ARG_F11,ARG_F12)
C
      CALL DIFF(II,ISTE,DKDX,DKDY)
      CALL DIFF(II,ISOM,DOMDX,DOMDY)
      CDKW1=(DKDX*DOMDX+DKDY*DOMDY)
      SIGMA_OMG2=(1./PR2_OM)
      CDKW1=CDKW1*2.*F(ISDEN+II)*SIGMA_OMG2/F(ISOM+II)
      CDKW=MAX(CDKW1,SMALL)

      ARG_F2=4.*F(ISDEN+II)*SIGMA_OMG2*F(ISTE+II)
      ARG_F2=ARG_F2/CDKW/(SMIN(II)+SMALL)/(SMIN(II)+SMALL)
      ARG1= MIN(ARG_F1,ARG_F2)
C
      F1_SST=TANH(ARG1**4)
C
      RETURN
      END
C
C***********************************************************************
C
      REAL FUNCTION F2_SST(II)
C
C***********************************************************************
C
      INCLUDE 'com2d'
C
      ARG_F21=2.*SQRT(abs(F(ISTE+II)))/0.09/(F(ISOM+II)+SMALL)/(SMIN(II)+SMALL)
      ARG_F22=500.*VISCOS/(SMIN(II)+SMALL)/(SMIN(II)+SMALL)/(F(ISOM+II)+SMALL)
      ARG_F2=MAX(ARG_F21,ARG_F22)
C
      F2_SST=TANH(ARG_F2**2)
C
      RETURN
      END
C
C************************************************************************
C                                                      
      SUBROUTINE SORCE8
C
C     CALCULATION OF SOURCE TERMS FOR V-SQUARE  EQUATION                      
C
C************************************************************************
C                                                      
      INCLUDE 'com2d'
C                                                
C    CROSS DERIVATIVE DIFFUSIVE FLUXES
C
      CALL CDFLUX(ISV2,8)
C
      TWOTHIRD=2./3.
C
      DO J=jsolv_b,jsolv_e
      DO I=isolv_b,isolv_e
         II=I+(J-1)*NI
         TERM=F(ISVOLP+II)*KBLK(II)
         SUII=F(ISTE+II)*F(ISF+II)*TERM
         F(ISSU+II)=F(ISSU+II)+max(0.,suii)
         EPSBYK=1./TFAC(I,J)
C        EPSBYK=F(ISED+II)/F(ISTE+II)
C        F(ISSP+II)=F(ISSP+II)-TERM*ENN/TFAC(I,J)
         GAMMA_V2=TERM*ENN*EPSBYK
         F(ISSP+II)=F(ISSP+II)-max(0.,gamma_v2)
C        
        END DO
      END DO
C
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
         su_unst=term*vsq(ii,iblk,1)
         sp_unst=term
        else
C
C SECOND ORDER TIME DISCRETISATION
C
         su_unst=term*(2*vsq(ii,iblk,2)-0.5*vsq(ii,iblk,1))
         sp_unst=1.5*term
        end if

        F(ISSU+ii)=f(issu+ii)+su_unst
        F(ISSp+ii)=f(issp+ii)-sp_unst

        end do
        end do

      RETURN                                                            
      END                                                               
C
C************************************************************************
C                                                      
      SUBROUTINE SORCE9
C
C     CALCULATION OF SOURCE TERMS FOR F- EQUATION                            
C
C************************************************************************
C                                                      
      INCLUDE 'com2d'
C                                                
      DO J=jsolv_b,jsolv_e
      DO I=isolv_b,isolv_e
         
         II=I+(J-1)*NI
         VTERM=F(ISVOLP+II)*KBLK(II)
         SUII=(C1F-ENN)*(F(ISV2+II)/(F(ISTE+II)+small))
         SUII=(SUII-(C1F-1.)*2./3.)/(TFAC(I,J)+SMALL)
         SUII=SUII-C2F*PROD(II)/(F(ISTE+II)+small)
         F(ISSU+II)=F(ISSU+II)-SUII*VTERM
         F(ISSP+II)=-VTERM
       
      END DO
      END DO

      RETURN                                                            
      END 
C
C************************************************************************
C                                                      
      SUBROUTINE SORCE10
C
C     CALCULATION OF SOURCE TERMS FOR NUTILDE EQUATION  
C
C************************************************************************
C                                                      
      INCLUDE 'com2d'
C                                                
C    CROSS DERIVATIVE DIFFUSIVE FLUXES
C
      CALL CDFLUX(ISNUT,10)
C                                                
      DO J=jsolv_b,jsolv_e
      DO I=isolv_b,isolv_e
        ii=i+(j-1)*ni

        VTERM=KBLK(II)*F(ISVOLP+II)
C
C EVALUATE THE GRADIENTS OF NUT,U AND V
C
        if (kblk(ii) .ne. 0) then
	CALL DIFF(II,ISNUT,DSDX,DSDY)
	CALL DIFF(II,ISU,DUDX,DUDY)
 	CALL DIFF(II,ISV,DVDX,DVDY)

        VORT=ABS(DUDY-DVDX)

         AMUTR=F(ISNUT+II)/VISCOS
         FV2=1.-AMUTR/(1.+AMUTR*FV1(AMUTR))
         ALENGTH=CAPPA*SMIN(II)
         STILDE=VORT+F(ISNUT+II)*FV2/ALENGTH/ALENGTH
         SUPROD=CB1*STILDE*F(ISNUT+II)
         SUPROD=SUPROD+CB2*(DSDX*DSDX+DSDY*DSDY)/PR(10)
         ARG=F(ISNUT+II)/(STILDE*ALENGTH*ALENGTH)
         SPII=CW1*FW(ARG)*F(ISNUT+II)*VTERM/SMIN(II)**2
         SUII=SUPROD*VTERM
         F(ISSU+II)=F(ISSU+II)+SUII
         IF(SPII.GE.0.) F(ISSP+II)=F(ISSP+II)-SPII                  
         IF(SPII.LT.0.) F(ISSU+II)=F(ISSU+II)-SPII*F(ISNUT+II)      
        end if ! end of kblk check            
        end do
        end do
C
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
         su_unst=term*anut(ii,iblk,1)
         sp_unst=term
        else
C
C SECOND ORDER TIME DISCRETISATION
C
         su_unst=term*(2*anut(ii,iblk,2)-0.5*anut(ii,iblk,1))
         sp_unst=1.5*term
        end if

        F(ISSU+ii)=f(issu+ii)+su_unst
        F(ISSP+ii)=f(issp+ii)-sp_unst
        end do
        end do

      RETURN                                                            
      END     
C     
C***********************************************************************
C     
      SUBROUTINE COEFF_F
C     
C     CALCULATE COEFFICIENTS FOR THE EQUATION OF ELLIPTIC RELAXATION
C     
C***********************************************************************
C     
       INCLUDE 'com2d'
C
      SQR(VALUE)=VALUE*VALUE                                            
C
      DO J=jsolv_b,jsolv_e
      DO I=isolv_b,isolv_e
           II=I+(J-1)*NI
           IP=II+1
              JP=II+NI
            FXIP=F(ISFX+IP)
            D11E=SQR(F(ISB11W+IP))+SQR(F(ISB12W+IP))
            VOLE=0.5*(F(ISVOLP+IP)+F(ISVOLP+II))
              RE=FXIP*F(ISR+IP)+(1.-FXIP)*F(ISR+II)
              DE=RE*RE*D11E/(VOLE+SMALL)
              F(ISAE+II)=DE*KBLK(II)
                 
               FXP=F(ISFX+II)
            D11W=SQR(F(ISB11W+II))+SQR(F(ISB12W+II))
            VOLW=0.5*(F(ISVOLP+II)+F(ISVOLP+II-1))
              RW=FXP*F(ISR+II)+(1.-FXP)*F(ISR+II-1)
              DW=RW*RW*D11W/(VOLW+SMALL)
              F(ISAW+II)=DW*KBLK(II)
              
               FYJP=F(ISFY+JP)
            D22N=SQR(F(ISB21S+JP))+SQR(F(ISB22S+JP))
            VOLN=0.5*(F(ISVOLP+JP)+F(ISVOLP+II))
              RN=FYJP*F(ISR+JP)+(1.-FYJP)*F(ISR+II)
              DN=RN*RN*D22N/(VOLN+SMALL)
              F(ISAN+II)=DN*KBLK(II)

               FYP=F(ISFY+II)
            D22S=SQR(F(ISB21S+II))+SQR(F(ISB22S+II))
            VOLS=0.5*(F(ISVOLP+II)+F(ISVOLP+II-NI))
              RS=FYP*F(ISR+II)+(1.-FYP)*F(ISR+II-NI)
              DS=RS*RS*D22S/(VOLS+SMALL)
              F(ISAS+II)=DS*KBLK(II)

           if (kblk(ii) .ne. 0)then
           ALSQR=ALFAC(I,J)**2
           end if
C
C     COMPUTE CROSS DERIVATIVE FLUXES FOR F
C
                 CALL CDFF(II,ALSQR)

           F(ISAW+II)=ALSQR*F(ISAW+II)
           F(ISAE+II)=ALSQR*F(ISAE+II)
           F(ISAN+II)=ALSQR*F(ISAN+II)
           F(ISAS+II)=ALSQR*F(ISAS+II)

         END DO
        END DO

      RETURN                               
      END                                                               
C
C************************************************************************
C                                                      
      SUBROUTINE CDFF(II,ALSQR)
C                                        
C     COMPUTES CROSS-DERIVATIVE FLUXES FOR VARIABLE F IN THE V2F MODEL
C                        
C***********************************************************************
C    
          INCLUDE 'com2d'                                                 
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
      UIJK =F(ISF+II)                                                 
      UIP  =F(ISF+II+1)                                               
      UIPJP=F(ISF+II+1+NI)                                            
      UIPJM=F(ISF+II+1-NI)                                            
      UIM  =F(ISF+II-1)                                               
      UIMJP=F(ISF+II-1+NI)                                            
      UIMJM=F(ISF+II-1-NI)                                            
      UJP  =F(ISF+II+NI)                                              
      UJM  =F(ISF+II-NI)                                              
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
        IP=II+1                                                         
      D12E=RE*RE*(F(ISB21W+IP)*F(ISB11W+IP)+F(ISB22W+IP)*F(ISB12W+IP))  
      D12W=RW*RW*(F(ISB21W+II)*F(ISB11W+II)+F(ISB22W+II)*F(ISB12W+II))  
        IP=II+NI                                                        
      D21N=RN*RN*(F(ISB21S+IP)*F(ISB11S+IP)+F(ISB22S+IP)*F(ISB12S+IP))  
      D21S=RS*RS*(F(ISB21S+II)*F(ISB11S+II)+F(ISB22S+II)*F(ISB12S+II))  
                                                                        
C    RECIPROCAL OF CELL VOLUMES AROUND DIFFERENT FACES
                 
      RVE=1./(0.5*(F(ISVOLP+II)+F(ISVOLP+II+1 ))+SMALL)                 
      RVW=1./(0.5*(F(ISVOLP+II)+F(ISVOLP+II-1 ))+SMALL)                 
      RVN=1./(0.5*(F(ISVOLP+II)+F(ISVOLP+II+NI))+SMALL)                 
      RVS=1./(0.5*(F(ISVOLP+II)+F(ISVOLP+II-NI))+SMALL)                 
                                                                        
C     CROSS-DERIVATIVE TERMS

        SUIJ=RVE*D12E*UNSE-RVW*D12W*UNSW  
     >      +RVN*D21N*UEWN-RVS*D21S*UEWS
        SUIJ=SUIJ*KBLK(II) 
        F(ISSU+II)=ALSQR*SUIJ       
C   
      RETURN
C                           
      END                                                               
C
C***********************************************************************
C
C   
      REAL FUNCTION TFAC(I,J)
C
C***********************************************************************
C
      INCLUDE 'com2d'
C
       II=I+(J-1)*NI

       TKEII=ABS(F(ISTE+II))
       EDII=ABS(F(ISED+II))
       T1=TKEII/(EDII+SMALL)
       T2=6.*(VISCOS/F(ISDEN+II)/(EDII+SMALL))**0.5
       TFAC=MAX(T1,T2)

       IF(.NOT.REAL_IND)   RETURN
C
C REALIZABLITY CONSTRAINTS
C
C       TR=0.6*F(ISTE+II)/SQRT(3.0)/F(ISV2+II)/CMU/STRN(II)
        TR=0.6*TKEII/SQRT(3.0)/F(ISV2+II)/CMU/STRN(II)
        TFAC=MIN(TFAC,TR)
        RETURN

      END
C
C***********************************************************************
C
      REAL FUNCTION ALFAC(I,J)
C
C***********************************************************************
C
      INCLUDE 'com2d'
C
           II=I+(J-1)*NI

           TKEII=ABS(F(ISTE+II))
           EDII=ABS(F(ISED+II))
           AL1=TKEII**1.5/(EDII+SMALL)
           AL2=CETA*VISCOS**0.75/((EDII+SMALL)**0.25)
          ALFAC=CL*MAX(AL1,AL2)
C
         IF(.NOT.REAL_IND)   RETURN

C REALIZABLITY CONSTRAINTS
C
        if (kblk(ii) .ne. 0) then
C        ALR=F(ISTE+II)**1.5/SQRT(3.0)/F(ISV2+II)/CMU/STRN(II)
         ALR=TKEII**1.5/SQRT(3.0)/F(ISV2+II)/CMU/STRN(II)
         AL=MIN(AL1,ALR)
         ALFAC=CL*MAX(AL,AL2)
      	end if
          RETURN
           END
C*********************************************************************
C   
         SUBROUTINE KEPSSTD(IP,IB,ISANOW,IFLAG,AN1,AN2,DELTA,AREA)
C   
C*********************************************************************
C   
           INCLUDE 'com2d'
C
          VR=AN2*F(ISU+IP)-AN1*F(ISV+IP)
          UWAL=F(ISU+IB)
          VWAL=F(ISV+IB)
          TWAL=F(ISS+IB)
          VPARL=AN2*UWAL-AN1*VWAL-VR+small
          SQRTK=SQRT(ABS(F(ISTE+IP)))
C
C  STANDARD WALL FUNCTION FORMULATION
C
          USTAR=CDQR*SQRTK
          YPLS=F(ISDEN+IP)*USTAR*DELTA/VISCOS+SMALL
          F01=MAX(0.,(YPLS-11.06)/(ABS(YPLS-11.06)+SMALL))
          TMULT=(1.-F01)*VISCOS/DELTA
     >          +F01*DENSIT*CDQR*SQRTK*CAPPA/LOG(ELOG*YPLS)
          TAUWAL(Ip,IFLAG)=TMULT*VPARL
          UTAU=SQRT(ABS(TMULT*VPARL)/DENSIT)
          YPLUS(IP,IFLAG)=(1.-F01)*DENSIT*UTAU*DELTA/VISCOS+F01*YPLS
C
          A11=1.-AN1*AN1
          A22=1.-AN2*AN2
          A12=AN1*AN2
          SHEAR=TMULT*AREA
           
          IF(NN.EQ.2) THEN

C      U-VELOCITY SOURCES
C
       F(ISSP+IP)=F(ISSP+IP)-SHEAR*A11
       F(ISSU+IP)=F(ISSU+IP)+SHEAR*((F(ISV+IP)-VWAL)*A12+UWAL*A11)
C      F(ISSU+IP)=F(ISSU+IP)+SHEAR*(F(ISV+IP)*A12+UWAL)
C
      ELSEIF(NN.EQ.3) THEN
C
C      V-VELOCITY SOURCES
C
       F(ISSP+IP)=F(ISSP+IP)-SHEAR*A22
       F(ISSU+IP)=F(ISSU+IP)+SHEAR*((F(ISU+IP)-UWAL)*A12+VWAL*A22)
C      F(ISSU+IP)=F(ISSU+IP)+SHEAR*(F(ISU+IP)*A12+VWAL)
C                                                                
           ELSEIF (NN.EQ.4) THEN 
C
C        WR-VELOCITY SOURCES
C
       RAT=F(ISR+IP)/(RN+SMALL)
       F(ISSP+IP)=F(ISSP+IP)-SHEAR/(RAT+SMALL)
       F(ISSU+IP)=F(ISSU+IP)+SHEAR*RAT*F(ISWR+IB)        


           ELSEIF (NN.EQ.5) THEN
           
           GENR=ABS(TMULT*VPARL*VPARL/DELTA)
           TERM=(1.-F01)*F(ISDEN+IP)*CDTQ*SQRTK*YPLS/DELTA
     >              +F01*F(ISDEN+IP)*CDTQ*SQRTK*LOG(ELOG*YPLS)/(CAPPA*DELTA)
           F(ISSP+IP)=BESPC(IP)-TERM*F(ISVOLP+IP)
           F(ISSU+IP)=BESUC(IP)+GENR*F(ISVOLP+IP)
           CONTINUE
C
           ELSEIF (NN.EQ.6) THEN
C
            TKEIP=F(ISTE+IP)
            EPSVIS=2.*TKEIP*VISCOS/DELTA/DELTA
            EPSLOG=  TKEIP**1.5*CDTQ/(CAPPA*DELTA)
            EPSNW=F01*EPSLOG+(1-F01)*EPSVIS
            F(ISSU+IP)=GREAT*EPSNW
            F(ISSP+IP)=-GREAT

         ELSEIF (NN.EQ.11) THEN

C TO INCLUDE SP AND SU FOR SCALAR FOR STD WALL

           END IF    ! END OF IF LOOP  FOR FLOW VARIABLES
C
C DELINK THE COEFFICIENTS CONNECTING THE WALL NODE
C
            F(ISANOW+IP)=0.
           
            RETURN
            END   
C
C*************************************************************************
C
         SUBROUTINE KEPS2LR(IP,IDD2,IB,ISANOW,IFLAG,AN1,AN2,DELTA,AREA)
C
C*************************************************************************
C
           INCLUDE 'com2d'
C
C  EVALUATION OF WALL SHEAR STRESS
C
          VR=AN2*F(ISU+IP)-AN1*F(ISV+IP)
          UWAL=F(ISU+IB)
          VWAL=F(ISV+IB)
	  TWAL=F(ISS+IB)
          VPARL=AN2*UWAL-AN1*VWAL-VR
          TAUWAL(IP,IFLAG)=VISCOS*VPARL/DELTA
          UTAU=SQRT(ABS(TAUWAL(IP,IFLAG))/F(ISDEN+IP))
          YPLUS(IP,IFLAG)=F(ISDEN+IP)*UTAU*DELTA/VISCOS
          TMULT=VISCOS/DELTA
C
C
          A11=1.-AN1*AN1
          A22=1.-AN2*AN2
          A12=AN1*AN2
          SHEAR=TMULT*AREA
           
          IF(NN.EQ.2) THEN
          
C      U-VELOCITY SOURCES
C
       F(ISSP+IP)=F(ISSP+IP)-SHEAR*A11
       F(ISSU+IP)=F(ISSU+IP)+SHEAR*((F(ISV+IP)-VWAL)*A12+UWAL*A11)
C      F(ISSU+IP)=F(ISSU+IP)+SHEAR*(F(ISV+IP)*A12+UWAL)
C
      ELSEIF(NN.EQ.3) THEN
C
C      V-VELOCITY SOURCES
C
       F(ISSP+IP)=F(ISSP+IP)-SHEAR*A22
       F(ISSU+IP)=F(ISSU+IP)+SHEAR*((F(ISU+IP)-UWAL)*A12+VWAL*A22)
C      F(ISSU+IP)=F(ISSU+IP)+SHEAR*(F(ISU+IP)*A12+VWAL)
C                                                                
           ELSEIF (NN.EQ.4) THEN 
C
C        WR-VELOCITY SOURCES
C
       RAT=F(ISR+IP)/(RN+SMALL)
       F(ISSP+IP)=F(ISSP+IP)-SHEAR/(RAT+SMALL)
       F(ISSU+IP)=F(ISSU+IP)+SHEAR*RAT*F(ISWR+IB)        
         
          ELSEIF (NN.EQ.5) THEN 
          
C
          ELSEIF (NN.EQ.6) THEN 
C
C    RODI'S TWO LAYER MODEL
C
             IBEG=IP
             INTV=-NI*MFLAG1+NI*MFLAG2-MFLAG3+MFLAG4
             JWALL=IP/NI+1
             IWALL=IP-(JWALL-1)*NI
             IFIN=IBEG+INTV*(MFLAG1*(JWALL-2)+MFLAG2*(NJ-1-JWALL)+
     >                         MFLAG3*(IWALL-2)+MFLAG4*(NI-1-IWALL))
C
             DO II=IBEG,IFIN,INTV
                 DIST = ABS((F(ISX+II)-F(ISCOX  +IDD2))*AN1
     >          +(F(ISY+II  )-F(ISCOY+IDD2))*AN2)+SMALL
                  DISTN(II)=min(dist,distn(ii))
                  RY=f(isden+ii)*sqrt(f(iste+ii))*distn(ii)/VISCOS
                  FACTM=1.-EXP(-RY*CDAMP*ADAMP)
                  ALMU(II)=CL1*DISTN(II)*FACTM
                  FACTE=1.+CEPS/CL1/RY
                  ALEPS(II)=CL1*DISTN(II)*FACTE
                  EPSII=F(ISTE+II)**1.5/(ALEPS(II)+small)
                  FACT=FACTM/FACTE
                  F(ISSP+II)=-GREAT        
                  F(ISSU+II)= GREAT*EPSII
                  YPLS=F(ISDEN+II)*UTAU*DISTN(II)/VISCOS
                  IF (FACT .GT. 0.95) GO TO 1001
              END DO
C
1001      CONTINUE
C                                    
         ELSEIF (NN.EQ.11) THEN
C
C        SCALAR  SOURCES
C
         COND=VISCOS*AREA/PRLAM/DELTA
         F(ISSP+IP)=F(ISSP+IP)-COND
         F(ISSU+IP)=F(ISSU+IP)+COND*TWAL

           END IF    ! END OF IF LOOP  FOR FLOW VARIABLES
C           
C DELINK THE COEFFICIENTS CONNECTING THE WALL NODE
C
             F(ISANOW+IP)=0.

            RETURN
            END
C
C*************************************************************************
C
          SUBROUTINE KEPSCHIEN (IP,IB,ISANOW,IFLAG,AN1,AN2,DELTA,AREA)
C
C*************************************************************************
C
           INCLUDE 'com2d'
C   
          VR=AN2*F(ISU+IP)-AN1*F(ISV+IP)
          UWAL=F(ISU+IB)
          VWAL=F(ISV+IB)
          TWAL=F(ISS+IB)
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
       F(ISSU+IP)=F(ISSU+IP)+SHEAR*((F(ISV+IP)-VWAL)*A12+UWAL*A11)
C      F(ISSU+IP)=F(ISSU+IP)+SHEAR*(F(ISV+IP)*A12+UWAL)
C
      ELSEIF(NN.EQ.3) THEN
C
C      V-VELOCITY SOURCES
C
       F(ISSP+IP)=F(ISSP+IP)-SHEAR*A22
       F(ISSU+IP)=F(ISSU+IP)+SHEAR*((F(ISU+IP)-UWAL)*A12+VWAL*A22)
C      F(ISSU+IP)=F(ISSU+IP)+SHEAR*(F(ISU+IP)*A12+VWAL)
C                                                                
           ELSEIF (NN.EQ.4) THEN 
C
C        WR-VELOCITY SOURCES
C
       RAT=F(ISR+IP)/(RN+SMALL)
       F(ISSP+IP)=F(ISSP+IP)-SHEAR/(RAT+SMALL)
      
           ELSEIF (NN.EQ.5) THEN 

           ELSEIF (NN.EQ.6) THEN 

          ELSEIF (NN.EQ.11) THEN
C
C        SCALAR  SOURCES
C
         COND=VISCOS*AREA/PRLAM/DELTA
         F(ISSP+IP)=F(ISSP+IP)-COND
         F(ISSU+IP)=F(ISSU+IP)+COND*TWAL

           END IF    ! END OF IF LOOP  FOR FLOW VARIABLES
C           
C DELINK THE COEFFICIENTS CONNECTING THE WALL NODE
C
             F(ISANOW+IP)=0.
 
          RETURN
          END  
C
C*************************************************************************
C
         SUBROUTINE KOMEGA(IP,IB,ISANOW,IFLAG,AN1,AN2,DELTA,AREA)
C
C*************************************************************************
C
          INCLUDE 'com2d'
C
          VR=AN2*F(ISU+IP)-AN1*F(ISV+IP)
          UWAL=F(ISU+IB)
          VWAL=F(ISV+IB)
	  TWAL=F(ISS+IB)
          VPARL=AN2*UWAL-AN1*VWAL-VR
          SQRTK=SQRT(ABS(F(ISTE+IP)))
C
          USTAR=CDQR*SQRTK
          YPLS=F(ISDEN+IP)*USTAR*DELTA/VISCOS+SMALL

          F01=MAX(0.,(YPLS-11.00)/(ABS(YPLS-11.00)+SMALL))

          TMULT=(1.-F01)*VISCOS/DELTA
     >          +F01*DENSIT*CDQR*SQRTK*CAPPA/LOG(ELOG*YPLS)
          TAUWAL(IP,IFLAG)=TMULT*VPARL
          UTAU=SQRT(ABS(TMULT*VPARL)/DENSIT)
          YPLUS(IP,IFLAG)=(1.-F01)*DENSIT*UTAU*DELTA/VISCOS+F01*YPLS

          A11=1.-AN1*AN1
          A22=1.-AN2*AN2
          A12=AN1*AN2
          SHEAR=TMULT*AREA
           
          IF(NN.EQ.2) THEN
          
C      U-VELOCITY SOURCES
C
       F(ISSP+IP)=F(ISSP+IP)-SHEAR*A11
       F(ISSU+IP)=F(ISSU+IP)+SHEAR*((F(ISV+IP)-VWAL)*A12+UWAL*A11)
C      F(ISSU+IP)=F(ISSU+IP)+SHEAR*(F(ISV+IP)*A12+UWAL)
C
      ELSEIF(NN.EQ.3) THEN
C
C      V-VELOCITY SOURCES
C
       F(ISSP+IP)=F(ISSP+IP)-SHEAR*A22
       F(ISSU+IP)=F(ISSU+IP)+SHEAR*((F(ISU+IP)-UWAL)*A12+VWAL*A22)
C      F(ISSU+IP)=F(ISSU+IP)+SHEAR*(F(ISU+IP)*A12+VWAL)
C                                                                
           ELSEIF (NN.EQ.4) THEN 
C
C        WR-VELOCITY SOURCES
C
       RAT=F(ISR+IP)/(RN+SMALL)
       F(ISSP+IP)=F(ISSP+IP)-SHEAR/(RAT+SMALL)
       F(ISSU+IP)=F(ISSU+IP)+SHEAR*RAT*F(ISWR+IB)        


           ELSEIF (NN.EQ.5) THEN
           IF(F01.EQ.1) THEN
           GENR=ABS(TMULT*VPARL*VPARL/DELTA)
            TERM=F(ISDEN+IP)*CDTQ*SQRTK*LOG(ELOG*YPLS)/(CAPPA*DELTA)

           F(ISSP+IP)=BESPC(IP)-TERM*F(ISVOLP+IP)
           F(ISSU+IP)=BESUC(IP)+GENR*F(ISVOLP+IP)
           END IF
           CONTINUE
C
           ELSEIF (NN.EQ.7) THEN

            AOMNW =6.*VISCOS/DELTA/DELTA/BETA 
            AOMW=10.*AOMNW

            F(ISSU+IP)=GREAT*AOMNW
            F(ISSP+IP)=-GREAT

C           F(ISOM+IB)=AOMW          
C           RETURN
C
           ELSEIF (NN.EQ.11) THEN
C
C        SCALAR  SOURCES
C
         COND=VISCOS*AREA/PRLAM/DELTA
         F(ISSP+IP)=F(ISSP+IP)-COND
         F(ISSU+IP)=F(ISSU+IP)+COND*TWAL

           END IF    ! END OF IF LOOP  FOR FLOW VARIABLES
C
C DELINK THE COEFFICIENTS CONNECTING THE WALL NODE
C
            F(ISANOW+IP)=0.
           
            RETURN
            END
C
C*************************************************************************
C
          SUBROUTINE V2F (IP,IB,ISANOW,IFLAG,AN1,AN2,DELTA,AREA)
C
C*************************************************************************
C
           INCLUDE 'com2d'
C
          VR=AN2*F(ISU+IP)-AN1*F(ISV+IP)
          UWAL=F(ISU+IB)
          VWAL=F(ISV+IB)
	  TWAL=F(ISS+IB)
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
       F(ISSU+IP)=F(ISSU+IP)+SHEAR*((F(ISV+IP)-VWAL)*A12+UWAL*A11)
C      F(ISSU+IP)=F(ISSU+IP)+SHEAR*(F(ISV+IP)*A12+UWAL)
C
      ELSEIF(NN.EQ.3) THEN
C
C      V-VELOCITY SOURCES
C
       F(ISSP+IP)=F(ISSP+IP)-SHEAR*A22
       F(ISSU+IP)=F(ISSU+IP)+SHEAR*((F(ISU+IP)-UWAL)*A12+VWAL*A22)
C      F(ISSU+IP)=F(ISSU+IP)+SHEAR*(F(ISU+IP)*A12+VWAL)
C                                                                
           ELSEIF (NN.EQ.4) THEN 
C
C        WR-VELOCITY SOURCES
C
       RAT=F(ISR+IP)/(RN+SMALL)
       F(ISSP+IP)=F(ISSP+IP)-SHEAR/(RAT+SMALL)
       F(ISSU+IP)=F(ISSU+IP)+SHEAR*RAT*F(ISWR+IB)        

          ELSEIF(NN.EQ.5) THEN
C   
          ELSEIF(NN.EQ.6) THEN

           F(ISSP+IP)=-GREAT
           F(ISSU+IP)= GREAT*2.*abs(F(ISTE+IP))*VISCOS/DELTA/DELTA     
C   
          ELSEIF(NN.EQ.8) THEN
C   
          ELSEIF(NN.EQ.9) THEN

	  RETURN

           ELSEIF (NN.EQ.11) THEN
C
C        SCALAR  SOURCES
C
         COND=VISCOS*AREA/PRLAM/DELTA
         F(ISSP+IP)=F(ISSP+IP)-COND
         F(ISSU+IP)=F(ISSU+IP)+COND*TWAL

         END IF    ! END OF IF LOOP  FOR FLOW VARIABLES
C
C DELINK THE COEFFICIENTS CONNECTING THE WALL NODE
C
            F(ISANOW+IP)=0.
C
          RETURN
          END  
C
C*************************************************************************
C
          SUBROUTINE SPL_ALM(IP,IB,ISANOW,IFLAG,AN1,AN2,DELTA,AREA)
C
C*************************************************************************
C
           INCLUDE 'com2d'
C
          VR=AN2*F(ISU+IP)-AN1*F(ISV+IP)
          UWAL=F(ISU+IB)
          VWAL=F(ISV+IB)
	  TWAL=F(ISS+IB)
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
       F(ISSU+IP)=F(ISSU+IP)+SHEAR*((F(ISV+IP)-VWAL)*A12+UWAL*A11)
C      F(ISSU+IP)=F(ISSU+IP)+SHEAR*(F(ISV+IP)*A12+UWAL)
C
      ELSEIF(NN.EQ.3) THEN
C
C      V-VELOCITY SOURCES
C
       F(ISSP+IP)=F(ISSP+IP)-SHEAR*A22
       F(ISSU+IP)=F(ISSU+IP)+SHEAR*((F(ISU+IP)-UWAL)*A12+VWAL*A22)
C      F(ISSU+IP)=F(ISSU+IP)+SHEAR*(F(ISU+IP)*A12+VWAL)
C                                                                
           ELSEIF (NN.EQ.4) THEN 
C
C        WR-VELOCITY SOURCES
C
       RAT=F(ISR+IP)/(RN+SMALL)
       F(ISSP+IP)=F(ISSP+IP)-SHEAR/(RAT+SMALL)
       F(ISSU+IP)=F(ISSU+IP)+SHEAR*RAT*F(ISWR+IB)        

          ELSEIF(NN.EQ.10) THEN

	    RETURN
C   
           ELSEIF (NN.EQ.11) THEN
C
C        SCALAR  SOURCES
C
         COND=VISCOS*AREA/PRLAM/DELTA
         F(ISSP+IP)=F(ISSP+IP)-COND
         F(ISSU+IP)=F(ISSU+IP)+COND*TWAL

         END IF    ! END OF IF LOOP  FOR FLOW VARIABLES
C
C DELINK THE COEFFICIENTS CONNECTING THE WALL NODE
C
            F(ISANOW+IP)=0.
C
          RETURN
          END  
C
C*************************************************************************
C
         SUBROUTINE SST(IP,IB,ISANOW,IFLAG,AN1,AN2,DELTA,AREA)
C
C*************************************************************************
C
          INCLUDE 'com2d'
C
          VR=AN2*F(ISU+IP)-AN1*F(ISV+IP)
          UWAL=F(ISU+IB)
          VWAL=F(ISV+IB)
	  TWAL=F(ISS+IB)
          VPARL=AN2*UWAL-AN1*VWAL-VR
          SQRTK=SQRT(ABS(F(ISTE+IP)))
C         
          USTAR=CDQR*SQRTK
          YPLS=F(ISDEN+IP)*USTAR*DELTA/VISCOS+SMALL
          F01=MAX(0.,(YPLS-11.06)/(ABS(YPLS-11.06)+SMALL))
          TMULT=(1.-F01)*VISCOS/DELTA
     >          +F01*DENSIT*CDQR*SQRTK*CAPPA/LOG(ELOG*YPLS)
          TAUWAL(IP,IFLAG)=TMULT*VPARL
          UTAU=SQRT(ABS(TMULT*VPARL)/DENSIT)
          YPLUS(IP,IFLAG)=(1.-F01)*DENSIT*UTAU*DELTA/VISCOS+F01*YPLS

          YPLS=YPLUS(IP,IFLAG)

          A11=1.-AN1*AN1
          A22=1.-AN2*AN2
          A12=AN1*AN2
          SHEAR=TMULT*AREA
           
          IF(NN.EQ.2) THEN
          
C      U-VELOCITY SOURCES
C
       F(ISSP+IP)=F(ISSP+IP)-SHEAR*A11
       F(ISSU+IP)=F(ISSU+IP)+SHEAR*((F(ISV+IP)-VWAL)*A12+UWAL*A11)
C      F(ISSU+IP)=F(ISSU+IP)+SHEAR*(F(ISV+IP)*A12+UWAL)
C
      ELSEIF(NN.EQ.3) THEN
C
C      V-VELOCITY SOURCES
C
       F(ISSP+IP)=F(ISSP+IP)-SHEAR*A22
       F(ISSU+IP)=F(ISSU+IP)+SHEAR*((F(ISU+IP)-UWAL)*A12+VWAL*A22)
C      F(ISSU+IP)=F(ISSU+IP)+SHEAR*(F(ISU+IP)*A12+VWAL)
C                                                                
           ELSEIF (NN.EQ.4) THEN 
C
C        WR-VELOCITY SOURCES
C
       RAT=F(ISR+IP)/(RN+SMALL)
       F(ISSP+IP)=F(ISSP+IP)-SHEAR/(RAT+SMALL)
       F(ISSU+IP)=F(ISSU+IP)+SHEAR*RAT*F(ISWR+IB)        


           ELSEIF (NN.EQ.5) THEN
           IF(F01.EQ.1) THEN
           GENR=ABS(TMULT*VPARL*VPARL/DELTA)
            TERM=F(ISDEN+IP)*CDTQ*SQRTK*LOG(ELOG*YPLS)/(CAPPA*DELTA)

           F(ISSP+IP)=BESPC(IP)-TERM*F(ISVOLP+IP)
           F(ISSU+IP)=BESUC(IP)+GENR*F(ISVOLP+IP)
           END IF
           CONTINUE
C
           ELSEIF (NN.EQ.7) THEN
            TKEIP=F(ISTE+IP)
            AOMVIS=6.*VISCOS/DELTA/DELTA/BETA1
            AOMLOG=  abs(UTAU)/(CDRT*CAPPA*DELTA)
            AOMNW=F01*AOMLOG+(1-F01)*AOMVIS
            F(ISSU+IP)=GREAT*AOMNW
            F(ISSP+IP)=-GREAT
C                                   
           ELSEIF (NN.EQ.11) THEN
C
C        SCALAR  SOURCES
C
         COND=VISCOS*AREA/PRLAM/DELTA
         F(ISSP+IP)=F(ISSP+IP)-COND
         F(ISSU+IP)=F(ISSU+IP)+COND*TWAL

           END IF    ! END OF IF LOOP  FOR FLOW VARIABLES
C
C DELINK THE COEFFICIENTS CONNECTING THE WALL NODE
C
            F(ISANOW+IP)=0.
           
            RETURN
            END
C                                                  
C***********************************************************************
C
      SUBROUTINE EDYVIS
C                                                  
C     UPDATE EFFECTIVE VISCOSITY AFTER 
C     EVALUATION OF K AND EPSILON FIELD
C                                                  
C***********************************************************************
C
      INCLUDE 'com2d'                                                
C
      RELAXP=RELAX(NVVIS-NVPP+1)
      RELAXM=1.-RELAXP

      DO J=1,NJ 
      DO I=1,NI 

       II=I+(J-1)*NI

      VISTOLD=F(ISVIS+II)-VISCOS
C
C EDDY VISCOSITY FOR STD K-EPSILON , 2-LAYER AND LOW RE CHIEN
C
       IF (TURBMOD.LE.3) THEN
      FACTOR=ALMU(II)/ALEPS(II)
      amu_t= FACTOR * F(ISDEN+II)
     >      *CMU*FMU(II)*F(ISTE+II)*F(ISTE+II)/(F(ISED+II)+SMALL)
        END IF
C
C EDDY VISCOSITY FOR  K-OMEGA
C
       IF (TURBMOD.EQ.4) THEN
      amu_t=  F(ISDEN+II)*F(ISTE+II)/(F(ISOM+II)+SMALL)
        END IF
C
C EDDY VISCOSITY FOR V2F 
C
       IF (TURBMOD.EQ.5) THEN
        amu_t=F(ISDEN+II)*CMU*F(ISV2+II)*TFAC(I,J)
        END IF
C
C EDDY VISCOSITY FOR SPALART-ALLMARAS MODEL
C
        IF(TURBMOD.EQ.6) THEN
C
C F(ISNUT+II) ARRAY IS USED FOR THE VARIABLE NUTILDE
C
      amutr=F(ISNUT+II)/VISCOS
      amu_t=F(ISNUT+II)*FV1(amutr)

        END IF
C
C EDDY VISCOSITY FOR SST MODEL
C
        IF(TURBMOD.EQ.7) THEN
C
	  if (kblk(ii) .ne. 0) then
          CALL DIFF(II,ISU,DUDX,DUDY)
          CALL DIFF(II,ISV,DVDX,DVDY)
C
          VORTIC(II)=ABS(DUDY-DVDX)
C
        AMU_T=A1*F(ISTE+II)/(MAX(A1*F(ISOM+II),VORTIC(II)*F2_SST(II))+SMALL)
	end if
C
        END IF
C
C  A LIMITER FOR MU_T IS APPLIED BASED ON THE SUGGEDTION BY
C  A.J.Lew et al, IJCFD,2001, VOl.14, PP. 201-200
C  MU_T IS LIMITED TO ENSURE POSITIVE VALUE. THIS NEEDS TO BE DONE SINCE
C  WE NOW ALLOW NEGATIVE VALUES OF K & EPSILON IN THE FIELD
C

       VISTNEW =max(amu_t,1.e-4*viscos)
C      VISTNEW =amu_t
      F(ISVIS +II)=RELAXP*VISTNEW + RELAXM*VISTOLD+VISCOS

        END DO          ! END OF I LOOP FOR COMPUTATION OF MUT
        END DO          ! END OF J LOOP FOR COMPUTATION OF MUT

      IF(SPECL(8)) CALL SPEBC8   
C
C COMPUTE CELL-FACE VISCOSITIES 
C
      DO  J=2,NJ                                                    
      DO  I=2,NI
      II=I+(J-1)*NI                                                    
      F(ISVISW+II)=F(ISFX+II)*F(ISVIS+II)+(1.-F(ISFX+II))*F(ISVIS+II-1)      
      F(ISVISS+II)=F(ISFY+II)*F(ISVIS+II)+(1.-F(ISFY+II))*F(ISVIS+II-NI)  
      END DO
      END DO
C                                                                        
C     MODIFY DIFFUSIVE COEFFICIENTS AT SOLID WALL                       
C
      IF(MWALN.NE.0) THEN                                               
        DO  I=1,MWALN                                                
     	   F(ISVISS+INDWN(I)+NI)=0.                                        
     	   F(ISVIS+INDWN(I)+NI)=0.                                        
        END DO 
      ENDIF                                                             

      IF(MWALS.NE.0) THEN                                               
        DO  I=1,MWALS                                                
     	   F(ISVISS+INDWS(I))=0.                                        
     	   F(ISVIS+INDWS(I)-NI)=0.                                        
        END DO 
      ENDIF                                                             

      IF(MWALE.NE.0) THEN                                               
        DO  I=1,MWALE                                                
     	   F(ISVISW+INDWE(I)+1)=0.                                        
     	   F(ISVIS+INDWE(I)+1)=0.                                        
        END DO 
      ENDIF                                                             

      IF(MWALW.NE.0) THEN                                               
        DO  I=1,MWALW                                                
     	   F(ISVISW+INDWW(I))=0.                                        
     	   F(ISVIS+INDWW(I)-1)=0.                                        
        END DO 
      ENDIF    
 
C
      RETURN                                                         

      END      
C
C***********************************************************************
C
        REAL FUNCTION FV1(R)
C
C    FUNCTION USED IN SA TURBULENCE MODEL
C
C***********************************************************************
C
        include 'com2d'
C
        FV1=R**3/(R**3+CV1**3)
C
        RETURN
        END
C
C******************************************************************
C
        REAL FUNCTION FW(R)
C
C    FUNCTION USED IN SA TURBULENCE MODEL
C
C******************************************************************
C
        include 'com2d'
C
        ONESIXTH=1./6.
C  AS R IS TOO LARGE FW REACHES AN ASYMPTOTIC VALUE
C  R IS LIMITED TO 10 OTHERWISE FOR VERY LARGE R
C  COMPUTATION GIVES 0. INSTEAD OF THE ASYMPTOTIC VALUE
C
        R=MIN(10.,R)
        G=R+CW2*(R**6-R)
        FW=G*((1.+CW3**6)/(G**6+CW3**6))**ONESIXTH

        RETURN
        END
C
C******************************************************************
C
      SUBROUTINE CHIEN_MOD     
C                              
C COMPUTES SOURCE TERMS FOR NEAR WALL EFFECTS IN ALL EQUATIONS
C                              
C***********************************************************************
C
      INCLUDE 'com2d'
C        
      SQR(VALUE)=VALUE*VALUE                                            
      VALUE=0.
      MWAL   =MWALN + MWALS + MWALE + MWALW

        if (mwal .eq. 0) return


C LOW RE VERSION OF CHIEN
C
C COMPUTE FMU,F2,D ETC. FOR THE K(TE) and EPSILON(ED) EQUATION
C
               DO J=JSOLV_B,JSOLV_E
               DO I=ISOLV_B,ISOLV_E
                  
                   II=I+(J-1)*NI

                IF(KBLK(II).EQ.1) THEN

                  dist_nor=SMIN(II)
                  IP=IWAL(II)
                  IFLAG=IWFLAG(II)
                  SHEAR=ABS(TAUWAL(IP,IFLAG))
                  utau=SQRT(SHEAR/DENSIT)  
	          ypls=densit*utau*dist_nor/viscos
                  REYT=DENSIT*F(ISTE+II)*F(ISTE+II)/VISCOS/F(ISED+II)
                  FMU(II)=1.-exp(-0.0115*YPLS)
                  F2(II)=1.-0.22*exp(-(REYT*REYT/36))
                  D(II)=-2.*VISCOS/dist_nor/dist_nor
                  E(II)=-2.*VISCOS*EXP(-0.5*YPLS)/dist_nor/dist_nor   
               END IF
    
	       END DO
	       END DO
      RETURN                                                            
       END
C
C***********************************************************************
C
      SUBROUTINE LOWRE_KOMEGA     
C                              
C COMPUTES SOURCE TERMS FOR NEAR WALL EFFECTS IN ALL EQUATIONS
C                              
C***********************************************************************
C                                                                        
      INCLUDE 'com2d'
C        
      SQR(VALUE)=VALUE*VALUE                                            
      VALUE=0.

C LOW RE VERSION OF KOMEGA MODEL 
C
C COMPUTE FMU,F1,F2 FOR THE K(TE) and OMEGA(OM) EQUATION
C
               DO J=JSOLV_B,JSOLV_E
               DO I=ISOLV_B,ISOLV_E
                  II=I+(J-1)*NI

                IF(KBLK(II).EQ.1) THEN

                  REYT=DENSIT*F(ISTE+II)/VISCOS/F(ISOM+II)
                  FMU(II)=(0.278+(REYT/8.)**4)/(1.+(REYT/8.)**4)
                  F1(II)=(0.025+REYT/6.)/(1+REYT/6.)
                  F2(II)=(0.1+REYT/2.7)/(1+REYT/2.7)/FMU(II)
               END IF
    
	       END DO
	       END DO
      RETURN                                                            
       END
C                              
C***********************************************************************
C                                                                       
      SUBROUTINE DISTFNC
C
C COMPUTES GENERALISED DISTANCE FUNCTION FOR ALL NODES
C                              
C***********************************************************************
C                                                                       
      INCLUDE 'com2d'                                               
C 

      DO I=1,NIJ
       SMIN(I)=1.E06
      END DO

        DO  nbound=1,nrbnda

        if(cbnd(nbound,2)(1:2).eq.'WA') then
 
        MFLAG1=0       
        MFLAG2=0       
        MFLAG3=0       
        MFLAG4=0 
      
        if(cbnd(nbound,1)(1:2).eq.'NO') then
           MFLAG1=1
           IFLAG=1
        end if
        if(cbnd(nbound,1)(1:2).eq.'SO') then 
           MFLAG2=1
           IFLAG=2
        end if
        if(cbnd(nbound,1)(1:2).eq.'EA') then 
           MFLAG3=1
           IFLAG=3
        end if
        if(cbnd(nbound,1)(1:2).eq.'WE') then 
           MFLAG4=1
           IFLAG=4
        end if

        do jw=nbnd(nbound,3),nbnd(nbound,4)
        do iw=nbnd(nbound,1),nbnd(nbound,2)

         IP=iw+(jw-1)*ni
         IB  =IP+( NI *  MFLAG1 -  NI* MFLAG2 +      MFLAG3 -     MFLAG4) 
         
         XWALL=F(ISX+IB)
         YWALL=F(ISY+IB)

         do j=2,njm
         do i=2,nim
         ii=i+(j-1)*ni
         xij=f(isx+ii)
         yij=f(isy+ii)   
         DX=xwall-xij       
         DY=ywall-yij         
         DIST=SQRT(DX*DX+DY*DY)

         IF (DIST.LE.SMIN(II)) THEN
            SMIN(II)=DIST
            IWAL(II)=IP
            IWFLAG(II)=IFLAG
         END IF

         END DO                     ! END OF NODE I SCAN LOOP OF THE FIELD
         END DO                     ! END OF NODE J SCAN LOOP OF THE FIELD

        END DO                      ! END OF IW LOOP FOR ALL WALLS
        END DO                      ! END OF JW LOOP FOR ALL WALLS

        END IF

        END DO                     ! END OF NBOUND LOOP FOR ALL BOUNDARIES
         

         RETURN
        END
