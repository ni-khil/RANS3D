
		SUBROUTINE SET_INIT_FLD 
C
C  SPECIFYING THE INITIAL FIELD FOR VELOCITY/TEMPERATURE AND TURBULENCE QUANTITIES
C
C********************************************************************
C
          INCLUDE 'com2d'      
C
        DIMENSION sinl(nmax),uinl(nmax),vinl(nmax),teinl(nmax),edinl(nmax),
     >           aominl(nmax),v2inl(nmax),finl(nmax),anutinl(nmax) 

C 
c   READING THE FLOW VARIABLE PROFILES FROM THE PIPE FLOW RESULTS
C
          OPEN(123,FILE='INPROF.DAT')
           READ (123,*)NPTS
           DO IPTS=1,NPTS
           READ(123,155)sinl(ipts),uinl(ipts),vinl(ipts),teinl(ipts),edinl(ipts)
     >                 ,aominl(ipts),v2inl(ipts),finl(ipts),anutinl(ipts)
           END DO
155       FORMAT(7E13.4) 
          CLOSE(123)
C
C    1D LINEAR INTERPOLATION TO COMPUTE THE GRID NODE VALUES
C    FOR THE JET IMPINGEMENT PROBLEM
C
           SCONS=0.                          
           DO J=NJ,NJ
           DO I=1,NI   
           II=I+(J-1)*NI
           IF(KN(II-NI).EQ.14) THEN
           SGIV=F(ISX+II)+SCONS
           CALL SEARCH(NPTS,SINL,SGIV,INTV)
           DY1=SINL(INTV+1)-SGIV
           DY2=SINL(INTV+1)-SINL(INTV)
           FACT=DY1/DY2
           F(ISU+II) =FACT*UINL(INTV)+(1.-FACT)*UINL(INTV+1)
           F(ISV+II) =FACT*VINL(INTV)+(1.-FACT)*VINL(INTV+1)
           F(ISTE+II)=FACT*TEINL(INTV)+(1.-FACT)*TEINL(INTV+1)
           F(ISED+II)=FACT*EDINL(INTV)+(1.-FACT)*EDINL(INTV+1)
           F(ISOM+II)=FACT*AOMINL(INTV)+(1.-FACT)*AOMINL(INTV+1)
           F(ISV2+II)=FACT*V2INL(INTV)+(1.-FACT)*V2INL(INTV+1)
           F(ISF+II)=FACT*FINL(INTV)+(1.-FACT)*FINL(INTV+1)
           F(ISNUT+II)=FACT*ANUTINL(INTV)+(1.-FACT)*ANUTINL(INTV+1)
c           write(119,*)i,sgiv,f(isu+ii),f(isv+ii),f(iste+ii),f(ised+ii),
c    >                  F(ISOM+II),F(ISV2+II),F(ISF+II),F(ISNUT+II)

           END IF 

           END DO
           END DO

           RETURN
           END
C
C*********************************************************************
C
          SUBROUTINE SEARCH(NDAT,PH,PHGIV,INTV)
C
C*********************************************************************
C
          INCLUDE 'com2d'      
C
          DIMENSION PH(NMAX)
C
          DO K=1,NDAT-1
          PH1=PH(K)-PHGIV
          PH2=PH(K+1)-PHGIV
          PRODUCT=PH1*PH2
          IF(PRODUCT.LE.0.) GO TO 100
          END DO
  100     INTV=K
          RETURN
          END
C
C****************************************************        
