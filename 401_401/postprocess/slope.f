C
C   CHECKING SLOPE FORMULA USING THREE POINT SCHEME
C
         DIMENSION F(364)
         
         NI = 11
         NJ = 11 
         ISX=  1
         ISY=ISX+NI*NJ
         ISVAR=ISY+NI*NJ
         WRITE(89,*) ISX,ISY,ISVAR

         DO J=1,NJ
         DO I=1,NI
         IJ=I+(J-1)*NI
         F(ISX+IJ)=FLOAT(I-1)/FLOAT(NI-1)
         F(ISY+IJ)=FLOAT(J-1)/FLOAT(NJ-1)
         F(ISVAR+IJ) = F(ISX+IJ)*F(ISX+IJ)+F(ISY+IJ)*F(ISY+IJ)
         WRITE(89,*)J,I,F(ISX+IJ),F(ISY+IJ),F(ISVAR+IJ)
         END DO
         END DO
         
         NDEG = 1
         
         DO J = 1,NJ
         IB = 1+(J-1)*NI
         INB = IB+1
         INNB = INB+1

         CALL INTERP(NDEG,IB,INB,INNB,ISX,ISY,ISVAR,F,GRAD)

         WRITE(978,*) IB,INB,INNB,ISVAR,GRAD

         END DO
         END
C
C*********************************************************************
        SUBROUTINE INTERP(NDEG,IB,INB,INNB,ISX,ISY,ISVAR,F,GRAD)
C
C*********************************************************************
C
C       THIS SUBROUTINE  FITS A STRAIGHT LINE OR A QUADRATIC FORM
C       GIVEN COORDINATES OF BOUNDARY AND NEAR BOUNDARY POINTS. 
C       IT RETURNS THE FIRST DERIVATIVE AT THE BOUNDARY NODE
C
C       INCLUDE 'postcom'

        DIMENSION F(364) 
        
        XB=F(ISX+IB) 
        XNB=F(ISX+INB)
        XNNB=F(ISX+INNB)
        YB=F(ISY+IB) 
        YNB=F(ISY+INB)
        YNNB=F(ISY+INNB)
        PHB=F(ISVAR+IB) 
        PHNB=F(ISVAR+INB)
        PHNNB=F(ISVAR+INNB)


C      NDEG = 1 IS FOR LINEAR INTERPOLATION BETWEEN IB AND INB NODES
C
        IF(NDEG.EQ.1) THEN
C               
C     COMPUTE THE SLOPE FROM THE VALUES AT WALL AND ONE NEARWALL POINT
C
        DIST = SQRT((XNB-XB)**2+(YNB-YB)**2)
        DIFVAR = PHB-PHNB          
        GRAD = DIFVAR/DIST
        WRITE(978,*) XB,XNB,YB,YNB,PHB,PHNB,GRAD
        ENDIF
        
        IF(NDEG.EQ.2) THEN
C               
C     COMPUTE THE SLOPE FROM THE VALUES AT WALL AND TWO OTHER NODES   
C
        DIF_1 = F(ISVAR+INB)-F(ISVAR+IB)
        DIF_2 = F(ISVAR+INNB)-F(ISVAR+IB)
        DX_1 = F(ISX+INB)-F(ISX+IB)
        DY_1 = F(ISY+INB)-F(ISY+IB)
        DS_1 = SQRT(DX_1*DX_1+DY_1*DY_1)
        DX_2 = F(ISX+INNB)-F(ISX+IB)
        DY_2 = F(ISY+INNB)-F(ISY+IB)
        DS_2 = SQRT(DX_2*DX_2+DY_2*DY_2)
        ANUM = DS_1*DS_1*DIF_2-DS_2*DS_2*DIF_1
        ADNM = DS_2*DS_1*DS_1-DS_1*DS_2*DS_2
        GRAD = ANUM/ADNM
        ENDIF
        RETURN
        END
