      SUBROUTINE OUT(II,IFLAG)       
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
                                                                        
      ISANOW = ISAN*MFLAG1 + ISAS*MFLAG2 + ISAE*MFLAG3 + ISAW*MFLAG4 
       
      ISPHI=ISPP+(NN-1)*NIJ
C     
C     DIRECTION 1 - ACROSS THE BOUNDARY
C     DIRECTION 2 - ALONG  THE BOUNDARY
C
          IF(TEMPSCH.EQ.1.OR.NTSTEP.EQ.1) GAMA=0.0
          IF(NTSTEP.GT.1.AND.TEMPSCH.EQ.2) GAMA=0.5

          DTF_1= (1.+GAMA)/UREF/DELT
          DTF_2= (1.+2.*GAMA)/UREF/DELT
          DTF_3= GAMA/UREF/DELT

c         write(99,*) ntstep,iblk,gamma, dtf_1,dtf_2,dtf_3

C        II  = INLEN(IM)*MFLAG1 + INLES(IM)*MFLAG2                      
C    >          + INLEE(IM)*MFLAG3 + INLEW(IM)*MFLAG4                      

         II_B = II + NI*MFLAG1 - NI*MFLAG2 + MFLAG3 - MFLAG4  
        
         IIM = II - NI*MFLAG1 + NI*MFLAG2 - MFLAG3 + MFLAG4          
          
         XB=F(ISX+II_B)
         XM=F(ISX+II)
         XMM=F(ISX+IIM)
         PHB=F(ISPHI+II_B)
         PHM=F(ISPHI+II)
         PHMM=F(ISPHI+IIM)
     
          IF(TEMPSCH.EQ.1.OR.NTSTEP.EQ.1) THEN  
       
           IF(NN.EQ.2) PHB_M=U(II_B,IBLK,1)
           IF(NN.EQ.3) PHB_M=V(II_B,IBLK,1)
           IF(NN.EQ.5) PHB_M=TE(II_B,IBLK,1)
           IF(NN.EQ.6) PHB_M=ED(II_B,IBLK,1)
           IF(NN.EQ.7) PHB_M=aom(II_B,IBLK,1)
           IF(NN.EQ.8) PHB_M=vsq(II_B,IBLK,1)
           IF(NN.EQ.10) PHB_M=anut(II_B,IBLK,1)
           IF(NN.EQ.11) PHB_M=scal(II_B,IBLK,1)
          END IF

         IF(TEMPSCH.EQ.2.AND.NTSTEP.GT.1) THEN
         IF(NN.EQ.2) THEN 
           PHB_M =U(II_B,IBLK,2)
           PHB_MM=U(II_B,IBLK,1)
         END IF
         IF(NN.EQ.3) THEN
           PHB_M =V(II_B,IBLK,2)
           PHB_MM=V(II_B,IBLK,1)
         END IF
         IF(NN.EQ.5) THEN
           PHB_M =TE(II_B,IBLK,2)
           PHB_MM=TE(II_B,IBLK,1)
         END IF
         IF(NN.EQ.6) THEN
           PHB_M =ED(II_B,IBLK,2)
           PHB_MM=ED(II_B,IBLK,1)
         END IF
         IF(NN.EQ.7) THEN
           PHB_M =AOM(II_B,IBLK,2)
           PHB_MM=AOM(II_B,IBLK,1)
         END IF
         IF(NN.EQ.8) THEN
           PHB_M =VSQ(II_B,IBLK,2)
           PHB_MM=VSQ(II_B,IBLK,1)
         END IF
         IF(NN.EQ.10) THEN
           PHB_M =ANUT(II_B,IBLK,2)
           PHB_MM=ANUT(II_B,IBLK,1)
         END IF
         IF(NN.EQ.11) THEN
           PHB_M =SCAL(II_B,IBLK,2)
           PHB_MM=SCAL(II_B,IBLK,1)
         END IF


         END IF
      
         IF(NN.EQ.9) THEN
         DTF_1=0
         DTF_2=0
         DTF_3=0
         END IF          
  
         A11=2.*XB + DTF_1*(XB*XB-XM*XM)
         A12=1.+DTF_1*(XB-XM)
         B1=PHB_M*DTF_2 -PHB_MM*DTF_3 - PHM*DTF_1
         
         A21=XM*XM-XMM*XMM
         A22=XM-XMM
         B2=PHM-PHMM

C        write(223,*) ntstep, iblk,nn,im,PHM,PHMM,PHB_M,PHB_MM
C        write(222,*) ntstep, iblk,nn,im,A11,A12,B1,A21,A22,B2
         
         A=(B1*A22-B2*A12)/((A11*A22-A21*A12)+small)
         B=(A11*B2-A21*B1)/((A11*A22-A21*A12)+small)
         C= PHM-A*XM*XM-B*XM

         PHB = A*XB*XB+B*XB+C
	 if (CONVECTIVE_BC) then
           F(ISPHI+II_B)=PHB
	 else 
           F(ISPHI+II_B)=F(ISPHI+II)
	 end if


         F(ISANOW+II)=0.
  
       RETURN                                                            
       END                                                               
 
