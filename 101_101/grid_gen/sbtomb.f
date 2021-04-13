C
C  THIS IS A ROUTINE TO CONVERT A SINGLE BLOCK GRID DATA TO MULTIBLOCK
C

       PARAMETER (NXS=  101,NYS = 101)
       PARAMETER (NXM=  102,NYM = 102)
C
       DIMENSION X(NXS,NYS),Y(NXS,NYS)
       DIMENSION XB(NXM,NYM),YB(NXM,NYM)
       DIMENSION IC(10),JC(10)
       LOGICAL PERIODI,PERIODJ

        OPEN(UNIT=8,FILE='sbgrid',FORM='UNFORMATTED')
        OPEN(UNIT=10,FILE='sbtomb.inp')
        OPEN(UNIT=11,ACCESS='DIRECT',FILE='grid',
     &  FORM='UNFORMATTED',recl=(2+(nxm*nym*2))*4)
C                                             
C      READ THE SINGLE BLOCK GRID
C
           READ(8) NI,NJ
           READ(8)((X(I,J),I=1,NI),J=1,NJ),
     >          ((Y(I,J),I=1,NI),J=1,NJ)

C                                             
C      READ INPUT DATA FOR MULTIBLOCK GENERATION
C
                                            
       READ (10,*)  
       READ (10,*) NSEGI,NSEGJ

       READ (10,*)  
       READ(10,*) (IC(I),I=1,NSEGI)
       READ (10,*)  
       READ(10,*) (JC(J),J=1,NSEGJ)

       READ (10,*)  
       READ (10,*) PERIODI,PERIODJ



       NBLK=(NSEGI-1)*(NSEGJ-1)

       DO JBND = 1, NSEGJ-1

         DO IBND = 1, NSEGI-1

         IBLK=IBND+(JBND-1)*(NSEGI-1)

         IEND=IC(IBND+1)-IC(IBND)+3  
         IF(.NOT.PERIODI.AND.IBND.EQ.1 ) IEND=IC(IBND+1)-IC(IBND)+2  
         IF(.NOT.PERIODI.AND.IBND.EQ.NSEGI-1) IEND=IC(IBND+1)-IC(IBND)+2  
         IF(NSEGI .EQ. 2) IEND = NI
         JEND=JC(JBND+1)-JC(JBND)+3  
         IF(.NOT.PERIODJ.AND.JBND.EQ.1 ) JEND=JC(JBND+1)-JC(JBND)+2  
         IF(.NOT.PERIODJ.AND.JBND.EQ.NSEGJ-1) JEND=JC(JBND+1)-JC(JBND)+2  
         IF(NSEGJ .EQ. 2) JEND = NJ
            
            DO J=1,JEND
            DO I=1,IEND  
            II=I+IC(IBND)-2
            JJ=J+JC(JBND)-2

            IF(PERIODI .AND. IBND.EQ.1 .AND. I.EQ.1) II=NI-1
            IF(PERIODJ .AND. JBND.EQ.1 .AND. J.EQ.1) JJ=NJ-1
            IF(PERIODI .AND. IBND.EQ.NSEGI-1 .AND. I.EQ.IEND) II=2   
            IF(PERIODJ .AND. JBND.EQ.NSEGJ-1 .AND. J.EQ.JEND) JJ=2 
            IF(.NOT. PERIODI .AND. IBND .EQ. 1)II = I 
            IF(.NOT. PERIODJ .AND. JBND .EQ. 1)JJ = J 
            XB(I,J)=X(II,JJ)
            YB(I,J)=Y(II,JJ)
            END DO
            END DO
            WRITE(11,rec=iblk) IEND,JEND,((XB(I,J),I=1,IEND),J=1,JEND),
     >            ((YB(I,J),I=1,IEND),J=1,JEND)
                                                  
         END DO

       END DO 

        STOP
        END
