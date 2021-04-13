C***********************************************************************
C
      SUBROUTINE BLOKCON
C
C     CORRECT EXIT VELOCITIES AND CONVECTIVE MASS FLUXES AT CELL 
C     FACES WITH OUTFLOW STATE TO ENSURE OVERALL MASS CONSERVATION
C     FOR THE BLOCK         
C
C***********************************************************************
C
      INCLUDE 'com2d'                                              
C
	if (iblk .eq. 1) then
            FLWOUT=0.                                                         
            FLWIN= 0.            
	end if
C
C  LOOP TO BE SET UP FOR ALL THE BOUNDARIES OF THE BLOCK
C
             do nbound=1,nrbnda

               if (cbnd(nbound,2)(1:2).eq. 'IN' .or.
     >                      cbnd(nbound,2)(1:2).eq. 'OU') then
               if (cbnd(nbound,1)(1:2).eq. 'WE') then
               
               is=nbnd(nbound,1)
               ie=nbnd(nbound,2)
               js=nbnd(nbound,3)
               je=nbnd(nbound,4)
              
                 DO j=js,je
                 DO i=is,ie 
                 II=I+(j-1)*NI
                 CWII=F(ISCW+II)
                 FLWIN=FLWIN+MAX(0.,CWII)
                 FLWOUT=FLWOUT-MIN(0.,CWII)
                 END DO
                 END DO
		end if

               if (cbnd(nbound,1)(1:2).eq. 'EA') then
               
               is=nbnd(nbound,1)+1
               ie=nbnd(nbound,2)+1
               js=nbnd(nbound,3)
               je=nbnd(nbound,4)
              
                 DO j=js,je
                 DO i=is,ie
                 II=I+(j-1)*NI
                 CWII=F(ISCW+II)
                 FLWIN=FLWIN-MIN(0.,CWII)
                 FLWOUT=FLWOUT+MAX(0.,CWII)
                 END DO
                 END DO
		end if
                
          if (cbnd(nbound,1)(1:2).eq. 'SO') then
               
               js=nbnd(nbound,3)
               je=nbnd(nbound,4)
               is=nbnd(nbound,1)
               ie=nbnd(nbound,2)

               DO j=js,je
               DO i=is,ie 
                 II=I+(j-1)*NI
                 CSII=F(ISCS+II)
                 FLWIN=FLWIN+MAX(0.,CSII)
                 FLWOUT=FLWOUT-MIN(0.,CSII)
                END DO
                END DO
		end if

            if (cbnd(nbound,1)(1:2).eq. 'NO') then
               js=nbnd(nbound,3)+1
               je=nbnd(nbound,4)+1
               is=nbnd(nbound,1)
               ie=nbnd(nbound,2)

               DO j=js,je
               DO i=is,ie
                 II=I+(j-1)*NI
                 CSII=F(ISCS+II)
                 FLWIN=FLWIN-MIN(0.,CSII)
                 FLWOUT=FLWOUT+MAX(0.,CSII)
                END DO
                END DO
		end if
		end if ! End of INF Check
	   end do !end of nbound loop
C
C EVALUATION OF THE CORRECTION FACTOR AS MULTIPLIER
C
	if (iblk .eq. nblock) then
        if (flwout .le. 0. .or. flwin .le. 0. )then
           FLORAT=1.
        else
           FLORAT=FLWIN/FLWOUT
        end if
	end if

       IF(MOD(NITER,1).EQ.0) THEN
       WRITE(6,*)'FLOWIN=',flwin,'FLOWOUT=',FLWOUT,'FACT=',FLORAT
       END IF
	if (ibswp .eq. 1 .and. ntstep .eq.1) florat=1.


        IF (.NOT. OVACOR) return

C
C --- CORRECT  CONVECTIVE FLUXES AT THE BOUNDARY FACES
C     BEHAVING AS OUT FLOW ONLY
C   

C
C  LOOP TO BE SET UP FOR ALL THE BOUNDARIES OF THE BLOCK
C
             do nbound=1,nrbnda

               if (cbnd(nbound,2)(1:2).eq. 'IN' .or.
     >                      cbnd(nbound,2)(1:2).eq. 'OU') then
               if (cbnd(nbound,1)(1:2).eq. 'WE') then
               
               is=nbnd(nbound,1)
               ie=nbnd(nbound,2)
               js=nbnd(nbound,3)
               je=nbnd(nbound,4)
              
                 DO j=js,je
                 DO i=is,ie 
                 II=I+(j-1)*NI
                 CWII=F(ISCW+II)
                 RAT=1.
                 IF(CWII.LE.0.) RAT=FLORAT
		 F(ISU+II-1)=F(ISU+II-1)*RAT
		 F(ISV+II-1)=F(ISV+II-1)*RAT
                 F(ISCW+II)=F(ISCW+II)*RAT
                 END DO
                 END DO
		end if
                
               if (cbnd(nbound,1)(1:2).eq. 'EA') then
               
               is=nbnd(nbound,1)+1
               ie=nbnd(nbound,2)+1
               js=nbnd(nbound,3)
               je=nbnd(nbound,4)
              
                 DO j=js,je
                 DO i=is,ie
                 II=I+(j-1)*NI
                 CWII=F(ISCW+II)
                 RAT=1.
                 IF(CWII.GE.0.) RAT=FLORAT
		 F(ISU+II)=F(ISU+II)*RAT
		 F(ISV+II)=F(ISV+II)*RAT
                 F(ISCW+II)=F(ISCW+II)*RAT
                 END DO
                 END DO
		end if

          if (cbnd(nbound,1)(1:2).eq. 'SO') then
               
               js=nbnd(nbound,3)
               je=nbnd(nbound,4)
               is=nbnd(nbound,1)
               ie=nbnd(nbound,2)

               DO j=js,je
               DO i=is,ie 
                 II=I+(j-1)*NI
                 CSII=F(ISCS+II)
                 RAT=1.
                 IF(CSII.LE.0.) RAT=FLORAT
		 F(ISU+II-NI)=F(ISU+II-NI)*RAT
		 F(ISV+II-NI)=F(ISV+II-NI)*RAT
                 F(ISCS+II)=F(ISCS+II)*RAT
                END DO
                END DO
		end if

            if (cbnd(nbound,1)(1:2).eq. 'NO') then
               js=nbnd(nbound,3)+1
               je=nbnd(nbound,4)+1
               is=nbnd(nbound,1)
               ie=nbnd(nbound,2)

               DO j=js,je
               DO i=is,ie
                 II=I+(j-1)*NI
                 CSII=F(ISCS+II)
                 RAT=1.
                 IF(CSII.GE.0.) RAT=FLORAT
		 F(ISU+II)=F(ISU+II)*RAT
		 F(ISV+II)=F(ISV+II)*RAT
                 F(ISCS+II)=F(ISCS+II)*RAT
                END DO
                END DO
		end if
		end if ! End of INF Check
	   end do !end of nbound loop

        RETURN                                                            
      END                                                               
