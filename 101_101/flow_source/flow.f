c**********************************************************************
C
C   THIS IS THE PROBLEM DEPENDENT PART OF THE PROGRAM RANS3D (2D Version)
C
C 
C***********************************************************************
        PROGRAM FLOW                                                      
C***********************************************************************
C                                            
      INCLUDE 'com2d'                                               
C                                            
C  FILE UNIT NUMBER 22 TO READ CONTROL INDICES
C  AND PHYSICAL INPUT DATA
C  INDEPENDENT OF THE BLOCKS

C
        CALL FILES(1)                                                     
        
       OPEN(22,FILE='CONREAD')

        CALL CONSTRD
c
        CALL SOLVE_INDEX
        
C 
C   THE TIME INTEGRATION LOOP STARTS HERE
C
        tinst=delt*(ntbeg-1)

        DO NTSTEP=NTBEG,NTEND

            OPEN(18,FILE='DTCONV',ACCESS='APPEND')
            OPEN(19,FILE='DTMONI',ACCESS='APPEND')
            OPEN(20,FILE='DTCOEF',ACCESS='APPEND')
            OPEN(21,FILE='DTCM',ACCESS='APPEND')
            OPEN(25,FILE='DTWAKEL',ACCESS='APPEND')
            OPEN(27,FILE='DTSEP',ACCESS='APPEND')
            OPEN(29,FILE='DTSTAG',ACCESS='APPEND')
            OPEN(24,FILE='DTERR',ACCESS='APPEND')
            OPEN(35,FILE='DTCLCD',ACCESS='APPEND')
C
C  THE ITERATIVE SWEEP LOOP STARTS HERE
C
	DO IBSWP=NBSWP_ST,NBSWP_END
C
C   ARRANGE OPENING OF APPROPRIATE FILES FOR STORAGE
C
	 IF (DOUBLE_PREC) THEN

C REQUIRED FILE SIZE WHEN CODE EXECUTED IN DOUBLE PRECISION

	    lenc=(15+idim1*18)*2*4
	    lenp=idim1*9*2*4
	    lenb=nmax*28*2*4
            lengv=(5+idim1*2)*2*4
	    lengg=idim1*4*2*4

	 ELSE

C REQUIRED FILE SIZE WHEN CODE EXECUTED IN SINGLE PRECISION

	    lenc=(15+idim1*18)*4
	    lenp=idim1*9*4
	    lenb=nmax*28*4
            lengv=(5+idim1*2)*4
	    lengg=idim1*4*4

	END IF ! END OF PRECISION CHECK
            
            OPEN(28,FILE='PREVRES',ACCESS='DIRECT',FORM='UNFORMATTED'
     >      ,recl=lenp)

            OPEN(8,FILE='CONTRES',ACCESS='DIRECT',FORM='UNFORMATTED'
     >      ,recl=lenc)

            OPEN(9,file='BLKTRAN',ACCESS='DIRECT',FORM='UNFORMATTED'
     >      ,recl=lenb)

            OPEN(116,FILE='FLOW_CUT',ACCESS='APPEND')
            OPEN(26,FILE='FLOMON',ACCESS='APPEND')
            OPEN(36,FILE='TURMON',ACCESS='APPEND')
C
C  INITIALISE THE GLOBAL RESIDUES
C
	   do ivar=1,11
	      globrs(ivar)=0.
	   end do

C  
C  REWINDING OF CERTAIN FILES TO FRESH READ/WRITE
C                                            
           IF(IBSWP.GT.1 .OR. NTSTEP.GT.1) CALL FILES(2)
C 
C  BLOCK ITERATION LOOP STARTS HERE
C                                            
	   DO IBLK=1,NBLOCK
	   
C
C READ BOUNDARY INDICES AND FLOW CONDITIONS FOR EACH BLOCK
C
	     CALL BLOKRD
C
C COMPUTE STARTING ADDRESS FOR EACH VARIABLE
C
             CALL ISTADR
C
C COMPUTE  APPROPRIATE CELL INDICES REQUIRED FOR BOUNDARY CONDITIONS
C
             CALL BINDEX 
C
C READING THE GRID COORDINATES                                         
C
             CALL GEOINP
C
C COMPUTE  GEOMETRICAL QUANTITIES RELATED TO THE GRID               
c
C

C GET THE GRID LOCATION AT THE PREVIOUS TIME STEP
C IN ORDER TO COMPUTE THE GRID VELOCITY NUMERICALLY
C
             CALL GEOMET  

C
C COMPUTE  GEOMETRICAL INTERPOLATION FACTORS                        
c
 	     CALL BINTER

C
C COMPUTE  MINIMUM DISTANCE OF EACH CELL FROM WALL BOUNDARY
c
             CALL DISTFNC
C
C COMPUTE  TURBULENCE MODEL CONSTANTS (CLOSURE COEFFICIENTS)                   
c
	     CALL TURBCONS
c
C COMPUTE OR ASSIGN INITIAL FIELD VALUES AND CORRESPONDING MASS FLUXES
C FOR THE START OF A SWEEP AND THE START OF THE TIME LOOP
C OTHERWISE(ELSE) FOR CONTINUATION RUN READ STORED RESULTS FROM DISC
C
             IF(IBSWP .EQ. 1 .AND. NTSTEP .EQ.1) THEN

	     IF (START_FIELD) CALL RDBLOK
C 
C ASSIGN INITIAL FIELD VALUES
C
               CALL INITIA

               CALL BNDRY_CONDN
c
c SUBROUTINES CWCS AND INFLUX COMPUTE THE MASS FLUXES
C SUBROUTINE EDYVIS COMPUTES THE EDDY VISCOSITY FIELD
C
C START_FIELD=.TRUE. : SOME SOLVED FIELD VALUES ARE GIVEN AS ITERATION STARTS
C            =.FALSE. OTHERWISE
C SAMEGRID_LEVEL =.TRUE.: THE STARTING FIELD VALUES ARE AT THE SAME GRID SIZE
C                =.FALSE. : OTHERWISE
C
	        IF (START_FIELD) THEN

                   IF (.NOT.SAMEGRID_LEVEL)THEN
                    
                    CALL CWCS

                    CALL EDYVIS

		   END IF ! END OF SAMEGRID_LEVEL CHECK

	        ELSE

	            CALL CWCS 

                    CALL EDYVIS

	        END IF ! END OF START_FIELD CHECK
c		   
C STORE FIELD VALUES AT FIRST TIME STEP
C
	        IF (DELT .LT. 1) THEN 
		     CALL STORVAR
                     IF (START_FIELD) CALL STORVAR_CWCS
		END IF

             ELSE !   ALTERNATIVE OF IF LOOP FOR CONTINUATION RUN
C
C  READ FIELD VALUES FROM THE FILE <CONTRES> FOR PREVIOUS TIME STEP / SWEEP
C
	          CALL RDBLOK

                  IF (IBSWP .EQ. 1) CALL BNDRY_CONDN
C
C CALL GET_STORVAR IF THE JOB HAS BEEN STOPPED AT A PARTICULAR TIME STEP(N)
C BEFORE CONVERGENCE (SAY DURING POWER FAILURE OR SYSTEM CRASH ETC.). 
C GET_STORVAR GETS THE FIELD DATA FOR THE (N-1) AND (N-2) TIME STEP VALUES
C NOTE : GET_STORVAR IS CALLED ONLY FOR THE UNSTEADY RUN
C
          IF (STOP_JOB .AND. DELT .LT. 1. .AND. IBSWP .EQ. NBSWP_ST) 
     >                  CALL GET_STORVAR

	          IF (IBSWP.EQ. 1 .AND. DELT .LT. 1) THEN 
C
C STORE FIELD VALUES AT THE PREVIOUS TIME STEP
C
                     CALL STORVAR

                     CALL STORVAR_CWCS
C
C ADDS THE FIELD VALUES AT EACH TIME STEP FOR TIME AVERAGING
C
	             CALL PHAVG

             END IF  ! END OF IBSWP.EQ.1 LOOP

             END IF       ! END OF IF LOOP FOR CONTINUATION RUN
c
C WRITE PROBLEM DETAILS WITH INPUT DATA ON MONITOR FILE
C CALLING TITLE AND INFO SUBROUTINES
C 
             IF (IBLK.EQ.1 ) CALL TITLE
c 
             CALL INFORM 

C
C READ FIELD DATA AT CUT PLANES STORED FROM THE NEIGHBOURING BLOCKS
C
             IF (IBSWP.NE.1) CALL CUTRD
C
C SEQUENTIAL SOLUTION OF THE RANS EQUATIONS STARTS HERE
C
             CALL RANSOL                                                     
c
c COMPUTE GLOBAL RESIDUE
C
             do ivar=1,11
	         globrs(ivar)=globrs(ivar)+rmom(ivar)*rmom(ivar)
	       end do
C
C  STORES FIELD VALUES AT MONITORING LOCATION
C
	        if(iblk.eq.ibmon) call moni_loc
C
C STORES FIELD VALUES IN THE FILE <CONTRES> 
C
             CALL WTBLOK                                          
C
C
C STORES FIELD VALUES IN THE FILE <POSTRES> AT SPECIFIED TIME STEPS
C
       if (ntstep .ge. ntsb .and. ntstep .le. ntse) CALL WTTIME
C
C WRITE FIELD DATA AT CUT PLANES FOR THE NEIGHBOURING BLOCKS
C
	       IF (NOCUT .NE. 0) CALL CUTWR
C
C
C WRITE FIELD DATA  ON MONITOR FILE AS DISPLAY               
C
C             IF(IDEBUG.EQ.1)   CALL FINAL

C TO STORE THE GRID POINT AT THE CURRENT INSTANT
C REQUIRED IF WE NEED TO COMPUTE THE GRID VELOCITY NUMERICALLY
C
C
	   END DO  !  END OF BLOCK LOOP  
C
C COMPUTING GLOBAL RESIDUE AND  MAXIMUM NORMALISED RESIDUE
C FROM THE RESIDUE OF U- V-MOMENTUM, CONTINUITY AND SCALAR EQUATION	
C   
           do ivar=1,11
	   globrs(ivar)=sqrt(globrs(ivar))/gnorm(ivar)
	   end do
	 
           resmax=max(globrs(1),globrs(2),globrs(3),globrs(11))
C
C WRITING THE RESIDUE OF ALL THE EQUATIONS ON FILE < RESID.DAT>

	  if (ntstep .le. ntend) then
	  if (ibswp .eq. 1) write(109,*)'time step= ',ntstep, 't =',tinst
	   write(109,*)ibswp,(globrs(ivar),ivar=1,11) 
	   if (resmax .le. eps .and. ibswp .gt. 1) go to 999
	  end if

	close(8)
	close(9)
        CLOSE(28)
        CLOSE(116)
        CLOSE(26)
        CLOSE(36)
C
        END DO         !END OF DO LOOP  FOR ITERATION SWEEP 
C                                                    
999     CONTINUE

         IF(IDEBUG.EQ.1)   CALL FINAL

C
C STOP_JOB IS MADE FALSE AND NBSWP_ST IS SET TO 1 AFTER CONVERGENCE OF 
C THAT TIME STEP (N) SO THAT THE NEXT TIME STEP (N+1) PROCEEDS IN 
C THE NORMAL WAY 
C
	IF (STOP_JOB) then
           STOP_JOB = .FALSE.
	   NBSWP_ST=1
	END IF

C
C  WRITE THE HEADINGS REQUIRED FOR PRINTED DATA
C
	if (ntstep .eq. 1) then
	 write(18,*)'ntstep,    tinst,   ibswp,  resmax'
	 write(19,*)'tinst,pitch,umon,    vmon,   pmon'
	 write(20,*)'tinst,  pitch     clp,     clf,    cdp,    cdf'
         write(21,*)'tinst, pitch,    cmz'
         write(25,*)'tinst,     wakel'
         write(27,*)'tinst,  no. sep pts    theta(i)   for block two'
         write(29,*)'tinst, no. stag pts     theta(i)   for block one'
	 write(35,*)'tinst,    cl,       cd'
         END IF 

	write(18,*)ntstep,tinst,ibswp,resmax
	write(19,*)tinst,pitch,umon,vmon,pmon
        write(20,88)tinst,pitch*180/pi,adfp_pr(2),adff_pr(2),adfp_pr(1),adff_pr(1)
        write(35,*)tinst,adfp_pr(2)+adff_pr(2),adfp_pr(1)+adff_pr(1)
88      format(6(e12.5,1x))
        write(21,*)tinst, pitch*180/pi,adm_pr(3)
        write(25,*)ntstep,tinst,wakel
        write(29,*)tinst,nsep(1),(theta(1,nsp),nsp=1,nsep(1))

        if (delt .lt. 1) then

       uerr_max=uerr(1)
       verr_max=verr(1)
       perr_max=perr(1)
       teerr_max=teerr(1)
       ederr_max=ederr(1)
       aomerr_max=aomerr(1)
       verr_max=vsqerr(1)
       anuterr_max=anuterr(1)
       scalerr_max=scalerr(1)

          DO nb=2,NBLOCK
           uerr_max=max(uerr_max,uerr(nb))
           scalerr_max=max(scalerr_max,scalerr(nb))
          END DO

        write(24,*)ntstep,uerr_max,verr_max,perr_max,teerr_max,ederr_max,
     >             aomerr_max,vsqerr_max,anuterr_max,scalerr_max
      
         err_max=max(uerr_max,verr_max,perr_max,scalerr_max)

         IF(err_max.le.EPS*DELT) GO TO 666        

         END IF  ! END OF UNSTEADY CHECK


C   TIME INSTANT IS INCREASED BY DELT
C
        tinst=tinst+delt

            CLOSE(8)
            CLOSE(18)
            CLOSE(19)
            CLOSE(20)
            CLOSE(24)
            CLOSE(25)
            CLOSE(27)
            CLOSE(29)
            CLOSE(35)

         END DO       ! END OF TIME LOOP 

666     continue

        IF(WRITE_EXIT) CALL WRITE_EXITPROF
C
C   CLOSE ALL FILES 
C
        CALL FILES(0)

         STOP
                                 
         END                                                               
C
C**********************************************************************        
C
        SUBROUTINE CONSTRD 
C
C 	THIS ROUTINE PROVIDES ALL CONTROL INDICES 
C       AND PHYSICAL INPUT DATA TO THE PROBLEM
C       BLOCK INDEPENDENT INPUT DATA FOR MULTIBLOCK VERSION
C
C***********************************************************************
C 
	INCLUDE 'com2d'
C
            PI=4.*atan(1.)

               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)NBLOCK
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)DOUBLE_PREC
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)EPS,NBSWP_ST,NBSWP_END,IDEBUG,STOP_JOB
C
C NBSWP_ST IS ENSURED TO BE GREATER THAN ONE IF THE JOB HAS BEEN STOPPED
C INTERMEDIATELY (I.E. IN THE MIDDLE OF A PARTICULAT TIME STEP) 
C BY AN EXTERNAL SOURCE 
C
          if (stop_job .and. NBSWP_ST .eq. 1) NBSWP_ST=2

               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)IMON,JMON,IBMON
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)KATO,AXISYM_X,AXISYM_Y,SWIRL,OVACOR
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)START_FIELD,SAMEGRID_LEVEL,NATCONV,INIT_FLD,BUOY_SORCE
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)CONVECTIVE_BC,WRITE_EXIT
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)(GNORM(I),I=1,11)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)(SOLVE(I),I=1,4),SOLVE(11)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)(NSWP(I),I=1,11)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)LINKPV
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)(RELAX(I),I=1,14)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)(NSOLV(I),I=1,11)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)(NDSCH(I),I=1,11)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)(FDEFER(I),I=1,11)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)DENSIT,VISIN_FACT,ALREF,UREF,VREF,TREF,TEREF,ALPHA
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)TURBMOD
               READ(22,*)
               READ(22,*)
               READ(22,*) 
               READ(22,*)
               READ(22,*)C_EPS,REAL_IND
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)REYNOLDS,PRLAM,RAYLEIGH
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)(PR(I),I=1,11)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)(SPFACT(I),I=1,11)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)DELT,NTBEG,NTEND,tbeg
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)NTSB,NTSE,NTSINT
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)NTAVG
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)TEMPSCH
               READ(22,*)
               READ(22,*)
               READ(22,*)
               READ(22,*)
        
               IF (TEMPSCH .EQ. 1) THEN
                  NLEVEL=1
               ELSE
                  NLEVEL=2
               END IF
C
C  EVALUATE VISCOS FROM GIVEN REYNOLDS NUMBER
C
        URES=SQRT(UREF*UREF+VREF*VREF)
C
        VISCOS =(URES*ALREF)/REYNOLDS


      if (turbmod.gt.7) then
       write(215,*)'Turbulence models for turbmod > 7 not installed '
       stop
      endif

                 alpha=alpha*pi/180.

              if (alpha .gt. 0) then
                 utemp=uref
                 uref=utemp*cos(alpha)
                 vref=utemp*sin(alpha)
              end if
C     
              RETURN
              END
C
C***************************************************************************
C
        SUBROUTINE BLOKRD
C
C       THIS ROUTINE PROVIDES ALL CONTROL INDICES 
C       AND PHYSICAL INPUT DATA TO THE PROBLEM
C       BLOCK DEPENDENT INPUT DATA FOR MULTIBLOCK VERSION
C
C***********************************************************************
C
       include 'com2d'
C
         if (iblk .eq. 1) then
         do iskp=1,6 
         read(23,*)
           end do
           end if
          read(23,*)
          read(23,*)
          read(23,*)
          read(23,*)
          read(23,*)iblkno,nocut
          read(23,*)
          read(23,*)
          read(23,*)
          read(23,*)
          read(23,*)ni,nj
          read(23,*)
          read(23,*)
          read(23,*)
          read(23,*)
          read(23,*)nrbnda
          do nbound=1,nrbnda
             read(23,*)
             read(23,*)(cbnd(nbound,j),j=1,2)
             if (cbnd(nbound,2)(1:2) .eq. 'CU') then
             read(23,*)
             read(23,*)nbblk(nbound),nbcutn(nbound),ndir(nbound)
             end if
             read(23,*)
             read(23,*)(nbnd(nbound,j),j=1,4)
             read(23,*) 
             read(23,*)uin(nbound),vin(nbound),win(nbound),tin(nbound)
          end do
          read(23,*)
          read(23,*)
          read(23,*)
          read(23,*)
          read(23,*)isolv_b,isolv_e,jsolv_b,jsolv_e
          read(23,*)
          read(23,*)
          read(23,*)
          read(23,*)
          read(23,*)nsyme,nsymw,nsymn,nsyms
          read(23,*)
          read(23,*)
          read(23,*)
          read(23,*)
          read(23,*)(specl(j),j=1,11)
          read(23,*)
          read(23,*)
          read(23,*)
C
C
        return
        end
C
C***************************************************************************
C
      SUBROUTINE TITLE                                                       
C
C     PRINTING TITLES AND OTHER IMPORTANT INFORMATION                   
C
C**********************************************************************
C
      INCLUDE 'com2d'   
C                  
      WRITE(6,20)                                                       
   20 FORMAT(/2X,'---------------------------------------------------'  
     >      ,/2X,'LAMINAR FLOW IN A LID DRIVEN SQUARE CAVITY '            
     >      ,/2X,'---------------------------------------------------'/)
C
       WRITE(6,*) ' REYNOLDS-NUMBER: ',REYNOLDS 
       WRITE(6,*) ' RAYLEIGH NUMBER: ',RAYLEIGH,' VISCOSITY-', VISCOS
       WRITE(6,*) ' PRANDTL NUMBER: ',PRLAM
  
C                           
      RETURN                                                            
      END
C
C***********************************************************************
C
      SUBROUTINE FILES(N)                                                       
C
C     CHANNEL 6 FOR STANDARD OUTPUT LISTING                             
C             7 TO READ GRID-COORDINATES                            
C             8 TO READ/WRITE COMPUTED RESULTS FOR EACH BLOCK    
C            17 TO WRITE RESIDUES AND MONITORING LOCATION PRESSURE
C                                   
C**************************************************************************
C
      INCLUDE 'com2d'
C
      IF(N.EQ.1) THEN
C
        IF (DOUBLE_PREC) THEN

C CODE EXECUTED IN DOUBLE PRECISION
C COMPUTES THE SIZE OF THE REQUIRED DATA     

            leng=(2+idim1*2)*2*4
            lent=(15+idim1*11)*2*4
            lenavg=(4+idim1*11)*2*4
          
        ELSE

C CODE EXECUTED IN SINGLE PRECISION
C COMPUTES THE SIZE OF THE REQUIRED DATA     
C
            leng=(2+idim1*2)*4           ! GRID COORDINATES
            lent=(15+idim1*11)*4         ! FULL FIELD DATA FOR POSTPROCESSING 
            lenavg=(4+idim1*11)*4!FULL FIELD DATA FOR AVERAGED FLOW VARIABLES  
 
        END IF ! END OF PRECISION CHECK

            OPEN(UNIT = 6,FILE='MONITOR')
            OPEN(26,FILE='FLOMON',ACCESS='APPEND')
            OPEN(36,FILE='TURMON',ACCESS='APPEND')
            OPEN(7,FILE='GRID',ACCESS='DIRECT',FORM='UNFORMATTED'
     >      ,recl=leng)
            OPEN(10,file='POSTRES',ACCESS='DIRECT',FORM='UNFORMATTED'
     >      ,recl=lent)
            OPEN(11,file='PHAVG',ACCESS='DIRECT',FORM='UNFORMATTED'
     >      ,recl=lenavg)
            OPEN(17,FILE='BLOKRESID')
            OPEN(23,FILE='BLKREAD')
            OPEN(109,FILE='RESID',ACCESS='APPEND')
            OPEN(116,FILE='FLOW_CUT',ACCESS='APPEND')
       
        END IF    ! IF LOOP ENDS FOR N=1
C                                 
        IF(N.EQ.2) THEN
      
        REWIND(23)
        REWIND(6)
        REWIND(17)
       
C       if (ibswp.le.1)REWIND(109)
C       if (ibswp.le.1)REWIND(116)
       
        END IF    ! IF LOOP ENDS FOR N=2
       
        IF(N.EQ.0) THEN

        CLOSE(6)                             
        CLOSE(7)
        CLOSE(8)
        CLOSE(9)
        CLOSE(17)                                                        
        CLOSE(22)
        CLOSE(23)
        CLOSE(24)
        CLOSE(26)
        CLOSE(36)
        CLOSE(109)                                                        
        CLOSE(116)                                                        
      
        END IF    ! IF LOOP ENDS FOR N=0

        RETURN                                                            
        END      
C
C************************************************************************
C                                                            
      SUBROUTINE INITIA 
C                                                     
C     SET UP INLET VALUES                                            
C     INITIALIZE VARIABLE FIELD 
C                                         
C************************************************************************
C                                                     
      INCLUDE 'com2d'
C
	visin=visin_fact*viscos

	IF (.NOT. START_FIELD) THEN

 	do i=1,nij
 	 f(isu+i)=uref 
 	 f(isv+i)=vref
         f(iss+i)=tref 

	 if (TURBMOD.GT.0) then
	    tein=teref*(uref*uref+vref*vref)
	    edin=cmu*densit*tein*tein/visin
	    f(iste+i)=tein
	    f(ised+i)=edin
	    f(isom+i)=edin/(tein+small)/CMU
            f(isv2+i)=(2./3.)*tein
            f(isf+i)=0.
            f(isnut+i)=visin
            f(isvis+i)=viscos+visin
	 end if

 	end do
       	END IF ! END OF START_FIELD CHECK

      RETURN                                                            
      END                                                                  
C                                         
C************************************************************************
C                                                     
      SUBROUTINE BNDRY_CONDN

C                                                     
C     SET UP INLET VALUES                                            
C     INITIALIZE VARIABLE FIELD 
C                                         
C************************************************************************
C                                                     
      INCLUDE 'com2d'
C
	visin=visin_fact*viscos

        IF (.NOT. START_FIELD .AND. NTSTEP .EQ.1 ) THEN

C  INITIALIZE ALL THE BOUNDARY CONDITION FOR THE FRESH RUN
C  WHEN NO INITIAL VALUES ARE KNOWN

	do nbound=1,nrbnda

	       if (alpha .gt. 0) then
                   utemp=uin(nbound)
                   uin(nbound)=utemp*cos(alpha)
                   vin(nbound)=utemp*sin(alpha)
	       end if

	  if (cbnd(nbound,1)(1:2).eq. 'WE') then
	       i=nbnd(nbound,1)-1
	       js=nbnd(nbound,3)
	       je=nbnd(nbound,4)
               DO j=js,je           
                 II=I+(j-1)*NI
                 F(ISU+II)=Uin(nbound)
                 F(ISV+II)=Vin(nbound)
		 f(isvis+ii)=viscos
                 f(iss+ii)=tin(nbound)                                         

	         if (TURBMOD.GT.0) then
	           tein=teref*(f(isu+ii)*f(isu+ii)+f(isv+ii)*f(isv+ii))
	           edin=cmu*densit*tein*tein/visin
	           f(iste+ii)=tein
	           f(ised+ii)=edin
	           f(isom+ii)=edin/(tein+small)/CMU
                   f(isv2+ii)=(2./3.)*tein
                   f(isf+ii)=0.
                   vistemp=cmu*densit*tein*tein/(edin+small)
                   f(isnut+ii)=vistemp
C                  f(isvis+ii)=vistemp+viscos
	         end if  ! END OF TURBMOD CHECK

                 if (cbnd(nbound,2)(1:2).eq. 'WA' )  then
                    F(ISU+II)=uin(nbound)
                    F(ISV+II)=vin(nbound)
		 end if

	       end do
	    end if ! END OF IF LOCK FOR WEST BOUNDARY

	    if (cbnd(nbound,1)(1:2).eq. 'EA') then
	       i=nbnd(nbound,1)+1
	       js=nbnd(nbound,3)
	       je=nbnd(nbound,4)
               DO j=js,je           
                 II=I+(j-1)*NI
                 F(ISU+II)=Uin(nbound)
                 F(ISV+II)=Vin(nbound)
		 f(isvis+ii)=viscos
                 f(iss+ii)=tin(nbound)                                         

	         if (TURBMOD.GT.0) then
	           tein=teref*(f(isu+ii)*f(isu+ii)+f(isv+ii)*f(isv+ii))
	           edin=cmu*densit*tein*tein/visin
	           f(iste+ii)=tein
	           f(ised+ii)=edin
	           f(isom+ii)=edin/(tein+small)/CMU
                   f(isv2+ii)=(2./3.)*tein
                   f(isf+ii)=0.
                   vistemp=cmu*densit*tein*tein/(edin+small)
                   f(isnut+ii)=vistemp
c	           f(isvis+ii)=vistemp+viscos
	         end if  ! END OF TURBMOD CHECK

                 if (cbnd(nbound,2)(1:2).eq. 'WA' )  then
                    F(ISU+II)=uin(nbound)
                    F(ISV+II)=vin(nbound)
		 end if

	       end do
	     end if ! END OF IF BLOCK FOR EAST BOUNDARY

	       if (cbnd(nbound,1)(1:2).eq. 'SO') then
	          js=1
	          je=nbnd(nbound,3)-1
	          is=nbnd(nbound,1)
	          ie=nbnd(nbound,2)
		  do j=js,je
                  DO i=is,ie           
                    II=I+(j-1)*NI
                 F(ISU+II)=Uin(nbound)
                 F(ISV+II)=Vin(nbound)
c		 f(isvis+ii)=viscos
                 f(iss+ii)=tin(nbound)                                         

	            if (TURBMOD.GT.0) then
	              tein=teref*(f(isu+ii)*f(isu+ii)+f(isv+ii)*f(isv+ii))
	              edin=cmu*densit*tein*tein/visin
	              f(iste+ii)=tein
	              f(ised+ii)=edin
	              f(isom+ii)=edin/(tein+small)/CMU
                      f(isv2+ii)=(2./3.)*tein
                      f(isf+ii)=0.
                      vistemp=cmu*densit*tein*tein/(edin+small)
                      f(isnut+ii)=vistemp
	              f(isvis+ii)=vistemp+viscos
	            end if  ! END OF TURBMOD CHECK

                 if (cbnd(nbound,2)(1:2).eq. 'WA' )  then
                    F(ISU+II)=uin(nbound)
                    F(ISV+II)=vin(nbound)
		 end if

	          end do ! end of i loop
	          end do ! end of j loop
	       end if ! END OF IF BLOCK FOR SOUTH  BOUNDARY

	       if (cbnd(nbound,1)(1:2).eq. 'NO') then
	           js=nbnd(nbound,3)+1
	           je=nj
                   is=nbnd(nbound,1)
	           ie=nbnd(nbound,2)
		  do j=js,je
                   DO i=is,ie           
                     II=I+(j-1)*NI
                 F(ISU+II)=Uin(nbound)
                 F(ISV+II)=Vin(nbound)
		 f(isvis+ii)=viscos
                 f(iss+ii)=tin(nbound)                                         

	             if (TURBMOD.GT.0) then
	               tein=teref*(f(isu+ii)*f(isu+ii)+f(isv+ii)*f(isv+ii))
	               edin=cmu*densit*tein*tein/visin
	               f(iste+ii)=tein
	               f(ised+ii)=edin
	               f(isom+ii)=edin/(tein+small)/CMU
                       f(isv2+ii)=(2./3.)*tein
                       f(isf+ii)=0.
                       vistemp=cmu*densit*tein*tein/(edin+small)
                       f(isnut+ii)=vistemp
c	               f(isvis+ii)=vistemp+viscos
	             end if  ! END OF TURBMOD CHECK

                 if (cbnd(nbound,2)(1:2).eq. 'WA' )  then
                    F(ISU+II)=uin(nbound)
                    F(ISV+II)=vin(nbound)
		 end if

	          end do ! end of i loop
	          end do ! end of j loop
	       end if ! END OF IF BLOCK FOR NORTH BOUNDARY

	    if (cbnd(nbound,1)(1:2).eq. 'BL') then
	       i=nbnd(nbound,1)+1
	       is=nbnd(nbound,1)
	       ie=nbnd(nbound,2)
	       js=nbnd(nbound,3)
	       je=nbnd(nbound,4)

               do j=js,je           
	        do i=is,ie
                 II=I+(j-1)*NI
                 F(ISU+II)=Uin(nbound)
                 F(ISV+II)=Vin(nbound)
		 f(isvis+ii)=viscos
                 f(iss+ii)=tin(nbound)                                         

	         if (TURBMOD.GT.0) then
	           tein=teref*(f(isu+ii)*f(isu+ii)+f(isv+ii)*f(isv+ii))
	           edin=cmu*densit*tein*tein/visin
	           f(iste+ii)=tein
	           f(ised+ii)=edin
	           f(isom+ii)=edin/(tein+small)/CMU
                   f(isv2+ii)=(2./3.)*tein
                   f(isf+ii)=0.
                   vistemp=cmu*densit*tein*tein/(edin+small)
                   f(isnut+ii)=vistemp
	           f(isvis+ii)=vistemp+viscos
	         end if  ! END OF TURBMOD CHECK

	         end do
	       end do
	    end if  ! END OF IF BLOCK FOR BLANKED REGION

	end do ! END OF NBOUND LOOP

        CALL BFLUX ! TO UPDATE THE BOUNDARY MASS FLUX

	ELSE ! IF START FIELD IS TRUE

C  
C  RE-INITIALIZE ONLY FOR INFLOW,WALL & BLOCK BOUNDARY CONDITIONS
C  FOR THE STARTING SOLUTION STARTING FROM A GIVEN FIELD VALUES 
C
	do nbound=1,nrbnda
          if (cbnd(nbound,2)(1:2).eq. 'WA' .or.
     >        cbnd(nbound,2)(1:2) .eq. 'IN' ) THEN 
	       if (alpha .gt. 0) then
                   utemp=uin(nbound)
                   uin(nbound)=utemp*cos(alpha)
                   vin(nbound)=utemp*sin(alpha)
	       end if

	  if (cbnd(nbound,1)(1:2).eq. 'WE') then
	       i=nbnd(nbound,1)-1
	       js=nbnd(nbound,3)
	       je=nbnd(nbound,4)
               DO j=js,je           
                 II=I+(j-1)*NI
                 F(ISU+II)=Uin(nbound)
                 F(ISV+II)=Vin(nbound)
		 f(isvis+ii)=viscos
                 f(iss+ii)=tin(nbound)                                         

	         if (TURBMOD.GT.0) then
	           tein=teref*(f(isu+ii)*f(isu+ii)+f(isv+ii)*f(isv+ii))
	           edin=cmu*densit*tein*tein/visin
	           f(iste+ii)=tein
	           f(ised+ii)=edin
	           f(isom+ii)=edin/(tein+small)/CMU
                   f(isv2+ii)=(2./3.)*tein
                   f(isf+ii)=0.
                   vistemp=cmu*densit*tein*tein/(edin+small)
                   f(isnut+ii)=vistemp
	           f(isvis+ii)=vistemp+viscos
	         end if  ! END OF TURBMOD CHECK

                 if (cbnd(nbound,2)(1:2).eq. 'WA' )  then
                    F(ISU+II)=uin(nbound)
                    F(ISV+II)=vin(nbound)
		 end if

	       end do
	    end if ! END OF IF LOCK FOR WEST BOUNDARY

	    if (cbnd(nbound,1)(1:2).eq. 'EA') then
	       i=nbnd(nbound,1)+1
	       js=nbnd(nbound,3)
	       je=nbnd(nbound,4)
               DO j=js,je           
                 II=I+(j-1)*NI
                 F(ISU+II)=Uin(nbound)
                 F(ISV+II)=Vin(nbound)
		 f(isvis+ii)=viscos
                 f(iss+ii)=tin(nbound)                                         

	         if (TURBMOD.GT.0) then
	           tein=teref*(f(isu+ii)*f(isu+ii)+f(isv+ii)*f(isv+ii))
	           edin=cmu*densit*tein*tein/visin
	           f(iste+ii)=tein
	           f(ised+ii)=edin
	           f(isom+ii)=edin/(tein+small)/CMU
                   f(isv2+ii)=(2./3.)*tein
                   f(isf+ii)=0.
                   vistemp=cmu*densit*tein*tein/(edin+small)
                   f(isnut+ii)=vistemp
	           f(isvis+ii)=vistemp+viscos
	         end if  ! END OF TURBMOD CHECK

                 if (cbnd(nbound,2)(1:2).eq. 'WA' )  then
                    F(ISU+II)=uin(nbound)
                    F(ISV+II)=vin(nbound)
		 end if

	       end do
	     end if ! END OF IF BLOCK FOR EAST BOUNDARY

	       if (cbnd(nbound,1)(1:2).eq. 'SO') then
	          js=1
	          je=nbnd(nbound,3)-1
	          is=nbnd(nbound,1)
	          ie=nbnd(nbound,2)

                  DO j=js,je           
                  DO i=is,ie           
                    II=I+(j-1)*NI
                 F(ISU+II)=Uin(nbound)
                 F(ISV+II)=Vin(nbound)
		 f(isvis+ii)=viscos
                 f(iss+ii)=tin(nbound)                                         

	            if (TURBMOD.GT.0) then
	              tein=teref*(f(isu+ii)*f(isu+ii)+f(isv+ii)*f(isv+ii))
	              edin=cmu*densit*tein*tein/visin
	              f(iste+ii)=tein
	              f(ised+ii)=edin
	              f(isom+ii)=edin/(tein+small)/CMU
                      f(isv2+ii)=(2./3.)*tein
                      f(isf+ii)=0.
                      vistemp=cmu*densit*tein*tein/(edin+small)
                      f(isnut+ii)=vistemp
	              f(isvis+ii)=vistemp+viscos
	            end if  ! END OF TURBMOD CHECK

                 if (cbnd(nbound,2)(1:2).eq. 'WA' )  then
                    F(ISU+II)=uin(nbound)
                    F(ISV+II)=vin(nbound)
		 end if

	          end do
	          end do
	       end if ! END OF IF BLOCK FOR SOUTH  BOUNDARY

	       if (cbnd(nbound,1)(1:2).eq. 'NO') then
	           js=nbnd(nbound,3)+1
	           je=nj
                   is=nbnd(nbound,1)
	           ie=nbnd(nbound,2)
                   DO j=js,je           
                   DO i=is,ie           
                     II=I+(j-1)*NI
	             cn=-f(iscs+ii)
	 	 if (cn .ge. 0) then !inflow check 
                 F(ISU+II)=Uin(nbound)
                 F(ISV+II)=Vin(nbound)
		 f(isvis+ii)=viscos
                 f(iss+ii)=tin(nbound)                                         

	             if (TURBMOD.GT.0) then
	               tein=teref*(f(isu+ii)*f(isu+ii)+f(isv+ii)*f(isv+ii))
	               edin=cmu*densit*tein*tein/visin
	               f(iste+ii)=tein
	               f(ised+ii)=edin
	               f(isom+ii)=edin/(tein+small)/CMU
                       f(isv2+ii)=(2./3.)*tein
                       f(isf+ii)=0.
                       vistemp=cmu*densit*tein*tein/(edin+small)
                       f(isnut+ii)=vistemp
	               f(isvis+ii)=vistemp+viscos
	             end if  ! END OF TURBMOD CHECK
		  end if ! inflow check

                 if (cbnd(nbound,2)(1:2).eq. 'WA' )  then
                    F(ISU+II)=uin(nbound)
                    F(ISV+II)=vin(nbound)
		 end if

	           end do
	           end do
	       end if ! END OF IF BLOCK FOR NORTH BOUNDARY
	       end if ! END OF WALL OR INFLOW BOUNDARY CHECK

	    if (cbnd(nbound,1)(1:2).eq. 'BL') then
	       i=nbnd(nbound,1)+1
	       is=nbnd(nbound,1)
	       ie=nbnd(nbound,2)
	       js=nbnd(nbound,3)
	       je=nbnd(nbound,4)

               do j=js,je           
	        do i=is,ie
                 II=I+(j-1)*NI
                 F(ISU+II)=Uin(nbound)
                 F(ISV+II)=Vin(nbound)
		 f(isvis+ii)=viscos
                 f(iss+ii)=tin(nbound)                                         

	         if (TURBMOD.GT.0) then
	           tein=teref*(f(isu+ii)*f(isu+ii)+f(isv+ii)*f(isv+ii))
	           edin=cmu*densit*tein*tein/visin
	           f(iste+ii)=tein
	           f(ised+ii)=edin
	           f(isom+ii)=edin/(tein+small)/CMU
                   f(isv2+ii)=(2./3.)*tein
                   f(isf+ii)=0.
                   vistemp=cmu*densit*tein*tein/(edin+small)
                   f(isnut+ii)=vistemp
	           f(isvis+ii)=vistemp+viscos
	         end if  ! END OF TURBMOD CHECK

	         end do
	       end do
	    end if  ! END OF IF BLOCK FOR BLANKED REGION

	end do ! END OF NBOUND LOOP

        CALL BFLUX ! TO UPDATE THE BOUNDARY MASS FLUX

        END IF ! END OF NOT_START_FIELD CHECK
c
C  SET USER SPECIFIED INITAL FLOW FIELD  
C
        IF(INIT_FLD) CALL SET_INIT_FLD

      RETURN                                                            
      END                                                                  
C
C*************************************************************************
C
C     THE FOLLOWING 8 SUBROUTINES ARE FOR SPECIAL          
C     BOUNDARY TREATMENTS TO BE SPECIFIED BY THE USER
C
C     SPEBC1--P', 2--U, 3--V, 4--WR ,5--TE, 6--ED, 7--S, 8--VIS         
C
C***********************************************************************
 
      SUBROUTINE SPEBC1
C                                                      
      INCLUDE 'com2d'
C              
      RETURN           
      END
                                                 
      SUBROUTINE SPEBC2                                                 
C
      INCLUDE 'com2d'                     
C
      RETURN
      END
                                                 
      SUBROUTINE SPEBC3                                                 
C
      INCLUDE 'com2d'                     
C
      RETURN                                                            
      END
                                                 
      SUBROUTINE SPEBC4                                                 
C
      INCLUDE 'com2d'                     
C
      RETURN                                                            
      END
                                                 
      SUBROUTINE SPEBC5                                                 
C
      INCLUDE 'com2d'                     
C
      RETURN                                                            
      END
                                                 
      SUBROUTINE SPEBC6                                                 
C
      INCLUDE 'com2d'                     
C
      RETURN                                                            
      END 
                                                
      SUBROUTINE SPEBC7                                                
C
      INCLUDE 'com2d'                     
C
      RETURN                                                            
      END
                                                 
      SUBROUTINE SPEBC8                                                 
C
      INCLUDE 'com2d'                     
C
      RETURN                                                            
      END
C                                                 
      SUBROUTINE SPEBC9
C                                                      
      INCLUDE 'com2d'
C
      RETURN           
      END
C                                                 
      SUBROUTINE SPEBC10
C                                                      
      INCLUDE 'com2d'
C              
      RETURN           
      END
C                                                 
      SUBROUTINE SPEBC11
C
C
      INCLUDE 'com2d'
C              
C     SPECIAL SOURCES TO OBTAIN ADIABATIC BOUNDARIES ON TOP & BOTTOM
C
       J = 2    ! FOR SOUTH BOUNDARY
       DO I=2,NI-1
       II = I+(J-1)*NI
       F(ISAS+II)=0.
       F(ISS+II-NI)=F(ISS+II)
C      WRITE(66,*) I, F(ISAS+II),F(ISS+II),F(ISS+II-NI)
       END DO
       J = NJ-1 ! FOR NORTH BOUNDARY
       DO I=2,NI-1
       II = I+(J-1)*NI
       F(ISAN+II)=0.
       F(ISS+II+NI)=F(ISS+II)
C      WRITE(66,*) I, F(ISAN+II),F(ISS+II),F(ISS+II+NI)
       END DO

      RETURN           
      END
C
C**********************************************************************
C                                                                       
      SUBROUTINE FINAL                                                  
C
C     PRINT FIELD VALUES FOR THE FINAL RESULTS                         
C
C***********************************************************************
C
      INCLUDE 'com2d'                     

      CALL PRINT(F(ISX  +1),' X-CENTRE ') 
      CALL PRINT(F(ISY  +1),' Y-CENTRE ') 
      CALL PRINT(F(ISU  +1),'U-VELOCITY') 
      CALL PRINT(F(ISV  +1),'V-VELOCITY')
      CALL PRINT(F(ISP  +1),' PRESSURE ')                               
      CALL PRINT(F(ISS  +1),' TEMPER')                               

	IF (TURBMOD .NE. 0) THEN
      CALL PRINT(F(ISTE  +1),' TURB TE ')                               
      CALL PRINT(F(ISED  +1),' TURB ED ')                               
      CALL PRINT(F(ISOM  +1),' OMEGA ')                               
      CALL PRINT(F(ISV2  +1),' V-SQ ')                               
      CALL PRINT(F(ISF  +1),' F ')                               
      CALL PRINT(F(ISVIS  +1),' VISCOSITY ')                               

        END IF
C             

      RETURN                                                            
      END  
C                                         
C*********************************************************************
C
		SUBROUTINE TURBCONS 
C
C********************************************************************
C
		INCLUDE 'com2d'
C          
C  CLOSURE COEFFICIENTS FOR K-EPSILON MODEL WITH WALL FUNCTION
C
             IF(TURBMOD.EQ.1) THEN
C
                C1 =1.44
                C2 =1.92
                CMU =0.09
                CDQR =0.548
                CDRT =0.3
                CDTQ =0.1643
                ELOG =9.793
                CAPPA =0.4187
                PR(5)=1.
                PR(6)=1.32

             END IF
C          
C  CLOSURE COEFFICIENTS FOR K-EPSILON MODEL WITH RODI'S TWO-LAYER MODIFICATION
C
             IF(TURBMOD.EQ.2) THEN
C
                C1 =1.44
                C2 =1.92
                CMU =0.09
                CDQR =0.548
                CDRT =0.3
                CDTQ =0.1643
                ELOG =9.793
                CAPPA =0.4187
                CL1=CAPPA/CDTQ
                CEPS=13.2
                CDAMP=0.495
                ADAMP=26
                PR(5)=1.
                PR(6)=1.32

             END IF
C          
C  CLOSURE COEFFICIENTS FOR K-EPSILON MODEL WITH CHIEN'S LOW-RE MODIFICATION
C
             IF(TURBMOD.EQ.3) THEN
C
                C1 =1.35
                C2 =1.8
                CMU =0.09
                CDQR =0.548
                CDRT =0.3
                CDTQ =0.1643
                ELOG =9.793
                CAPPA =0.4187
                PR(5)=1.
                PR(6)=1.3
                PR(11)=0.9

              END IF

C          
C  CLOSURE COEFFICIENTS FOR K-OMEGA MODEL
C
             IF(TURBMOD.EQ.4) THEN

                GAMMA=5./9.
                BETA=3./40.
                BETA_STR=0.09
                CMU  =0.09
                CDQR =0.548
                CDRT =0.3
                ELOG =9.793
                CAPPA =0.4187
                PR(5)=2.
                PR(7)=2.

             END IF
C          
C  CLOSURE COEFFICIENTS FOR V2F MODEL
C
             IF(TURBMOD.EQ.5) THEN

                ENN=6.
C
C CONSTANT USED FOR PIPE/CHANNEL FLOW WHICH GIVES THE RIGHT CF
C

C               CL=0.17
C               CETA=70
C               C1=1.0 
C               C2=1.92
C               C1F=1.4
C               C2F=0.3
C               CMU=0.22
C               PR(5)=1.
C               PR(6)=1.3
C               PR(8)=1.
C               PR(9)=1.

C
C
C LIEN AND DURBIN (CTR-1996) CONSTANTS
C
                CL=0.17
                CETA=70
                C1=1.0 
                C2=1.9
                C1F=1.4
                C2F=0.3
                CMU=0.19
                PR(5)=1.
                PR(6)=1.5
                PR(8)=1.
                PR(9)=1.
C
C
C FLUENT CONSTANTS
C
C                CL=0.23
c                CETA=70
c                C1=1.0 
c                C2=1.9
c                C1F=1.4
c                C2F=0.3
c                CMU=0.22
c                PR(5)=1.
c                PR(6)=1.3
c                PR(8)=1.
c                PR(9)=1.
C
C               CL=0.25
C               CETA=85
C               C1=1.0 
C               C2=1.9
C               C1F=1.4
C               C2F=0.3
C               CMU=0.22
C               PR(5)=1.
C               PR(6)=1.3
C               PR(8)=1.
C               PR(9)=1.

             END IF     ! END OF IF (TURBMOD.EQ.5 ) THEN LOOP
C
C  CLOSURE COEFFICIENTS FOR SPALART-ALLMARAS MODEL
C
             IF(TURBMOD.EQ.6) THEN

              CMU =0.09
              CAPPA=0.41
              CB1=0.1355
              CB2=0.622
              CV1=7.1
              PR(10)=2./3.
              CW1=CB1/CAPPA/CAPPA+(1.+CB2)/PR(10)
              CW2=0.3
              CW3=2.

             END IF
C
C  CLOSURE COEFFICIENTS FOR SST MODEL
C
             IF(TURBMOD.EQ.7) THEN
C
             BETA1=0.075
             BETA2=0.0828
             BETA_STR=0.09
             GAMA1=0.555
             GAMA2=0.44
             A1=0.31
C
             CMU=0.09
             CDRT =0.3
             CDQR =0.548
C            CAPPA =0.4187
C            ELOG =9.793
             CAPPA =0.41
             ELOG =9.
	
C
             PR(5)=2.
             PR(7)=2.
C            PR_KSST=1.
C            PR_EPSST=1.3
C            SIGMA_OMG2=1/1.3
             PR1_TE=2.
             PR1_OM=2.
             PR2_TE=1.
             PR2_OM=1.3

C
             END IF
C
C     INITIALISE THE DAMPING FUNCTION ARRAYS
C
      DO  II=1,NIJ
		F1(II)=1.
		F2(II)=1.
		FMU(II)=1.
		D(II)=0.
		E(II)=0.
		ALMU(II)=1.
		ALEPS(II)=1.
      	END DO
   
		return
		end 
C
C**********************************************************************
C
	SUBROUTINE STORVAR
C
        INCLUDE 'com2d'
C
C**********************************************************************
C
	IF (TEMPSCH .EQ. 1 .OR. NTSTEP .EQ.1) THEN
C
C FOR FIRST ORDER DISCRETIZATION SCHEME(TEMPSCH=1 & NLEVEL=1) 
C ONLY PREVIOUS(N-1) TIME STEP VALUE IS STORED AT EACH TIME STEP(N)
C
C AND FOR SECOND ORDER (TEMPSCH=2 & NLEVEL=2) AT FIRST TIME STEP(N=1)
C THEN (N-1) AND (N-2) HAVE THE INITIAL VALUES 
C
	   DO IL=1,NLEVEL
	    DO J=1,NJ
	     DO I=1,NI
	       II=I+(J-1)*NI
	       U(ii,iblk,il)=F(ISU+II)
	       V(ii,iblk,il)=F(ISV+II)
	       PRESS(ii,iblk,il)=F(ISP+II)
	       SCAL(ii,iblk,il)=F(ISS+II)
	       te(ii,iblk,il)=f(iste+ii)
	       ed(ii,iblk,il)=f(ised+ii)
	       aom(ii,iblk,il)=f(isom+ii)
	       vsq(ii,iblk,il)=f(isv2+ii)
	       anut(ii,iblk,il)=f(isnut+ii)
	      END DO
	    END DO
	   END DO

C
C STORING THE N-2 TIME STEP FIELD VALUES (BASICALLY NEEDED FOR 
C CONTINUATION RUN ONLY) 
C
	 if (tempsch .eq. 1)  then
C
C STORE N-1 TIME STEP VALUES -- FOR 1st ORDER TEMPORAL DISCRETISATION
C -- REQUIRED ONLY FOR STOPED JOBS
C
         write(28,rec=iblk)(u(ii,iblk,1),v(ii,iblk,1),press(ii,iblk,1),
     >                     scal(ii,iblk,1),te(ii,iblk,1),
     >    ed(ii,iblk,1),aom(ii,iblk,1),vsq(ii,iblk,1),anut(ii,iblk,1),ii=1,nij)

	 else
C
C FOR 2nd ORDER TEMPORAL DISCRETISATIO
C STORE N-1 TIME STEP VALUES
C
         write(28,rec=iblk)(u(ii,iblk,2),v(ii,iblk,2),press(ii,iblk,2),
     >                     scal(ii,iblk,2),te(ii,iblk,2),
     >     ed(ii,iblk,2),aom(ii,iblk,2),vsq(ii,iblk,2),anut(ii,iblk,2),ii=1,nij)
C
C STORE N-2 TIME STEP VALUES -- REQUIRED ONLY FOR STOPED JOBS
C
         write(28,rec=nblock+iblk)(u(ii,iblk,1),v(ii,iblk,1),press(ii,iblk,1),
     >                     scal(ii,iblk,1),te(ii,iblk,1),
     >     ed(ii,iblk,1),aom(ii,iblk,1),vsq(ii,iblk,1),anut(ii,iblk,1),ii=1,nij)
	end if
	ELSE
C
C TO STORE THE (N-2) & (N-1) AT EVERY TIME STEP (N>1) FOR SECOND 
C ORDER TEMPORAL DISCRETISATION
C
C
C READING THE N-2 TIME STEP FIELD VALUES FOR THE FILE (BASICALLY NEEDED FOR 
C CONTINUATION RUN ONLY) 
C
         read(28,rec=iblk)(u(ii,iblk,2),v(ii,iblk,2),press(ii,iblk,2),
     >                     scal(ii,iblk,2),te(ii,iblk,2),
     >     ed(ii,iblk,2),aom(ii,iblk,2),vsq(ii,iblk,2),anut(ii,iblk,2),ii=1,nij)
	  DO J=1,NJ
	   DO I=1,NI
	      II=I+(J-1)*NI
	      U(ii,iblk,1)=U(ii,iblk,2)
	      U(ii,iblk,2)=F(ISU+II)
	      V(ii,iblk,1)=V(ii,iblk,2)
	      V(ii,iblk,2)=F(ISV+II)
	      PRESS(ii,iblk,1)=PRESS(ii,iblk,2)
	      PRESS(ii,iblk,2)=F(ISP+II)
	      SCAL(ii,iblk,1)=SCAL(ii,iblk,2)
	      SCAL(ii,iblk,2)=F(ISS+II)
	      TE(ii,iblk,1)=TE(ii,iblk,2)
	      te(ii,iblk,2)=f(iste+ii)
	      ed(ii,iblk,1)=ed(ii,iblk,2)
	      ed(ii,iblk,2)=f(ised+ii)
	      aom(ii,iblk,1)=aom(ii,iblk,2)
	      aom(ii,iblk,2)=f(isom+ii)
	      vsq(ii,iblk,1)=vsq(ii,iblk,2)
	      vsq(ii,iblk,2)=f(isv2+ii)
	      anut(ii,iblk,1)=anut(ii,iblk,2)
	      anut(ii,iblk,2)=f(isnut+ii)
	    END DO
	  END DO

C
C STORING THE N-2 TIME STEP FIELD VALUES (BASICALLY NEEDED FOR 
C CONTINUATION RUN ONLY) 
C
         write(28,rec=iblk)(u(ii,iblk,2),v(ii,iblk,2),press(ii,iblk,2),
     >                     scal(ii,iblk,2),te(ii,iblk,2),
     >   ed(ii,iblk,2),aom(ii,iblk,2),vsq(ii,iblk,2),anut(ii,iblk,2),ii=1,nij)
C
C STORE N-2 TIME STEP VALUES -- REQUIRED ONLY FOR STOPED JOBS
C
         write(28,rec=nblock+iblk)(u(ii,iblk,1),v(ii,iblk,1),press(ii,iblk,1),
     >                     scal(ii,iblk,1),te(ii,iblk,1),
     >    ed(ii,iblk,1),aom(ii,iblk,1),vsq(ii,iblk,1),anut(ii,iblk,1),ii=1,nij)
	END IF

	return
	end 
C
C**********************************************************************
C
	SUBROUTINE STORVAR_CWCS
C
        INCLUDE 'com2d'
C
C**********************************************************************
C
C
C
	    DO J=1,NJ
	     DO I=1,NI
	       II=I+(J-1)*NI
	       CW_OLD(ii,iblk)=F(ISCW+II)
	       CS_OLD(ii,iblk)=F(ISCS+II)
	      END DO
	    END DO

	return
	end 
C
C**********************************************************************
C
	SUBROUTINE DELPHI_MAX
C
C  COMPUTES THE MAXIMUM CHANGE OF EACH VARIABLE DURING ONE ITERATION
C 
C**********************************************************************
C
         INCLUDE 'com2d'
C
	uerr(iblk)=0.
	verr(iblk)=0.
	perr(iblk)=0.
	teerr(iblk)=0.
	ederr(iblk)=0.
	aomerr(iblk)=0.
	vsqerr(iblk)=0.
	ffnerr(iblk)=0.
	anuterr(iblk)=0.
	scalerr(iblk)=0.
C
	il=nlevel

c	DO J=1,NJ
c	DO I=1,NI
      DO J=JSOLV_B,JSOLV_E
      DO I=ISOLV_B,ISOLV_E

	II=I+(J-1)*NI

	uerr(iblk)=max(abs(f(isu+ii)-u(ii,iblk,il)),uerr(iblk))
	verr(iblk)=max(abs(f(isv+ii)-v(ii,iblk,il)),verr(iblk))
	perr(iblk)=max(abs(f(isp+ii)-press(ii,iblk,il)),perr(iblk))
	teerr(iblk)=max(abs(f(iste+ii)-te(ii,iblk,il)),teerr(iblk))
	ederr(iblk)=max(abs(f(ised+ii)-ed(ii,iblk,il)),ederr(iblk))
	aomerr(iblk)=max(abs(f(isom+ii)-aom(ii,iblk,il)),aomerr(iblk))
	vsqerr(iblk)=max(abs(f(isv2+ii)-vsq(ii,iblk,il)),vsqerr(iblk))
	ffnerr(iblk)=max(abs(f(isf+ii)-ffn(ii,iblk,il)),ffnerr(iblk))
	anuterr(iblk)=max(abs(f(isnut+ii)-anut(ii,iblk,il)),anuterr(iblk))
	scalerr(iblk)=max(abs(f(iss+ii)-scal(ii,iblk,il)),scalerr(iblk))

	END DO
	end DO

	return
	end 
C
C**********************************************************************
C
	SUBROUTINE MONI_LOC 
C
C**********************************************************************
C
	include 'com2d'
C
	iimon=imon+(jmon-1)*ni
	umon=f(isu+iimon)		
	vmon=f(isv+iimon)		
	pmon=f(isp+iimon)		
	phimon=f(iss+iimon)		

	return
	end   
C
C**********************************************************************
C
	SUBROUTINE PHAVG
C
C**********************************************************************
C
        INCLUDE 'com2d'
C
	IF (NTSTEP .LT. NTAVG+1) THEN
C
C PHSUM IS DONE ONLY AFTER NTAVG TIME STEPS HAVE ELAPSED IN THE TIME LOOP
C
	    NTCONT=0
	    DO J=1,NJ
	    DO I=1,NI
	      II=I+(J-1)*NI
	        pavg(ii)=0.
	        uavg(ii)=0.
	        vavg(ii)=0.
	        wavg(ii)=0.
	       edavg(ii)=0.
	       teavg(ii)=0.
	      aomavg(ii)=0.
	      vsqavg(ii)=0.
	      ffnavg(ii)=0.
	     anutavg(ii)=0.
	     scalavg(ii)=0.
	    END DO
	    END DO
	ELSE
	   IF (NTSTEP .GT.NTAVG+1) THEN
C
C  READING THE PHSUM OF PREVIOUS TIME STEP FOR EACH BLOCK
C 	    
            READ(11,rec=iblk) NTPREV,NTCONT,ni,nj,
     >                        (pavg(ii),uavg(ii),vavg(ii),
     >                         wavg(ii),edavg(ii),teavg(ii),
     >                         aomavg(ii),vsqavg(ii),ffnavg(ii),
     >                         anutavg(ii),scalavg(ii),ii=1,ni*nj)

	    END IF
C
C PH VALUE OF THE PREVIOUS ITERATION IS ADDED TO PHSUM
C
	    NTCONT=NTCONT+1
	    DO J=1,NJ
	    DO I=1,NI
   	       II=I+(J-1)*NI
	       pavg(ii)=pavg(ii)+F(ISP+II)
	       uavg(ii)=uavg(ii)+F(ISU+II)
	       vavg(ii)=vavg(ii)+F(ISV+II)
	       wavg(ii)=wavg(ii)+F(ISW+II)
	       edavg(ii)=edavg(ii)+F(ISED+II)
	       teavg(ii)=teavg(ii)+F(ISTE+II)
	       aomavg(ii)=aomavg(ii)+F(ISOM+II)
	       vsqavg(ii)=vsqavg(ii)+F(ISV2+II)
	       ffnavg(ii)=ffnavg(ii)+F(ISF+II)
	       anutavg(ii)=anutavg(ii)+F(ISNUT+II)
	       scalavg(ii)=scalavg(ii)+F(ISS+II)
	     END DO
	     END DO
                WRITE(11,rec=iblk) NTSTEP-1,NTCONT,ni,nj,
     >                        (pavg(ii),uavg(ii),vavg(ii),
     >                         wavg(ii),edavg(ii),teavg(ii),
     >                         aomavg(ii),vsqavg(ii),ffnavg(ii),
     >                         anutavg(ii),scalavg(ii),ii=1,ni*nj)
C
	END IF

	return
	end 
C
C**********************************************************************
C
        SUBROUTINE GET_STORVAR
        INCLUDE 'com2d'

          
        IF (TEMPSCH .EQ. 1) THEN
         read(28,rec=iblk)(u(ii,iblk,1),v(ii,iblk,1),press(ii,iblk,1),
     >                     scal(ii,iblk,1),te(ii,iblk,1),
     >                     ed(ii,iblk,1),aom(ii,iblk,1),
     >                     vsq(ii,iblk,1),anut(ii,iblk,1),ii=1,nij)
          else

         read(28,rec=iblk)(u(ii,iblk,2),v(ii,iblk,2),press(ii,iblk,2),
     >                     scal(ii,iblk,2),te(ii,iblk,2),
     >                     ed(ii,iblk,2),aom(ii,iblk,2),
     >                     vsq(ii,iblk,2),anut(ii,iblk,2),ii=1,nij)
C
C STORE N-2 TIME STEP VALUES -- REQUIRED ONLY FOR STOPED JOBS
C
         read(28,rec=nblock+iblk)(u(ii,iblk,1),v(ii,iblk,1),press(ii,iblk,1),
     >                     scal(ii,iblk,1),te(ii,iblk,1),
     >                     ed(ii,iblk,1),aom(ii,iblk,1),
     >                     vsq(ii,iblk,1),anut(ii,iblk,1),ii=1,nij)
        end if

        return
        end

C                                         
C*********************************************************************
C
		SUBROUTINE SOLVE_INDEX 
C
C  SOLVING INDICES FOR DIFFERENT TURBULENCE MODELS
C
C********************************************************************
C
		INCLUDE 'com2d'
C
C  LAMINAR FLOW
C
              IF (TURBMOD.EQ.0) THEN
                   SOLVE(5)=.FALSE.
                   SOLVE(6)=.FALSE.
                   SOLVE(7)=.FALSE.
                   SOLVE(8)=.FALSE.
                   SOLVE(9)=.FALSE.
              END IF
C
C  STD K-EPSILON MODEL 
C
              IF (TURBMOD.EQ.1) THEN
                   SOLVE(5)=.TRUE.
                   SOLVE(6)=.TRUE.
                   SOLVE(7)=.FALSE.
                   SOLVE(8)=.FALSE.
                   SOLVE(5)=.FALSE.
                   SOLVE(6)=.FALSE.
                   SOLVE(7)=.FALSE.
                   SOLVE(8)=.FALSE.
                   SOLVE(9)=.FALSE.
              END IF
C
C  STD K-EPSILON MODEL 
C
              IF (TURBMOD.EQ.1) THEN
                   SOLVE(5)=.TRUE.
                   SOLVE(6)=.TRUE.
                   SOLVE(7)=.FALSE.
                   SOLVE(8)=.FALSE.
                   SOLVE(9)=.FALSE.
              END IF
C
C   K-EPSILON TWO LAYER MODEL 
C
              IF (TURBMOD.EQ.2) THEN
                   SOLVE(5)=.TRUE.
                   SOLVE(6)=.TRUE.
                   SOLVE(7)=.FALSE.
                   SOLVE(8)=.FALSE.
                   SOLVE(9)=.FALSE.
              END IF
C
C   K-EPSILON LOW RE CHIEN MODEL 
C
              IF (TURBMOD.EQ.3) THEN
                   SOLVE(5)=.TRUE.
                   SOLVE(6)=.TRUE.
                   SOLVE(7)=.FALSE.
                   SOLVE(8)=.FALSE.
                   SOLVE(9)=.FALSE.
              END IF
C
C   K-OMEGA MODEL 
C
              IF (TURBMOD.EQ.4) THEN
                   SOLVE(5)=.TRUE.
                   SOLVE(6)=.FALSE.
                   SOLVE(7)=.TRUE.
                   SOLVE(8)=.FALSE.
                   SOLVE(9)=.FALSE.
              END IF
C
C   K-EPSILON V2F  MODEL 
C
              IF (TURBMOD.EQ.5) THEN
                   SOLVE(5)=.TRUE.
                   SOLVE(6)=.TRUE.
                   SOLVE(7)=.FALSE.
                   SOLVE(8)=.TRUE.
                   SOLVE(9)=.TRUE.
              END IF
C
C   SPALART-ALLMARAS  MODEL
C
              IF (TURBMOD.EQ.6) THEN
                   SOLVE(5)=.FALSE.
                   SOLVE(6)=.FALSE.
                   SOLVE(7)=.FALSE.
                   SOLVE(8)=.FALSE.
                   SOLVE(9)=.FALSE.
                   SOLVE(10)=.TRUE.
              END IF
C
C   SST  MODEL  ( SHEAR STRESS TRANSPORT MODEL )
C
              IF (TURBMOD.EQ.7) THEN
                   SOLVE(5)=.TRUE.
                   SOLVE(6)=.FALSE.
                   SOLVE(7)=.TRUE.
                   SOLVE(8)=.FALSE.
                   SOLVE(9)=.FALSE.
                   SOLVE(10)=.FALSE.
              END IF

         RETURN
         END
C
C**********************************************************************

		SUBROUTINE WRITE_EXITPROF 
C
C  WRITING THE EXIT PROFILE FOR CHANNEL FLOW RESULTS
C
C********************************************************************
C
          INCLUDE 'com2d'      
C
c   WRITING THE FLOW VARIABLE PROFILES FROM THE CHANNEL FLOW RESULTS
C
          OPEN(123,FILE='INPROF.DAT')
           WRITE(123,*)NJ   
           DO J   =1,NJ  
            II=NI+(J-1)*NI
            SINL=F(ISY+II)
            UINL=F(ISU+II)
            VINL=F(ISV+II)
            TEINL=F(ISTE+II)
            EDINL=F(ISED+II)
            AOMINL=F(ISOM+II)
            V2INL=F(ISV2+II) 
            FINL=F(ISF+II)
            ANUTINL=F(ISNUT+II)  

         WRITE(123,155)sinl,uinl,vinl,teinl,edinl
     >                 ,aominl,v2inl,finl,anutinl
           END DO
155       FORMAT(9E10.3) 
          CLOSE(123)
           RETURN
           END
C
C****************************************************        
