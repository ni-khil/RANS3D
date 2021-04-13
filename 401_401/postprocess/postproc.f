C
C       ***** POSTPROC.F   *******
C
C     THIS IS A POST-PROCESSING ROUTINE FOR PLOTTING 
C         1) PRESSURE DISTRIBUTION
C         2) SKIN FRICTION ALONG THE BODY
C         3) CONTOUR /GRID/VELOCITY PLOT ON TECPLOT
C
      
      INCLUDE 'postcom'
      DIMENSION NV(NVTOT),RNORM(11),IND_FLAG(4)
      DIMENSION PH(NX,NY),VOR(NX,NY),ab(nmax)
	dimension xg(max_stat),yg(max_stat)
        dimension iblk_xgiv(max_stat),iblk_ygiv(max_stat)
       DIMENSION YM(NY),ANUSLT(NY)
       
       LOGICAL DOUBLE_PREC

C
C  FILE <contres> IS FLOW FIELD FOR STEADY FLOW
C  FILE <postres> IS FLOW FIELD AT DIFFERENT TIMESTEPS FOR UNSTEADY FLOW
C  FILE <postdat> IS INPUT DATA INFORMATION OF THIS PROGRAM 
C  FILE <cpcf.dat> GIVES PRESSURE AND SKIN FRICTION DISTRIBUTION 
C                  TO BE PLOTTED USING TECPLOT SOFTWARE
C  FILE <uvpt_fld.dat> GIVES THE PRESSURE AND VELOCITY FIELD IN TECPLOT FORMAT
C  FILE <turb_fld.dat> GIVES THE TURBULENCE SCALAR FIELD IN TECPLOT FORMAT
C  FILE <grid.dat> GIVES THE GRID DATA INFORMATION IN THE TECPLOT FORMAT
C  FILE <xprofile.dat> &FILE <YPROFILE.DAT>GIVES THE U,V,K FIELD DATA INFORMATION IN THE TECPLOT FORMAT AT DIFFERENT STATIONS ALONG X & Y
C  FILE <wall_dist.datT> GIVES THE WALL DISTANCES                             
C
              PI=4.*atan(1.) 
C
C    OPEN STATEMENTS FOR SEQUENTIAL FILES     
C

      OPEN(UNIT=9,FILE='postdat')
      OPEN(UNIT=50,FILE='xprofile.dat')
      OPEN(UNIT=51,FILE='yprofile.dat')
      OPEN(UNIT=2,FILE='uvpt_fld.dat')
      OPEN(UNIT=4,FILE='turb_fld.dat')
      OPEN(UNIT=3,FILE='grid.dat') 
      OPEN(UNIT=12,FILE='cpcf.dat')  

      READ(9,*)
      READ(9,*)DOUBLE_PREC

         IF (DOUBLE_PREC) THEN
C
C  CALCULATION OF RECORD LENGTH FOR DIFFERENT FILES
C
	    leng=(2+idim1*2)*2*4
	    lent=(15+idim1*11)*2*4
            lenc=(15+idim1*18)*2*4

C CODE EXECUTED IN DOUBLE PRECISION

	ELSE

C CODE EXECUTED IN SINGLE PRECISION
C
C  CALCULATION OF RECORD LENGTH FOR DIFFERENT FILES
C
	    leng=(2+idim1*2)*4
	    lent=(15+idim1*11)*4
            lenc=(15+idim1*18)*4

	END IF    ! END OF IF LOOP FOR DOUBLE PRECISION CHECK
C
C    OPEN STATEMENTS FOR DIRECT ACCESS FILES 
C
            OPEN(7,FILE='GRID',ACCESS='DIRECT',FORM='UNFORMATTED'
     >,recl=leng)
          OPEN(8,FILE='CONTRES',ACCESS='DIRECT',FORM='UNFORMATTED'
     >,recl=lenc)
             OPEN(20,file='POSTRES',ACCESS='DIRECT',FORM='UNFORMATTED'
     >,recl=lent)

C                       
C *************************************
C READ INFORMATION FROM <postdat> FILE
C *************************************
C 
C GIVE THE VALUE OF REYNOLDS NO., UIN, DENSITY
C
      READ(9,*)
      READ(9,*)RE,UREF,DENSIT
      VISCOS=1./RE
C
C READ THE NUMBER OF BLOCKS
C
      READ(9,*)     
      READ(9,*)nblock

C READ THE I J AND BLOCK NUMBER FOR THE PRESSURE REFERENCE POINT
C
      READ(9,*)
      READ(9,*)ipref,jpref,ibref_pr
C
C READ  THE  VALUE OF NFILE = 1 ( STEADY) >1 (UNSTEADY)          
C no_xstat = nos. of x stations where profile is to be plotted 
C XG(i), iblk_xgiv(i) = X value and Block number for the stations
C no_ystat = nos. of y stations where profile is to be plotted 
C YG(i), iblk_ygiv(i) = Y value and Block number for the stations
C
	read(9,*)
	read(9,*)NFILE
	read(9,*)
	read(9,*)no_xstat
	read(9,*)
	read(9,*)(XG(i),iblk_xgiv(i),i=1,no_xstat)
	read(9,*)
	read(9,*)no_ystat
	read(9,*)
	read(9,*)(YG(i),iblk_ygiv(i),i=1,no_ystat)
C
C READ  THE RELEVANT RECORD NO TO BE READ FROM THE FILE <POSTRES>
C NTREC_BEG : STARTING RECORD NUMBER
C NTREC_END : ENDING RECORD NUMBER
C NTREC_STP : NO OF RECODS TO SKIP 
C
      read(9,*)
      read(9,*)ntrec_beg,ntrec_end,ntrec_stp
      read(9,*)
      read(9,*)
	
        if (NFILE .EQ.1 ) then
	   ntrec_beg=1
	   ntrec_end=1
	   ntrec_end=1
	end if

	write(*,*)'no_xstat,no_ystat'
	write(*,*)no_xstat,no_ystat
c
c     write headings for file='uvpt_fld.dat'	
c
      WRITE(2,*)'TITLE = "uvpt_fld"'
      WRITE(2,*)'VARIABLES = XC, YC, U, V, P,T'
c
c     write headings for file='turb_fld.dat'	
c
      WRITE(4,*)'TITLE = "turb_fld"'
      WRITE(4,*)'VARIABLES = XC, YC, k, epsilon,  MUT'
c
c     write headings for file='grid.dat'	
c
	WRITE(3,*)'TITLE = "grid"'
      WRITE(3,*)'VARIABLES = X, Y'
c
c     write headings for file='xprofile.dat'
c
      WRITE(50,*)'TITLE = "xprofile"'
      WRITE(50,*)'VARIABLES = Y,U,V,T'
c
c     write headings for file='yprofile.dat'
c
      WRITE(51,*)'TITLE = "yprofile"'
      WRITE(51,*)'VARIABLES = X,U,V,T'
c
c
c     write headings for file='cpcf.dat' 
c
         WRITE(12,*)'TITLE = "CP"'
      WRITE(12,*)'VARIABLES = i, XC, CP, CF, YPLUS'

C
C TIME STEP LOOP STARTS (FOR TIME INSTANTS STORED IN POSTRES)
C
	do nt=ntrec_beg,ntrec_end,ntrec_stp
	   ntsp_rec=nt

	if (nt .gt. ntrec_beg) then
	     rewind(9)
	     do nl=1,14
	       read(9,*)
	     end do
	end if

c
c       BLOCK LOOP STARTS
c
      do iblk=nblock,1,-1	
c
      	READ(9,*)
      	READ(9,*)
      	READ(9,*)ni,nj
	    
      	READ(9,*)
      	READ(9,*)isolve_b,isolve_e,jsolve_b,jsolve_e
      	READ(9,*)
      	READ(9,*)
      	READ(9,*)
C
C   IFLAG = 1 FOR WALL AT NORTH  IFLAG = 2 FOR WALL AT SOUTH
C   IFLAG = 3 FOR WALL AT EAST   IFLAG = 4 FOR WALL AT WEST
C
	READ(9,*)NWALS,(ind_flag(nw),nw=1,nwals)

      	NIM=NI-1
      	NJM=NJ-1
      	NIJ=NI*NJ
C
      	NV(1)=1
      	DO IV=2,NVTOT
    	      NV(IV)=NV(IV-1)+NIJ
	    end do
C
C     CALCULATE STARTING ADDRESSES
C  

        ISU   =NV( 1)-1
        ISV   =NV( 2)-1
        ISWR  =NV( 3)-1
        ISTE  =NV( 4)-1
        ISED  =NV( 5)-1
        ISOM  =NV( 6)-1
        ISV2  =NV( 7)-1
        ISF   =NV( 8)-1
        ISNUT =NV( 9)-1
        ISS   =NV( 10)-1
        ISP   =NV( 11)-1
        ISDEN =NV( 12)-1
        ISVIS =NV( 13)-1
        ISVISW=NV( 14)-1
        ISVISS=NV( 15)-1
        ISW   =NV( 16)-1
        ISCW  =NV( 17)-1
        ISCS  =NV( 18)-1
        ISCOX =NV( 19)-1
        ISCOY =NV( 20)-1
        ISX   =NV( 21)-1
        ISY   =NV( 22)-1
  
C       
c      READ THE GRID CO-ORDINATES AND CALCULATE CELL CENTRE COORDS.
c 
	CALL GEOINP

        CALL GEOMET

	write(*,*)iblk,ni,nj

C 
C  READ THE FIELD VALUES FROM THE FILE CONTRES 
C  FOR STEADY FLOW CASE
C

       if (NFILE .EQ.1) then
        IPHLST=18*NIJ
      read(8,rec=iblk) IBSWP,NTSTEP,RNORM,delt,tinst,
     >                  (f(iphi),iphi=1,iphlst)

        ELSE
C 
C  READ THE FIELD VALUES FROM THE FILE CONTRES 
C  FOR UNSTEADY FLOW CASE
C
        ipb=1
        ipe=10*ni*nj
        irnum=iblk+nblock*(nt-1)
         read(20,rec=irnum) IBSWP,NTSTEP,RNORM,delt,tinst,
     >      (f(ip),ip=ipb,ipe)

C
        end if

c
c      COMPUTATION OF BOUNDARY INDICES
c
       CALL BINDEX
c
C
C WRITE THE FIELD AND GRID VALUES TO BE PLOTTED ON TECPLOT
C
          WRITE(2,*)'ZONE T = " BL=',iblk,' ",I= ',NI,',J = ',NJ
          WRITE(3,*)'ZONE T = " BL=',iblk,' ",I= ',NI,',J = ',NJ
          WRITE(4,*)'ZONE T = " BL=',iblk,' ",I= ',NI,',J = ',NJ

                             
	    DO J =  1,nj
            DO I = 1,ni
              II = I+(J-1)*NI
C		
              YCENT=F(ISY+II)
              XCENT=F(ISX+II)
	        XCOR=F(ISCOX+II)
	        YCOR=F(ISCOY+II)
 			U = F(ISU + II)
 			V = F(ISV + II)
			P = F(ISP + II)
                        T = F(ISS + II)
			TE = F(ISTE + II)
			ED = F(ISED + II)
		        AMUT=F(ISVIS+II)-VISCOS
			WRITE(2,112)XCENT,YCENT,U,V,P,T
			WRITE(3,*)XCOR,YCOR
			WRITE(4,112)XCENT,YCENT,TE,ED,AMUT
 
		   END DO

           END DO

112	format(6(E14.6,1x))
c****************************************************************************
C
c        TRANSVERSE PROFILES OF VARIABLE AT GIVEN X 
C
C**************************************************************************


	  do nxs=1,no_xstat

		XGIV=XG(nxs)
		if (iblk.eq.iblk_xgiv(nxs)) then


	  njtot=jsolve_e-jsolve_b+3
		
          WRITE(50,*)
     >    'ZONE T =" xg=',xgiv,', ib=',iblk, '",I=',1,',J=',NJTOT

C
C   SEARCH TO FIND OUT BETWEEN WHICH TWO I=CONST PLANES THE GIVEN PLANE LIES
C
                do i=1,ni
                  ii=i+(1-1)*ni
                  ab(i)=f(isx+ii)
                end do

          CALL SEARCH(ni,ab,XGIV,FIND,IL)
C
          WRITE(*,*)'IL=',IL

          DO j=jsolve_b-1,jsolve_e+1
       
          II = IL + (J-1)*NI

          xl=f(isx+II)
          xr=f(isx+II+1)

           y=f(isy+II)
          ul=f(isu+II)/uref
          ur=f(isu+II+1)/uref

          CALL LINT(xl,xr,ul,ur,XGIV,u)

          vl=f(isv+II)/uref
          vr=f(isv+II+1)/uref

          CALL LINT(xl,xr,vl,vr,XGIV,v)

          Tl=f(ISS+II)
          Tr=f(ISS+II+1)

          CALL LINT(xl,xr,Tl,Tr,XGIV,T)


	  WRITE(50,234)Y,U,V,T
          END DO
          end if
	  end do ! end x= constant station loop
C
C      EVALUATION OF LOCAL AND AVERAGE NUSSELT NOS.
C
C      AS PER THREE POINT FORMULA
C
              DO J=1,NJ
              II1=1+(J-1)*NI
              II2=II1+1
              II3=II2+1
              TL=F(ISS+II1)
              TR=F(ISS+II2)
              TRR=F(ISS+II3)
              XL=F(ISX+II1)
              XR=F(ISX+II2)
              XRR=F(ISX+II3)
              YM(J)=F(ISY+II1)
              DYDX1=-TL*(XRR+XR)/XR/XRR
              DYDX2=(XRR*XRR*TR - XR*XR*TRR)/(XR*XRR*(XRR-XR))

              ANUSLT(J)=DYDX1+DYDX2 
              ANU_TP = (TL-TR)/(XL-XR)     

              WRITE(55,*)J,YM(J),TL,TR,TRR, ANU_TP,ANUSLT(J)
              END DO

234       FORMAT (4E20.8)

C      EVALUATION OF  AVERAGE NUSSELT NOS.
C
	SUM=0.
	ALT=0.
        I= 1
	DO J=1,NJ
	IJB=I+(J-1)*NI
	IJT=IJB+NI
	DY=F(ISY+IJT)-F(ISY+IJB)
	DX=F(ISX+IJT)-F(ISX+IJB)
        DS=SQRT(DX*DX+DY*DY)
        SUM=SUM+ANUSLT(J)*DS
        ALT=ALT+DS
        WRITE(56,*) J,DX,DY,DS,ANUSLT(J),SUM,ALT
        END DO

        AVG_NUSLT=ABS(SUM)
C       WRITE(55,*)'AVERAGE NU =',AVG_NUSLT

	
C
C***********************************************************************
c        TRANSVERSE PROFILES OF VARIABLE AT GIVEN Y 
C**************************************************************************

	  do nys=1,no_ystat

		ygiv=yg(nys)

		if (iblk.eq.iblk_ygiv(nys)) then


	  nitot=isolve_e-isolve_b+3

          WRITE(51,*)
     >     'ZONE T =" yg=',ygiv,', ib=',iblk, '",I=',NITOT,',J=',1
C
C     LONGITUDINAL PROFILES OF VARIABLE AT GIVEN Y (i,e X vs V )
C
                do j=1,nj
                  ii=1+(j-1)*ni
                  ab(j)=f(isy+ii)
                end do

            CALL SEARCH(nj,ab,ygiv,FIND,JL)
C
          WRITE(*,*)'JL=',JL

          DO  I=isolve_b-1,isolve_e+1
              II = I + (JL-1)*ni

            x=f(isx+II)  
            yl=f(isy+II)
            yr=f(isy+II+NI)
            ul=f(isu+II)/uref
            ur=f(isu+II+NI)/uref

            CALL LINT(yl,yr,ul,ur,YGIV,U)

          vl=f(isv+II)/uref
          vr=f(isv+II+NI)/uref

          CALL LINT(yl,yr,vl,vr,YGIV,v)

          Tl=f(ISS+II)
          Tr=f(ISS+II+NI)

          CALL LINT(yl,yr,Tl,Tr,yGIV,T)


            WRITE(51,234)x,u,v,t
            ENDDO
		end if

	end do ! end of y=constant loop

C************************************************************************8888

c          IDENTIFYING THE REFERENCE POINT PRESSURE.
c
           if (iblk .eq. ibref_pr) then
              iipref=ipref+(jpref-1)*ni
              pref=f(isp+iipref)
          end if
              write(*,*)iblk,'PREF   ',ipref,jpref,pref

c
c          PLOTTING THE SURFACE PRESSURE AND SKIN FRICTION.
c

	do nw=1,nwals ! nesting over no of walls in each block

              IFLAG=IND_FLAG(NW)

C   IFLAG = 1 FOR WALL AT NORTH  IFLAG = 2 FOR WALL AT SOUTH
C   IFLAG = 3 FOR WALL AT EAST   IFLAG = 4 FOR WALL AT WEST

	   if (IFLAG .EQ. 1 .OR. IFLAG .EQ. 2) THEN

	   NPTS=ISOLVE_E-ISOLVE_B+1
            WRITE(12,*)'ZONE T = " BL=',iblk,',flg=',iflag,' ",I=',NPTS,',J=',1

C Computation of cp, cf & yplus for NORTH/SOUTH WALL

	     if (iflag .eq. 1) then
C North wall
		jnum=jsolve_e
		ndel=ni
	      else
C South wall
                JNUM=jsolve_b
		ndel=-ni
	      end if

C
C COMPUTATION OF CP, CF AND YPLUS
C
             do i=ISOLVE_B,ISOLVE_E
             II = I + (JNUM-1)*NI 
             IP=II
             CP=(F(ISP+II)-PREF)/(0.5*UREF*UREF*DENSIT)
C
C  AVOID REFERENCING FOR THE BLANKED REGION
C
             CP=CP*(1-KBLK(II))
             YSURF=F(ISY+II)
             XSURF=F(ISX+II)
c
        CALL FCOMP(IP,IFLAG,PDRAG,FDRAG,PLIFT,FLIFT,tau,yplus,PARL)
c
              cf=tau/(0.5*uref*uref*densit)
             CF=CF*(1-KBLK(II))
              YPLUS=YPLUS*(1-KBLK(II))
              if (i .eq.isolve_e) then 
              utsq=abs(tau)/densit
              tau_w=tau
              end if 

	     if (NFILE .EQ. 1) then
             WRITE(12,5)i,xsurf,cp,-cf,yplus
             END IF

              END DO ! end of i loop 

	   else


C Computation of cp, cf & yplus for EAST/WEST WALL
C   IFLAG = 3 FOR WALL AT EAST   IFLAG = 4 FOR WALL AT WEST

	   NPTS=JSOLVE_E-JSOLVE_B+1
          WRITE(12,*)'ZONE T = " BL=',iblk,',flg=',iflag,' ",I=',NPTS,',J=',1
             if (iflag .eq. 3) then
C East  wall
                inum=isolve_e
                ndel=1
              else
C West  wall
                INUM=isolve_b
                ndel=-1
              end if

C COMPUTATION OF CP, CF AND YPLUS
C
             do j=JSOLVE_B,JSOLVE_E
              II = INUM + (J-1)*NI
              IP=II
              CP=(F(ISP+II)-PREF)/(0.5*UREF*UREF*DENSIT)
              YSURF=F(ISY+II)-YREF
              XSURF=F(ISX+II)
c
        CALL FCOMP(IP,IFLAG,PDRAG,FDRAG,PLIFT,FLIFT,tau,yplus,PARL)
c
              cf=tau/(0.5*uref*uref*densit)

             if (NFILE .EQ. 1) then
             WRITE(12,5)J,ysurf,cp,-cf,yplus
             end if

              END DO ! end of j loop

	    end if! end of iflag check

	end do ! end of no. of walls loop
C
C  EVALUATION OF THE VORTICITY FIELD DUE TO NATURAL CONVECTION
C
c      COMPUTATION OF BOUNDARY INDICES
c
       CALL BINDEX
c
C    EVALUATION OF DUDX AND DUDY AT ALL CELL CENTRES
C    
        
        DO J=1,NJ
        DO I=1,NI
         II=I+(J-1)*NI
          PH(I,J)=F(ISU+II)
        END DO
        END DO
        
        DO J = 2,NJ-1
        DO I = 2,NI-1
          II=I+(J-1)*NI
          
             CALL DIFF(I,J,PH,DPDX,DPDY)
             
              DUDX(I,J)=DPDX*(1-KBLK(II))
              DUDY(I,J)=DPDY*(1-KBLK(II))
                                                                   
          END DO
          END DO
c
C    EVALUATION OF DVDX AND DVDY AT ALL CELL CENTRES
C    
          
        DO J=1,NJ
        DO I=1,NI
         II=I+(J-1)*NI
          PH(I,J)=F(ISV+II)
        END DO
        END DO
        
        DO J = 2,NJ-1
        DO I = 2,NI-1
          II=I+(J-1)*NI
          
             CALL DIFF(I,J,PH,DPDX,DPDY)
             
              DVDX(I,J)=DPDX*(1-KBLK(II))
              DVDY(I,J)=DPDY*(1-KBLK(II))
          END DO
          END DO
c
C       EVALUATION OF BOUNDARY DERVATIVES
C       
        CALL BNDDER
C

	DO J=1,NJ
	DO I=1,NI
	VOR(I,J)=DUDY(I,J)-DVDX(I,J)
C       divij=dudx(i,j)+dvdy(i,j)
c	write(15,*)j,i,vor(i,j),divij
	end do
	end do	
C         
	  WRITE(222,*)'ZONE T = "',ntstep,' ",I = ',NI,',J = ',NJ
c
	    DO J = 1,nj
            DO I = 1,ni
              II = I+(J-1)*NI
              YCENT=F(ISY+II)
              XCENT=F(ISX+II)
			WRITE(222,*)XCENT,YCENT,vor(i,j)
		   END DO
            END DO

c
c   BLOCK LOOP ENDS HERE
c

C
        END DO    ! END OF BLOCK LOOP 
C

5       FORMAT(i3,4(F14.8,1X))
6       FORMAT(5(F10.5,1X))


	end do
999	continue
	CLOSE(20)
             
     
      STOP
      END
c
c*******************************
        SUBROUTINE DIFF(IGIV,JGIV,PH,DPDX,DPDY)
C*****************************************************************
C      
C      EVALUATION OF SPATIAL GRADIENT OF A VARIABLE
C      IN A FIELD PROBLEM USING GREEN'S THEOREM
C     
      INCLUDE 'postcom'
      
      COMMON  PHC(5)
      DIMENSION PH(NX,NY)

C       
        CALL CORNPHI(IGIV,JGIV,PH)
C       
        
        CALL DIFFPHI(IGIV,JGIV,DPDX,DPDY)

        RETURN
        END
        
C**********************************************************************
           SUBROUTINE CORNPHI(IGIV,JGIV,PH)
C**********************************************************************
C  THIS SUBROUTINE CALCULATES ALL THE FOUR NEIGHBOURING CORNER PHI VALUES
C FOR A GIVEN CENTRAL POINT IN A CONTROL VOLUME SPECIFIED WITH I,J VALUES
C***********************************************************************

      INCLUDE 'postcom'
      COMMON  PHC(5)
      DIMENSION PH(NX,NY)
      
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
        
        SUBROUTINE DIFFPHI(IGIV,JGIV,DPDX,DPDY)
C***********************************************************************
C SUBROUTINE DIFFPHI CALCULATES DIFFERENTIAL OF PHI W.R.T X OR Y
C************************************************************************

      INCLUDE 'postcom'
      
      COMMON  PHC(5)
      DIMENSION XX(5),YY(5)
C     
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
        
                                                                          
        DPDX=0.
        DPDY=0.
        DXAB=XX(1)-XX(3)
        DYAB=YY(1)-YY(3)
        DXCD=XX(2)-XX(4)
        DYCD=YY(2)-YY(4)
        
         AREA=ABS(0.5*(DXAB*DYCD-DXCD*DYAB))
         
        DO ISIDE = 1,4
        PHAVG=0.5*(PHC(ISIDE)+PHC(ISIDE+1))
        AN1=-(YY(ISIDE)-YY(ISIDE+1))
        AN2=(XX(ISIDE)-XX(ISIDE+1))
        DPDX=DPDX+PHAVG*AN1
        DPDY=DPDY+PHAVG*AN2
C	write(16,*)area,an1,an2,phavg	
        END DO
        DPDX=DPDX/AREA
        DPDY=DPDY/AREA
C	write(16,*)dpdx,dpdy	
        
        RETURN
        END
        
        SUBROUTINE BNDDER
                                                                          
      INCLUDE 'postcom'
      
        I=1
        
        DO J=1,NJ
          DUDX(I,J)=DUDX(I+1,J)
          DUDY(I,J)=DUDY(I+1,J)
          DVDX(I,J)=DVDX(I+1,J)
          DVDY(I,J)=DVDY(I+1,J)
        END DO
        
        I=NI
        
        DO J=1,NJ
          DUDX(I,J)=DUDX(I-1,J)
                                                                          
          DUDY(I,J)=DUDY(I-1,J)
          DVDX(I,J)=DVDX(I-1,J)
          DVDY(I,J)=DVDY(I-1,J)
        END DO
        
        J=1
        
        DO I=1,NI
          DUDX(I,J)=DUDX(I,J+1)
          DUDY(I,J)=DUDY(I,J+1)
          DVDX(I,J)=DVDX(I,J+1)
          DVDY(I,J)=DVDY(I,J+1)
        END DO
        
        J=NJ
        DO I=1,NI
          DUDX(I,J)=DUDX(I,J-1)
          DUDY(I,J)=DUDY(I,J-1)
          DVDX(I,J)=DVDX(I,J-1)
          DVDY(I,J)=DVDY(I,J-1)
        END DO
        
        RETURN
                                                                          
        END
C                                                                               
C***********************************************************************        
C                                                                               
      SUBROUTINE FCOMP(IP,IFLAG,PDRAG,FDRAG,PLIFT,FLIFT,tau,yplus,PARL)
C                                                                               
C     CALCULATES THE LIFT DRAG ETC.                                             
C***********************************************************************        
C                                                                               
       INCLUDE 'postcom'

	cdqr=0.548
	elog=9.793
	cappa=0.4187
                                                                        
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

      IB  =IP+( NI *  MFLAG1 -  NI* MFLAG2 +      MFLAG3 -     MFLAG4) 
      ID1 =IP+((NI+1)*MFLAG1 +      MFLAG2+(NI+1)*MFLAG3 +  NI*MFLAG4)  
      ID2 =IP+( NI *  MFLAG1 +   0* MFLAG2 +      MFLAG3 +  0 *MFLAG4)  
C                                                                       
C     NORMAL VECTOR, AREA AND NORMAL DISTANCE   
C                        
             DX=(F(ISCOX+ID1)-F(ISCOX+ID2))
             DY=(F(ISCOY+ID1)-F(ISCOY+ID2)) 
	     DN1=-2*DY
	     DN2=2*DX
	     AREA=0.5*SQRT(DN1*DN1+DN2*DN2)
	     A1=0.5*DN1/(AREA+SMALL)
	     A2=0.5*DN2/(AREA+SMALL)
	     XSURF=0.5*(F(ISCOX+ID1)+F(ISCOX+ID2))
	     YSURF=0.5*(F(ISCOY+ID1)+F(ISCOY+ID2))
	     DELTA=ABS((F(ISX+IP)-XSURF)*A1 + (F(ISY+IP)-YSURF)*A2)
c    
C  SHEAR STRESS AT THE WALL
C                                          
	     A11=1.-A1*A1
	     A22=1.-A2*A2
	     A12=A1*A2

             UWAL=F(ISU+IB)
             VWAL=F(ISV+IB)

	     VXCOMP=A11*F(ISU+IP)-A12*F(ISV+IP)
	     VYCOMP=A22*F(ISV+IP)-A12*F(ISU+IP)
	     VP=SQRT(VXCOMP*VXCOMP+VYCOMP*VYCOMP)
	     PARL=(DX*VXCOMP+DY*VYCOMP)
   
	  sqrtk=sqrt(f(iste+ip))
          XYPLUS=F(ISDEN+IP)*SQRTK*CDQR*DELTA/VISCOS+SMALL       
       
          F01=MAX(0.,(XYPLUS-11.63)/(ABS(XYPLUS-11.63)+SMALL))        
          TMULT=(1.-F01)*VISCOS/DELTA                                   
     >          +F01*F(ISDEN+IP)*CDQR*SQRTK*CAPPA/LOG(ELOG*XYPLUS)

        tau=tmult*(-a2*(f(isu+ip)-uwal)+a1*(f(isv+ip)-vwal))

          utau= sqrt(abs(tau)/densit)
          yplus = delta*utau*densit/VISCOS
	APROJX=AREA*A1
	APROJY=AREA*A2


 	FDRAG=TMULT*AREA*VXCOMP
        PDRAG=-((F(ISP+IP)-pref))*APROJX

        PLIFT=-(F(ISP+IP)-pref)*APROJY
 	FLIFT=TMULT*AREA*VYCOMP

	
       RETURN                                                                   
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

       INCLUDE 'postcom'
C
      READ(7,rec=iblk) nim,njm,((f(iscox+i+(j-1)*ni),I=2,ni),j=2,nj),
     >                         ((f(iscoy+i+(j-1)*ni),I=2,NI),j=2,nj)


      DO  J=1,NJ                                                      
	ii=1+(j-1)*ni
	f(iscox+ii)=f(iscox+ii+1)
	f(iscoy+ii)=f(iscoy+ii+1)
      end do
C 
      DO I=1,NI                                                      
	ii=i+(1-1)*ni
	f(iscox+ii)=f(iscox+ii+ni)
	f(iscoy+ii)=f(iscoy+ii+ni)
      end do	
	return

      END 

	SUBROUTINE GEOMET

       INCLUDE 'postcom'

	do j=1,nj
	  do i=1,ni
            ISKP=1                                                            
            JSKP=NI                                                           
            IF(I.EQ.NI) ISKP=0                                                
            IF(J.EQ.NJ) JSKP=0                                                
                                                                        
	   ii=i+(j-1)*ni
           xp=0.25*(F(IScox+II)+F(IScox+II+JSKP)                             
     >            +F(IScox+II+JSKP+ISKP)+F(IScox+II+ISKP))                  
           yp=0.25*(F(IScoy+II)+F(IScoy+II+JSKP)                             
     >            +F(IScoy+II+JSKP+ISKP)+F(IScoy+II+ISKP))                  

           F(ISX+II)=XP                                                      
           F(ISY+II)=YP                                                      
	  end do
	end do
	return

      END 
c
        SUBROUTINE BINDEX
        
      INCLUDE 'postcom'
      
        DO J=1,NJ
        DO I=1,NI
        
          II=I+(J-1)*NI
          
        KW(II)=0
        KE(II)=0
        KS(II)=0
                                                                          
        KN(II)=0
        KBLK(II)=0
        
        END DO
        END DO
        
        I=2
        DO J=2,NJ-1
          II=I+(J-1)*NI
         KW(II)=1
        END DO
        
        
        I=NI-1
C       
        DO J=2,NJ-1
          II=I+(J-1)*NI
        KE(II)=1
        END DO
        
        
        
        J=2
        
        DO I=2,NI-1
          II=I+(J-1)*NI
        KS(II)=1
        
        END DO
        
        J=NJ-1
        DO I=2,NI-1
          II=I+(J-1)*NI
                                                                         
        KN(II)=1
        END DO

        
        RETURN
        
        END
C
C  *****************************************************************
           SUBROUTINE SEARCH(NDAT,PH1,PHGIV,FIND,INTV)
C
C   *****************************************************************

         INCLUDE 'postcom'
          DIMENSION PH1(nmax)
C
          KFIN=NDAT-1
          FIND=0
          DO K=1,KFIN
           xl=ph1(k)
           xr=ph1(k+1)
          if (phgiv .ge. xl .and. phgiv .le. xr) then
             intv=k
c             write(99,*)phgiv,xl,xr,intv
              FIND =1
             return

          end if
             END DO
           intv=k
	           find=0

          RETURN
          END
C
C*********************************************************************
        SUBROUTINE LINT(XL,XR,YL,YR,XOUT,YOUT)
C
C*********************************************************************
C
C       THIS SUBROUTINE  FITS A STRAIGHT LINE FORM
C       GIVEN END POINTS. IT RETURNS THE ORDINATE(YOUT)
C       AT A GIVEN POINT X=XOUT.
C
        YOUT=YL+(XOUT-XL)*(YR-YL)/(XR-XL)

        return
        end
C
C*********************************************************************
