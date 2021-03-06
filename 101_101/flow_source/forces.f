***********************************************************************        
      SUBROUTINE FORCES
C      
C     CALCULATES THE FORCES AND MOMENTS                                             
C************************************************************************        
      INCLUDE 'com2d'                                                           
C                                                                               
	dimension adfp(3),adff(3),adm(3)
	dimension adfp_w(3),adff_w(3),adm_w(3)
      dimension vnwal(nmax)
    
	utemp=sqrt(uref*uref+vref*vref)
	PLANFORM=1
	if (iblk .eq. 1) then
	do i=1,3
	  adfp_pr(i) = 0.
	  adff_pr(i) = 0.
	  adm_pr(i) = 0. 
        end do
	end if

	do i=1,3
	  adfp(i)=0.
	  adff(i)=0.
	  adm(i)=0.
	end do

 	do ifl=1,4

C   This change is valid only for problems having only wall at south


         call bforce(ifl,adfp_w,adff_w,adm_w,vnwal)                                                          

	 do i=1,3
	   adfp(i)=adfp_w(i)+adfp(i)
	   adff(i)=adff_w(i)+adff(i)
	   adm(i)=adm(i)+adm_w(i)
	 end do
	
	end do
C                       
C EVALUATION OF THE FORCES AND MOMEMTS ALONG THE FLOW DIRECTION
C
C FORCE COMPONENTS DUE TO PRESSURE FORCE ONLY 
C
       adf_flow_1 = adfp(1)*cos(alpha)+ adfp(2)*sin(alpha)
	adfp_pr(1)=adfp_pr(1)+adf_flow_1

       adf_flow_2 = -adfp(1)*sin(alpha)+adfp(2)*cos(alpha)
	adfp_pr(2)=adfp_pr(2)+adf_flow_2
C
C FORCE COMPONENTS DUE TO FRICTION FORCE ONLY 
C
       adf_flow_1 = adff(1)*cos(alpha)+adff(2)*sin(alpha)
	adff_pr(1)=adff_pr(1)+adf_flow_1

       adf_flow_2 = -adff(1)*sin(alpha)+adff(2)*cos(alpha)
	adff_pr(2)=adff_pr(2)+adf_flow_2

C
C MOMENT COMPONENTS DUE TO FRICTION AND PRESSUE FORCES
C
        adm_flow_3 = -adm(1)*sin(alpha)*cos(beta)+adm(3)*cos(alpha)
     >            -adm(2)*sin(alpha)*sin(beta)
	adm_pr(3)=adm_pr(3)+adm_flow_3

	if (iblk .eq. nblock) then
	do i=1,3
	  adfp_pr(i) = adfp_pr(i)/(0.5*Utemp*Utemp*DENSIT)/(PLANFORM)
	  adff_pr(i) = adff_pr(i)/(0.5*Utemp*Utemp*DENSIT)/(PLANFORM)
	  adm_pr(i) = adm_pr(i)/(0.5*Utemp*Utemp*DENSIT)/(PLANFORM)
        end do
	end if
C
C  Evaluation of recirculation length at any instant
C
          if (iblk.eq.2) call rec_bl

C
C  Evaluation of the angles for stagnation, separation and reattachment etc.
C            
         call sep_angle(2,vnwal) 

        RETURN                                                                  
        END                                                                     
C                                                                               
C***********************************************************************        
      SUBROUTINE BFORCE(IFLAG,ADFP_W,ADFF_w,ADM_W, vnwal)                                                  
C     CALCULATES THE FORCES AND MOMENTS ABOUT THE CARTESIAN AXES
C     DUE THE PRESENCE OF WALL AT ANY FOUR BOUNDARIES OF THE DOMAIN
C     IFLAG = 1(NORTH), 2(SOUTH), 3(EAST), 4(WEST)
C***********************************************************************        
      INCLUDE 'com2d'                                                           
      dimension adfp_w(3),adff_w(3),adm_w(3),vnwal(nmax)
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
C                                                                               
      MWAL = MWALN*MFLAG1 + MWALS*MFLAG2 + MWALE*MFLAG3 + MWALW*MFLAG4                        
      ISANOW =  ISAN*MFLAG1 +  ISAS*MFLAG2 +  ISAE*MFLAG3 + ISAW*MFLAG4                        
	do i=1,3
	  adfp_w(i)=0.
	  adff_w(i)=0.
	  adm_w(i)=0.
	end do

      DO  IM=1,MWAL                                                         
          IIP= INDWN(IM)*MFLAG1 + INDWS(IM)*MFLAG2 + INDWE(IM)*MFLAG3           
     >       + INDWW(IM)*MFLAG4            
	jwall=iip/NI+1
	iwall=iip-(jwall-1)*ni
      IIB=IIP+ NI*MFLAG1 - NI*MFLAG2 + MFLAG3 -MFLAG4            
      ID1=IIP+(1+NI)*MFLAG1+MFLAG2+(1+NI)*MFLAG3+NI*MFLAG4            
      ID2=IIP+NI*MFLAG1 +0*MFLAG2 + MFLAG3+0*MFLAG4           
	idir=-1*mflag1+1*mflag2+1*mflag3-1*mflag4

C     NORMAL VECTOR, AREA AND NORMAL DISTANCE                                   
             DAX=(F(ISCOX+ID1)-F(ISCOX+ID2))*idir                            
             DAY=(F(ISCOY+ID1)-F(ISCOY+ID2))*idir                              
         DN1=-2*DAY                                                    
         DN2=2*DAX                                                    
        AREA=0.5*SQRT(DN1*DN1+DN2*DN2)                                  
          A1=0.5*DN1/(AREA+SMALL)                                               
          A2=0.5*DN2/(AREA+SMALL)                                               
      XSURF=0.5*(F(ISCOX+ID1)+F(ISCOX+ID2))          
      YSURF=0.5*(F(ISCOY+ID1)+F(ISCOY+ID2))          
       DELTA=ABS((F(ISX+IIP)-XSURF)*A1+(F(ISY+IIP)-YSURF)*A2)+SMALL                                   
C     SHEAR STRESS AT THE WALL                                                  
            A11=1.-A1*A1                                                        
            A22=1.-A2*A2                                                        
            A12=A1*A2                                                           
          SQRTK=SQRT(ABS(F(ISTE+IIP)))                                    

          YPLS=F(ISDEN+IIP)*SQRTK*CDQR*DELTA/VISCOS                            
          TMULT=VISCOS/(DELTA+SMALL)  
         IF(YPLS.GT.11.6)                                                      
     >    TMULT=F(ISDEN+IIP)*CDQR*SQRTK*CAPPA/LOG(ELOG*YPLS)                  
          FACTOR=TMULT*AREA                                                     

          UPUB=F(ISU+IIP)-F(ISU+IIB)                                            
          VPVB=F(ISV+IIP)-F(ISV+IIB)                                            

          VXCOMP=UPUB*A11-VPVB*A12
          VYCOMP=VPVB*A22-UPUB*A12
  
          SX=-FACTOR*VXCOMP*mysolv(iip)
          SY=-FACTOR*VYCOMP*mysolv(iip)

          VNWAL(IM)=DAX*VXCOMP+DAY*VYCOMP
          
          
          PSURF=F(ISP+IIP)

          PX=-PSURF*AREA*A1*mysolv(iip)
          PY=-PSURF*AREA*A2*mysolv(iip)

          adfp_w(1)=adfp_w(1)+PX
          adfp_w(2)=adfp_w(2)+PY

          adff_w(1)=adff_w(1)-SX
          adff_w(2)=adff_w(2)-SY

          adm_w(3)=adm_w(3)+(PX-SX)*F(ISY+IIP)-(PY-SY)*F(ISX+IIP)

	END DO

       RETURN                                                                   
       END                                                                      

C***********************************************************************        
      SUBROUTINE rec_bl
C                                                         
C     CALCULATES THE LENGTH OF THE RECIRCULATION ZONE AT ANY INSTANT
C                                             
C************************************************************************        
      INCLUDE 'com2d' 
      DIMENSION UAX(NMAX),XAX(NMAX)                                                          
C      
C  Calculate the wake centerline location ( between I= IFIND and I=IFIND-1) 
C        
      YGIV = 0.
             
      DO  J=JSOLV_B,JSOLV_E

       DO I=2,NI
       II=I+(J-1)*NI
        YT=F(ISY+II)-YGIV
        YB=F(ISY+II-1)-YGIV
        IF(YT*YB.LE.0) THEN
         IFIND=I
         GO TO 999
        END IF
       END DO   !   END OF I LOOP

 999  CONTINUE

       IIT=IFIND+(J-1)*NI
       IIB=IIT-1
       YT=F(ISY+IIT)
       YB=F(ISY+IIB)
       XT=F(ISX+IIT)
       XB=F(ISX+IIB)
       UT=F(ISU+IIT)
       UB=F(ISU+IIB)
       UAX(J)=UB+(UB-UT)*(YGIV-YB)/(YB-YT)
       XAX(J)=XB+(XB-XT)*(YGIV-YB)/(YB-YT)

      END DO    !  END OF J LOOP
       
      UGIV=0.
      WAKEL=0.
                   
      DO J= JSOLV_B,JSOLV_E
      P1=UGIV-UAX(J)
      P2=UGIV-UAX(J+1)
      PRO=P1*P2
      IF(PRO.LE.0) THEN
      WAKEL=XAX(J)+(XAX(J)-XAX(J+1))*(UGIV-UAX(J))/(UAX(J)-UAX(J+1))
      GO TO 1000
      END IF
      END DO    ! END OF J LOOP
 
 1000 CONTINUE  
C
C  WAKEL IS TO BE MEASURED FROM THE REAR STAGNATION POINT OF THE CYLINDER
C
      WAKEL= WAKEL- 0.5

      RETURN
      END 

C***********************************************************************        
      SUBROUTINE sep_angle(iflag,vnwal)
C                                                         
C CALCULATES THE LOCATION OF THE SEPARATION,STAGNATION AND REATTACHMENT POINTS
C                                             
C************************************************************************        
      INCLUDE 'com2d' 
      DIMENSION vnwal(nmax)

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
C                                                                               
        MWAL = MWALN*MFLAG1 + MWALS*MFLAG2 + MWALE*MFLAG3 + MWALW*MFLAG4                        
C
C Loop along the wall boundary 
C To check the sign of tangential velocity at near wall node     
C
          PI=4.*ATAN(1.)

          isep=0

          DO IM=2,MWAL-1
 
         IIPM= INDWN(IM)*MFLAG1 + INDWS(IM)*MFLAG2 + INDWE(IM)*MFLAG3           
     >                          + INDWW(IM)*MFLAG4 

         IIBM=IIPM+ NI*MFLAG1 - NI*MFLAG2 + MFLAG3 -MFLAG4                       
         XSM=F(ISX+IIBM)
         YSM=F(ISY+IIBM)
         IIPP= INDWN(IM+1)*MFLAG1 + INDWS(IM+1)*MFLAG2 + INDWE(IM+1)*MFLAG3           
     >                           + INDWW(IM+1)*MFLAG4 

         IIBP=IIPP+ NI*MFLAG1 - NI*MFLAG2 + MFLAG3 -MFLAG4                       
         XSP=F(ISX+IIBP)
         YSP=F(ISY+IIBP)

         pro=vnwal(im)*vnwal(im+1)
         pro=pro/abs(pro)
         if(pro .lt. 0.) then
          isep=isep+1 
         FACT=vnwal(im)/(vnwal(im)-vnwal(im+1))   
         XSEP=FACT*XSP+(1.-FACT)*XSM
         YSEP=FACT*YSP+(1.-FACT)*YSM
         if(iblk.eq.1) THETA(iblk,isep)=ATAN2(-XSEP,-YSEP)-PI/2.            
         if(iblk.eq.2) THETA(iblk,isep)=ATAN2(XSEP,YSEP)+PI/2. 
         THETA(iblk,isep)=THETA(iblk,isep)*180./pi           
          end if 

         end do        ! end of IM loop  
          nsep(iblk)=isep

         return
         end   

