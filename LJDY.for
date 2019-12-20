	
!     THIS PROGRAM IS FOR CALCULATING THE DENSITIES OF ONE-COMPONENT FLUID
!     NEAR A SPHERICAL PARTICULATE.                   
!                                 ---EDITED BY YE YI,  6,11,2018
         
	PROGRAM HARDWALLS_HS_LJ

	IMPLICIT NONE
      include 'omp_lib.h'
	INTEGER NZ,I,K,NS,NW,NG,dim,RI,nc,IC,JJ,nr,id,j,KK,II,imax
	DOUBLE PRECISION DZ,PI,F,ETA,Y,SADD,SUMP,AMU,MUATT,Rc0,S,kmct
	PARAMETER (NZ=1024*8,dim=2)
	DOUBLE PRECISION RHOD(2*NZ),RO(NZ),RHOs(NZ),RD(NZ)
	DOUBLE PRECISION FFs(NZ),FFD(NZ),RH(dim,0:NZ),ROU
	DOUBLE PRECISION SIGMA,Z_R,TEMP,RCI,RHOBC,TEN,r,RHOW(NZ)
	DOUBLE PRECISION RC,RS,A,B,X,DATT,DT,JT,dk,sum,sita,DEFFN
	
	DOUBLE PRECISION S1,SD(2),FATT,SUM3,attb(NZ),mx
	PARAMETER (PI=3.141592654)
	DOUBLE PRECISION RHO2(NZ),RO2(NZ),RHO12(NZ),RD2(NZ)
	DOUBLE PRECISION FF2(NZ),RH2(NZ),RHOW2(NZ),hsb(NZ),rmax
	DOUBLE PRECISION RHOB,TF,DEFF,H,Z0,ZH,TEMP0(3000),FEX(2),MSD,MSD2
	DOUBLE PRECISION RHOC,RHOBX(0:4),MUBX(0:4),PHIX(0:4),ADS,SUM1,SUM01
	DOUBLE PRECISION r_vari(0:nz),k_vari(0:nz),sk(nz),mct1,mct
      DOUBLE PRECISION dff(dim,NZ),roudff(dim,0:NZ),droudff(dim,0:NZ)
      DOUBLE PRECISION RHT(dim,NZ),df1(Dim,NZ),FF(dim,0:NZ),zz1
      DOUBLE PRECISION PR(NZ),PII(NZ),FR(NZ),FI(NZ),FROU(NZ)
      DOUBLE PRECISION rr,cr(dim,dim,NZ),etan,cstif(dim,dim,NZ)
       DOUBLE PRECISION z1,phi,crstiff(nz),RHON(NZ),FFN(NZ)
       DOUBLE PRECISION AA,BB,rr1,dye,DXT2(0:NZ),shear(dim)
      DOUBLE PRECISION at,bt,at2,bt2,GR(0:NZ),fre
      DOUBLE PRECISION temprr1,tempGz(0:NZ)
      CHARACTER(LEN=100) :: FILENAME,FORM
	COMMON /MUEX/ attb,hsb
	COMMON /HARDWALL/ Z0,ZH,TF,DEFF,DEFFN,RHOB,S,SIGMA
      COMMON /angle/ sita,a,b
      COMMON /dyna/ dt

	COMMON /step/ dz

      SIGMA=1.0D0
      
	TF=4.  !400./46.
      
      shear(1)=5.0
      shear(2)=0.0
      fre=10.0
      
      etan=2.0
      
      

	Rc0=5.0D0
	F=0.1D0
      

	DEFF=SIGMA*(1+0.2977D0*TF)/(1.D0+0.33163D0*TF+0.00104771D0*TF**2)
      DEFFN=4.0*deff
	DZ=0.05D0*deff
      dk=2.*pi/dz/nz/2.
      ZH=DZ*NZ*2.
      kmct=int(2.*pi/DEFF/dk)
      
      do i=0,nz
	r_vari(i)=i*dz
	k_vari(i)=i*dk
	if(i.eq.0)then
	r_vari(i)=1.0e-10
	k_vari(i)=1.0e-10
      endif
      enddo
     
91    FORMAT(1X,2000g15.6)
      
	
		    
        
	OPEN(3,FILE='roud0.dat') 
	DO I=1,NZ
	read(3,*) RHOD(i)
     
      enddo
	CLOSE(3)    
      
     
     
            
      
      DO I=1,NZ
          Y=abs(DZ*DFLOAT(I-nz/2))
	RHOs(i)=(100./pi)**0.5*exp(-100.*y**2.)/deffn**2.
c      if(RHOs(i).le.1.0e-20) RHOs(i)=1.0e-20
      ENDDO
c      RHOs(1)=RHOs(2)*2.-RHOs(3)
      
      
      sum01=0.0
      DO I=201,nz-200
          Y=DZ*DFLOAT(I)
	sum01=sum01+RHOs(i)*DZ
      ENDDO
      
      
      DO I=1,NZ
      zz1=dz*float(i-1)
      if(zz1.le.2.**(1./6.)*deff)then
          cr(2,2,i)=-1./tf
          elseif(zz1.le.3.*deff)then
	cr(2,2,i)=4./tf*((DEFF/zz1)**12.-(DEFF/zz1)**6.)
          else
              cr(2,2,i)=0.0
          endif
          
      IF(zz1.le.0.5*(deff+deffn)) THEN
     	cr(2,1,i)=0.0
      ELSE   
	cr(2,1,i)=-etan/Tf*exp(-(zz1-0.5*(deff+deffn))/deff)
      ENDIF
      cr(1,2,i)=0.0   !cr(1,2,i)  
      cr(1,1,i)=0.0

      ENDDO
      
      
c      GOTO 1142
      
      do ii=1,dim
      do jj=1,dim
   	do kk=1,NZ
	zz1=kk*DZ
	rmax=10.d0-zz1
	imax=int(rmax/dz)
      if (rmax.LE.1.0e-4) then
      cstif(ii,jj,kk)=0.0
	else
	sum=0.0
	do i=1,imax
	rr=sqrt(zz1**2.+(dz*float(i))**2.)
	sum=sum+cr(ii,jj,int(rr/dz))*(dz*float(i))*dz
	enddo 
      cstif(ii,jj,kk)=sum*2.*pi
      endif
      enddo
      enddo
      enddo
	
1142  CONTINUE
	
	
	
      
     
	
      
      
     
	

!  calculate the recurrence functions and the effective external potential 
	NG=1
      dt=0.00001
       JT=0.0
       OPEN(191,FILE='mct.dat')
	do ii=1,10000000
      if(ii.gt.10000) dt=0.00005
      if(ii.gt.100000) dt=0.0001
      if(ii.gt.1000000) dt=0.0005

          CALL FFLR(1,NZ,DZ,shear,cstif,FFs,RHOs,RHOD,dim)

          CALL FFLR(2,NZ,DZ,shear,cstif,FFD,RHOs,RHOD,dim)
      
      rht=0.0
      do j=201,nz-200
              rht(1,J)=FFs(j)
              rht(2,J)=FFd(j)
             
      ENDDO
      
            
      
      
      sum1=0.0
      DO I=201,nz-200
          Y=DZ*DFLOAT(I)
	sum1=sum1+rht(1,I)*DZ
      ENDDO
      
      do i=201,nz-200
	rht(1,i)= rht(1,i)/SUM1*SUM01 
      enddo
          
      
      
      if(modulo(log10(dble(ii)),1.0).eq.0.0)then
              write(FORM,*) ii
              write(FILENAME,*) "roun",TRIM(FORM),".dat"
              
      open(unit=13,file=FILENAME)
      DO I=1,NZ
	Y=DZ*FLOAT(I-1)
	WRITE(13,*) Y,(RHt(j,i),j=1,dim)
      enddo
      close(13)
      endif
      
      
      if(modulo(log10(dble(ii)/2.),1.0).eq.0.0)then
              write(FORM,*) ii
              write(FILENAME,*) "roun",TRIM(FORM),".dat"
              
      open(unit=13,file=FILENAME)
      DO I=1,NZ
	Y=DZ*FLOAT(I-1)
	WRITE(13,*) Y,(RHt(j,i),j=1,dim)
      enddo
      close(13)
      endif
      
      
      if(modulo(log10(dble(ii)/3.),1.0).eq.0.0)then
              write(FORM,*) ii
              write(FILENAME,*) "roun",TRIM(FORM),".dat"
              
      open(unit=13,file=FILENAME)
      DO I=1,NZ
	Y=DZ*FLOAT(I-1)
	WRITE(13,*) Y,(RHt(j,i),j=1,dim)
      enddo
      close(13)
      endif
      
      
      if(modulo(log10(dble(ii)/4.),1.0).eq.0.0)then
              write(FORM,*) ii
              write(FILENAME,*) "roun",TRIM(FORM),".dat"
              
      open(unit=13,file=FILENAME)
      DO I=1,NZ
	Y=DZ*FLOAT(I-1)
	WRITE(13,*) Y,(RHt(j,i),j=1,dim)
      enddo
      close(13)
      endif
      
      if(modulo(log10(dble(ii)/5.),1.0).eq.0.0)then
              write(FORM,*) ii
              write(FILENAME,*) "roun",TRIM(FORM),".dat"
              
      open(unit=13,file=FILENAME)
      DO I=1,NZ
	Y=DZ*FLOAT(I-1)
	WRITE(13,*) Y,(RHt(j,i),j=1,dim)
      enddo
      close(13)
      endif
      
      
      if(modulo(log10(dble(ii)/6.),1.0).eq.0.0)then
              write(FORM,*) ii
              write(FILENAME,*) "roun",TRIM(FORM),".dat"
              
      open(unit=13,file=FILENAME)
      DO I=1,NZ
	Y=DZ*FLOAT(I-1)
	WRITE(13,*) Y,(RHt(j,i),j=1,dim)
      enddo
      close(13)
      endif
      
      
      if(modulo(log10(dble(ii)/7.),1.0).eq.0.0)then
              write(FORM,*) ii
              write(FILENAME,*) "roun",TRIM(FORM),".dat"
              
      open(unit=13,file=FILENAME)
      DO I=1,NZ
	Y=DZ*FLOAT(I-1)
	WRITE(13,*) Y,(RHt(j,i),j=1,dim)
      enddo
      close(13)
      endif
      
      
      if(modulo(log10(dble(ii)/8.),1.0).eq.0.0)then
              write(FORM,*) ii
              write(FILENAME,*) "roun",TRIM(FORM),".dat"
              
      open(unit=13,file=FILENAME)
      DO I=1,NZ
	Y=DZ*FLOAT(I-1)
	WRITE(13,*) Y,(RHt(j,i),j=1,dim)
      enddo
      close(13)
      endif
      
      if(modulo(log10(dble(ii)/9.),1.0).eq.0.0)then
              write(FORM,*) ii
              write(FILENAME,*) "roun",TRIM(FORM),".dat"
              
      open(unit=13,file=FILENAME)
      DO I=1,NZ
	Y=DZ*FLOAT(I-1)
	WRITE(13,*) Y,(RHt(j,i),j=1,dim)
      enddo
      close(13)
      endif
      
      
          
      JT=JT+dt
      write(*,*)'jt=', jt
    
 
      
      if(modulo(NG,100).EQ.1)THEN
          
      OPEN(3,FILE='rout.dat') 
	DO I=1,NZ
	Y=DZ*FLOAT(I-1)
	WRITE(3,91) Y,(RHt(j,i),j=1,dim)
      enddo
	CLOSE(3)
      
           
      
      if(modulo(NG,100).EQ.1)THEN
      MSD=0.0
          MSD2=0.0
          DO i=201,nz-200
          Y=abs(DZ*DFLOAT(I-nz/2))
          MSD=MSD+RHt(1,i)*Y**2*DZ
          MSD2=MSD2+RHt(1,i)*Y**4*DZ
          ENDDO
          
      
      write(191,91) JT, MSD,MSD2
      ENDIF
      endif
      
      NG=NG+1
      
c      sum1=0.0
c      do j=1,nz
c       sum1=sum1+RHt(1,J)*dz   
c       sum2=sum2+RHt(1,J)*dz     
c      ENDDO

              
       do j=201,nz-200
              RHOs(J)=RHt(1,J)
              
c              if(j.le.300.and.RHt(2,J).gt.RHt(2,J+1))then
c                  RHOd(J)=RHOd(J)
c                  else
              RHOd(J)=RHt(2,J)
c                  endif
                  
c              if(j.le.500.and.RHt(3,J).gt.RHt(3,J+1))then
c                  RHOn(J)=RHOn(J)
c                  else
           
c                  endif
                  
      ENDDO
       

      
      enddo
      CLOSE(191)
      

	END






	SUBROUTINE FFLR(NN,NZ,DZ,shear,cstif,rht,RHO1,RHO2,dim)
	IMPLICIT NONE
      include 'omp_lib.h'
	INTEGER NZ,I,K,IR,J,KIN,KIP,nk1,NN,Dim
	DOUBLE PRECISION DZ,AMU,Z1,SUM1,S,DATT,y,RHOBC1,S1
	DOUBLE PRECISION RHO(NZ),PHI,FF(NZ),MUATT,rhobc
	DOUBLE PRECISION SIGMA,RHOB,TF,DEFF,H,Z0,ZH,sum,sumc
	DOUBLE PRECISION R1,R2,RHOW(NZ),ROU,HD,ROU1,ROU2,ASEXPR1
	DOUBLE PRECISION MU(NZ),AFREE1,ASFREE1,AFREE2,ASFREE2,def(dim)
	DOUBLE PRECISION df(NZ),f0(NZ),rhodf(NZ),RDFW,DEFFN,cstif(dim,dim,nz)
      DOUBLE PRECISION miuchain(NZ),miuangle(NZ),DEFFI(dim)
      DOUBLE PRECISION sita,a,b,RHO1(NZ),RHO2(NZ),DFN(6,NZ)
      DOUBLE PRECISION con,RHO3(NZ),RH(dim,NZ),shear(dim),rht(NZ)
      DOUBLE PRECISION df1(NZ),dff(NZ),roudff(NZ),droudff(NZ),dt
      DOUBLE PRECISION ratio
	COMMON /HARDWALL/ Z0,ZH,TF,DEFF,DEFFN,RHOB,S,SIGMA
      COMMON /angle/ sita,a,b
      COMMON /dyna/ dt
      
91    FORMAT(1X,2000g15.6)
       
      
 
       
      DEF(1)=deffn
      DEF(2)=deff
      
      
      DO I=1,NZ
          RH(1,I)=RHO1(I)
          RH(2,I)=RHO2(I)
          
       ENDDO
           
      
       
      if(nn.eq.2) then
	call SFM1Dchain(DZ,NZ,RH,def,miuchain)
      endif
       
      !$OMP PARALLEL NUM_THREADS(12)
      
      !$OMP do private(z1,sum1,datt)
            
	DO 10 I=190,NZ-190
	Z1=DZ*FLOAT(I-1)
	CALL SFM1D(NN,DEF,DZ,NZ,Z1,SUM1,RH,DIM)
      CALL convol(nn,Z1,DZ,NZ,RH,cstif,DATT,dim)
       FF(I)=SUM1+DATT
      if(nn.eq.2)then
      FF(I)=FF(I)+miuchain(i)
      endif       
10    CONTINUE
      
      !$OMP END do
      
          
      
      !$OMP do
      do j=190,NZ-190
      df1(J)=(RH(NN,J+1)-RH(NN,J-1))/dz/2.0
      ENDDO
c      df1(I,NZ)=df1(I,NZ-1)
c      df1(I,1)=df1(I,2)*2.-df1(I,3)
      !$OMP END do

      !$OMP do
          do j=190,NZ-190
	dff(j)=(ff(j+1)-ff(j-1))/dz/2.       
      enddo 
c          dff(i,NZ)=dff(i,NZ-1)
c         dff(i,1)=dff(i,2)*2.-dff(i,3)
       !$OMP END do
      
      !$OMP do
          do j=190,NZ-190            
	roudff(j)=df1(J)+dff(j)*RH(nn,j)  !*cos(JT*pi*fre)
    
          enddo  
           !$OMP END do
c          roudff(i,0)=roudff(i,2)*2.-roudff(i,1)
      
       
          
      !$OMP do
          do j=190,NZ-190             
	droudff(j)=(roudff(j+1)-roudff(j-1))/dz/2.-shear(nn)*RH(nn,j)
      enddo 
      !$OMP END do
c          droudff(i,1)=droudff(i,2)*2.-droudff(i,3)
c          droudff(i,nz)=droudff(i,nz-1)
      !$OMP do private(z1,ratio)
      do j=190,NZ-190
          
      Z1=DZ*FLOAT(I-1)
      
	CALL RATDENS(nn,DEFF,DZ,NZ,Z1,ratio,RH,dim)
      
c      if(abs(droudff(1,j)).ge.8.)then
c          droudff(2,j)=8.*droudff(1,j)/abs(droudff(1,j))    
c          endif
	rht(j)=RH(nn,j)+ratio*dt*(droudff(j))/DEF(nn)*deff
      IF(rht(j).LE.0.0) rht(j)=(rht(j-1)+rht(j+1))/2.
   
           
      enddo
      !$OMP END do
      !$OMP END PARALLEL
      

      
      
   
      
	RETURN
	END



      SUBROUTINE SFM1D0(J,SIGMA,DZ,NZ,Z,SUM1,DFN,Dim)   !!!!!!!!!!!
	IMPLICIT NONE
	INTEGER NZ,I,NP,J,NZ1,nn,Dim
	DOUBLE PRECISION PI,DZ,Z1,SUM1,SF21,SF31,SIGMA(Dim),SF3X,temp
	DOUBLE PRECISION SF21VI,Z,DF21,DF31,DF21VI,ZX,zd,zdd
	DOUBLE PRECISION DFN(6,NZ),DEFF,DF3X
	DOUBLE PRECISION k1,k2
	PARAMETER (PI=3.141592654)
	DOUBLE PRECISION SF01,SF11,SF11VI,df01,df11,df11vi




      IF(ABS(Z).LT.1.e-7)THEN
      Z=1.D-7
	ENDIF

      DEFF=SIGMA(J)
	   
	NP=INT(DEFF/2/DZ+1.d-10)


	SF21=0.0
	SF31=0.0
	SF21VI=0.0
	SF01=0.0
	SF11=0.0
	SF11VI=0.0


c	if(z.gt.64)then
c	SF21=0.0
c	endif

	K1=ABS(Z-DZ*DFLOAT(NP))
	K2=Z+DZ*DFLOAT(NP)+1.D-15

	IF(Z.LT.1.D-10)THEN
	K1=0.5*DEFF
	ENDIF



	DO 10 Z1=K1,K2,DZ
c	DO 10 I=1,2*NP+1
c	Z1=Z-NP*DZ+DZ*FLOAT(I-1)  

c	if(z.ge.k2)then
c	SF21=0.0
c	endif

      temp=ABS(Z1+1.D-20)

	NZ1=NINT(temp/DZ+1.D-5)
      

c	CALL FMN1D(SIGMA,DZ,NZ,Z1,DFN,RHO,Dim,sigf)

c	DF21=DFN(1)*DEFF                              !! from a wall 
c	DF31=DFN(2)*(0.25*DEFF**2-(Z-Z1)**2)                
c	DF21VI=-DFN(3)*(Z-Z1)    !!! + - sig difference             

c	df21=dfn(1,NZ1)*DEFF*z1/z
c	df31=dfn(2,NZ1)*(0.25*DEFF**2-(z-z1)**2)*z1/z            !! around a particle
c	df21vi=dfn(3,NZ1)*(z1**2-z**2+DEFF**2/4)/(z)/2.0   !!!!*z2/(z1**2.)
c	df01=dfn(4,NZ1)*DEFF*z1/z/PI/DEFF**2
c	df11=dfn(5,NZ1)*DEFF*z1/z/PI/DEFF/2.0            !! around a particle
c	df11vi=dfn(6,NZ1)*(z1**2-z**2+DEFF**2/4)/(z)/2.0/PI/DEFF/2.0   !!!!



	DF21=DFN(1,NZ1)*DEFF                                       !! from a wall 
	DF31=DFN(2,NZ1)*(0.25*DEFF**2-(Z-Z1)**2)
	DF21VI=-DFN(3,NZ1)*(Z-Z1)
	DF01=DFN(4,NZ1)*DEFF/PI/DEFF**2
	DF11=DFN(5,NZ1)*DEFF/2./PI/DEFF
	DF11VI=-DFN(6,NZ1)*(Z-Z1)/2./PI/DEFF


      IF(Z1.LE.1.D-20)THEN
	df31=0.D0
      df21vi=0.D0
      df21=0.D0
      df01=0.0
	df11=0.D0          !! around a particle
	df11vi=0.D0
	ENDIF

	IF(Z1.EQ.K1.OR.ABS(Z1-K2).LT.1.D-5)THEN
	DF21=DF21*0.5
	DF31=DF31*0.5
	DF21VI=DF21VI*0.5
      DF01=DF01*0.5
	DF11=DF11*0.5
	DF11VI=DF11VI*0.5
	ENDIF

c	IF(I.EQ.1.OR.I.EQ.2*NP+1)THEN
c	DF21=DF21*0.5
c	DF31=DF31*0.5
c	DF21VI=DF21VI*0.5
c     DF01=DF01*0.5
c	DF11=DF11*0.5
c	DF11VI=DF11VI*0.5
c	ENDIF


	SF21=SF21+DF21
	SF31=SF31+DF31
	SF21VI=SF21VI+DF21VI
      SF01=SF01+DF01
	SF11=SF11+DF11
	SF11VI=SF11VI+DF11VI

	
10	CONTINUE

c	SF21=PI*SF21*DZ                       ! from a wall 
c	SF31=PI*SF31*DZ
c	SF21VI=2.*PI*SF21VI*DZ

	SF21=PI*SF21*DZ                       ! around a particle
	SF31=PI*SF31*DZ
	SF21VI=2.*PI*SF21VI*DZ                !!!!!!!!!!!
	SF01=PI*SF01*DZ                       ! around a particle
	SF11=PI*SF11*DZ
	SF11VI=2.*PI*SF11VI*DZ            !!!!!!!!!!!

c      IF(ABS(Z-DEFF*0.49).LT.1.D-5)THEN
c	SF3X=0.D0
c	ENDIF


c      SF3X=0.0D0
c      IF((Z-DEFF/2).LE.0.D0)THEN
c      ZX=DEFF/2-Z+1.D-10

c      DO ZD=0.D0,ZX,DZ
c          NZ1=NINT(ZD/DZ) 
c	ZDD=ZD
c	IF(ZD.LT.1.D-6)THEN
c	ZDD=0.D0
c	ENDIF
c	DF3X=4.*PI*ZDD**2*dfn(2,NZ1)*DZ
c	IF(ABS(ZD-ZX).LT.1.D-5)THEN
c	DF3X=DF3X*0.5
c	ENDIF
c      SF3X=SF3X+DF3X
c      ENDDO
c      ELSE 
c      SF3X=0.0D0
c	ENDIF

c	SUM1=(SF21+SF31+SF21VI)
      SUM1=(SF21+SF31+SF21VI+SF01+SF11+SF11VI)  !+SF3X

	RETURN
	END	

c calculate the derivatives of the free energy density wrt weighted densities
	SUBROUTINE FMN1D0(SIGMA,DZ,NZ,DFN,RHO,Dim)  !!!!!!!!!!!!
	IMPLICIT NONE
	INTEGER NZ,I,NP,IZ2,Dim,J,ii,jj,n,nn,NZ1
	DOUBLE PRECISION PI,PI2,PI4,DZ,Z1,Z2,XI1,XI2,DFN3,temp
	DOUBLE PRECISION AN2,AN3,ANV2I,F21,F31,F21VI,DEFF,F
	DOUBLE PRECISION AN31,DN2F1,DN2F2,DN2F3,DN2F4,SIGMA(Dim)
	DOUBLE PRECISION DN3F1,DN3F2,DN3F3,DN3F4,DN3F5
	DOUBLE PRECISION DNVF1,DNVF2,DNVF3,AN31S,ANV22
	PARAMETER (PI=3.141592654,PI2=2.*PI,PI4=4.*PI)
	DOUBLE PRECISION DFN(6,NZ),RHO(Dim,NZ),AN3X,zx,zd,zdd,af3x
	DOUBLE PRECISION SIGMAII,k1,k2,N3OUT(NZ),N2OUT(NZ)

	DOUBLE PRECISION AN0,AN1,ANV1I,f01,f11vi,f11



c      IF(ABS(Z1).LT.0.001)THEN
c      Z1=1.D-25
c	ENDIF

91    FORMAT(1X,2000g15.6)



      do n=1,NZ

	Z1=DZ*float(n)   !60.0+

c      IF(ABS(Z1).LT.1.e-10)THEN
c      Z1=1.D-25
c	ENDIF


	AN2=0.0
	AN3=0.0
	ANV2I=0.0
      AN0=0.0
	AN1=0.0
	ANV1I=0.0

	
      do ii=1,Dim

      SIGMAII=SIGMA(ii)
	NP=INT(SIGMAII/2/DZ+1.d-10)

	
	K1=abs(Z1-DZ*DFLOAT(NP))
	K2=Z1+DZ*DFLOAT(NP)+1.D-10
      
	DO Z2=K1,K2,DZ
c      I=INT((Z2-K1)/DZ+1.D-10)+1
c	DO I=1,2*NP+1
c	Z2=Z1-NP*DZ+DZ*FLOAT(I-1)

      temp=ABS(Z2+1.D-20)

	IZ2=INT(temp/DZ) 
c      IF(IZ2.LE.1) THEN
c	IZ2=1  !!!
c	ELSE
c	IZ2=INT(Z2/DZ+(1E-8))+1  !!!    
c	ENDIFIF(IZ2.LT.0) THEN
c      IF(IZ2.LT.0) THEN
c	IZ2=INT(Z2/DZ-(1E-8))+1
c	ELSE
c	IZ2=INT(Z2/DZ+(1E-8))+1
c	ENDIF
      IF (IZ2.Lt.1) THEN
	F=0.0  !RHO(ii,abs(IZ2)+1)
	ELSEIF (IZ2.LT.NZ) THEN
   	F=RHO(ii,IZ2)
	ELSE
	F=RHO(ii,2*NZ-IZ2)
	ENDIF

CC      WRITE(*,*) Z2, IZ2
	IF(Z2.EQ.K1.OR.ABS(Z2-K2).LT.1.D-6)THEN
	F=F*0.5
	ENDIF


c	IF(I.EQ.1.OR.I.EQ.2*NP+1)THEN
c	F=F*0.5
c	ENDIF



c      f21=f*SIGMAII*z2/z1
c	f31=f*((0.5*SIGMAII)**2.-(z1-z2)**2.)*z2/z1          !! around a particle
c	f21vi=f*(z1**2.-z2**2.+SIGMAII**2./4.0)*z2/(z1**2.)/2.0
c	f01=f21/PI/SIGMAII**2.      !! around a particle
c      f11=f21/PI/SIGMAII/2.0
c	f11vi=f21vi/PI/SIGMAII/2.0

      F21=F*SIGMAII                                     !! from a wall 
	F31=F*(0.25*SIGMAII**2-(Z1-Z2)**2)
	F21VI=F*(Z1-Z2)
	F01=F21/SIGMAII**2/PI
	F11=F21/2./SIGMAII/PI
      F11VI=F21VI/2./SIGMAII/PI



	AN2=AN2+F21
	AN3=AN3+F31
      ANV2I=ANV2I+F21VI
      
      AN0=AN0+F01
	AN1=AN1+F11
      ANV1I=ANV1I+F11VI

	
	enddo

c      AN3X=0.0D0
c      IF((Z1-SIGMAII/2).LE.0.D0)THEN
c      ZX=SIGMAII/2-Z1+1.D-15

c      DO ZD=0.D0,ZX,DZ
c      I=INT(ZD/DZ+1.D-10)+1
c	ZDD=ZD
c	IF(ZD.LT.1.D-6)THEN
c	ZDD=0.D0
c	ENDIF
c   	F21=RHO(ii,I)
c	IF(ABS(ZD-ZX).LT.1.D-5)THEN
c	F21=F21*0.5
c	ENDIF
c	AF3X=4.*ZDD**2*F21
c      AN3X=AN3X+AF3X
c	ENDDO

c      AN3=AN3+AN3X
c	ENDIF

      
       enddo

      

c	AN2=PI*DZ*AN2 ! 4-33b n2      !  from a wall
c	AN3=PI*DZ*AN3 ! 4-33a n3
c	ANV2I=2.*PI*DZ*ANV2I ! 4-33c nv2	

	AN2=PI*AN2*DZ ! 4-33b n2      ! around a particle
	AN3=PI*AN3*DZ ! 4-33a n3
	ANV2I=2.*PI*ANV2I*DZ ! 4-33c nv2	!!!!!!!!!!!

	AN0=PI*AN0*DZ ! 4-33b n2      ! around a particle
	AN1=PI*AN1*DZ ! 4-33a n3
	ANV1I=2.*PI*ANV1I*DZ ! 4-33c nv2	!!!!!!!!!!!

     
      if(AN3.ge.1./sqrt(2.)) AN3=1./sqrt(2.)
c      AN3X=0.0D0


c      IF((Z1-SIGMAII/2.0-1.d-6).LE.0.D0)THEN

c      ZX=SIGMAII/2-Z1+1.D-15

c      DO ZD=0.D0,ZX,DZ
c      I=INT(ZD/DZ+1.D-20)+1
c	ZDD=ZD
c	IF(ZD.LT.1.D-10)THEN
c	ZDD=0.D0
c	ENDIF
c c  	F=RHO(i)
c	IF(ABS(ZD-ZX).LT.1.D-8)THEN
c	F=F*0.5d0
c	ENDIF
c	AF3X=4.*ZDD**2*F
c      AN3X=AN3X+AF3X
c	ENDDO


c      at=0.0
c	bt=SIGMAII/2-Z1+1.D-15
c c     do i=1,nn     
c      ZDD=ksi(i)*(bt-at)/2.0+(bt+at)/2.0
c	NZ1=INT(ZDD/DZ)

c	AF3X=4.*PI*ZDD**2*RHO(NZ1)*dksi(i)*(bt-at)/2.0
c	IF(ABS(ZDD-bt).LT.1.D-10)THEN
c	AF3X=AF3X*0.5
c	ENDIF
c      AN3X=AN3X+AF3X
c	enddo

c      AN3=AN3+AN3X
c	endif



    
	IF(AN3.LE.1.E-9)THEN
	DFN(1,n)=0.0d0  !  n2
	DFN(2,n)=0.0d0  !  n3
	DFN(3,n)=0.0d0  !  nv2
	DFN(4,n)=0.0d0  !  n0
	DFN(5,n)=0.0d0  !  n1
	DFN(6,n)=0.0d0  !  nv1

      else



c      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     rosenfeld Fai



c	ANV22=ANV2I**2

c      AN31=1.0-AN3
c      AN31S=AN31**2

c	DN2F2=AN1/AN31
c	DN2F3=(3.0*AN2**2-ANV22*3.0)/PI4/6.0/AN31S
c	DFN(1,n)=DN2F2+DN2F3    ! n2


c	DN3F1=AN0/AN31
c	DN3F2=(AN1*AN2-ANV1I*ANV2I)/AN31S
c	DN3F3=(AN2**3-3.0*AN2*ANV22)/PI4/3.0/AN31**3
c	DFN(2,n)=DN3F1+DN3F2+DN3F3 ! n3


c	DNVF1=-ANV1I/AN31 
c	DNVF2=-AN2*ANV2I/PI4/AN31**2
c	DFN(3,n)=DNVF1+DNVF2 ! nv2

c	DFN(4,n)=-DLOG(ABS(AN31)) ! n0
c	DFN(5,n)=AN2/AN31 ! n1
c	DFN(6,n)=-ANV2I/AN31 ! nv1

c      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      

	ANV22=ANV2I**2

     
      AN31=abs(1.0-AN3)
      AN31S=AN31**2

	DN2F2=AN1/AN31
	DN2F3=(3.0*AN2**2-ANV22*3.0)/PI4/9.0/AN3**2*DLOG(ABS(AN31))
	DN2F4=(3.0*AN2**2-ANV22*3.0)/PI4/9.0/AN3/AN31S
	DFN(1,n)=DN2F2+DN2F3+DN2F4    ! n2


	DN3F1=AN0/AN31
	DN3F2=(AN1*AN2-ANV1I*ANV2I)/AN31S
	DN3F3=(AN2**3-3.0*AN2*ANV22)/PI4/9.0*(-2.0*DLOG(ABS(AN31))/AN3**3)
	DN3F4=(AN2**3-3.0*AN2*ANV22)/PI4/9.0*(-1.0/AN31/AN3**2)
	DN3F5=(AN2**3-3.0*AN2*ANV22)/PI4/9.0*((3.0*AN3**2-AN3)/
     &	  (AN31*AN3)**3)
	DFN(2,n)=DN3F1+DN3F2+DN3F3+DN3F4+DN3F5 ! n3


	DNVF1=-ANV1I/AN31 
	DNVF2=-2.0*AN2*ANV2I*3.0/(PI4*9.0*AN3**2)*DLOG(ABS(AN31))
	DNVF3=-2.0*AN2*ANV2I*3.0/(PI4*9.0*AN3*AN31S)
	DFN(3,n)=DNVF1+DNVF2+DNVF3 ! nv2

	DFN(4,n)=-DLOG(ABS(AN31)) ! n0
	DFN(5,n)=AN2/AN31 ! n1
	DFN(6,n)=-ANV2I/AN31
      ENDIF
	
	
      enddo

      

	RETURN
      END
      
      SUBROUTINE FMN1Ds0(SIGMA,DZ,NZ,DFN,RHO,Dim)  !!!!!!!!!!!!
	IMPLICIT NONE
	INTEGER NZ,I,NP,IZ2,Dim,J,ii,jj,n,nn,NZ1
	DOUBLE PRECISION PI,PI2,PI4,DZ,Z1,Z2,XI1,XI2,DFN3,temp
	DOUBLE PRECISION AN2,AN3,ANV2I,F21,F31,F21VI,DEFF,F
	DOUBLE PRECISION AN31,DN2F1,DN2F2,DN2F3,DN2F4,SIGMA(Dim)
	DOUBLE PRECISION DN3F1,DN3F2,DN3F3,DN3F4,DN3F5
	DOUBLE PRECISION DNVF1,DNVF2,DNVF3,AN31S,ANV22
	PARAMETER (PI=3.141592654,PI2=2.*PI,PI4=4.*PI)
	DOUBLE PRECISION DFN(6,NZ),RHO(Dim,NZ),AN3X,zx,zd,zdd,af3x
	DOUBLE PRECISION SIGMAII,k1,k2,N3OUT(NZ),N2OUT(NZ)

	DOUBLE PRECISION AN0,AN1,ANV1I,f01,f11vi,f11



c      IF(ABS(Z1).LT.0.001)THEN
c      Z1=1.D-25
c	ENDIF

91    FORMAT(1X,2000g15.6)



      do n=1,NZ

	Z1=DZ*float(n)   !60.0+

c      IF(ABS(Z1).LT.1.e-10)THEN
c      Z1=1.D-25
c	ENDIF


	AN2=0.0
	AN3=0.0
	ANV2I=0.0
      AN0=0.0
	AN1=0.0
	ANV1I=0.0

	
      do ii=2,Dim

      SIGMAII=SIGMA(ii)
	NP=INT(SIGMAII/2/DZ+1.d-10)

	
	K1=abs(Z1-DZ*DFLOAT(NP))
	K2=Z1+DZ*DFLOAT(NP)+1.D-10
      
	DO Z2=K1,K2,DZ
c      I=INT((Z2-K1)/DZ+1.D-10)+1
c	DO I=1,2*NP+1
c	Z2=Z1-NP*DZ+DZ*FLOAT(I-1)

      temp=ABS(Z2+1.D-20)

	IZ2=INT(temp/DZ) 
c      IF(IZ2.LE.1) THEN
c	IZ2=1  !!!
c	ELSE
c	IZ2=INT(Z2/DZ+(1E-8))+1  !!!    
c	ENDIFIF(IZ2.LT.0) THEN
c      IF(IZ2.LT.0) THEN
c	IZ2=INT(Z2/DZ-(1E-8))+1
c	ELSE
c	IZ2=INT(Z2/DZ+(1E-8))+1
c	ENDIF
      IF (IZ2.Lt.1) THEN
	F=0.0  !RHO(ii,abs(IZ2)+1)
	ELSEIF (IZ2.LT.NZ) THEN
   	F=RHO(ii,IZ2)
	ELSE
	F=RHO(ii,2*NZ-IZ2)
	ENDIF

CC      WRITE(*,*) Z2, IZ2
	IF(Z2.EQ.K1.OR.ABS(Z2-K2).LT.1.D-6)THEN
	F=F*0.5
	ENDIF


c	IF(I.EQ.1.OR.I.EQ.2*NP+1)THEN
c	F=F*0.5
c	ENDIF



c      f21=f*SIGMAII*z2/z1
c	f31=f*((0.5*SIGMAII)**2.-(z1-z2)**2.)*z2/z1          !! around a particle
c	f21vi=f*(z1**2.-z2**2.+SIGMAII**2./4.0)*z2/(z1**2.)/2.0
c	f01=f21/PI/SIGMAII**2.      !! around a particle
c      f11=f21/PI/SIGMAII/2.0
c	f11vi=f21vi/PI/SIGMAII/2.0

      F21=F*SIGMAII                                     !! from a wall 
	F31=F*(0.25*SIGMAII**2-(Z1-Z2)**2)
	F21VI=F*(Z1-Z2)
	F01=F21/SIGMAII**2/PI
	F11=F21/2./SIGMAII/PI
      F11VI=F21VI/2./SIGMAII/PI



	AN2=AN2+F21
	AN3=AN3+F31
      ANV2I=ANV2I+F21VI
      
      AN0=AN0+F01
	AN1=AN1+F11
      ANV1I=ANV1I+F11VI

	
	enddo

c      AN3X=0.0D0
c      IF((Z1-SIGMAII/2).LE.0.D0)THEN
c      ZX=SIGMAII/2-Z1+1.D-15

c      DO ZD=0.D0,ZX,DZ
c      I=INT(ZD/DZ+1.D-10)+1
c	ZDD=ZD
c	IF(ZD.LT.1.D-6)THEN
c	ZDD=0.D0
c	ENDIF
c   	F21=RHO(ii,I)
c	IF(ABS(ZD-ZX).LT.1.D-5)THEN
c	F21=F21*0.5
c	ENDIF
c	AF3X=4.*ZDD**2*F21
c      AN3X=AN3X+AF3X
c	ENDDO

c      AN3=AN3+AN3X
c	ENDIF

      
       enddo

      

c	AN2=PI*DZ*AN2 ! 4-33b n2      !  from a wall
c	AN3=PI*DZ*AN3 ! 4-33a n3
c	ANV2I=2.*PI*DZ*ANV2I ! 4-33c nv2	

	AN2=PI*AN2*DZ ! 4-33b n2      ! around a particle
	AN3=PI*AN3*DZ ! 4-33a n3
	ANV2I=2.*PI*ANV2I*DZ ! 4-33c nv2	!!!!!!!!!!!

	AN0=PI*AN0*DZ ! 4-33b n2      ! around a particle
	AN1=PI*AN1*DZ ! 4-33a n3
	ANV1I=2.*PI*ANV1I*DZ ! 4-33c nv2	!!!!!!!!!!!

      N3OUT(n)=AN3
      N2OUT(n)=AN2

c      AN3X=0.0D0


c      IF((Z1-SIGMAII/2.0-1.d-6).LE.0.D0)THEN

c      ZX=SIGMAII/2-Z1+1.D-15

c      DO ZD=0.D0,ZX,DZ
c      I=INT(ZD/DZ+1.D-20)+1
c	ZDD=ZD
c	IF(ZD.LT.1.D-10)THEN
c	ZDD=0.D0
c	ENDIF
c c  	F=RHO(i)
c	IF(ABS(ZD-ZX).LT.1.D-8)THEN
c	F=F*0.5d0
c	ENDIF
c	AF3X=4.*ZDD**2*F
c      AN3X=AN3X+AF3X
c	ENDDO


c      at=0.0
c	bt=SIGMAII/2-Z1+1.D-15
c c     do i=1,nn     
c      ZDD=ksi(i)*(bt-at)/2.0+(bt+at)/2.0
c	NZ1=INT(ZDD/DZ)

c	AF3X=4.*PI*ZDD**2*RHO(NZ1)*dksi(i)*(bt-at)/2.0
c	IF(ABS(ZDD-bt).LT.1.D-10)THEN
c	AF3X=AF3X*0.5
c	ENDIF
c      AN3X=AN3X+AF3X
c	enddo

c      AN3=AN3+AN3X
c	endif



    
	IF(AN3.LE.1.E-9)THEN
	DFN(1,n)=0.0d0  !  n2
	DFN(2,n)=0.0d0  !  n3
	DFN(3,n)=0.0d0  !  nv2
	DFN(4,n)=0.0d0  !  n0
	DFN(5,n)=0.0d0  !  n1
	DFN(6,n)=0.0d0  !  nv1

      else



c      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     rosenfeld Fai



c	ANV22=ANV2I**2

c      AN31=1.0-AN3
c      AN31S=AN31**2

c	DN2F2=AN1/AN31
c	DN2F3=(3.0*AN2**2-ANV22*3.0)/PI4/6.0/AN31S
c	DFN(1,n)=DN2F2+DN2F3    ! n2


c	DN3F1=AN0/AN31
c	DN3F2=(AN1*AN2-ANV1I*ANV2I)/AN31S
c	DN3F3=(AN2**3-3.0*AN2*ANV22)/PI4/3.0/AN31**3
c	DFN(2,n)=DN3F1+DN3F2+DN3F3 ! n3


c	DNVF1=-ANV1I/AN31 
c	DNVF2=-AN2*ANV2I/PI4/AN31**2
c	DFN(3,n)=DNVF1+DNVF2 ! nv2

c	DFN(4,n)=-DLOG(ABS(AN31)) ! n0
c	DFN(5,n)=AN2/AN31 ! n1
c	DFN(6,n)=-ANV2I/AN31 ! nv1

c      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      

	ANV22=ANV2I**2

     
      AN31=abs(1.0-AN3)
      AN31S=AN31**2

	DN2F2=AN1/AN31
	DN2F3=(3.0*AN2**2-ANV22*3.0)/PI4/9.0/AN3**2*DLOG(ABS(AN31))
	DN2F4=(3.0*AN2**2-ANV22*3.0)/PI4/9.0/AN3/AN31S
	DFN(1,n)=DN2F2+DN2F3+DN2F4    ! n2


	DN3F1=AN0/AN31
	DN3F2=(AN1*AN2-ANV1I*ANV2I)/AN31S
	DN3F3=(AN2**3-3.0*AN2*ANV22)/PI4/9.0*(-2.0*DLOG(ABS(AN31))/AN3**3)
	DN3F4=(AN2**3-3.0*AN2*ANV22)/PI4/9.0*(-1.0/AN31/AN3**2)
	DN3F5=(AN2**3-3.0*AN2*ANV22)/PI4/9.0*((3.0*AN3**2-AN3)/
     &	  (AN31*AN3)**3)
	DFN(2,n)=DN3F1+DN3F2+DN3F3+DN3F4+DN3F5 ! n3


	DNVF1=-ANV1I/AN31 
	DNVF2=-2.0*AN2*ANV2I*3.0/(PI4*9.0*AN3**2)*DLOG(ABS(AN31))
	DNVF3=-2.0*AN2*ANV2I*3.0/(PI4*9.0*AN3*AN31S)
	DFN(3,n)=DNVF1+DNVF2+DNVF3 ! nv2

	DFN(4,n)=-DLOG(ABS(AN31)) ! n0
	DFN(5,n)=AN2/AN31 ! n1
	DFN(6,n)=-ANV2I/AN31
      ENDIF
	
	
      enddo

      
      
      
c      stop 
      


	RETURN
      END
      
      
      SUBROUTINE SFM1D(nn,DEFF,DZ,NZ,Z,SUM1,RHO,dim)
	IMPLICIT NONE
	INTEGER NZ,I,NP,iz1,dim,nn
	DOUBLE PRECISION PI,DZ,Z1,SUM1,SF21,SF31,DX,K1,K2
	DOUBLE PRECISION SF21VI,Z,DF21,DF31,DF21VI,ZD,ZX,DF3X,SF3X
	DOUBLE PRECISION DFN(3),DFA(3),RHO(dim,NZ),DEFF(dim),SUM2,ZDD
	PARAMETER (PI=3.141592654)
      COMMON /SS/ SUM2 
c      DZ=0.01D0*DEFF
c      DZ=DX
      IF(ABS(Z).LT.0.001)THEN
      Z=1.D-7
	ENDIF

	NP=INT(DEFF(nn)/2/DZ+1.d-10)

	SF21=0.0
	SF31=0.0
	SF21VI=0.0

      K1=Z-DZ*DFLOAT(NP)
      IF((Z-DEFF(nn)/2).LT.0.D0)THEN
	K1=ABS(Z-DZ*DFLOAT(NP))
	ENDIF
	K2=Z+DZ*DFLOAT(NP)+1.D-10

	IF(Z.LT.1.D-4)THEN
	K1=0.5D0
	ENDIF

C	DO 10 I=1,2*NP+1
C	Z1=Z-NP*DZ+DZ*FLOAT(I-1)

	DO 10 I=1,2*NP+1   !Z1=K1,K2,DZ
      Z1=Z-NP*DZ+DZ*FLOAT(I-1)   !I=INT((Z1-K1)/DZ+1.D-10)+1
          
	CALL FMN1D(nn,DEFF,DZ,NZ,Z1,DFN,RHO,dim)

	
      
      DF21=DFN(1)*DEFF(nn)                                       !! from a wall 
	DF31=DFN(2)*(0.25*DEFF(nn)**2-(Z-Z1)**2)
	DF21VI=-DFN(3)*(Z-Z1)

      IF(Z1.LT.1.D-6)THEN
	df31=0.D0
      df21vi=0.D0
      df21=0.D0
	ENDIF


	IF(I.EQ.1.OR.I.EQ.2*NP+1)THEN

C      IF((Z-DEFF/2).GT.0.D0)THEN

C	IF(Z1.EQ.K1.OR.I.EQ.2*NP+1)THEN
c	IF(Z1.EQ.K1.OR.ABS(Z1-K2).LT.1.D-6)THEN
	DF21=DF21*0.5
	DF31=DF31*0.5
	DF21VI=DF21VI*0.5
	ENDIF
C	ENDIF

	SF21=SF21+DF21
	SF31=SF31+DF31
	SF21VI=SF21VI+DF21VI
10	CONTINUE

	SF21=PI*SF21*DZ
	SF31=PI*SF31*DZ
	SF21VI=2.*PI*SF21VI*DZ

c      IF(ABS(Z-0.49D0).LT.1.D-5)THEN
c	SF3X=0.D0
c	ENDIF


c      SF3X=0.0D0
c      IF((Z-DEFF/2).LE.0.D0)THEN
c      ZX=DEFF/2-Z+1.D-15

c      DO ZD=0.D0,ZX,DZ
c      CALL FMN1D(DEFF,DZ,NZ,ZD,DFN,RHO)
c	ZDD=ZD
c	IF(ZD.LT.1.D-6)THEN
c	ZDD=0.D0
c	ENDIF
c	DF3X=4.*PI*ZDD**2*dfn(2)*DZ
c	IF(ABS(ZD-ZX).LT.1.D-5)THEN
c	DF3X=DF3X*0.5
c	ENDIF
c      SF3X=SF3X+DF3X
c	ENDDO
c	ENDIF

c	SUM2=(SF21+SF31+SF21VI)+SF3X
      SUM1=(SF21+SF31+SF21VI)
	RETURN
	END

	


! calculate the derivatives of the free energy density wrt weighted densities
	SUBROUTINE FMN1D(nn,DEFF,DZ,NZ,Z1,DFN,RHO,dim)
	IMPLICIT NONE
	INTEGER NZ,I,NP,IZ2,X,nn,ii,dim,ni
	DOUBLE PRECISION PI,PI2,PI4,DZ,Z1,Z2,XI1,XI2,DFN3,K1,K2
	DOUBLE PRECISION AN2,AN3,ANV2I,F21,F31,F21VI,DEFF(dim)
      DOUBLE PRECISION AN0,AN1,ANV1I
	DOUBLE PRECISION AN31,DN2F1,DN2F2,DN2F3,DN2F4,DX,ZDD
	DOUBLE PRECISION DN3F1,DN3F2,DN3F3,DN3F4,DN3F5,ZD,ZX,AF3X,AN3X
	DOUBLE PRECISION DNVF1,DNVF2,DNVF3,AN31S,ANV22,f
	PARAMETER (PI=3.141592654,PI2=2.*PI,PI4=4.*PI)
	DOUBLE PRECISION DFN(3),RHO(dim,NZ)

      IF(ABS(Z1).LT.0.001)THEN
      Z1=1.D-25
	ENDIF

C      DZ=DX
      

	AN2=0.0
	AN3=0.0
	ANV2I=0.0
      AN0=0.0
      AN1=0.0
      ANV1I=0.0
      
      if(nn.eq.1) then
          ni=2
      else
          ni=1
      endif
      
      do ii=ni,dim

	NP=INT(DEFF(ii)/2/DZ+1.d-10)

      K1=Z1-DZ*DFLOAT(NP)
      IF((Z1-DEFF(ii)/2).LE.0.D0)THEN
	K1=ABS(Z1-DZ*DFLOAT(NP))
	ENDIF
	K2=Z1+DZ*DFLOAT(NP)+1.D-10

	DO I=1,2*NP+1
	Z2=Z1-NP*DZ+DZ*FLOAT(I-1)

c	DO Z2=K1,K2,DZ
c      I=INT((Z2-K1)/DZ+1.D-10)+1

	IZ2=INT(Z2/DZ)+1
      IF(IZ2.LT.0) THEN
	IZ2=INT(Z2/DZ-(1E-8)) +1
	ELSE
	IZ2=INT(Z2/DZ+(1E-8)) +1
	ENDIF
      IF (IZ2.LT.1) THEN
	F=0.0
	ELSEIF (IZ2.LT.NZ) THEN
   	F=RHO(ii,IZ2)
	ELSE
	F=RHO(ii,NZ*2-IZ2)
	ENDIF

	IF(I.EQ.1.OR.I.EQ.2*NP+1)THEN

C      IF((Z1-DEFF/2).GT.0.D0)THEN

C	IF(Z2.EQ.K1.OR.I.EQ.2*NP+1)THEN
c	IF(Z2.EQ.K1.OR.ABS(Z2-K2).LT.1.D-6)THEN
	F=F*0.5
	ENDIF
C	ENDIF

c	f31=f21*((0.5*deff)**2-(z1-z2)**2)*z2/z1
c      if(z2.lt.z1)then
c	f21vi=f21*(z1**2-z2**2+deff**2/4)*z2/(z1**2)/2
c	else
c	f21vi=f21*(z1**2-z2**2+deff**2/4)*z2/(z1**2)/2
c	endif
c      f21=f21*deff*z2/z1
      
      F21=F*deff(ii)                                     !! from a wall 
	F31=F*(0.25*deff(ii)**2-(Z1-Z2)**2)
	F21VI=F*(Z1-Z2)

	AN2=AN2+F21
	AN3=AN3+F31
      ANV2I=ANV2I+F21VI 
      AN0=AN0+F21/pi/deff(ii)**2.
      AN1=AN1+F21/2./pi/deff(ii)
      ANV1I=ANV1I+F21VI/2./pi/deff(ii)
      
      ENDDO
      
      enddo
	AN2=PI*DZ*AN2
	AN3=PI*DZ*AN3
	ANV2I=2.*PI*DZ*ANV2I
	AN0=PI*DZ*AN0
      AN1=PI*DZ*AN1
      ANV1I=2.*PI*DZ*ANV1I
c      AN3X=0.0D0
c      IF((Z1-DEFF/2).LE.0.D0)THEN
c      ZX=DEFF/2-Z1+1.D-15

c      DO ZD=0.D0,ZX,DZ
c      I=INT(ZD/DZ+1.D-10)+1
c	ZDD=ZD
c	IF(ZD.LT.1.D-6)THEN
c	ZDD=0.D0
c	ENDIF
c   	F21=RHO(I)
c	IF(ABS(ZD-ZX).LT.1.D-5)THEN
c	F21=F21*0.5
c	ENDIF
c	AF3X=4.*PI*ZDD**2*F21*DZ
c      AN3X=AN3X+AF3X
c	ENDDO

c      AN3=AN3+AN3X
c	ENDIF

c      if(AN3.ge.1./sqrt(2.)) AN3=1./sqrt(2.)

	IF(AN3.LE.1.E-9)THEN
	DFN(1)=0.0D0
	DFN(2)=0.0D0
	DFN(3)=0.0D0
	RETURN
	ENDIF
      
	ANV22=ANV2I**2
	XI1=1.-ANV22/AN2**2
	XI2=1.-3.0*ANV22/AN2**2

      AN31=1.0-AN3
      AN31S=AN31**2
	DN2F1=-DLOG(ABS(AN31))/PI/DEFF(nn)**2
	DN2F2=AN2/2./PI/DEFF(nn)/AN31+AN1/AN31
	DN2F3=(3.0*AN2**2-ANV22*3.0)/PI4/9.0/AN3**2*DLOG(ABS(AN31))
	DN2F4=(3.0*AN2**2-ANV22*3.0)/PI4/9.0/AN3/AN31S
	DFN(1)=DN2F1+DN2F2+DN2F3+DN2F4

	DN3F1=AN0/AN31
	DN3F2=(AN1*AN2-ANV1I*ANV2I)/AN31S
	DN3F3=-(AN2**3-3.*AN2*ANV22)/(PI4*9.0*AN3**2*AN31)
	DN3F4=-2.0*(AN2**3-3.*AN2*ANV22)/(PI4*9.0*AN3**3)*DLOG(ABS(AN31))    
	DN3F5=-AN2**3*XI2*(1.0-3.0*AN3)/(PI4*9.0*AN31**3*AN3**2)     !
	DFN(2)=DN3F1+DN3F2+DN3F3+DN3F4+DN3F5

	DNVF1=-ANV2I/2./PI/DEFF(nn)/AN31-ANV1I/AN31
	DNVF2=-2.0*AN2*ANV2I*3.0/(PI4*9.0*AN3**2)*DLOG(ABS(AN31))
	DNVF3=-2.0*AN2*ANV2I*3.0/(PI4*9.0*AN3*AN31S)
	DFN(3)=DNVF1+DNVF2+DNVF3


c	DN2F1=-DLOG(ABS(AN31))/PI/DEFF**2
c	DN2F2=AN2/PI/DEFF/AN31
c	DN2F3=(3.0*AN2**2-ANV22*3.0)/PI4/6.0/AN31S
c	DFN(1)=DN2F1+DN2F2+DN2F3


c	DN3F1=AN2/PI/DEFF**2/AN31
c	DN3F2=AN2**2*XI1/PI/2.0/DEFF/AN31S
c	DN3F3=AN2**3*XI2/(PI4*3.0*AN31*AN31S)
c	DFN(2)=DN3F1+DN3F2+DN3F3

c	DNVF1=-ANV2I/PI/DEFF/AN31
c	DNVF2=-AN2*ANV2I/(PI4*AN31S)
c	DFN(3)=DNVF1+DNVF2


	RETURN
      END

      
      SUBROUTINE RATDENS(nn,DEFF,DZ,NZ,Z1,ratio,RHO,dim)
	IMPLICIT NONE
	INTEGER NZ,I,NP,IZ2,X,nn,ii,dim,ni
	DOUBLE PRECISION PI,PI2,PI4,DZ,Z1,Z2,XI1,XI2,DFN3,K1,K2
	DOUBLE PRECISION AN2,AN3,ANV2I,F21,F31,F21VI,DEFF(dim)
      DOUBLE PRECISION AN0,AN1,ANV1I
	DOUBLE PRECISION AN31,DN2F1,DN2F2,DN2F3,DN2F4,DX,ZDD
	DOUBLE PRECISION DN3F1,DN3F2,DN3F3,DN3F4,DN3F5,ZD,ZX,AF3X,AN3X
	DOUBLE PRECISION DNVF1,DNVF2,DNVF3,AN31S,ANV22,f
	PARAMETER (PI=3.141592654,PI2=2.*PI,PI4=4.*PI)
	DOUBLE PRECISION DFN(3),RHO(dim,NZ),ratio

      IF(ABS(Z1).LT.0.001)THEN
      Z1=1.D-25
	ENDIF

C      DZ=DX
      

	AN3=0.0
	      
      if(nn.eq.1) then
          ni=2
      else
          ni=1
      endif
      
      do ii=ni,dim

	NP=INT(DEFF(ii)/2/DZ+1.d-10)

      K1=Z1-DZ*DFLOAT(NP)
      IF((Z1-DEFF(ii)/2).LE.0.D0)THEN
	K1=ABS(Z1-DZ*DFLOAT(NP))
	ENDIF
	K2=Z1+DZ*DFLOAT(NP)+1.D-10

	DO I=1,2*NP+1
	Z2=Z1-NP*DZ+DZ*FLOAT(I-1)

c	DO Z2=K1,K2,DZ
c      I=INT((Z2-K1)/DZ+1.D-10)+1

	IZ2=INT(Z2/DZ)+1
      IF(IZ2.LT.0) THEN
	IZ2=INT(Z2/DZ-(1E-8)) +1
	ELSE
	IZ2=INT(Z2/DZ+(1E-8)) +1
	ENDIF
      IF (IZ2.LT.1) THEN
	F=0.0
	ELSEIF (IZ2.LT.NZ) THEN
   	F=RHO(ii,IZ2)
	ELSE
	F=RHO(ii,NZ*2-IZ2)
	ENDIF

	IF(I.EQ.1.OR.I.EQ.2*NP+1)THEN

C      IF((Z1-DEFF/2).GT.0.D0)THEN

C	IF(Z2.EQ.K1.OR.I.EQ.2*NP+1)THEN
c	IF(Z2.EQ.K1.OR.ABS(Z2-K2).LT.1.D-6)THEN
	F=F*0.5
	ENDIF
C	ENDIF

c	f31=f21*((0.5*deff)**2-(z1-z2)**2)*z2/z1
c      if(z2.lt.z1)then
c	f21vi=f21*(z1**2-z2**2+deff**2/4)*z2/(z1**2)/2
c	else
c	f21vi=f21*(z1**2-z2**2+deff**2/4)*z2/(z1**2)/2
c	endif
c      f21=f21*deff*z2/z1
      
                                      !! from a wall 
	F31=F*(0.25*deff(ii)**2-(Z1-Z2)**2)
	
	AN3=AN3+F31
      
      
      ENDDO
      
      enddo
	
	AN3=PI*DZ*AN3
      
c      IF(AN3.GE.0.99) AN3=0.99
      ratio=(1.-AN3)**3./(1.+2.*AN3+1.492*AN3*(1.-AN3)**3.)

	RETURN
	END

     

! Calculate the excess local chemical potential at point Z due to two-Yukawa attraction
	SUBROUTINE convol(nn,Z,DZ,NZ,f2,f1,DATT,dim)
	IMPLICIT NONE
	INTEGER NZ,K1,K2,i,j,dim,i2,i1,A,B,k,nn,ii
      real(kind=8) PAI,DZ,DATT,Z,Z1
      real(kind=8) RCI,ZZ1,temp,zz2,dr
      real(kind=8) f1(dim,dim,NZ)
      real(kind=8) f2(dim,NZ)
	PARAMETER (PAI=3.141592654)

      A=int((Z-10.D0)/DZ)
	IF(Z-10.D0.LE.0.D0)THEN
      A=1
	ENDIF
	B=int((Z+10.D0)/DZ)
	IF(Z+10.D0.GE.DZ*float(NZ))THEN
      B=int((Z+10.D0)/DZ)
	ENDIF
	DATT=0.0
      do ii=1,dim
	DO 10 K1=A,B
	Z1=FLOAT(K1-1)*DZ
	ZZ1=ABS(Z-Z1)
      k=Nint((zz1)/dz)
	DATT=DATT+f2(II,K1)*f1(nn,ii,k)
10    CONTINUE
      enddo
      DATT=DATT*DZ  !*2.*PAI
	RETURN
	END

	


C Calcuulate the integration of rc(r) from |z-z'| to infinite		
	FUNCTION RCI(ZZ1)
      IMPLICIT NONE	
	INTEGER J	
	DOUBLE PRECISION RCI,ZZ1,zz2,ZY(2),X(2,7),TEMP1(2)     
     &                ,cgt10,cgt11,cgt12,cgt13,Rc,XX	
      COMMON /C2DYUK/ ZY,X,cgt10,cgt11,cgt12,cgt13,Rc	
	DOUBLE PRECISION SIGMA,RHOB,TF,DEFF,H,Z0,ZH,TY(2),S  
	DOUBLE PRECISION RCC,RS,A,B,DPF,DEFFN	 	
	COMMON /HARDWALL/ Z0,ZH,TF,DEFF,DEFFN,RHOB,S,SIGMA

	IF (ZZ1.GT.Rc) THEN	   	
	  RCI=0.0	
	ELSEIF (ZZ1.GE.SIGMA) 	THEN 	

        DO J=1,2
        TY(J)=-X(J,1)/ZY(J)*DEFF**2
     &	*(EXP(-ZY(J)*(Rc/DEFF-1.D0))-EXP(-ZY(J)*(ZZ1/DEFF-1.D0)))
    
        ENDDO
        RCI=TY(1)-TY(2)                            !  ¡Òrc(r)dr, zz1 to rc

	ELSEIF (ZZ1.GE.DEFF) THEN 		
	  ZZ2=SIGMA	

        DO J=1,2
        TY(J)=-X(J,1)/ZY(J)*DEFF**2
     &	*(EXP(-ZY(J)*(Rc/DEFF-1.D0))-EXP(-ZY(J)*(ZZ2/DEFF-1.D0)))
          
        ENDDO
        RCI=TY(1)-TY(2)	                        !  ¡Òrc(r)dr, ¦Ò to rc
	  	  	
	ELSE	
	
	  DO J=1,2		
	  TEMP1(J)=-(X(J,1)+X(J,2))*(1.-EXP(ZY(J)*(1.-ZZ1/DEFF)))/ZY(J)     
     &          +X(J,3)*(1.-EXP(ZY(J)*(ZZ1/DEFF-1.0)))/ZY(J)
     &          +X(J,4)*(1.-(ZZ1/DEFF)**5)/5.0
     &          +X(J,5)*(1.-(ZZ1/DEFF)**3)/3.0            ! ¡Òrc(r)dr, zz1 to d
     &          +X(J,6)*(1.-(ZZ1/DEFF)**2)/2.0
     &          +X(J,7)*(1.-(ZZ1/DEFF))
     	  ENDDO		
        ZZ2=SIGMA		

        DO J=1,2
        TY(J)=-X(J,1)/ZY(J)*DEFF**2
    
     &	*(EXP(-ZY(J)*(Rc/DEFF-1.D0))-EXP(-ZY(J)*(ZZ2/DEFF-1.D0)))
        
        ENDDO
        RCI=TY(1)-TY(2)	 

	  RCI=(TEMP1(1)-TEMP1(2))*DEFF**2+RCI
	ENDIF	
	END



C calculate the coefficients of the direct correlation funciton 
	



      SUBROUTINE EXTP(Z,PHI)
      IMPLICIT NONE
      DOUBLE PRECISION Z0,ZH,Z,Z1,DELTA,SIG,RHOB,SIGMA,
     &         TF,SIGW,EPSW,TEMP1,TEMP2,PHI,DEFF,S
      DOUBLE PRECISION RS,RC,DPF,A,B,EPS,Y,RCI,DEFFN
      Parameter (DELTA=0.8044,SIGW=0.903,epsw=12.96)
	COMMON /HARDWALL/ Z0,ZH,TF,DEFF,DEFFN,RHOB,S,SIGMA
	COMMON /EPSPHI/EPS

      
c	DPF=1.D-10
	
c      RS=(26.D0/7.D0)**(1.D0/6.D0)*SIGMA
c      RC=(67.D0/48.D0)*RS
c	A=-(24192.D0/3211.D0)/RS**2
c	B=-(387072.D0/61009.D0)/RS**3


c	IF (Z.LE.DPF) THEN
c	    PHI=1.D10
c      ELSEIF (Z.LE.(DPF+RS)) THEN
c	    PHI=4.0D0*EPS*((SIGMA/(Z-DPF))**12-(SIGMA/(Z-DPF))**6)
c      ELSEIF (Z.LE.(DPF+RC)) THEN
c	    PHI=A*(Z-DPF-RC)**2*EPS+B*(Z-DPF-RC)**3*EPS
c	ELSE
c	    PHI=0.0D0
c	ENDIF

c	PHI=PHI/TF !  1/TF=¦Â*¦Å 
      Y=Z  !zh/2.-Z+1.e-8
      
      IF(y.lE.3.*deff) THEN
C          PHI=1.E8
C          ELSE
     	PHI=4./Tf*((deff/(Y))**12.
     & -(deff/(Y))**6.-(1./3.)**12.+(1./3.)**6.)	
	ELSE
	PHI=0.0
      ENDIF

      END


      subroutine KKFFT(PR,PI,N,K,FR,FI,L,IL)
	dimension PR(N),PI(N),FR(N),FI(N)
	real(kind=8) PR,PI,FR,FI,P,Q,S,VR,VI,PODDR,PODDI
	do 20 IT=0,N-1
	  M=IT
	  IS=0
	  do 10 I=0,K-1
	    J=M/2
	    IS=2*IS+(M-2*J)
	    M=J
10	  continue
	  FR(IT+1)=PR(IS+1)
	  FI(IT+1)=PI(IS+1)
20	continue
	PR(1)=1.0
	PI(1)=0.0
	PR(2)=COS(6.283185306/N)
	PI(2)=-SIN(6.283185306/N)
	if (L.NE.0) PI(2)=-PI(2)
	do 30 I=3,N
	  P=PR(I-1)*PR(2)
	  Q=PI(I-1)*PI(2)
	  S=(PR(I-1)+PI(I-1))*(PR(2)+PI(2))
	  PR(I)=P-Q
	  PI(I)=S-P-Q
30	continue
	do 40 IT=0,N-2,2
	  VR=FR(IT+1)
	  VI=FI(IT+1)
	  FR(IT+1)=VR+FR(IT+2)
	  FI(IT+1)=VI+FI(IT+2)
	  FR(IT+2)=VR-FR(IT+2)
	  FI(IT+2)=VI-FI(IT+2)
40	continue
	M=N/2
	NV=2
	do 70 L0=K-2,0,-1
	  M=M/2
	  NV=2*NV
	  do 60 IT=0,(M-1)*NV,NV
	  do 60 J=0,(NV/2)-1
	    P=PR(M*J+1)*FR(IT+J+1+NV/2)
	    Q=PI(M*J+1)*FI(IT+J+1+NV/2)
	    S=PR(M*J+1)+PI(M*J+1)
	    S=S*(FR(IT+J+1+NV/2)+FI(IT+J+1+NV/2))
	    PODDR=P-Q
	    PODDI=S-P-Q
	    FR(IT+J+1+NV/2)=FR(IT+J+1)-PODDR
	    FI(IT+J+1+NV/2)=FI(IT+J+1)-PODDI
	    FR(IT+J+1)=FR(IT+J+1)+PODDR
	    FI(IT+J+1)=FI(IT+J+1)+PODDI
60	  continue
70	continue
	if (L.NE.0) then
	  do 80 I=1,N
	    FR(I)=FR(I)/N
	    FI(I)=FI(I)/N
80	  continue
	end if
	if (IL.NE.0) then
	  do 90 I=1,N
	    PR(I)=SQRT(FR(I)*FR(I)+FI(I)*FI(I))
	    PI(I)=ATAN(FI(I)/FR(I))*360.0/6.283185306
90	  continue
	end if
	return
      end  

   
 
      
   
      SUBROUTINE SFM1Dchain(DZ,NZ,RO,sigf,miuchain)  !!!!!!!!!!!!
	IMPLICIT NONE
	INTEGER NZ,I,NP,IZ2,Dim,J,ii,jj,NN,n,Dim21
	DOUBLE PRECISION PI,DZ,Z1,Z2
	DOUBLE PRECISION F,h,t,yab1,yab2,g_SEXP,roz
	PARAMETER (PI=3.141592654)
	DOUBLE PRECISION k1,k2,xtemp
	DOUBLE PRECISION RO(2,NZ),sigf(2),SIGMAII
      DOUBLE PRECISION Nab(NZ)
      DOUBLE PRECISION Tab(NZ),pTpNab(NZ)
      DOUBLE PRECISION ksi2(NZ),ksi3(NZ),Yab(NZ)
      DOUBLE PRECISION pYp2ab(NZ),pYp3ab(NZ)
      DOUBLE PRECISION Zab(NZ),pZp2ab(NZ)
      DOUBLE PRECISION pZp3ab(NZ),miuchain21(NZ)
      DOUBLE PRECISION Temp1(NZ),Temp2(NZ),Temp3(NZ)
      DOUBLE PRECISION Tempa(NZ),Tempb(NZ),Temp4(NZ)
      DOUBLE PRECISION temp,miuchain(NZ),sigD,y,Temp5(NZ)



      NN=Dim21
      
91    FORMAT(1X,200g15.6)
	
      do n=1,NZ

	Z1=float(n-1)*DZ

      Nab(n)=0.0
      SIGMAII=sigf(2)
	NP=INT(SIGMAII/2/DZ+1.d-10)

	
	DO I=1,2*NP+1
	Z2=Z1-NP*DZ+DZ*FLOAT(I-1)
	
      IZ2=INT(Z2/DZ)+1
      IF(IZ2.LT.0) THEN
	IZ2=INT(Z2/DZ-(1E-8)) +1
	ELSE
	IZ2=INT(Z2/DZ+(1E-8)) +1
	ENDIF
      IF (IZ2.LT.1) THEN
	F=0.0
	ELSEIF (IZ2.LT.NZ) THEN

   	F=RO(2,IZ2)
	ELSE

	F=RO(2,2*NZ-IZ2)
	ENDIF

	IF(Z2.EQ.K1.OR.ABS(Z2-K2).LT.1.D-6)THEN
	F=F*0.5
	ENDIF     
      Nab(n)=Nab(n)+f   !*Z2/Z1
      enddo
      Nab(n)=Nab(n)*DZ/2.0/SIGMAII+1.0e-20


c      write(*,*) i,j,n,log(Nab(i,j,n))

	Tab(n)=log((Nab(n)))
	pTpNab(n)=1.0/abs(Nab(n))
      enddo
	



      do n=1,NZ

	Z1=float(n-1)*DZ

      ksi2(n)=0.0   
      ksi3(n)=0.0   	  	
	do ii=1,2
      SIGMAII=sigf(ii)
	NP=INT(SIGMAII/2/DZ+1.d-10)

	
	DO I=1,2*NP+1
	Z2=Z1-NP*DZ+DZ*FLOAT(I-1)
	
      IZ2=INT(Z2/DZ)+1
      IF(IZ2.LT.0) THEN
	IZ2=INT(Z2/DZ-(1E-8)) +1
	ELSE
	IZ2=INT(Z2/DZ+(1E-8)) +1
	ENDIF
	
      IF (IZ2.LT.1) THEN
	F=0.0
	ELSEIF (IZ2.LT.NZ) THEN
   	F=RO(ii,IZ2)
	ELSE
	F=RO(ii,2*NZ-IZ2)
	ENDIF
	IF(Z2.EQ.K1.OR.ABS(Z2-K2).LT.1.D-6)THEN
	F=F*0.5
      ENDIF  
      
      ksi2(n)=ksi2(n)+f*(SIGMAII**2-(Z2-Z1)**2)*PI*SIGMAII**2/
     &	    8.0/SIGMAII**3  !*Z2/Z1
      
      ksi3(n)=ksi3(n)+f*(SIGMAII**2-(Z2-Z1)**2)*PI*SIGMAII**3/
     &	    8.0/SIGMAII**3   !*Z2/Z1
                    
      enddo
      
      ksi2(n)=ksi2(n)*DZ
      ksi3(n)=ksi3(n)*DZ
      enddo
      enddo

      h=1.0e-5
      
      do n=1,NZ
          
          Yab(n)=1.0/(1.0-ksi3(n))+3.0*sigf(2)*sigf(2)/
     &	(sigf(2)+sigf(2))*ksi2(n)/(1.0-ksi3(n))**2+2.0*(sigf(2)*
     &	sigf(2)/(sigf(2)+sigf(2)))**2*ksi2(n)**2/(1.0-ksi3(n))**3

      pYp2ab(n)=3.0*sigf(2)*sigf(2)/
     &	(sigf(2)+sigf(2))/(1.0-ksi3(n))**2+4.0*(sigf(2)*
     &	sigf(2)/(sigf(2)+sigf(2)))**2*ksi2(n)/(1.0-ksi3(n))**3

      pYp3ab(n)=1.0/(1.0-ksi3(n))**2+6.0*sigf(2)*sigf(2)/
     &	(sigf(2)+sigf(2))*ksi2(n)/(1.0-ksi3(n))**3+6.0*(sigf(2)*
     &	sigf(2)/(sigf(2)+sigf(2)))**2*ksi2(n)**2/(1.0-ksi3(n))**4

      Zab(n)=log(abs(Yab(n)))
      pZp2ab(n)=1.0/Yab(n)*pYp2ab(n)
      pZp3ab(n)=1.0/Yab(n)*pYp3ab(n)

     
	      
      enddo     
     
	
      do n=1,NZ
      Tempa(n)=RO(2,n)*pZp2ab(n)
      Tempb(n)=RO(2,n)*pZp3ab(n)

      enddo
     

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
      
     
      do n=1,NZ
	Z1=float(n-1)*DZ

      Temp1(n)=0.0
      SIGMAII=sigf(2)
	NP=INT(SIGMAII/2/DZ+1.d-10)

	
	DO I=1,2*NP+1
	Z2=Z1-NP*DZ+DZ*FLOAT(I-1)
	
      IZ2=INT(Z2/DZ)+1
      IF(IZ2.LT.0) THEN
	IZ2=INT(Z2/DZ-(1E-8)) +1
	ELSE
	IZ2=INT(Z2/DZ+(1E-8)) +1
	ENDIF

      IF (IZ2.LT.1) THEN
	F=0.0
	ELSEIF (IZ2.LE.NZ) THEN
   	F=RO(2,IZ2)*pTpNab(IZ2)
	ELSE
	F=RO(2,2*NZ-IZ2)*pTpNab(2*NZ-IZ2)
	ENDIF
	IF(Z2.EQ.K1.OR.ABS(Z2-K2).LT.1.D-6)THEN
	F=F*0.5
	ENDIF     
      Temp1(n)=Temp1(n)+f !*Z2/Z1
      enddo
      Temp1(n)=Temp1(n)*DZ/2.0/SIGMAII
      
      Temp3(n)=0.0
      SIGMAII=sigf(2)
	NP=INT(SIGMAII/2/DZ+1.d-10)

	
	DO I=1,2*NP+1
	Z2=Z1-NP*DZ+DZ*FLOAT(I-1)
	
      IZ2=INT(Z2/DZ)+1
      IF(IZ2.LT.0) THEN
	IZ2=INT(Z2/DZ-(1E-8)) +1
	ELSE
	IZ2=INT(Z2/DZ+(1E-8)) +1
	ENDIF

      IF (IZ2.LT.1) THEN
	F=0.0
	ELSEIF (IZ2.LT.NZ) THEN
   	F=Tempa(IZ2)
	ELSE
	F=Tempa(2*NZ-IZ2)
	ENDIF
	IF(Z2.EQ.K1.OR.ABS(Z2-K2).LT.1.D-6)THEN
	F=F*0.5
	ENDIF  
      Temp3(n)=Temp3(n)+f*(SIGMAII**2-(Z2-Z1)**2)*PI*SIGMAII**2/
     &	    8.0/SIGMAII**3 !*Z2/Z1

      enddo

      Temp3(n)=Temp3(n)*DZ

      Temp4(n)=0.0
      SIGMAII=sigf(2)
	NP=INT(SIGMAII/2/DZ+1.d-10)

	
	DO I=1,2*NP+1
	Z2=Z1-NP*DZ+DZ*FLOAT(I-1)
	
      IZ2=INT(Z2/DZ)+1
      IF(IZ2.LT.0) THEN
	IZ2=INT(Z2/DZ-(1E-8)) +1
	ELSE
	IZ2=INT(Z2/DZ+(1E-8)) +1
	ENDIF

      IF (IZ2.LT.1) THEN
	F=0.0
	ELSEIF (IZ2.LT.NZ) THEN
   	F=Tempb(IZ2)
	ELSE
	F=Tempb(2*NZ-IZ2)
	ENDIF
	IF(Z2.EQ.K1.OR.ABS(Z2-K2).LT.1.D-6)THEN
	F=F*0.5
	ENDIF  
      Temp4(n)=Temp4(n)+f*(SIGMAII**2-(Z2-Z1)**2)*PI*SIGMAII**3/
     &	    8.0/SIGMAII**3 !*Z2/Z1

      enddo
      Temp4(n)=Temp4(n)*DZ


      miuchain21(n)=( (1.0-Zab(n)-Tab(n))
     &	-Temp1(n)-Temp3(n)-Temp4(n) )  !
	      
      enddo


      do n=1,NZ
      miuchain(n)=miuchain21(n)   !*(nn-1.)/nn
   			      
      enddo




	RETURN
      END
      
      SUBROUTINE CPS(M,R,ROU,AMU)
      INTEGER M
	DOUBLE PRECISION  R(M),ROU(M),AMU(M),amu1(M)
	DOUBLE PRECISION  ATT(2),PI,PIA,A1,A2,A3,A4
	DOUBLE PRECISION  AKX0,AKX1,AKX2,AKX3,F2,DF3,PHS
	DOUBLE PRECISION  AKX31,AKX31S,AKX31C
	PARAMETER (PI=3.141592654)
	AKX0=0.0
	AKX1=0.0
	AKX2=0.0
	AKX3=0.0
	PIA=PI/6.
	DO 10 I=1,M
	AKX0=AKX0+ROU(I)
	AKX1=AKX1+ROU(I)*R(I)/2           ! R is diameter 
	AKX2=AKX2+ROU(I)*R(I)*R(I)*PI
	AKX3=AKX3+ROU(I)*R(I)**3*PIA
10	CONTINUE

	AKX31=1.0-AKX3
	AKX31S=AKX31*AKX31
	AKX31C=AKX31S*AKX31

	PHS=AKX0/AKX31+AKX1*AKX2/AKX31S+(3.-AKX3)*AKX2**3/AKX31C/36./PI

	F2=DLOG(AKX31)+AKX3/AKX31-0.5*(AKX3/AKX31)**2
	F3=2.*DLOG(AKX31)+AKX3*(1.+AKX31)/AKX31

	DO 20 I=1,M
	A1=-DLOG(AKX31)  +   AKX2/AKX31*R(I)/2.

	A2=AKX1/AKX31*PI*R(I)**2
     &	+3*AKX2**2*(DLOG(AKX31)*AKX31S+AKX3)/36./PI/AKX3**2/AKX31S
     &    *PI*R(I)**2

	A3=(AKX0/AKX31+AKX1*AKX2/AKX31S)*R(I)**3*PIA

	A4=(-AKX2**3/36./PI/AKX3**2/AKX31
     &    -2.*AKX2**3/36./PI/AKX3**3*DLOG(AKX31)
     &    -AKX2**3*(1.-3.*AKX3)/36./PI/AKX31**3/AKX3**2)*R(I)**3*PIA

	AMU(I)=A1+A2+A3+A4

20	CONTINUE
	RETURN
	END
