c     Author      : Sumanta Kundu
c     Date        : 19 May 2023
c     Reference   : Inequality of avalanche sizes in models of fracture.
c                   [https://doi.org/10.48550/arXiv.2303.10168]
c     Description : Avalanche dynamics of Fuse Model uing CG Method.
c                   "Exit criterion" - Checking of conductivity of the network.
c------------------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      parameter(L=128,vmin=1.0d-06,nconf=500)
      parameter(Lp1=L+1,Lm1=L-1,L2=L*L,LLpL=L2+L,L3=L2-L,LL2=L2*2)
      parameter(beta=1.00d00)
 
      dimension v(L,L),p(L2),r(L2),ap(L2),v1(L,L)
      dimension n(2,L3),threshold(2,L3),Ds(LL2),Dt(LL2)
      integer p1(L),m1(L),pbc(0:Lp1)

      character*30 filename

      common /b2/ n  
      common /b3/ p1,m1,pbc 

      call randinit(%val(82353669))

      Ds=0.0d00
      Dt=0.0d00
      count=0.0d00

      do i=1,L
      p1(i)=i+1
      m1(i)=i-1
      pbc(i)=i
      enddo
      pbc(Lp1)=1
      pbc(0)=L

      vtop=1.0d00
      do i=1,L
      v1(i,L)=vtop
      v1(i,1)=0.0d00
      do j=2,Lm1
      v1(i,j)=(j-1)*(vtop/Lm1)
      enddo
      enddo

      do 2000 iconf=1,nconf

      write(filename,108)"Avalb1.00L",L,"_",iconf,".dat"
      open(unit=3,file=trim(filename))
      write(3,*)"L, beta, iconf",L,beta,iconf
      write(3,*)" "
      write(3,*)"         Naval,      size"

      threshmin=1.0d20
      do i=1,L3
      do idir=1,2
      n(idir,i)=1
      R1=rann()*2.0d00-1.0d00
      threshold(idir,i)=10.0d00**(beta*R1)
      if(threshold(idir,i).lt.threshmin)then
      threshmin=threshold(idir,i)
      isite=i
      id=idir
      endif
      enddo
      enddo

      vtop=(L-1)*threshmin
      do i=1,L
      v(i,1)=0.0d00
      v(i,L)=vtop
      do j=2,Lm1
      v(i,j)=(j-1)*threshmin
      enddo
      enddo

      n(id,isite)=0
      ncutbond=1

C-----Avalanche begins here--------

      naval=0
 100  naval=naval+1

      itime=0
      isize=0
 200  itime=itime+1  
      isize=isize+ncutbond

      call volt_dist(v)

      ncutbond=0
      do j=1,Lm1
      do i=1,L
      k=(j-1)*L+i
      vdiff=dabs(v(i,p1(j))-v(i,j))
      if(n(2,k).eq.1.and.threshold(2,K).lt.vdiff)then
      n(2,k)=0
      ncutbond=ncutbond+1
      endif

      vdiff=dabs(v(pbc(m1(i)),p1(j))-v(i,j))
      if(n(1,k).eq.1.and.threshold(1,K).lt.vdiff)then
      n(1,k)=0
      ncutbond=ncutbond+1
      endif
      enddo
      enddo

      if(ncutbond.gt.0)go to 200

c-----Checking of Connectivity--------------------

      call volt_dist(v1)

      conduct=0.0d00
      do j=1,Lm1
      do i=1,L
      k=(j-1)*L+i
      conduct=conduct+n(1,k)*(v1(i,j)-v1(pbc(m1(i)),p1(j)))**2
     &               +n(2,k)*(v1(i,j)-v1(i,p1(j)))**2
      enddo
      enddo

      if(conduct.lt.(0.1d00/L)) go to 1000

c-----Starting of a new avalanche---------------------

      count=count+1.0d00
      Ds(isize)=Ds(isize)+1.0d00
      Dt(itime)=Dt(itime)+1.0d00

      write(3,*)naval,isize

      isite=0
      factormin=1.0d10
      do j=1,Lm1
      do i=1,L
      k=(j-1)*L+i
      vdiff=dabs(v(i,p1(j))-v(i,j))
c     if(n(2,k).eq.1.and.vdiff.gt.vmin)then
      if(n(2,k).eq.1)then
      factor=threshold(2,k)/vdiff
      if(factor.lt.factormin)then
      factormin=factor
      isite=k
      idir=2
      endif
      endif

      vdiff=dabs(v(pbc(m1(i)),p1(j))-v(i,j))
c     if(n(1,k).eq.1.and.vdiff.gt.vmin)then
      if(n(1,k).eq.1)then
      factor=threshold(1,k)/vdiff
      if(factor.lt.factormin)then
      factormin=factor
      isite=k
      idir=1
      endif
      endif
      enddo
      enddo

      if(isite.gt.0)then
      vtop=vtop*factormin
      do j=2,L
      do i=1,L
      v(i,j)=v(i,j)*factormin
      enddo
      enddo
      n(idir,isite)=0
      ncutbond=1
      go to 100
      endif

 1000 continue

      if(iconf/10*10.eq.iconf)then
      open(1,file='aval_size_dist_L128.dat')
      open(2,file='aval_time_dist_L128.dat')
      do i=1,LL2
      if(Ds(i).gt.0.0d00)write(1,*)i,Ds(i)/count
      if(Dt(i).gt.0.0d00)write(2,*)i,Dt(i)/count
      enddo
      close(1)
      close(2)
      endif

      close(unit=3)
 2000 continue
 
 108  format(a10,I3.3,a1,I3.3,a4)
      end




       subroutine volt_dist(v)
c------To know the Potential at each node by solving Kirchhoff eq.
c------Using Conjugate Gradient Method:
 
      implicit double precision (a-h,o-z)
      parameter(L=128,eps=1.0d-12)
      parameter(Lp1=L+1,Lm1=L-1,L2=L*L,LLpL=L2+L,L3=L2-L,LL2=L2*2)

      dimension v(L,L),p(L2),r(L2),ap(L2)
      dimension n(2,L3)
      integer p1(L),m1(L),pbc(0:Lp1)

      common /b2/ n  
      common /b3/ p1,m1,pbc 

c     do ij=1,llpl
c     v(ij)=float((ij-1)/l)/l
c     enddo


      do i=1,L
      do j=2,Lm1
      k=i+(j-1)*L
      p(k)=
     &   -n(1,k                 )*(v(pbc(m1(i)),p1(j))-v(i,j))
     &   -n(2,k                 )*(v(i,p1(j)         )-v(i,j))
     &   -n(1,pbc(p1(i))+(j-2)*L)*(v(pbc(p1(i)),m1(j))-v(i,j))
     &   -n(2,i+(j-2)*L         )*(v(i,m1(j)         )-v(i,j))
      r(k)=p(k)
      enddo
      enddo

      do i=1,L
      p(i   )=0.0d00
      p(i+L3)=0.0d00
      enddo

      do iter=1,LL2

      rTr=0.0d00
      do i=Lp1,L3
      rTr=rTr+r(i)*r(i)
      enddo
      if(rTr.le.eps) goto 5000

      do i=Lp1,L3
      i1=i-1+L
      if(i1/L*L.eq.i1)i1=i1+L
      i2=i+L
      i3=i-L+1
      if((i3-1)/L*L.eq.(i3-1))i3=i3-L
      i4=i-L
      pi=p(i)
      ap(i)=n(1,i )*(p(i1)-pi)+n(2,i )*(p(i2)-pi)
     &      +n(1,i3)*(p(i3)-pi)+n(2,i4)*(p(i4)-pi)
      enddo

      pTap=0.0d00
      do  i=Lp1,L3
      pTap=pTap+p(i)*ap(i) 
      enddo

      alpha=rTr/pTap

      do i=1,L
      do j=2,Lm1
      k=i+(j-1)*L
      v(i,j)=v(i,j)+alpha*p(k)
      enddo
      enddo

      do i=Lp1,L3
      r(i)=r(i)-alpha*ap(i)
      enddo

      rnTrn=0.0d00

      do i=Lp1,L3
      rnTrn=rnTrn+r(i)*r(i)
      enddo

      bita=rnTrn/rTr

      do i=Lp1,L3
      p(i)=r(i)+bita*p(i)
      enddo
c End of iteration

      enddo

5000  continue

      return
      end 






*-----RANDOM NUMBER GENERATOR 1 ------------------------------------------

      double precision function ranf1(iran1)
      iran1=iran1*1566083941
      if(iran1.lt.0)iran1=iran1+2147483647+1
      iran1=iran1*1566083941
      if(iran1.lt.0)iran1=iran1+2147483647+1
      ranf1=iran1*4.6566128752458D-10
      return
      end

*-----RANDOM NUMBER GENERATOR 2 ------------------------------------------

      double precision function ranf2(iran2)
      iran2=iran2*1664525
      if(iran2.lt.0)iran2=iran2+2147483647+1
      iran2=iran2*1664525
      if(iran2.lt.0)iran2=iran2+2147483647+1
      ranf2=iran2*4.6566128752458D-10
      return
      end

*-----RANDOM NUMBER GENERATOR 3 ------------------------------------------

      double precision function ranf3(iran3)
      iran3=iran3*16807
      if(iran3.lt.0)iran3=iran3+2147483647+1
      iran3=iran3*16807
      if(iran3.lt.0)iran3=iran3+2147483647+1
      ranf3=iran3*4.6566128752458D-10
      return
      end

