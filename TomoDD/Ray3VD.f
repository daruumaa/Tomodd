c----This file is extracted from Simul2000 (Clifford Thurber and 
c    Donna Eberhart-Phillips)
c    The main purpose of this file is to read 3D starting velocity model,
c    and set up the grid necessary to conduct pseudo-bending ray tracing. For 
c    parameter changes, check file ray_common.inc. 
c    --------------------- modified by Haijiang Zhang, 
c                          hjzhang@geology.wisc.edu
c-----------------------------------------------------------------------------
cz 5/14/2004 Add iusep,iuseq,invdel,ifixl in the subroutine input3

      subroutine path(xe, ye, ze, xr, yr, zr, isp, dlta, ttime)
c  this routine determines the minimum-time ray path
c  in two steps:  first, an approximate path is
c  determined using approximate ray tracing
c  then the approximate path is used as a starting point in
c  shooting ray tracing routine to determine the path in 3-d
c  ***  note - current version does not contain full 3-d ray
c  ***  tracing - routines are under development
c
c  declaration statements:
c
c  common block variables:
      implicit none
      include 'ray_common.inc'
      integer jflag, ncrold, ndpold
      real xe, ye, ze, xr, yr, zr
      real dlta,ttime, fstime, tdif
      integer isp, nitpbu, no, jpb

c  determine 3-d path using circular raypaths.
c  determine approximate path
  120 jflag=jfl
      ncrold=nco
      ndpold=ndo
c
      call rayweb(isp,xe,ye,ze,xr,yr,zr,
     * fstime,jflag,ncrold,ndpold)
      ttime=fstime
      ttc=ttime
      nco=ncrold
      ndo=ndpold
c
c  do psuedo-bending if i3d>0
      if (i3d.eq.0) return
c
  125 continue
c  number of pb iter depends on distance
      nitpbu=nitpb(1)
      if(dlta.gt.delt1) nitpbu=nitpb(2)
      call minima(no,isp,ttime,nitpbu,jpb)
      if(jpb.lt.nitpbu) goto 130
cz        write(26,2601) no,jpb,nitpbu
 2601   format(' Minima: no=',i4,2x,a4,', used maximum ',
     2  'number PB iter.: j=',i3,', nitpb=',i3)
  130 continue
      ttc=ttime
c  write out travel-time differences to file 19
      tdif=fstime-ttime
cz      write(19,1900) ne,dlta,fstime,ttime,tdif
 1900 format(i4,1x,a4,f7.2,3f8.4)
c
      return
c***** end of subroutine path *****
      end
c

      subroutine rayweb(isp,xe,ye,ze,xr,yr,zr,
     *            fstime,jflag,ncrold,ndpold)
c  approximate ray tracing package art2
c    with fast ray tracing code
c     by Cliff Thurber (from his simul3l version)
c
c  common block variables:
      implicit none
      include 'ray_common.inc'
      integer isp,jflag,ncrold,ndpold
      integer ns1,ii,ndip1,ndip2,ncr0,ncr1,iz0,npt2,npt3      
      real sn1, xstep, ystep, zstep, z0

c
c  declaration statements:
c  parameters
      real xe,ye,ze,xr,yr,zr,fstime
c  local variables
      real delx,dely,delz,sep,delsep,pthsep(530),strpth(1590),
     * fstpth(1590),dipvec(3,9),disvec(1590,9),trpath(1590,9),
     * trtime(9),tmin,tt,trpth1(1590)
      integer nd,ns,npt,ncr,i,ic,ndp,nc,n1,n2,n3,nn,np,ndpfst
c
c  compute source-receiver separation
      delx=xr-xe
      dely=yr-ye
      delz=zr-ze
      sep=sqrt(delx*delx+dely*dely+delz*delz)
c  determine integer parameters for set of curves to be constructed
      call setup(sep,scale1,scale2,nd,ns,npt,ncr)
c
c  set up pthsep array for straight-line travel time calculation
      sn1=1.0/ns
      delsep=sep*sn1
      do 20 i=1,ns
         pthsep(i)=delsep
   20 continue
c
c  determine points along straight-line path
      xstep=delx*sn1
      ystep=dely*sn1
      zstep=delz*sn1
c
      ic=0
      ns1=ns+1
      do 25 ii=1,ns1
c
         i=ii-1
c
         ic=ic+1
         strpth(ic)=xe+xstep*i
         fstpth(ic)=strpth(ic)
         ic=ic+1
         strpth(ic)=ye+ystep*i
         fstpth(ic)=strpth(ic)
         ic=ic+1
         strpth(ic)=ze+zstep*i
         fstpth(ic)=strpth(ic)
   25 continue
c
c  compute travel time along straight-line path
      call ttimed(isp,ns,npt,strpth,pthsep,fstime)
c
      if (ncr.eq.1) go to 65
c
c  compute the dip vectors of length scale2
      call cmpdpv(xe,ye,ze,xr,yr,zr,scale2,ndip,dipvec)
c
c  compute the basic set of displacement vectors
      call cmpdsv(ndip,iskip,ns,dipvec,disvec)
c
c  set first and last points of all trial paths to source and receiver
      n1=3*npt-2
      n2=n1+1
      n3=n1+2
      ndip1=1+iskip
      ndip2=ndip-iskip
c
c  fast ray tracing code
c
      nz1=nz-1
      ncr0=1
      ncr1=ncr-1
      if (jflag.eq.0) go to 28
      ncr0=ncrold-1
      if (ncr0.lt.1) ncr0=1
      ncr1=ncrold+1
      if (ncr1.gt.ncr) ncr1=ncr
c  ndip was changed (in main) for ihomo iterations., make ndpold=vertical plane.
c     if(jfl.eq.1) ndpold=(ndip+1)/2     ! 21-feb-86, this is done in strt, so unnecessary here
      ndip1=ndpold-1
      if (ndip1.lt.1+iskip) ndip1=1+iskip
      ndip2=ndpold+1
      if (ndip2.gt.ndip-iskip) ndip2=ndip-iskip
   28 continue
c  set "old" values to straight line
      ncrold=0
      ndpold=(ndip+1)/2
c
c
      do 30 ndp=ndip1,ndip2
         trpath(1,ndp)=xe
         trpath(2,ndp)=ye
         trpath(3,ndp)=ze
         trpath(n1,ndp)=xr
         trpath(n2,ndp)=yr
         trpath(n3,ndp)=zr
   30 continue
      trpth1(1)=xe
      trpth1(2)=ye
      trpth1(3)=ze
      trpth1(n1)=xr
      trpth1(n2)=yr
      trpth1(n3)=zr
c
c  loop over the curve sets
      do 40 nc=ncr0,ncr1
         iz0=0
c
c  loop over different dips for one set of curves
         do 42 ndp=ndip1,ndip2
c
            npt2=npt-2
c  loop to determine points along one path
            do 44 np=1,npt2
               n1=3*np+1
               n3=n1+2
               do 43 nn=n1,n3
                  trpath(nn,ndp)=nc*disvec(nn,ndp)+strpth(nn)
                  trpth1(nn)=trpath(nn,ndp)
   43          continue
               if((i3d.lt.3).or.(ze.gt.zn(nz1))) goto 44
c  Use less curvature if below "moho"
               if(trpth1(n3).gt.zn(nz1)) then
                  if(iz0.eq.0) then
                     z0=trpth1(n3)-disvec(n3,ndp)
                     iz0=1
                  endif
                  trpath(n3,ndp)=disvec(n3,ndp)+z0
c                 write(6,600) trpth1(n3),zn(nz1),trpath(n3,ndp)
  600             format('trpth1(n3)',f7.2,',zn(nz1)',f7.2,
     2               'new trpth1(n3)',f7.2)
                  trpth1(n3)=trpath(n3,ndp)
               endif
   44       continue
c
c  set up pthsep array for travel time calculations
            if (ndp.eq.ndip1) call cmpsep(trpth1,pthsep,ns)
c  compute travel time along one path
            call ttimed(isp,ns,npt,trpth1,pthsep,tt)
            trtime(ndp)=tt
   42    continue
c
c  sort through trtime to find fastest path from current set
         tmin=1.0e15
         do 50 ndp=ndip1,ndip2
            if (trtime(ndp).gt.tmin) go to 50
            tmin=trtime(ndp)
            ndpfst=ndp
   50    continue
c
c  compare fastest trtime to current value of fstime
c  replace fstime and fstpth if needed
         if (tmin.ge.fstime) go to 40
         fstime=tmin
c  reset "old" values
         ncrold=nc
         ndpold=ndpfst
c
         npt3=3*npt
         do 52 np=1,npt3
            fstpth(np)=trpath(np,ndpfst)
   52    continue
c
   40 continue
c
c  put fstpth into rp array
   65 continue
      do 60 np=1,npt
         n3=3*np
         n1=n3-2
         n2=n1+1
         rp(1,np)=fstpth(n1)
         rp(2,np)=fstpth(n2)
         rp(3,np)=fstpth(n3)
   60 continue
      nrp=npt
c
c***** end of subroutine rayweb *****
      return
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -

      subroutine setup(sep,scale1,scale2,nd,ns,npt,ncr)
      implicit none
c
c  parameters
      real sep,scale1,scale2
c
      integer nd,ns,npt,ncr
c
c  determine the number of path divisions - use scale1
      nd=1+nint(3.32193*log10(sep/scale1))
      if (nd.gt.7) nd=7
c  number of segments along path
      ns=2**nd
c  number of points on path
      npt=ns+1
c
c  determine the number of curves - use scale2
      ncr=1+0.5*sep/scale2
c
      if (sep.gt.scale1) return
      nd=0
      ns=1
      npt=2
      ncr=1
c
c***** end of subroutine setup *****
      return
      end


      subroutine ttimed(isp,ns,npt,pathr,pthsep,tt)

      implicit none
      integer isp, np1
c
c  travel time along path via trapezoidal rule integration
c
c  parameters
      real pathr(1590),pthsep(530),tt
c
      integer ns,npt
c  local variables
      real vpt(530),x,y,z,v
c
      integer ip,np
c
c  loop over points along path to determine velocities
      ip=0
      do 10 np=1,npt
         ip=ip+1
         x=pathr(ip)
         ip=ip+1
         y=pathr(ip)
         ip=ip+1
         z=pathr(ip)
         call vel3(isp,x,y,z,v)
   10 vpt(np)=v
c
      tt=0.0
c  sum up travel time - use trapezoidal rule
      do 20 np=1,ns
         np1=np+1
c  Check for value outside defined area
         if((vpt(np).le.0.0).or.(vpt(np1).le.0.0)) goto 99
         tt=tt+pthsep(np)/(vpt(np)+vpt(np1))
   20 continue
      tt=2.0*tt
c
      return
   99 write(16,100) np,vpt(np),np1,vpt(np1)

  100 format(' **** ERROR IN TTIME ****, velocity le 0.0',
     2 /,'    np   vpt(np)   np1   vpt(np1)',/,1x,2(i5,f10.3))
      stop
c***** end of subroutine ttime *****
      end


      subroutine vel3(isp,x,y,z,v)
      implicit none
c  This routine is Cliff Thurber's
c  common block variables:
      include 'ray_common.inc'
      real wv(8)
      integer ip, jp, kp, kpg
      common/weight/ wv,ip,jp,kp,kpg

      integer isp,kp1,ip1,jp1
      real x,y,z,v,xf,yf,zf,xf1,yf1,zf1
c
c  use Prothero's intmap here
      call intmap(x,y,z,ip,jp,kp)
c
      ip1=ip+1
      jp1=jp+1
      kp1=kp+1
c
c	write(16,100)x,xl,ip,y,yl,jp,z,zl,kp
c100	format(3(2f7.3,i3))
      xf=(x-xn(ip))/(xn(ip1)-xn(ip))
      yf=(y-yn(jp))/(yn(jp1)-yn(jp))
      zf=(z-zn(kp))/(zn(kp1)-zn(kp))
      xf1=1.0-xf
      yf1=1.0-yf
      zf1=1.0-zf
c
      wv(1)=xf1*yf1*zf1
      wv(2)=xf*yf1*zf1
      wv(3)=xf1*yf*zf1
      wv(4)=xf*yf*zf1
      wv(5)=xf1*yf1*zf
      wv(6)=xf*yf1*zf
      wv(7)=xf1*yf*zf
      wv(8)=xf*yf*zf
c  calculate velocity
c  S-velocity is stored after P-velocity
c  (or V*Q if iuseq=1)
      kpg=kp
      if(isp.eq.1) kp=kp+nz
      kp1=kp+1
      v=wv(1)*vel(ip,jp,kp)+wv(2)*vel(ip1,jp,kp)
     2 +wv(3)*vel(ip,jp1,kp)+wv(4)*vel(ip1,jp1,kp)
     * +wv(5)*vel(ip,jp,kp1)+wv(6)*vel(ip1,jp,kp1)
     * +wv(7)*vel(ip,jp1,kp1)+wv(8)*vel(ip1,jp1,kp1)
      return
c***** end of subroutine vel3 *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine intmap_old(x,y,z,ip,jp,kp)
      implicit none
c  Modified by W. Prothero so a single call can get the indices
c  common block variables:
      include 'ray_common.inc'
      integer ip, jp, kp
      real x,y,z
c
      ip=int((x+xl)/bld)
      ip=ixloc(ip)
      jp=int((yl+y)/bld)
      jp=iyloc(jp)
      kp=int((z+zl)/bld)
      kp=izloc(kp)
c  If an array element=0, the position is off the map.
      return
c***** end of subroutine intmap *****
      end

      subroutine intmap(x,y,z,ip,jp,kp)
      implicit none
c  Modified by W. Prothero so a single call can get the indices
c  common block variables:
      include 'ray_common.inc'
      integer ip, jp, kp
      real x,y,z
c change nint to int
      ip=int((x+xl)/bld)
      ip=ixloc(ip)
      jp=int((yl+y)/bld)
      jp=iyloc(jp)
      kp=int((z+zl)/bld)
      kp=izloc(kp)
c  If an array element=0, the position is off the map.
      return
c***** end of subroutine intmap *****
      end

      subroutine cmpdpv(xe,ye,ze,xr,yr,zr,scale2,ndip,dipvec)
      implicit none
c
c  parameters
      integer n1, n2
      real xe,ye,ze,xr,yr,zr,scale2,dipvec(3,9)
c
      integer ndip
c  local variables
      real dx,dy,dz,xh1,yh1,zh1,xh2,yh2,zh2,size,xv,yv,zv,
     *     rescal,x451,y451,z451,x452,y452,z452
c
      integer nv
c
      dx=xr-xe
      dy=yr-ye
      dz=zr-ze
c
c  near-vertical vector
      xv=-dx*dz
      yv=-dy*dz
      zv=dx*dx+dy*dy
c  rescale vector to length scale2
      size=sqrt(xv*xv+yv*yv+zv*zv)
      rescal=scale2/size
c
      xv=xv*rescal
      yv=yv*rescal
      zv=zv*rescal
c
c  store this vector
      nv=(ndip+1)/2
      dipvec(1,nv)=xv
      dipvec(2,nv)=yv
      dipvec(3,nv)=zv
c
      if (ndip.eq.1) return
c
c  horizontal vectors
      xh1=dy
      yh1=-dx
      zh1=0.0
      xh2=-dy
      yh2=dx
      zh2=0.0
c  rescale the vectors to length scale2
      size=sqrt(xh1*xh1+yh1*yh1)
      rescal=scale2/size
c
      xh1=xh1*rescal
      yh1=yh1*rescal
      xh2=xh2*rescal
      yh2=yh2*rescal
c
c  store these two vectors
      dipvec(1,1)=xh1
      dipvec(2,1)=yh1
      dipvec(3,1)=zh1
c
      dipvec(1,ndip)=xh2
      dipvec(2,ndip)=yh2
      dipvec(3,ndip)=zh2
c
      if (ndip.eq.3) return
c
c  determine two 45 degree dip vectors
      rescal=0.7071068
c
      n1=(1+nv)/2
      n2=(nv+ndip)/2
c
      x451=(xh1+xv)*rescal
      y451=(yh1+yv)*rescal
      z451=(zh1+zv)*rescal
c
      x452=(xh2+xv)*rescal
      y452=(yh2+yv)*rescal
      z452=(zh2+zv)*rescal
c
      dipvec(1,n1)=x451
      dipvec(2,n1)=y451
      dipvec(3,n1)=z451
c
      dipvec(1,n2)=x452
      dipvec(2,n2)=y452
      dipvec(3,n2)=z452
c
      if (ndip.eq.5) return
c
c  determine four 22.5 degree dip vectors
      rescal=0.5411961
c
      dipvec(1,2)=(xh1+x451)*rescal
      dipvec(2,2)=(yh1+y451)*rescal
      dipvec(3,2)=(zh1+z451)*rescal
c
      dipvec(1,4)=(x451+xv)*rescal
      dipvec(2,4)=(y451+yv)*rescal
      dipvec(3,4)=(z451+zv)*rescal
c
      dipvec(1,6)=(xv+x452)*rescal
      dipvec(2,6)=(yv+y452)*rescal
      dipvec(3,6)=(zv+z452)*rescal
c
      dipvec(1,8)=(x452+xh2)*rescal
      dipvec(2,8)=(y452+yh2)*rescal
      dipvec(3,8)=(z452+zh2)*rescal
c
c***** end of subroutine cmpdpv *****
      return
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine cmpdsv(ndip,iskip,ns,dipvec,disvec)
      implicit none
c
c  parameters
      integer ndip1, ndip2
      real dipvec(3,9),disvec(1590,9)
c
      integer ndip,iskip,ns
c  local variables
      real darc(129)
c
      integer inc,narc,ndp,np,n,nd
c  coefficients of standard arc
c
       data darc/0.,0.0342346,0.0676860,0.1003707,0.1323047,0.1635029,
     *.1939796,.2237485,.2528226,.2812141,.3089350,.3359963,
     *.3624091,.3881833,.4133289,.4378553,.4617713,.4850857,
     *.5078065,.5299417,.5514988,.5724850,.5929070,.6127717,
     *.6320853,.6508538,.6690831,.6867787,.7039459,.7205898,
     *.7367154,.7523272,.7674298,.7820274,.7961241,.8097238,
     *.8228301,.8354468,.8475771,.8592244,.8703916,.8810817,
     *.8912975,.9010416,.9103164,.9191245,.9274679,.9353487,
     *.9427691,.9497307,.9562353,.9622845,.9678797,.9730224,
     *.9777138,.9819550,.9857470,.9890908,.9919872,.9944367,.9964401,
     *.9979979,.9991102,.9997776,1.0000000,.9997776,.9991102,.9979979,
     *.9964401,.9944367,.9919872,.9890908,.9857470,.9819550,.9777138,
     *.9730224,.9678797,.9622845,.9562353,.9497307,.9427691,.9353487,
     *.9274679,.9191245,.9103164,.9010416,.8912975,.8810817,.8703916,
     *.8592244,.8475771,.8354468,.8228301,.8097238,.7961241,.7820274,
     *.7674298,.7523272,.7367154,.7205898,.7039459,.6867787,.6690831,
     *.6508538,.6320853,.6127717,.5929070,.5724850,.5514988,.5299417,
     *.5078065,.4850857,.4617713,.4378553,.4133289,.3881833,.3624091,
     *.3359963,.3089350,.2812141,.2528226,.2237485,.1939796,.1635029,
     *.1323047,.1003707,.0676860,.0342346,0.0/
      inc=128/ns
        ndip1=1+iskip
      ndip2=ndip-iskip
c
c  loop over dips
      do 30 ndp=ndip1,ndip2
         narc=1
         nd=3
c  loop over points on the path (skip first and last)
         do 20 np=2,ns
            narc=narc+inc
c
c  loop over x,y,z
            do 10 n=1,3
               nd=nd+1
               disvec(nd,ndp)=darc(narc)*dipvec(n,ndp)
c
   10       continue
   20    continue
   30 continue
c
c***** end of subroutine cmpdsv *****
      return
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c

      subroutine ttmder(isp,ttime, dth1, stepl)
      implicit none
c
c  declaration statements:
      include 'ray_common.inc'
      real wv(8)
      integer ip,jp,kp,kpg
      common/weight/ wv,ip,jp,kp,kpg
      integer in(8)     
      real xe, ye,ze,v,us,dx,dy,dz,ds,uds,tt,half
      real rx,ry,rz,sl,fnsegi,ssl,dxs,dys,dzs,xp,yp,zp
      real dt
      integer i,nrp1,i1,nseg,is,izg,kk,kk1,jj,jj1,ii,ii1,ijk,ini,ne,no
      integer j, inp
c
c  common block variables:
      real dth1(4), ttime, stepl
      integer isp

c
c
c   This subroutine has been modified for Q inversion 13-March-1998.
c   Following Andreas Rietbrock's version.
c   Then (when iuseq=1), the second part of the model array
c   is used to store Vp*Qp (not Vp/Vs).
c   The solution is only for Q, not for Vp though.  
c   Thus the solution arrays that would otherwise be related
c   to the first part of the model array (Vp) are used.
c   Therefore in ttmder, vel3 is called with "1" to retrieve
c   the V*Q value; and vel3 is called with "0" to get the 
c   indices for the solution arrays.
c
c
c  resegment ray path
c  calculate travel time derivatives with respect to hypocentral
c  parameters - use coords of first two points on raypath
c  determine slowness at source
      xe=rp(1,1)
      ye=rp(2,1)
      ze=rp(3,1)              
      call vel3(isp,xe,ye,ze,v)       
      us=1.0/v

c**
c  determine cartesian derivatives from direction cosines
      dx=rp(1,2)-rp(1,1)
      dy=rp(2,2)-rp(2,1)
      dz=rp(3,2)-rp(3,1)
      ds=sqrt(dx*dx+dy*dy+dz*dz)
      uds=-us/ds
c  hypocentral derivatives
      dth1(2)=uds*dx
      dth1(3)=uds*dy
      dth1(4)=uds*dz
c  origin time derivative
      dth1(1)=1.0

c  skip next section if all velocity nodes are fixed (nparvi=0)
      if(nparvi.eq.0) goto 88
c  skip next section if only doing earthquake location
c     if(kopt.eq.0) go to 88
c
c  travel time and velocity partial derivatives
      tt=0.0
      half=0.5
cz  initialize the velocity derivative
      do i=1, npari
         dtm(i)=0.0
      enddo
c  loop over segments comprising the ray path
c     write(*,*)'Number of segments and stepl:',nrp, stepl
      nrp1=nrp-1
      pl=0.0
      do 50 i=1,nrp1
         i1=i+1
         rx=rp(1,i)
         ry=rp(2,i)
         rz=rp(3,i)
         dx=rp(1,i1)-rx
         dy=rp(2,i1)-ry
         dz=rp(3,i1)-rz
c  compute segment length
         sl=sqrt(dx*dx+dy*dy+dz*dz)
         pl=pl+sl
c  decide on number of subsegments and compute length
         nseg=nint(sl/stepl)+1
         fnsegi=1.0/float(nseg)
         ssl=sl*fnsegi
         dxs=dx*fnsegi
         dys=dy*fnsegi
         dzs=dz*fnsegi
c
         xp=rx-half*dxs
         yp=ry-half*dys
         zp=rz-half*dzs
c  loop over subsegments
         do 55 is=1,nseg
            xp=xp+dxs
            yp=yp+dys
            zp=zp+dzs
c
          call vel3(isp, xp, yp, zp, v)
          dt=ssl/v
c
c
c  The next section is a change from 'block' to 'linear'
c   partial derivatives, by C.Thurber,may10,1983.
c  Nodes with non-zero weight
            in(1)=ip-1+nx2*(jp-2)+nxy2*(kp-2)-nxy2*(2*isp)
            in(2)=in(1)+1
            in(3)=in(1)+nx2
            in(4)=in(3)+1
            in(5)=in(1)+nxy2
            in(6)=in(5)+1
            in(7)=in(5)+nx2
            in(8)=in(7)+1
c
c  Assign zero weight to boundary nodes (these nodes are not 
c  included in the inversion, but are in the velocity array,
c  thus we want to avoid writing to negative or incorrect 
c  elements of the partial derivative matrix)
c
            if(ip.eq.1) then
c              write(16,1610) xp,yp,zp,v,ip,jp,kpg
               wv(1)=0.0
               wv(3)=0.0
               wv(5)=0.0
               wv(7)=0.0
            else
               if(ip.eq.nx1) then
c                 write(16,1610) xp,yp,zp,v,ip,jp,kpg
                  wv(2)=0.0
                  wv(4)=0.0
                  wv(6)=0.0
                  wv(8)=0.0
               end if
            endif
c
            if(jp.eq.1) then
c              write(16,1610) xp,yp,zp,v,ip,jp,kpg
               wv(1)=0.0
               wv(2)=0.0
               wv(5)=0.0
               wv(6)=0.0
            else
               if(jp.eq.ny1) then
c                 write(16,1610) xp,yp,zp,v,ip,jp,kpg
                  wv(3)=0.0
                  wv(4)=0.0
                  wv(7)=0.0
                  wv(8)=0.0
               endif
            endif
c
            if((kpg.eq.1).or.(kpg.eq.(nz1+1))) then
c              write(16,1610) xp,yp,zp,v,ip,jp,kpg
               do 30 izg=1,4
                  wv(izg)=0.0
   30          continue
            else
               if((kpg.eq.nz1).or.(kpg.eq.(2*nz1))) then
c                 write(16,1610) xp,yp,zp,v,ip,jp,kpg
                  do 35 izg=5,8
                     wv(izg)=0.0
   35             continue
               endif
            endif
 1610       format(' ASSIGNING ZERO WEIGHTS IN TTMDER',
     2         'xp=',f7.2,',yp=',f7.2,',zp=',
     3         f7.2,',v=',f5.3,',ip=',i2,',jp=',i2,',kpg=',i2)
c
c  Accumulate model partial derivatives
            do 48 kk=1,2
               kk1=kk-1
               do 47 jj=1,2
                  jj1=jj-1
                  do 46 ii=1,2
                     ii1=ii-1
                     ijk=ii+2*jj1+4*kk1
c skip boundary nodes
                     if(wv(ijk).lt.0.05) goto 46
c DEP
c write out DWS for all nodes (including fixed and linked) when nitmax=1
c (useful for planning fixed and linked)
c Include weight factor like for inversion
c Note that for shots on nit=0, combined weight, including residual
c wt'g is not yet calculated so use wtsht.
cz--- Keep it unused temperarily ( by H. Zhang)
cz             if(nitmax.eq.1) then
cz               if((ne.gt.(neqs+nbls)).and.(nit.eq.0)) then
cz                 hitall(in(ijk))=hitall(in(ijk))+wv(ijk)*wtsht
cz               else
cz                 hitall(in(ijk))=hitall(in(ijk))+wv(ijk)*wtcomb(no,ne)
cz               endif
cz             endif
c  skip fixed nodes
                     if(nfix(in(ijk)).eq.1) goto 46
                     ini=ndexfx(in(ijk))
c
c  start cht 1998
      if (imerge(in(ijk)).eq.1) ini=ndexfx(jequal(in(ijk)))
c  end cht 1998
c
c check for writing to an array element that is outside of the inversion
c  solution array
                     if((ini.lt.1).or.(ini.gt.nparvi)) then
                        write(16,1606) ini,ijk
 1606                   format(' *** Error in TTMDER, accessing',
     2                     ' gridpoint outside of velocity inversion',
     3                     ' gridpoints, ini=',i5,', ijk=',i5,/,
     4                     22x,'Probably boundary gridpoints are',
     5                     ' too close (TTMDER tries to write DTM',
     6                     ' elements with wv >= 0.05)')
                        write(16,1603) ne,no,xp,yp,zp,v,ip,jp,kp,kpg
 1603                   format(' ne=',i5,', no=',i5,', xp=',f8.2,
     2                     ', yp=',f8.2,', zp=',f8.2,', v=',f8.3,/,
     3                     21x,'ip=',i6,',   jp=',i6,',   kp=',i6,
     4                     ',   kpg=',i6)
                        write(16,1607) (j,in(j),j,wv(j),j=1,8)
 1607                   format(' in(',i1,')=',i6,' wv(',i1,')=',e15.5)
                        write(16,1608)
 1608                   format(' * * * * STOP * * * * (to avoid',
     2                     ' writing outside of defined DTM array)')
                        stop
                     end if
cz                   inp=ini+(no-1)*npari
                     inp=ini
cDEP Now include weight factor for hit(DWS)
cz                     if((ne.gt.(neqs+nbls)).and.(nit.eq.0)) then
cz                       hit(in(ijk))=hit(in(ijk))+wv(ijk)*wtsht
cz                     else
cz                       hit(in(ijk))=hit(in(ijk))+wv(ijk)*wtcomb(no,ne)
cz                     endif
cz                     hit(in(ijk))=hit(in(ijk))+wv(ijk)*TimeWeight

c     set up the model derivative matrix
                       dtm(inp)=dtm(inp)+wv(ijk)*ssl
c***
   46             continue
   47          continue
   48       continue
c
   55    continue
   50 continue
  88    continue
c**
c
      return
c***** end of subroutine ttmder *****
      end
c

c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c

      subroutine bldmap
      implicit none
c  common block variables:
      include 'ray_common.inc'
      integer ixmax, iymax, izmax,ix,i,ix1,iy,iy1,iz,iz1
      real xnow, ynow, znow
      
c
c     array size limits
c     ixkms=iykms=izkms=1500
c
c     write(6,400)
c 400 format(' subroutine bldmap')
      xl=bld-xn(1)
      ixmax=(xn(nx)+xl)/bld
      yl=bld-yn(1)
      iymax=(yn(ny)+yl)/bld
      zl=bld-zn(1)
      izmax=(zn(nz)+zl)/bld
c     write(6,402)ixmax,iymax,izmax
c 402 format(' array sizes: ',3i5)
c
c  Check for array size overflow
      if(ixmax.gt.ixkms.or.iymax.gt.iykms.or.izmax.gt.izkms)goto 330
      ix=1
      do 10 i=1,ixmax
c
         ix1=ix+1
c
         xnow=float(i)*bld-xl
         if (xnow.ge.xn(ix1)) ix=ix1
c
         ixloc(i)=ix
   10 continue
c  Fill remainder of array with zeroes.
      do 12 i=ixmax,ixkms
         ixloc(i)=0
   12 continue
c
c
      iy=1
      do 15 i=1,iymax
c
         iy1=iy+1
c
         ynow=float(i)*bld-yl
         if (ynow.ge.yn(iy1)) iy=iy1
c
         iyloc(i)=iy
   15 continue
c
c  Fill rest of array with zeroes.
      do 17 i=iymax,iykms
         iyloc(i)=0
 17   continue
c
      iz=1
      do 20 i=1,izmax
c
         iz1=iz+1
c
         znow=float(i)*bld-zl
         if (znow.ge.zn(iz1)) iz=iz1
c
         izloc(i)=iz
   20 continue
c
c  Fill remainder of array with zeroes.
      do 22 i=izmax,izkms
         izloc(i)=0
  22  continue
      return
 330   continue
      write(16,331)ixkms,iykms,izkms
 331  format(' ***** error in array size in common/locate/',/,
     *' maximum map dimensions (km)=',/,' x=',i5,' y=',i5,' z=',i5)
      write(16,332)ixmax,iymax,izmax
  332 format(' Actual map size (km): ',/,' x=',i5,' y=',i5,' z=',i5)
      stop
c***** end of subroutine bldmap *****
      end

      subroutine cmpsep(path,pthsep,ns)
      implicit none
      integer n
c
c  parameters
      real path(1590),pthsep(530)
c
      integer ns
c  local variables
      integer nx,ny,nz,nx1,ny1,nz1
c
      nx=-2
c  loop over pairs of points in one set of stored vectors
      do 10 n=1,ns
         nx=nx+3
         ny=nx+1
         nz=nx+2
         nx1=nx+3
         ny1=nx+4
         nz1=nx+5
c
         pthsep(n)=sqrt((path(nx1)-path(nx))**2
     *              +(path(ny1)-path(ny))**2
     *              +(path(nz1)-path(nz))**2)
c
   10 continue
c
c***** end of subroutine cmpsep *****
      return
      end

      subroutine minima(no,isp,fstime,npbmax,jpb)
      implicit none

      include 'ray_common.inc'
c
c*****this routine finds the minimum path using pseudo-bending
c
c  common block variables:
      integer no, isp, npbmax, jpb
      real fstime
      real x(530),y(530),z(530),v(530),vq(530),tra,qtra
      integer n, nn
      common/pathm/x,y,z,v,vq,tra,qtra,n,nn
      real xtemp(530),ytemp(530),ztemp(530),rtemp(530),ttemp(530)
      common/temp/xtemp,ytemp,ztemp,rtemp,ttemp
      integer nt, j
      real xt(530),yt(530),zt(530),rt(530),rm,rs,rtm,rts
      common/pat/xt,yt,zt,rt,rm,rs,rtm,rts,nt,j
      real upbtot,nupb
      common/upb/ upbtot,nupb
      integer i
      real xp, yp, zp, vp,vpp, ta, deltat
      
c
      tra=fstime
      n=nrp
c
      do 20 i=1,n
         x(i)=rp(1,i)
         y(i)=rp(2,i)
         z(i)=rp(3,i)
   20 continue
c
      do 21 i=1,n
         xp=x(i)
         yp=y(i)
         zp=z(i)
         xtemp(i)=xp
         ytemp(i)=yp
         ztemp(i)=zp
         call vel3(isp,xp,yp,zp,vp)
         v(i)=vp
         call vel3(1,xp,yp,zp,vpp)
         vq(i)=vpp
   21 continue
      call travel
c
      nn=n-1
      nupb=nupb+1
c     write(6,6000) nupb
 6000 format(' nubp=',i8)
      do 100 j=1,npbmax
         ta=tra
         call bend(isp,xfac)
         call travel
         deltat=ta-tra
c
         if (deltat.lt.0.0) go to 102
         do 22 i=1,n
            x(i)=xtemp(i)
            y(i)=ytemp(i)
            z(i)=ztemp(i)
   22    continue
         if(deltat.le.tlim) go to 102
100   continue
c
c
  102 continue
      if(iuseq .eq. 1) then
        call qtravel
      endif


      if(j.gt.npbmax) j=npbmax
      upbtot=upbtot + float(j)
      jpb=j
  105 do 300 i=1,n
         rp(1,i)=x(i)
         rp(2,i)=y(i)
         rp(3,i)=z(i)
300   continue
      if(iuseq .eq. 0) then
        fstime=tra
      else
      fstime=qtra
      endif
c
      return
c ***** end of subroutine minima *****
      end

      subroutine travel
      implicit none
      real xtemp(530),ytemp(530),ztemp(530),rtemp(530),ttemp(530)
c
      common/temp/xtemp,ytemp,ztemp,rtemp,ttemp
      
c
      real x(530),y(530),z(530),v(530),vq(530),tra,qtra
      integer n, nn
      common/pathm/x,y,z,v,vq,tra,qtra,n,nn
      integer i1, i
      real xd, yd, zd, ds,tv
c
      tra=0
      do 60 i=2,n
         i1=i-1
         xd=xtemp(i)-xtemp(i1)
         yd=ytemp(i)-ytemp(i1)
         zd=ztemp(i)-ztemp(i1)
         ds=sqrt(xd*xd+yd*yd+zd*zd)
         tv=ds*(1.0/v(i)+1.0/v(i1))
         tra=tra+tv
  60  continue
      tra=0.5*tra
c
      return
c ***** end of subroutine travel *****
      end

      subroutine qtravel
      implicit none
c
      real xtemp(530),ytemp(530),ztemp(530),rtemp(530),ttemp(530)
      common/temp/xtemp,ytemp,ztemp,rtemp,ttemp
c
      real x(530),y(530),z(530),v(530),vq(530),tra,qtra
      integer n, nn
      common/pathm/x,y,z,v,vq,tra,qtra,n,nn
      integer i, i1
      real xd, yd, zd,ds, tv
c
      qtra=0
      do 60 i=2,n
         i1=i-1
         xd=xtemp(i)-xtemp(i1)
         yd=ytemp(i)-ytemp(i1)
         zd=ztemp(i)-ztemp(i1)
         ds=sqrt(xd*xd+yd*yd+zd*zd)
         tv=ds*(1.0/vq(i)+1.0/vq(i1))
         qtra=qtra+tv
  60  continue
      qtra=0.5*qtra
c
      return
c ***** end of subroutine qtravel *****
      end

      subroutine bend(isp,xfac)
      implicit none
      integer isp
      real xfac
c*****this routine perturbs the initial path in the direction
c      of the normal to the ray path tangent at each point
c      by the optimal distance r
c
      real x(530),y(530),z(530),v(530),vq(530),tra,qtra
      integer n,nn
      common/pathm/x,y,z,v,vq,tra,qtra,n,nn
      real xtemp(530),ytemp(530),ztemp(530),rtemp(530),ttemp(530)
      common/temp/xtemp,ytemp,ztemp,rtemp,ttemp
      integer k, kk, kkk
      real dx, dy, dz, dn, ddn, rdx, rdy, rdz
      real xk, yk, zk, vk, vx, vy, vz, vrd, rvx, rvy, rvz, rvs,rcur
      real xxk, yyk, zzk,vkk
c
c ***
      xtemp(1)=x(1)
      ytemp(1)=y(1)
      ztemp(1)=z(1)
c ***
      do 200 k=2,nn
c
         kk=k-1
         kkk=k+1
c
c*****compute the normal direction of maximum gradient of velocity
c
         dx=x(kkk)-xtemp(kk)
         dy=y(kkk)-ytemp(kk)
         dz=z(kkk)-ztemp(kk)
         dn=dx*dx+dy*dy+dz*dz
         ddn=sqrt(dn)
         rdx=dx/ddn
         rdy=dy/ddn
         rdz=dz/ddn
c
         xk=0.5*dx+xtemp(kk)
         yk=0.5*dy+ytemp(kk)
         zk=0.5*dz+ztemp(kk)
c ***
         call vel3(isp,xk,yk,zk,vk)
         call veld(isp,xk,yk,zk,vx,vy,vz)
c
c ***
         vrd=vx*rdx+vy*rdy+vz*rdz
         rvx=vx-vrd*rdx
         rvy=vy-vrd*rdy
         rvz=vz-vrd*rdz
c
         rvs=sqrt(rvx*rvx+rvy*rvy+rvz*rvz)
         if(rvs.eq.0.0) goto 200
         rvx=rvx/rvs
         rvy=rvy/rvs
         rvz=rvz/rvs
c
c*****compute the optimal distance r
         rcur=vk/rvs

         if(rcur*rcur.gt.0.25*dn) then
            rtemp(k)=rcur-sqrt(rcur*rcur-0.25*dn)
         else
            rtemp(k)=rcur
         endif
c
c*****compute the new points and distance of perturbations
c
         xxk=xk+rvx*rtemp(k)
         yyk=yk+rvy*rtemp(k)
         zzk=zk+rvz*rtemp(k)
c
c  convergence enhancement
         xxk=xfac*(xxk-x(k))+x(k)
         yyk=xfac*(yyk-y(k))+y(k)
         zzk=xfac*(zzk-z(k))+z(k)
c
         ttemp(k)=sqrt((x(k)-xxk)**2+(y(k)-yyk)**2+(z(k)-zzk)**2)
         xtemp(k)=xxk
         ytemp(k)=yyk
         ztemp(k)=zzk
         call vel3(isp,xxk,yyk,zzk,vk)
         v(k)=vk
         call vel3(1,xxk,yyk,zzk,vkk)
         vq(k)=vkk
c ***
200   continue
c
      return
c ***** end of subroutine bend *****
      end

      subroutine veld(isp,xx,yy,zz,vx,vy,vz)
      implicit none
      include 'ray_common.inc'
      integer isp
      real xx, yy, zz, vx, vy, vz
c

c*****this routine computes the derivatives of velocity
c     in x, y, and z directions
c
c  common block variables:
      integer ip, jp, kp,ip1,jp1, kp1
      real xd,yd,zd,xf,yf,zf,xf1,yf1,zf1
c
c  use Prothero's intmap here
      call intmap(xx,yy,zz,ip,jp,kp)
c
      ip1=ip+1
      jp1=jp+1
      kp1=kp+1
c
      xd=xn(ip1)-xn(ip)
      yd=yn(jp1)-yn(jp)
      zd=zn(kp1)-zn(kp)
c
      xf=(xx-xn(ip))/xd
      yf=(yy-yn(jp))/yd
      zf=(zz-zn(kp))/zd
c
      xf1=1.0-xf
      yf1=1.0-yf
      zf1=1.0-zf
c
c  S-velocity is stored in the 2nd half of the velocity array
      if(isp.eq.1) kp=kp+nz
      kp1=kp+1
c
c*****calculate derivatives of velocity
c
      vx=(yf1*zf1*(vel(ip1,jp,kp)-vel(ip,jp,kp))
     *+yf*zf1*(vel(ip1,jp1,kp)-vel(ip,jp1,kp))
     *+yf1*zf*(vel(ip1,jp,kp1)-vel(ip,jp,kp1))
     *+yf*zf*(vel(ip1,jp1,kp1)-vel(ip,jp1,kp1)))/xd
c
      vy=(xf1*zf1*(vel(ip,jp1,kp)-vel(ip,jp,kp))
     *+xf*zf1*(vel(ip1,jp1,kp)-vel(ip1,jp,kp))
     *+xf1*zf*(vel(ip,jp1,kp1)-vel(ip,jp,kp1))
     *+xf*zf*(vel(ip1,jp1,kp1)-vel(ip1,jp,kp1)))/yd
c
      vz=(xf1*yf1*(vel(ip,jp,kp1)-vel(ip,jp,kp))
     *+xf*yf1*(vel(ip1,jp,kp1)-vel(ip1,jp,kp))
     *+xf1*yf*(vel(ip,jp1,kp1)-vel(ip,jp1,kp))
     *+xf*yf*(vel(ip1,jp1,kp1)-vel(ip1,jp1,kp)))/zd
c
      return
c ***** end of subroutine veld *****
      end

      subroutine input3
      implicit none
c  this routine reads in the initial velocity model in the
c  form of velocity specified on a uniform but not
c  necessarily evenly spaced grid of points
c  (reads from file03 )
c
c  common block variables:
      include 'ray_common.inc'
c
c  declaration statements:
      integer n, ierror,k,i, kv,k2,j,ks,iizf,mnode,kk,n1,n2
      integer in,nma,km,jm,im,inf2,i1, is
      real atemp 
      integer ixf(maxpar),iyf(maxpar),izf(maxpar)
c
c  start cht 1998

c  end cht 1998
      character*1 vtype(2)
      real zero
      integer izero
      parameter(zero=0.0,izero=0)
c
      
      iusep=1
      iuseq=0
      invdel=1
      ifixl=0

      do n=1,maxpar
        cnode(n)='0'
      enddo

      ierror=0
      vtype(1)='P'
      vtype(2)='S'
c
c  for this version the gridpoints can be unevenly spaced
c  the origin of the coordinate system is at (x,y,z)=(0,0,0)
c  which will not in general correspond to the point
c  xn(1),yn(1),zn(1).
c  xn,yn,zn should be factors of bld (ie: a.0 for bld=1.0 or a.b for bld=0.1)
c
c input the number of gridpoints in x, y and z directions
c  and bld factor (1.0 or 0.1 km) used to set up velocity interpolation grid
      read(3,3002) bld,nx,ny,nz
 3002 format(f4.1,3i3)
      if((bld.ne.1.0).and.(bld.ne.0.1) .and.(bld.ne.0.01)) then
        write(16,1625) bld
 1625   format(/, '******** STOP *********, bld must be 1.0,0.1,0.01,
     2   not ',f6.2)
      endif
        atemp=iuses*(nx-2)*(ny-2)*(nz-2)
        if(atemp.le.maxpar)goto 40
        write(16,42)
 42     format('0Too many nodes for program array sizes.')
        stop
 40     continue
c
c  start cht 1998
      do 123 k=1,atemp
      imerge(k)=0
      jequal(k)=0
 123  continue
c
c  end cht 1998
c
c  input the x grid, y grid, and z grid
cfh read in free format (makes life easier...)
c        read(3,3004) (xn(i),i=1,nx)
c        read(3,3004) (yn(i),i=1,ny)
c        read(3,3004) (zn(i),i=1,nz)
        read(3,*) (xn(i),i=1,nx)
        read(3,*) (yn(i),i=1,ny)
        read(3,*) (zn(i),i=1,nz)
 3003 format(3i3)
 3004 format(20f6.1)
c
      write(16,3005) bld,nx,ny,nz
 3005 format(//,' velocity grid size:',/,
     * 'bld =',f5.2,5x,' nx =',i3,5x,'ny =',i3,5x,'nz =',i3)
c
cfh give all these numbers the same format
      write(16,3006) (xn(i),i=1,nx)
cfh 3006 format(/,' xgrid',/,3x,12f7.1,8f6.1)
 3006 format(/,' xgrid',/,3x,20f8.3)
      write(16,3007) (yn(i),i=1,ny)
cfh 3007 format(/,' ygrid',/,3x,12f7.1,8f6.1)
 3007 format(/,' ygrid',/,3x,20f8.3)
      write(16,3008) (zn(i),i=1,nz)
cfh 3008 format(/,' zgrid',/,3x,8f6.1,12f7.1/)
 3008 format(/,' zgrid',/,3x,20f8.3/)
c
c  set all nodes to have fixed velocity - cht 2002
cz this set is removed by HZ
      inf=0


c
c  start cht 1998
c  lines moved followed by new code
c  compute total number of gridpoints (nodes)
      nodes=nx*ny*nz
      nxy=nx*ny
      nx2=nx-2             ! number non-edge nodes in row
      nxy2=nx2*(ny-2)      ! number non-edge nodes in layer
      nz2=nz-2
      nodes2=nz2*nxy2
c  peripheral nodes
      nx1=nx-1
      ny1=ny-1
      nz1=nz-1
c
c
c  end cht 1998
c
c  now read in the velocity values
   65 write(16,3101)
c     do 38 kv=1,iuses
         kv=1
         do 37 k=1,nz
            k2=k + (kv-1)*nz
            write(16,3015) k,vtype(kv),zn(k)
            do 36 j=1,ny
cfh               read(3,3011) (vel(i,j,k2),i=1,nx)
               read(3,*) (vel(i,j,k2),i=1,nx)
               write(16,3013) (vel(i,j,k2),i=1,nx)
   36       continue
   37    continue
c  38 continue
c CHANGE FOR VP/VS INVERSION
      if((iuses.eq.2).or.(iuseq.eq.1)) then
        do 100 k=1,nz
          if(iuseq .eq. 0) then
              write(16,3016) k,zn(k)
          else
              write(16,3017) k,zn(k)
          endif
           do 99 j=1,ny
             if(iuseq .eq. 0) then
cfh                 read(3,3011) (vpvs(i,j,k),i=1,nx)
                 read(3,*) (vpvs(i,j,k),i=1,nx)
                 write(16,3013) (vpvs(i,j,k),i=1,nx)
             else
               read(3,3014) (qval(i,j,k),i=1,nx)
               write(16,3014) (qval(i,j,k),i=1,nx)
             endif
   99     continue
  100  continue
c  compute Vs from Vp and Vp/Vs or compute 1/tstar
        kv=2
       if(iuseq .eq. 0) then
           do 120 k=1,nz
              ks=k+nz
              write(16,3015) k,vtype(kv),zn(k)  
              do 115 j=1,ny
                 do 110 i=1,nx
                    vel(i,j,ks)=vel(i,j,k)/vpvs(i,j,k)
  110           continue
                 write(16,3013) (vel(i,j,ks),i=1,nx)
  115        continue
  120     continue
       else
         do 140 k=1,nz
            ks=k+nz
            write(16,3018) k,zn(k)
            do 135 j=1,ny
               do 130 i=1,nx
                  vel(i,j,ks)=vel(i,j,k)*qval(i,j,k)
  130           continue
                 write(16,3014) (vel(i,j,ks),i=1,nx)
  135        continue
  140     continue
       endif
      endif
c
 3013 format(20f7.3)
 3014 format(20f8.3)
 3015 format(/,' layer',i3,5x,a1,' velocity',10x,'z =',f8.3)
 3016 format(/,' layer',i3,5x,'Vp/Vs',10x,'z =',f8.3)
 3017 format(/,' layer',i3,5x,'Q',10x,'z =',f8.3)
 3018 format(/,' layer',i3,5x,'Q * Vp',10x,'z =',f8.3)
 3011 format(20f5.3)
 3101 format(//,' velocity values on three-dimensional grid')
c  Number of medium parameters to invert for
      npar=nodes2*iuses
      nparv=npar
cz      if(invdel.ne.0)npar=(npar+nsts*iuses)
c
c  Check to see whether medium parameters fit within array sizes
      if(nparv.gt.maxpar) goto 980
c
cfhdmep
c get number of Vp and Vp/Vs nodes that are free in the inversion
c    nodes reduced by fixednodes
      nparpi=nodes2
      nparsi=nodes2
c  fix specified nodes by setting nfix(k)=1, else=0
      if(inf.eq.0) goto 496
        do 70 i=1,inf
           iizf=izf(i)-2
c  if s velocity node
cfhdmep           if(izf(i).gt.nz) iizf=izf(i)-4
           if(izf(i).gt.nz) then
              iizf=izf(i)-4
              nparsi=nparsi-1
           else
              nparpi=nparpi-1
           endif
c
           k=iizf*nxy2 + (iyf(i)-2)*nx2 + (ixf(i)-1)
           nfix(k)=1
           if(cnode(k).ne.'0') then
             write(16,1686) cnode(mnode)
 1686   format(' *-*-* ERROR velocity input.  This node has',
     2  ' already been ',/,' *-*-* assigned cnode= ',a1)
             ierror=1
             write(16,1681) ixf(i),iyf(i),izf(i),k
 1681        format('input3 fixed node',i5,' ixf,iyf,izf:',3i3,
     2       ' node number:',i5)
           endif
           cnode(k)='F'
   70   continue
c
       write(16,1611)
 1610 format(/,' velocity FIXED at the following nodes(1):')
 1611 format(/,' VELOCITY INVERSION GRID    0=free, F=Fixed',/,
     2  '   M=Master, C=Constant Pert. Link, ',
     3  'L=Linear Pert. Link')
  311 do 495 kv=1,iuses
         nz1=nz-1
         ny1=ny-1
         do 320 k=2,nz1
          if(iuseq.eq.0) then
            if(kv.eq.1) write(16,1009) k,vtype(kv),zn(k)
 1009       format(/,' layer',i3,5x,a1,'-velocity nodes',
     2         10x,'z =',f7.1)
            if(kv.eq.2) write(16,3016) k,zn(k)
          else
            write(16,3017) k,zn(k)
          endif
            kk=k+(kv-1)*nz2

            do 310 j=2,ny1
               n1=(kk-2)*nxy2+(j-2)*nx2+1
               n2=n1+nx2-1
               write(16,1006) (cnode(i),i=n1,n2)
c               write(16,1005) (nfix(i),i=n1,n2)
  310       continue
  320    continue
c 1005 format('    ',18i6)
 1006 format('  ',40(2x,a1))
  495 continue
  496 continue
c
c  ndexfx: index from full nodes to nodes reduced by fixed (invert nodes)
c  mdexfx: index from inversion solution nodes to full velocity nodes
      in=0
c
c  start cht 1998
      infl=inf+ilink
c
      write(16,6161) inf,ilink,infl
 6161 format(/,' number of fixed, linked, fixed+linked nodes: ',3i5)
c
      do 80 n=1,nparv
c
c  remove fixed and linked nodes from inversion solution node set
c
         if(nfix(n).eq.1) goto 80
c  start of imerge if-then-else
         if(imerge(n).eq.0) go to 888
c
c  calculate x-y-z indices of linked velocity grid
         k=(n-1)/nxy2+2
         j=2+(n-1+(2-k)*nxy2)/nx2
         i=1+n+nx2*(2-j)+nxy2*(2-k)
         if(k.ge.nz) k=k+2      ! if s velocity node
c
c  calculate x and z indices of "master" velocity grid
         nma=jequal(n)
c
         km=(nma-1)/nxy2+2
         jm=2+(nma-1+(2-km)*nxy2)/nx2
         im=1+nma+nx2*(2-jm)+nxy2*(2-km)
         if(km.ge.nz) km=km+2      ! if s velocity node
c
cDEP It seems better to have "constant" link for velocity derivatives
cDEP  to result in constant perturbation between master and linked,
cDEP  rather than just constant velocity.  So do not change initial
cDEP  model here. Of course a constant velocity initial model could
cDEP  be input if that was desired.
cDEPc  start of constant/linear link if-then-else
cDEP      if (ltype(nma).eq.1) then
cDEP      write(16,2626) i,j,k,n,im,jm,km,nma
cDEP 2626 format('  setting value for constant-type linked node: ',
cDEP     &  4i4,2x,4i4)
cDEPc
cDEP       if(iuseq.eq.0) then
cDEP           vel(i,j,k)=vel(im,jm,km)
cDEP           if (k.gt.nz)
cDEP     &     vpvs(i,j,k-nz)=vpvs(im,jm,km-nz)
cDEPc
cDEP       else
cDEP           vel(i,j,k+nz)=vel(im,jm,km+nz)
cDEP       endif
cDEPc
cDEP      else
cDEPc   continuing constant/linear link if-then-else
cDEP       il=i
cDEP       jl=j
cDEP       ktemp=k
cDEP       kmtemp=km
cDEP      if (iuseq.ne.0) then
cDEP       k=k+nz
cDEP       km=km+nz
cDEP      endif
cDEP       kl=k
cDEPc
cDEP       if (i.ne.im) then
cDEP       il=2*i-im
cDEP       vel(i,j,k)=vel(im,jm,km)+(vel(il,jl,kl)-vel(im,jm,km))*
cDEP     & (xn(i)-xn(im))/(xn(il)-xn(im))
cDEP       endif
cDEPc
cDEP       if (j.ne.jm) then
cDEP       jl=2*j-jm
cDEP       vel(i,j,k)=vel(im,jm,km)+(vel(il,jl,kl)-vel(im,jm,km))*
cDEP     & (yn(i)-yn(im))/(yn(il)-yn(im))
cDEP       endif
cDEPc
cDEP       if (k.ne.km) then
cDEP       kl=2*k-km
cDEP       vel(i,j,k)=vel(im,jm,km)+(vel(il,jl,kl)-vel(im,jm,km))*
cDEP     & (zn(i)-zn(im))/(zn(il)-zn(im))
cDEP       endif
cDEPc
cDEP      write(16,2627) i,j,k,n,im,jm,km,nma,il,jl,kl,vel(i,j,k)
cDEP 2627 format('  setting value for linear-type linked nodes: ',/,
cDEP     &  4i4,2x,4i4,2x,3i4,f7.2)
cDEPc
cDEP      if (k.gt.nz) then
cDEP       if (i.ne.im) then
cDEP       il=2*i-im
cDEP       vpvs(i,j,k)=vpvs(im,jm,km)+(vpvs(il,jl,kl)-vpvs(im,jm,km))*
cDEP     & (xn(i)-xn(im))/(xn(il)-xn(im))
cDEP       endif
cDEPc
cDEP       if (j.ne.jm) then
cDEP       jl=2*j-jm
cDEP       vpvs(i,j,k)=vpvs(im,jm,km)+(vpvs(il,jl,kl)-vpvs(im,jm,km))*
cDEP     & (yn(i)-yn(im))/(yn(il)-yn(im))
cDEP       endif
cDEPc
cDEP       if (k.ne.km) then
cDEP       kl=2*k-km
cDEP       vpvs(i,j,k)=vpvs(im,jm,km)+(vpvs(il,jl,kl)-vpvs(im,jm,km))*
cDEP     & (zn(i)-zn(im))/(zn(il)-zn(im))
cDEP       endif
cDEP      endif
cDEPc
cDEP      endif
c
      go to 80
c
c  end of linked node section
c
  888 continue
c  add to index if not linked or fixed
c
         in=in+1
         ndexfx(n)=in
         mdexfx(in)=n
c         write(16,1888) in,n
c 1888    format(' adding inversion node number ',
c     &   i5,' corresponding to gridpoint number ',i5)
c
   80 continue
      inf2=nparv-in
      if(inf2.eq.infl) goto 85
c
c  end cht 1998
c
        write(16,1615) infl,nparv,in,inf2
 1615   format(/,' **** number of fixed and linked nodes input,',i4,
     2  ', does not equal velocity nodes,',i4,', minus invert',
     3  ' nodes,',i4,'.  Continue with inf=',i5,' ****',/)
        infl=inf2
   85 continue
      nparvi=nparv-infl
      npari=nparvi
cz      if(invdel.ne.0) npari=nparvi+nstsi*iuses
      if(invdel.eq.0) goto 95
c  also set indices if station delays are included in inversion
      i1=nparv+1
      do 90 i=i1,npar
         is=i-nparv
c s-delay
cz         if(is.gt.nsts) is=is-nsts
         if(nfixst(is).eq.1) goto 90
         in=in+1
         ndexfx(i)=in
         mdexfx(in)=i
   90 continue
      npari=in
   95 continue
      write(16,1620) npar,nparv,npari,nparvi
 1620 format(' INPUT3:npar,nparv,npari,nparvi',4i6)
c  Check to see whether medium parameters fit within array sizes
      if(npari.gt.mxpari) goto 990
      if(npar.gt.maxpar) goto 995
c  Stop if there was an error reading in fixed and linked nodes
      if(ierror.ne.0) then
        write(16,1683)
        write(6,1683)
 1683   format(/,'STOP SIMUL2000! , Error in velocity input file ')
        stop
      endif
c
c  Set up an array which is used to point to node indices, for any x,y,z
      call bldmap
c
      return
c
  980 continue
      write(6,1683)
      write(16,1698) nparv,maxpar
 1698 format(/,'  ****** STOP ******',/,i8,' velocity nodes, program',
     2 ' arrays only allow',i6)
      stop
  990 continue
      write(6,1683)
      write(16,1699) npari,mxpari
 1699 format(/,'  ****** STOP ******',/,i8,' parameters to invert for',
     2 ', program arrays only allow',i6)
      stop
  995 continue
      write(6,1683)
      write(16,1695) npar,maxpar
 1695 format(/,'  ****** STOP ******',/,i8,' parameters, arrays are ',
     2 'for ',i8)
      stop
c
c***** end of subroutine input3 *****
      end
