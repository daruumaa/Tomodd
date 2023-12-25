c   mxpari: maximum number of parameters to invert for.  This 
c           includes Vp and Vp/Vs nodes and P and S station 
c           corrections.  mxpari will generally be much lower
c           than maxpar because velocity nodes and station
c           corrections can be fixed. (Station corrections
c           are generally fixed to zero in the region of the
c           velocity inversion.
	integer mxpari
	parameter(mxpari=900000)
c
c   maxpar: maximum number of potential parameters that could
c           be included in the inversion.  This includes all 
c           velocity nodes (includes both Vp and Vp/Vs) that
c           define the 3-D velocity model (nparv). Does not 
c           include edge nodes.
c           If station corrections are included in the inversion,
c           maxpar also includes stations (nparv+nsts or 
c           nparv + 2*nsts if P and S-P are used).
c           maxpar is set much larger than mxpari to allow
c           for assigning regional velocity (from prior
c           coarse inversions) values (fixed) in peripheral
c           areas where there is nil resolution.
	integer maxpar
      	parameter(maxpar=900000)
c
c   maxnx: maximum number of nodes in x direction
	integer maxnx
      	parameter(maxnx=90)
c   maxny: maximum number of nodes in y direction
	integer maxny
      	parameter(maxny=90)     
c   maxnz: maximum number of nodes in z direction
	integer maxnz
      	parameter(maxnz=90)     
c
c
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c  common block variables:
c
c   control parameters
c
	integer i3d, ndip, iskip, iusep, iuses, invdel, ifixl, nitpb(2),
     *          iuseq
	real	delt1, xfac, scale1, scale2, tlim
      	common/contrl/ i3d,delt1,ndip,iskip,scale1,scale2,iusep,
     *          iuses,invdel,ifixl,xfac,tlim,nitpb,iuseq

c
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c   event observation and station arrays
	integer netemp,nbtemp,neb
      	common/ntemp/ netemp,nbtemp,neb

	integer jfl, nco, ndo
      	common/fastrc/ jfl,nco,ndo
c


c
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c   velocity arrays
c
c   /vmod3d/
c    vel(nx,ny,2*nz)     ! 1st half is P velocity, 2nd half is S velocity
	integer maxnz2
      	parameter(maxnz2=maxnz*2)
c
	real bld,xn(maxnx),yn(maxny),zn(maxnz),vel(maxnx,maxny,maxnz2),
     *   vpvs(maxnx,maxny,maxnz), qval(maxnx,maxny,maxnz)

	integer nx, ny, nz
        
      	common/vmod3d/bld,nx,ny,nz,xn,yn,zn,
     *  vel,vpvs,
     *  qval
c
c   /locate/
c   ixkms: size of map in x direction.  The distance between 
c            furthest nodes.  This will be in units of km if bld=1.0 
c            or in units of 0.1km if bld=0.1
	integer ixkms
      	parameter(ixkms=100000)
c   iykms=size of map in y dir.
	integer iykms
      	parameter(iykms=100000)
c   izkms=size of map in z dir.
	integer izkms
      	parameter(izkms=100000)
c
	real xl, yl, zl
	integer ixloc(ixkms), iyloc(iykms), izloc(izkms)
      	common/locate/ xl,yl,zl,ixloc,iyloc,izloc
c
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  inversion arrays
c
c
cfhdmep nparpi, nparsi, nrowp,nrows,nrowst added here

	integer nz2, nx1, ny1, nz1,nparvi,nstsi,npari,inf
        integer nfix(maxpar),ndexfx(maxpar), mdexfx(mxpari),
     *          nparv,npar,nodes, nxy,khit(maxpar),nx2,nxy2,
     *          nodes2,nparpi,nparsi,nrowp,nrows,nrowst,
     *          imerge(maxpar),jequal(maxpar),ltype(maxpar),
     *          infl,imast,ilink, nfixst(1000)
	real dtm(mxpari), hit(maxpar),hitall(maxpar),vadj(maxpar),
     *       sumhit(maxnz2),qadj(maxpar)
	character cnode(maxpar)*1

      common/modinv/ nz2,nx1,dtm,
     *  ny1,nz1,nparvi,nstsi,npari,inf,nfix,ndexfx,
     *  mdexfx,nparv,npar,nodes,nxy,hit,hitall,
     *  vadj,khit,nx2,nxy2,nodes2,sumhit,
     *  qadj,nparpi,nparsi,nrowp,nrows,nrowst,cnode,
     *  imerge,jequal,ltype,infl,imast,ilink,
     *  nfixst
c  above line added cht 1998
c
c
	real rp(3,530), ttc, xmov, pl
	integer nrp
      	common/rpath/ rp,nrp,ttc,xmov, pl



