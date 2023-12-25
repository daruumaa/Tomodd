c-----------------------------------------------------------
c Modified by Haijiang Zhang from getinp.f released with
c HhypoDD. Some parameters to control 3D ray tracing are 
c read. 
c----------------------------------------------------
cz 8/31/2003, add the parameter (air_dep) to determine if an event is
cz airquake

cz 5/14/2004, add output Vp and Vs file names
cz            remove iusep, iuseq,invdel,ifixl from parameters
cz            add them instead in subroutine input3

c Get input parameters

	subroutine getinp_tomoDD (MAXEV,log,fn_inp,
     &	fn_cc, fn_ct, fn_sta, fn_eve, fn_vel, fn_abs,
     &	fn_loc, fn_reloc, fn_res, fn_stares, fn_srcpar,
     &	fn_vp, fn_vs, idata, iphase,
     &	minobs_cc, minobs_ct,
     &	amaxres_cross, amaxres_net, amaxdcc, amaxdct,
     &	noisef_dt, maxdist,
     &	awt_ccp, awt_ccs, awt_ctp, awt_cts, awt_ctd, adamp,
     &	istart, maxiter, isolv, niter, aiter, ajoint, threshold,
     &	iclust, ncusp, icusp, 
     &  lat_Orig, lon_Orig, dep_Orig, iorig,
     &  rota, stepl, CC_format,
     &  weight1, weight2, weight3, air_dep)

	implicit none

	include'ray_common.inc'

c	Parameters:
	integer*8	maxev		! Array dimension
	integer		log		! Log-file identifier
	character	fn_inp*80	! File of control info.
	character	fn_cc*80	! File of cross-corr. times
	character	fn_ct*80	! File of catalog times
	character	fn_sta*80	! Station file
	character	fn_eve*80	! Event file
	character	fn_loc*80	! Output file of original locs.
	character	fn_reloc*80	! Output file of final locs.
	character	fn_res*80	! Output residual file
	character	fn_stares*80	! Output station file
	character	fn_srcpar*80	! Output source-parameter file
	character	fn_abs*80	! Absolute travel time data file
	character       fn_vel*80       ! The velocity from each inversion
	character       fn_vp*80        ! Output the final Vp model
        character       fn_vs*80        ! Output the final Vs model
	integer		idata		! 0: Synthetics
					! 1: Cross-correlation
					! 2: Catalog
					! 3: Both
	integer		iphase		! 1: P; 2: S; 3: Both
	integer		minobs_cc	! Min. obs./pair for ccor. data
	integer		minobs_ct	! Min. obs./pair for cat. data
	real		amaxres_cross(20)! [1..niter] Ccor. res. thresh.
	real		amaxres_net(20)	! [1..niter] Cat. res. thresh.
	real		amaxdcc(20)	! [1..niter] Ccor. link-dist. limit
	real		amaxdct(20)	! [1..niter] Cat. link-dist. limit
	real		noisef_dt	! Synthetic noise
	real		maxdist		! Max. cluster-station distance
	real            threshold(20)   ! Determine the threshold
	real            awt_ctd(20)     ! relative weighting between abs 
                                        ! and diff cat
	real		awt_ccp(20)	! [1..niter] Wts. for ccor. P
	real		awt_ccs(20)	! [1..niter] Wts. for ccor. S
	real		awt_ctp(20)	! [1..niter] Wts. for cat. P
	real		awt_cts(20)	! [1..niter] Wts. for cat. S
	real		adamp(20)	! [1..niter] Damping (lsqr only)
	integer         ajoint(20)      ! Joint inversion or not
	integer		istart		! 1: From single source
					! 2: From network sources
	integer		maxiter
	integer		isolv		! 1: SVD; 2: LSQR
	integer		niter		! No. of iteration sets
	integer		aiter(0:20)	! [1..niter] Iterations/set
	integer		iclust		! Cluster to relocate (0: all).
	integer		ncusp		! No. of event keys in icusp[]
	integer		icusp(maxev)	! [1..ncusp] Events to relocate
	real            stepl
	real            lat_Orig        ! Coordinate center latitude
 	real		lon_Orig        ! Coordinate center longitude
	real            dep_Orig        ! Coordinate center depth
	real		rota            ! Rotation angle 
                                        ! -: anti-clockwise
                                        ! +: clockwise
	integer         iorig           ! 1: choose (lat_Orig,lon_Orig,dep_Orig)
                                        ! as the center, otherwise use the hypocenter	
	integer         CC_format       ! CC data format: 1 hypoDD format, 
                                        !                 2 tomoDD format
	real		weight1         ! Weight at X
	real		weight2         ! Weight at Y
	real		weight3         ! Weight at Z
	real            air_dep         ! Events shallower than this depth
                                        ! are removed from the inversion.

c	Local variables:
	integer		fu_inp
	integer		i
	integer		ii
	integer		l
	character	line*80
	integer		trimlen


c--- newest format: 083000 with iteration step dependent weighting
c-- open input file:
      call freeunit (fu_inp)
      open (fu_inp,status='unknown',file=fn_inp,err=998)
      ncusp= 0
      niter= 0  ! number of iteration blocks
      l = 1
      ii= 1

c-- Loop to read each parameter lines, skipping comments
 210  read (fu_inp,'(a)',end=220) line
      if (line(1:1).eq.'*' .or. line(2:2).eq.'*') goto 210
      if (l.eq.1) read (line,'(a)',err=999) fn_cc
      if (l.eq.2) read (line,'(a)',err=999) fn_ct
      if (l.eq.3) read (line,'(a)',err=999) fn_abs
      if (l.eq.4) read (line,'(a)',err=999) fn_eve
      if (l.eq.5) read (line,'(a)',err=999) fn_sta
      if (l.eq.6) read (line,'(a)',err=999) fn_loc
      if (l.eq.7) read (line,'(a)',err=999) fn_reloc
      if (l.eq.8) read (line,'(a)',err=999) fn_stares
      if (l.eq.9) read (line,'(a)',err=999) fn_res
      if (l.eq.10) read (line,'(a)',err=999) fn_srcpar
      if (l.eq.11) read (line,'(a)',err=999) fn_vel
      if (l.eq.12) read (line,'(a)',err=999) fn_vp
      if (l.eq.13) read (line,'(a)',err=999) fn_vs
      if (l.eq.14) read (line,*,err=999) idata, iphase, maxdist
      if (l.eq.15) read (line,*,err=999) minobs_cc,minobs_ct, CC_format
      if (l.eq.16) then 
          read (line,*,err=999) istart, isolv, niter, 
     &        weight1, weight2, weight3, air_dep
      end if
      if (l.eq.17) then ! parameters for pseudo-bending ray tracing
               read (line,*,err=999) i3d, delt1, ndip, iskip, scale1, 
     &               scale2, iuses 
      end if
      if (l.eq.18) then
               read (line,*,err=999) xfac, tlim, 
     &               nitpb(1), nitpb(2), stepl
      end if
      if (l.eq.19) then 
	 read (line,*,err=999) lat_Orig, lon_Orig, dep_Orig, iorig, rota
      end if

       

c--Read iteration instructions
      if (l.ge.20 .and. l.le.19+niter) then
         i=l-19
         read (line,*,err=999) aiter(i),
     &  awt_ccp(i), awt_ccs(i), amaxres_cross(i), amaxdcc(i),
     &  awt_ctp(i), awt_cts(i), amaxres_net(i), amaxdct(i), awt_ctd(i),
     &  adamp(i),ajoint(i),threshold(i)
      endif

c--Read specific clusters/events to relocate
      if (l.eq.20+niter) read (line,*,err=999) iclust
      if (l.ge.21+niter) then
         read (line,*,err=999,end=230) (icusp(i),i=ii,ii+7)
230      ii= i
      endif
      l= l+1
      goto 210
220   close (fu_inp)
      ncusp= ii-1

c- rearrange aiter:
      do i=2,niter
        aiter(i)= aiter(i-1)+aiter(i)
      enddo

c- check files
      call exist (fn_eve)
      call exist (fn_sta)
      call exist (fn_abs)

      if ((idata.eq.1 .or.idata.eq.3).and.trimlen(fn_cc).gt.1)
     & call exist(fn_cc)
      if ((idata.eq.2 .or.idata.eq.3).and.trimlen(fn_ct).gt.1)
     & call exist (fn_ct)

      maxiter= aiter(niter)
c synthetic noise:
      noisef_dt= 0.002

c write log output: of newest format
600   if (trimlen(fn_loc).lt.2) fn_loc= 'tomoDD.loc'
      if (trimlen(fn_reloc).lt.2) fn_reloc= 'tomoDD.reloc'
      write (6,'("INPUT FILES:",/,
     &"cross dtime data: ",a,/,"catalog dtime data: ",a,/,
     &"events: ",a,/,"stations: ",a,/,"OUTPUT FILES:",/,
     &"initial locations: ",a,/,"relocated events: ",a,/,
     &"event pair residuals: ",a,/,"station residuals: ",a,/,
     &"source parameters: ",a)')
     &fn_cc(1:trimlen(fn_cc)),
     &fn_ct(1:trimlen(fn_ct)),
     &fn_eve(1:trimlen(fn_eve)),
     &fn_sta(1:trimlen(fn_sta)),fn_loc(1:trimlen(fn_loc)),
     &fn_reloc(1:trimlen(fn_reloc)),fn_res(1:trimlen(fn_res)),
     &fn_stares(1:trimlen(fn_stares)),fn_srcpar(1:trimlen(fn_srcpar))

      write (log,'("Input parameters: (from ",a,")",/,
     &"  cross dtime file: ",a,/,"  catalog dtime file: ",a,/,
     &"  station file: ",a,/,"  event file: ",a,/,
     &"  initial locations: ",a,/,"  absoute data: ",a)')
     &fn_inp(1:trimlen(fn_inp)),
     &fn_cc(1:trimlen(fn_cc)),
     &fn_ct(1:trimlen(fn_ct)),
     &fn_sta(1:trimlen(fn_sta)),
     &fn_eve(1:trimlen(fn_eve)),fn_loc(1:trimlen(fn_loc)),
     &fn_reloc(1:trimlen(fn_abs))

      write (log,'(
     &"  event pair file: ",a,/,"  station residual file: ",a,/,
     &"  source parameter file: ",a,/,
     &"  IDATA= ",i2,2X,"IPHASE= ",i2,2x,"MAXDIST= ",f5.0,/,
     &"  MINOBS_CC= ",i3,2x,"MINOBS_CT= ",i3,/,"  ISTART= ",i1,2x,
     &"ISOLV= ",i1,2x)')
     &fn_res(1:trimlen(fn_res)),
     &fn_stares(1:trimlen(fn_stares)),fn_srcpar(1:trimlen(fn_srcpar)),
     &idata,iphase,maxdist,minobs_cc,minobs_ct,istart,isolv

      aiter(0)=0
      write (log, '("  ITER ",i2,"-",i2,
     &": DAMP= "f5.1,/,"    WT_CCP= ",f7.4,2X,"WT_CCS= ",f7.4,2x,
     &"MAXR_CC= ",f7.4,2X,"MAXD_CC= ",f7.2,2X,/,
     &"    WT_CTP= ",f7.4,2x,"WT_CTS= ",f7.4,2x,"MAXR_CT= ",f7.4,2x,
     &"MAXD_CT= ",f7.2, 2x, "JOINT= ", i3)')
     &(aiter(i-1)+1,aiter(i),adamp(i),awt_ccp(i),awt_ccs(i),
     & amaxres_cross(i),
     & amaxdcc(i),awt_ctp(i),awt_cts(i), amaxres_net(i), amaxdct(i),
     & ajoint(i),
     & i=1,niter)

c--- write the pseudo-bending parameters
      write(log,*)'i3d=',i3d
      write(log,*)'delt1=',delt1
      write(log,*)'ndip=',ndip
      write(log,*)'iskip=',iskip
      write(log,*)'scale1=',scale1
      write(log,*)'scale2=',scale2
      write(log,*)'iuses=',iuses
c--- smoothing constraint
      write(log,*)'Smoothing applied.....'
      write(log,*)'Weight1=',weight1
      write(log,*)'Weight2=',weight2
      write(log,*)'Weight3=',weight3
c--- CC data format
      if(CC_format.eq.1) then
	write(log,*)'hypoDD CC format is used'
      elseif(CC_format.eq.2) then
	write(log,*)'Another CC format is used! It has the following format:'
	write(log,*)'EveID1 EveID2 Station Diff_time CC_coef CC_phase'
      else
	write(log,*)'Neither of the format is chosen! Stop the program!'
	write(log,*)'CC_format must be 1 or 2!'
	stop
      endif	

c--Repeat number of clusters, events to relocate
      if (iclust.eq.0) then
        write (*,*) 'Relocate all clusters'
        write (log,*) 'Relocate all clusters'
      else
        write (*,*) 'Relocate cluster number ',iclust
        write (log,*) 'Relocate cluster number ',iclust
      end if

      if (ncusp.eq.0) then
        write (*,*) 'Relocate all events'
        write (log,*) 'Relocate all events'
      else
        write (*,*) 'Relocate ',ncusp,' events'
        write (log,*) 'Relocate ',ncusp,' events'
      end if
      return

c--Input error handling
998   write(*,*)'>>> ERROR OPENING CONTROL PARAMETER FILE'
      goto 1000

999   write (*,*)'>>> ERROR READING CONTROL PARAMETERS IN LINE ',l
      write (*,*) line
1000  stop 'Program run aborted.'
      end  ! of subroutine getinp

