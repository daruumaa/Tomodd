C------------------------------------------------------
c  Modified from dtres.f released with hypoDD
c  This file is modified so that the absolute arrival times can
c  also included in the system.
c-------------------- By Haijiang Zhang----------------------

	subroutine dtres_tomoDD(log, ndt, stdim, nsrc, nsta, 
     &	dt_dt, dt_idx,
     &	dt_ista, dt_ic1, dt_ic2,
     &	src_cusp, src_t, tmp_ttp, tmp_tts,
     &	dt_cal, dt_res, eve_sta)

        use tomoDD
	implicit none

c	Parameters:
	integer		log		! Log-file identifier
	integer		ndt		! No. of data
	integer		stdim		! Column dimenson of arrays tmp_tt[ps]
	integer		nsrc		! No. of sources
	integer         nsta
	real		dt_dt(MAXDATA)	! [1..ndt] Observed time differences
	integer		dt_idx(MAXDATA)	! [1..ndt]
	integer		dt_ista(MAXDATA)! [1..ndt] Station indices
	integer		dt_ic1(MAXDATA)	! [1..ndt] Event indices
	integer		dt_ic2(MAXDATA)	! [1..ndt] Event indices
	integer		src_cusp(MAXEVE)! [1..nsrc] Event keys
	real		src_t(MAXEVE)	! [1..nsrc] Event times
	real		tmp_ttp(MAXOBS,MAXEVE)! [1.., 1..nsrc]
	real		tmp_tts(MAXOBS,MAXEVE)! [1.., 1..nsrc]
	real		dt_cal(MAXDATA)	! [1..ndt] Theoretical time differences
	real		dt_res(MAXDATA)	! [1..ndt] Time-difference residuals
        integer         eve_sta(MAXEVE,MAXOBS+1)

c	Local variables:
	integer		i, k1, k2
	real		tt1, tt2


      write(log,'("~ getting residual vector...")')

	if (nsrc.eq.1) then
c       Single source
	   do i=1,ndt
	      dt_res(i) = dt_dt(i)
	   enddo
	else
c       Mulitple sources
	   tt1 = 0.0
	   tt2 = 0.0
	   do i=1,ndt
	      if (dt_idx(i).eq.1 .or. dt_idx(i).eq.3) then
c       P phase
                 call find_id(eve_sta,dt_ista(i),dt_ic1(i),k1)
	         call find_id(eve_sta,dt_ista(i),dt_ic2(i),k2)
		 tt1 = tmp_ttp(k1,dt_ic1(i)) -
     &    	      src_t(dt_ic1(i))/1000
		 tt2 = tmp_ttp(k2,dt_ic2(i)) -
     &  	      src_t(dt_ic2(i))/1000
c		 if (tt1.eq.0 .or. tt2.eq.0) then
c		    write(*,*)' P Phase'
c                   write(*,*)'k1,k2',k1,k2,eve_sta(dt_ic1(i),1),eve_sta(dt_ic2(i),1)
c		    write(*,*)  tmp_ttp(k1,dt_ic1(i))
c                   write(*,*)  tmp_ttp(k2,dt_ic2(i))
c		    write(*,*)  dt_ic1(i), dt_ic2(i)
c		    write(*,*)  dt_ista(i)
c		    write(*,'("FATAL ERROR (theor tt)")')
c		    stop
c		 endif
	      elseif (dt_idx(i).eq.2 .or. dt_idx(i).eq.4) then
c       S phase
                 call find_id(eve_sta,dt_ista(i),dt_ic1(i),k1)
                 call find_id(eve_sta,dt_ista(i),dt_ic2(i),k2)
		 tt1 = tmp_tts(k1,dt_ic1(i)) -
     &              src_t(dt_ic1(i))/1000
		 tt2 = tmp_tts(k2,dt_ic2(i)) -
     &              src_t(dt_ic2(i))/1000
c		 if (tt1.eq.0 .or. tt2.eq.0) then
c		    write(*,*)' S Phase'
c		    write(*,*) tmp_tts(k1,dt_ic1(i))
c		    write(*,*) tmp_tts(k2,dt_ic2(i))
c		    write(*,*)  dt_ic1(i), dt_ic2(i)
c		    write(*,*)  dt_ista(i)
c		    write(*,'("FATAL ERROR (theor tt)")')
c		    stop 	       
c		 endif
		 
	      endif
	      if (dt_ic1(i) .ne. dt_ic2(i) ) then ! difference time
		 dt_cal(i) = tt1 - tt2
		 dt_res(i) = dt_dt(i) - dt_cal(i)
	      else
		 dt_res(i) = dt_dt(i) - tt1 ! absolute catalog time
	      endif 
	   enddo
	endif
	
	end			!of subroutine dtres_tomoDD



