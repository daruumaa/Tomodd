        module tomoDD
c hypoDD.inc: Stores parameters that define array dimensions in hypoDD.
c             Modify to fit size of problem and available computer memory.
c Parameter Description:
c MAXEVE:   Max number of events (must be at least the size of the number 
c           of events listed in the event file)
c MAXDATA:  Max number of observations (must be at least the size of the 
c           number of observations).  
c MAXEVE0:  Max number of events used for SVD. If only LSQR is used, 
c           MAXEVE0 can be set to 2 to free up memory. 
c MAXDATA0: Max number of observations used for SVD. If only LSQR is used, 
c           MAXDATA0 can be set to 1 to free up memory. 
c MAXLAY:   Max number of model layers.
c MAXSTA:   Max number of stations.
c MAXCL:    Max number of clusters allowed.
c MAXND:    Control MAXND*MAXDATA
c MAXNODE:    Max number of inversion nodes for each ray 
	integer*8 MAXEVE, MAXLAY, MAXDATA, MAXSTA, MAXEVE0, MAXDATA0, MAXCL
	integer*8 MAXND, MAXNODE, MAXOBS

c Define variables for size calculations
        integer*8 sizeof
        integer*8 s_tp, s_tpi, s_ts, s_tsi

        parameter(MAXEVE=   2000,    
     &            MAXDATA=  1000000,
     &            MAXEVE0=  2,	
     &            MAXDATA0= 1,	
     &            MAXLAY=   50,		
     &            MAXSTA=   100,		
     &            MAXCL=    50,
     &            MAXND=    50,
     &            MAXOBS=   100,
     &            MAXNODE=  180)		

        real,allocatable ::            tmp_vp(:,:,:)
        integer,allocatable ::         tmp_vp_index(:,:,:)
        real,allocatable ::            tmp_vs(:,:,:)
        integer,allocatable ::         tmp_vs_index(:,:,:)

        contains

c       allocate the arrays
        subroutine alloc_DD()

        allocate(tmp_vp(MAXOBS,MAXEVE,MAXNODE))
        if (.not.allocated(tmp_vp)) print *, "Memory allocation ERROR for tmp_vp"
        s_tp=4
        s_tp=(s_tp*MAXOBS*MAXEVE*MAXNODE)
        print *, "Size of tmp_vp=",s_tp

        allocate(tmp_vp_index(MAXOBS,MAXEVE,MAXNODE+1))
        if (.not.allocated(tmp_vp_index)) print *, "Memory allocation ERROR for tmp_vp_index"
        s_tpi=4
        s_tpi=(s_tpi*MAXOBS*MAXEVE*(MAXNODE+1))
        print *, "Size of tmp_vp_index=",s_tpi

        allocate(tmp_vs(MAXOBS,MAXEVE,MAXNODE))
        if (.not.allocated(tmp_vs)) print *, "Memory allocation ERROR for tmp_vs"
        s_ts=4
        s_ts=(s_ts*MAXOBS*MAXEVE*MAXNODE)
        print *, "Size of tmp_vs=",s_ts

        allocate(tmp_vs_index(MAXOBS,MAXEVE,MAXNODE+1))
        if (.not.allocated(tmp_vs_index)) print *, "Memory allocation ERROR for tmp_vs_index"
        s_tsi=4
        s_tsi=(s_tsi*MAXOBS*MAXEVE*(MAXNODE+1))
        print *, "Size of tmp_vs_index=",s_tsi

        end subroutine

        end module
