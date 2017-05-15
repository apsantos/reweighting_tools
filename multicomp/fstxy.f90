
MODULE hdata
	integer min_part1,max_part1,min_part2,max_part2
	integer nfiles ! number of files being read
	integer tnum	! temperature counter
	character*20, dimension(:), allocatable :: filename
	real*8 betap,mu1p,mu2p		! trial parameters
	real*8 avgnum1,avgnum2	! average number of each particle type
	real*8 nliq1,ngas1,nliq2,ngas2 ! number of particles on coexistence
	real*8 logZ,logZgas,logZliq	! logs of partition functions, exponent of calculation
	real*8 nl1,ng1,nl2,ng2	! initial guess for number of particles on coexistence
	real*8, dimension(:), allocatable :: beta
	real*8, dimension(:), allocatable :: mu1
	real*8, dimension(:), allocatable :: mu2
	real*8, dimension(:), allocatable :: nentry
	real*8, dimension(:), allocatable :: weight	
	real*8, dimension(:), allocatable :: oldweight
	real*8, dimension(:), allocatable :: dens1 ! n1 distribution
	real*8, dimension(:), allocatable :: dens2 ! n2 distribution
	real*8, dimension(:,:), allocatable :: dens ! n1-n2 distribution

	!define data type - this is to required to reduce memory usage
	type column1
		integer, pointer ::flength(:) ! Number of entries in the file
	end type
	type column2
		real*8, pointer :: flength2(:) ! Number of entries in the file
	end type
	type(column1), dimension(:), allocatable :: n1,n2
	type(column2), dimension(:), allocatable :: e

END MODULE hdata

program phlst

! v. 2.21  6/24/08 - formats (from Antti-Pekka's version)
! version 2.2 - 10th Sep 2004 (APH), histogram has now one initial text line
! version 2.1 - 7/13/99 (AZP)
! version 1.1 - 11/24/97  This program iterates to find the weights for
! patching the individual histograms together.  The histogram data is read
! in the form of N1, N2, E for each configuration.  Reads the weights in
! from the convergence program 'nwfsp1.f90'.

! modifications
! 11/24/97 - added code to calculate the average energy of each phase

! 3/9/98 - added a scheme to calculate the midpoint of the energy 
! distribution.  This allows the histograms to be added together in
! any order.

! 5/28/98 - modified storage of histogram lists to minimize memory
!           usage.

! 12/16/98 - changed iteration method from bisection to Newton's: V4.
USE hdata
implicit none
! LIST FORMAT VARIABLES
integer, parameter :: mentry = 1e6	! maximum number of observations
integer iostat				! iostat becomes -1 on end-of-file
integer ipart,i,idiff	! local counters
integer min_part1f,max_part1f,min_part2f,max_part2f		! min. and max. number of particles in file
integer ifile,jfile,kfile,iter
integer, target :: n1_temp(mentry)	! n1 temp array
integer, target :: n2_temp(mentry)	! n2 temp array
real*8 e_temp(mentry)			! energy temp array
real*8 eng_min,eng_max 	!min and max of energies stored in the lists
real*8 prob,maxd,val,y,ylog,tmp1	
real*8 xdim,ydim,zdim	! linear dimensions of box on which data were taken
real*8 xdim1,ydim1,zdim1	! linear dimensions of box for current file


!PHASE COEXISTENCE VARIABLES
integer, parameter :: maxtemp = 10	! maximun number of temperatures 
integer, parameter :: maxiter=500	! maximum iterations for phase equil.
integer ntemp	! number of temperatures to calc phase coexistence
real*8 eng_liq,eng_gas 		! ave energy of the liquid and gas phases
real*8 dmu1p
real*8 nbelow, nbelow1	 ! for Newton's method
real*8 mu2min(maxtemp),mu2max(maxtemp)
real*8 mu2_incr(maxtemp)	! for calculation of phase diagram
real*8 t_new(maxtemp)
character*20 ftemp,ftemp2,fname,fname2	!internal read/write for file naming
real*8 mu1_old,mu1incr	! for continuation

! VARIABLES FOR PHASE COEXISTENCE PRINTOUTS
integer nstart,j,nend,m,num
character*20 file,file1,file2,file3
real*8 beta_tmp,mu1_tmp,mu2_tmp ! temporary variables

write (*,*) '--------------------------------------------------------------------------'
write (*,*) ' FS Txy phase behavior, v. 2.21, July 08 (AZP) [run fspatch first]'
write (*,*) '--------------------------------------------------------------------------'

open(unit=2,file='inp_fs2.dat',status='old')
read(2,*) nfiles

allocate(n1(nfiles))
allocate(n2(nfiles))
allocate(e(nfiles))
allocate(filename(nfiles))
allocate(beta(nfiles))
allocate(mu1(nfiles))
allocate(mu2(nfiles))
allocate(nentry(nfiles))
allocate(weight(nfiles))
allocate(oldweight(nfiles))

do ifile = 1,nfiles
  read(2,*) filename(ifile),nentry(ifile),oldweight(ifile),beta(ifile),&
            mu1(ifile),mu2(ifile)
  beta(ifile) = 1./beta(ifile)
enddo
close(unit=2)

min_part1 = 1e5 + 1	! so that these are later set to correct values
max_part1 = -1	
min_part2 = 1e5 + 1	
max_part2 = -1	
eng_min = 1e5
eng_max = -1e5

write(*,'(A)') '         T        mu1       mu2    minp1   minp2  maxp1  maxp2  npoints'
										
do ifile= 1,nfiles
  
  min_part1f = 1e5 + 1	! so that these are later set to correct values
  max_part1f = -1	
  min_part2f = 1e5 + 1	
  max_part2f = -1
   		
    open (23,file='his'//trim(filename(ifile))//'.dat')
	read (23,*)
	read (23,*) beta_tmp,mu1_tmp,mu2_tmp,xdim1,ydim1,zdim1
	beta_tmp = 1./beta_tmp
	if (abs(beta_tmp-beta(ifile))+abs(mu1_tmp-mu1(ifile))+abs(mu2_tmp-mu2(ifile)) > 1.e-5) then
		write(*,*) 'Histogram file ',trim(filename(ifile)),' inconsistent with inp_fs2.dat'
		stop
	endif
    if (ifile==1) then
      xdim = xdim1
      ydim = ydim1
      zdim = zdim1
    else
      if (xdim1/=xdim.or.ydim1/=ydim.or.zdim1/=zdim) then
	write (*,*) 'System size or temperature in file ',trim(filename(ifile)),' inconsistent with first file'
	stop
      endif
    endif
    iostat = 0 
    i = 0
    do while (iostat/=-1) 
      i = i + 1    
      
      if(i+1.gt.mentry) then
        write(*,'(2f12.1)') 'maximum number of entries exceeded',&
			    nentry(ifile),mentry
        stop
      endif
        
      read (23,*,iostat=iostat) n1_temp(i),n2_temp(i),e_temp(i)     
      
      if(iostat.ne.-1) then
        eng_min = min(eng_min,e_temp(i))
        eng_max = max(eng_max,e_temp(i))
        min_part1 = min(n1_temp(i),min_part1)
        max_part1 = max(n1_temp(i),max_part1)
        min_part2 = min(n2_temp(i),min_part2) 
        max_part2 = max(n2_temp(i),max_part2)
      
        min_part1f = min(n1_temp(i),min_part1f)
        max_part1f = max(n1_temp(i),max_part1f)
        min_part2f = min(n2_temp(i),min_part2f) 
        max_part2f = max(n2_temp(i),max_part2f)
      endif
            
    enddo
    nentry(ifile) = i - 1
    
    if(nentry(ifile).gt.mentry) then
      write(*,*) 'maximum number of entries exceeded', nentry(ifile),mentry
      stop
    endif
    
    allocate(n1(ifile)%flength(1:i-1))
    allocate(n2(ifile)%flength(1:i-1))
    allocate(e(ifile)%flength2(1:i-1))
    
    n1(ifile)%flength(1:i-1) = n1_temp(1:nentry(ifile))
    n2(ifile)%flength(1:i-1) = n2_temp(1:nentry(ifile))
    e(ifile)%flength2(1:i-1) = e_temp(1:nentry(ifile))
    write(*,'(2x,A,f8.3,2x,f8.3,2x,f8.3,2x,i5,2x,i5,2x,i5,2x,i5,2x,f10.0)') &
      trim(filename(ifile)), &
      1./beta(ifile),mu1(ifile),mu2(ifile),min_part1f,&
      min_part2f,max_part1f,max_part2f,nentry(ifile)
       
enddo

dmu1p = mu1(1)/10000.

allocate(dens(min_part1:max_part1,min_part2:max_part2))
allocate(dens1(min_part1:max_part1))
allocate(dens2(min_part2:max_part2))

close(unit=23)
iter = 0 
maxd = 100

! start phase coexistence calculation

open(unit=1,file='inp_txy.dat',status='old')
read(1,*)
read(1,*) mu1p, mu1incr
read(1,*)
iostat = 0
i = 1
do while(iostat.ne.-1)  
  read (1,*,iostat=iostat) t_new(i),mu2min(i),mu2max(i),&
                           mu2_incr(i),nl1,nl2,ng1,ng2
  if(iostat.ne.-1) then
     i = i + 1
  endif
enddo
close(1)            

ntemp = i-1
  
open(unit=25, file='phase.dat', status = 'unknown')
write (25,'(A,20A5/3h /*,50A/3h /*,50A)') '/* Files = '&
	,(trim(filename(ifile)),ifile=1,nfiles)
write (25,'(A)') '/* T         mu1       mu2   N_liq1  N_liq2  N_gas1  Ngas2  ln(Zliq) ln(Zgas)'

write (*,'(A)') ' ----    T       mu1    mu2     N_liq1  N_liq2  N_gas1   Ngas2 ---------'
do tnum=1,ntemp
	betap = 1./t_new(tnum)  
	num = 0
	if(mu2_incr(tnum).gt.0) then
		mu2p = mu2min(tnum)
	else
		mu2p = mu2max(tnum)
	end if

	do while (mu2p.lt.mu2max(tnum)+.00001.and.mu2p.gt.mu2min(tnum)-.00001)    
  	
		iter = 0      
		num = num + 1
		nbelow = 0.0
		do while(abs(0.5-nbelow) > 1.e-4.and.iter.le.maxiter)      	
      
			if(iter > 0) then
				mu1p = mu1p + dmu1p
				call Zcalc(nbelow1)
				! set ln(x) - ln(1-x) to zero as this behaves better in convergence
				mu1p = mu1p - (log(nbelow1)-log(1.-nbelow1))*dmu1p/ &
						 (log(nbelow1)-log(nbelow) - log(1.-nbelow1) + log(1.-nbelow))
			endif
		
			call Zcalc(nbelow)
		
			iter = iter + 1

!			write (*,'(A,i4,g12.4,f10.4)') 'Iteration, nbelow, mu1 =',iter,nbelow,mu1p  ! debug

		enddo
	 	
		ngas1 = ngas1/nbelow
		ngas2 = ngas2/nbelow
		nliq1 = nliq1/(1.-nbelow)
		nliq2 = nliq2/(1.-nbelow)

		nl1 = nliq1
		nl2 = nliq2
		ng1 = ngas1
		ng2 = ngas2
            
		write (25,'(f7.3,2f10.3,f8.3,3f8.3,2f8.3,2x,f6.2,2x,f6.3)') &
			1./betap,mu1p,mu2p, nliq1,nliq2,ngas1,ngas2,logZliq,logZgas
		  
            
		write (*,'(5x,2f8.3,f7.2,2x,4f8.3)') 1./betap,mu1p,mu2p,nl1,nl2,ng1,ng2

		mu1_old = mu1p
		if (tnum == 1) then
			mu1p = mu1p + mu1incr
		else
			mu1p = 2*mu1p - mu1_old	! continuation
		endif

		! use internal read/write to turn integers into characters    
		write(ftemp,*) tnum
		write(ftemp2,*) num
		read(ftemp,*) fname
		read(ftemp2,*) fname2
    
!		file ='denseq'//trim(fname)//'t'//trim(fname2)//'a.dat'
!		file2 = 'd1n'//trim(fname)//'t'//trim(fname2)//'a.dat'
!		file3 = 'd2n'//trim(fname)//'t'//trim(fname2)//'a.dat'

!		open(unit=30, file = file, status = 'unknown')
!		write(30,*) '  Tstar       Mu1	Mu2','     Xdim     Ydim     Zdim'
!		write(30,'(f8.3,2x,f10.2,2x,f10.2,3i5)') 1./betap,mu1p,mu2p,int(xdim),int(ydim),int(zdim)
      
		do ipart = min_part1, max_part1
      		nstart = -1
			j = min_part2
			do while(nstart.eq.-1.and.j.lt.max_part2)
				if(dens(ipart,j).ne.0) then
					nstart = j
				else
 					j = j + 1
				end if 
			enddo
	    
			if(j.lt.max_part2) then
				nend = -1
				m = max_part2
				do while(nend.eq.-1)
					if(dens(ipart,m).ne.0) then
						nend = m
					else
						m = m -1
					end if 
				enddo
!				write(30,'(3i5)') ipart,nend-nstart+1,nstart
!				write(30,'(5g15.5e3)') (exp(dens(ipart,m)-logZ),m=nstart,nend)
			end if 
		enddo
	  
!		open(unit=31, file = file2, status ='unknown')
!		write(31,'(A,f8.3,2x,f10.2,2x,f10.2,3i5)') '/*', 1./betap,mu1p,mu2p,int(xdim),int(ydim),int(zdim)
!		do ipart = min_part1,max_part1
!			write(31,'(i6,g15.4e3)') ipart, exp(dens1(ipart)-logZ)
!		enddo
!		close(unit=31)

!		open(unit=32, file = file3, status ='unknown')
!		write(32,'(A,f8.3,2x,f10.2,2x,f10.2,3i5)') '/*',1./betap,mu1p,mu2p,int(xdim),int(ydim),int(zdim)
!		do ipart = min_part2,max_part2
!			write(32,'(i6,g15.4e3)') ipart,exp(dens2(ipart)-logZ)
!		enddo
!		close(unit=32)
    
		mu2p = mu2p + mu2_incr(tnum)
	enddo ! end of chemical potential loop
enddo ! end of temperature loop
 
close(unit=25)
end

!=========================================================================   
Subroutine Zcalc(nbelow)
USE hdata

integer ipart,i,idiff	! local counters
real*8 prob,maxd,val,y,ylog,tmp1
real*8 nbelow	
real dist_liq,dist_gas	! distances from liquid and gas phase average N values

logZ = -1.e9
logZgas = -1.e9
logZliq = -1.e9
dens = -1.e9
dens1 = -1.e9
dens2 = -1.e9
avgnum1 = 0.0
avgnum2 = 0.0
nbelow = 0.
ngas1 = 0.0
nliq1 = 0.0
ngas2 = 0.0
nliq2 = 0.0
		
do ifile = 1,nfiles
	do i=1,nentry(ifile)
		y = 0.0
		ylog = -1e9  ! exp(ylog) = y = 0.0
		do jfile=1,nfiles
			prob =(-(beta(jfile)-betap)*e(ifile)%flength2(i) +&
			(beta(jfile)*mu1(jfile)-betap*mu1p)*n1(ifile)%flength(i) +&
			(beta(jfile)*mu2(jfile)-betap*mu2p)*n2(ifile)%flength(i))

!           y = y+nentry(jfile)*exp(prob-log(weight(jfile)))
                   
! The following is a code rearrangement to avoid overflow problems
! This trick was compliments of Peter-Lawrence Montgomery 2/12/98
! Replace y with log(y)
! Define log_expsum(y,z) = log(exp(y)+exp(z))
! Evaluate log_expsum(y,z) = max(y,z) + log(1+exp(-abs(y-z)))

			tmp1 = prob + log(nentry(jfile)) - oldweight(jfile)
			ylog = max(ylog,tmp1) + log(1.+exp(-abs(ylog-tmp1)))
		enddo
		logZ = max(logZ,-ylog) + log(1.+exp(-abs(logZ+ylog)))
		dist_liq = (n1(ifile)%flength(i)-nl1)**2+(n2(ifile)%flength(i)-nl2)**2
		dist_gas = (n1(ifile)%flength(i)-ng1)**2+(n2(ifile)%flength(i)-ng2)**2
		if(dist_liq>dist_gas) then
  			logZgas = max(logZgas,-ylog) + log(1.+exp(-abs(logZgas+ylog)))
		else 
  			logZliq = max(logZliq,-ylog) + log(1.+exp(-abs(logZliq+ylog)))
		end if
		
		dens(n1(ifile)%flength(i),n2(ifile)%flength(i)) = &
			max(dens(n1(ifile)%flength(i),n2(ifile)%flength(i)),-ylog) + &
			log(1.+exp(-abs(dens(n1(ifile)%flength(i),n2(ifile)%flength(i))+ylog)))
		dens1(n1(ifile)%flength(i)) = &
			max(dens1(n1(ifile)%flength(i)),-ylog) + &
			log(1.+exp(-abs(dens1(n1(ifile)%flength(i))+ylog)))
		dens2(n2(ifile)%flength(i)) = &
			max(dens2(n2(ifile)%flength(i)),-ylog) + &
			log(1.+exp(-abs(dens2(n2(ifile)%flength(i))+ylog)))
	enddo  
enddo
		
		      
do ipart = min_part1,max_part1
	do i = min_part2,max_part2
		avgnum1 = avgnum1 + ipart * exp(dens(ipart,i)-logZ)
		avgnum2 = avgnum2 + i * exp(dens(ipart,i)-logZ)
		dist_liq = (ipart-nl1)**2+(i-nl2)**2
		dist_gas = (ipart-ng1)**2+(i-ng2)**2
		if (dist_liq > dist_gas) then
  			nbelow = nbelow + exp(dens(ipart,i)-logZ)
  			ngas1 = ngas1 + ipart * exp(dens(ipart,i)-logZ)
  			ngas2 = ngas2 + i * exp(dens(ipart,i)-logZ)
		else
  			nliq1 = nliq1 + ipart * exp(dens(ipart,i)-logZ)
  			nliq2 = nliq2 + i * exp(dens(ipart,i)-logZ)
		endif     
	enddo
enddo


	
end subroutine Zcalc
