program entropy2
! version 2.39 - April '14: fixed +1 indexing error for energy, added arbitrary precision (NAM)
! version 2.38 - Oct. '05: printing of unconverged files
! version 2.37 - June '04: directory read from entr2_par.dat file
! version 2.36 - May '04: small change in output format for MATLAB compatibility
! version 2.35 - Jul '03: new obj. function w/ more weight to high points
! version 2.34 - Jun '03: fixed bug on 2D capability
! version 2.33 - Aug '02: 2D capability
! version 2.3 - Jan '00: specexp 
! version 2.2 - Apr '99: allocatable maxpart, npoints
! version 2.1 - May '97: ising.dat at fixed location; produces pvt.dat
! version 2.0 - February '97 (AZP)
! based on version 1.1 - July '96 (AZP)
! program takes as input observations of histograms of N and E
! and calculates the entropy S(N,E) = ln(Omega) by the
! Ferrenberg and Swedsen algorithm (PRL, 1989)

implicit none

integer maxpart,npoints 
real tolerance ! tolerance for convergence of weights
real accfactor ! acceleration factor
logical success
real, allocatable :: entropy(:,:)   ! contains logs of counts for a certain energy+number
									! entropy (0,maxpart) is the starting (lowest) value of the energy
									! entropy (-1,maxpart) is the highest active value of the energy
real, parameter :: min_ene = -999999.	! minimum energy
real  auxmin	! auxilliary minimum energy
real, allocatable :: omega(:,:)   ! temporary holding variable, reads contr. to entropy from
								! each file being read
real energy_start

real, allocatable :: counts(:)		! counts for every file being read per (E,N) point
real, allocatable :: counts_file(:,:)	! total counts for each file
real, allocatable :: g_n(:,:)		! related to "correlation time" for each file - g_n in S&W paper
								! if g_n=1 every configuration is considered independent, 
								! greater values mean correlated configurations.
integer ncount	! counts how many observations there are for a certain number of particles
integer min_part,max_part,ipart,i,idiff	! local counters
integer min_part1,max_part1		! min. and max. number of particles in file
real, parameter :: incr_width = 0.1	! width for energy/N in counts file
real, allocatable :: beta(:,:),mu(:,:) ! temperature, chemical potential
real width				! width of energy distributions in each histogram file
real beta1,mu1			! temperature and chem. pot. of interest for calcs.
real smix					! s field mixing parameter
integer, parameter :: maxcrit	= 10000 ! discretization of field mixing parameter
real p_l(0:maxcrit)			! ordering operator distribution
real x						! ordering operator distribution x axis
real max_x, min_x, avg_x	! limits of x, average x
integer index				! index for x	
real ncrit				! number of points for p_l
real variance				! variance of p_l
character(20), allocatable :: filename(:)	! histogram file name - if the form his?????.dat; only the ????? part is stored
real, allocatable :: weight(:,:)		! f_n in Ferr. + Swed. paper
real, allocatable :: oldweight(:,:)  	! old weight
real ln_num_P				! Log of Numerator for P(S,K) in eq. (3) of Ferr. + Swed. paper
real ln_probab				! probability of microstate at (E,N) for given T,mu
integer nfiles,ifile,jfile,jsuffix	! how many files being read
integer iostat				! iostat becomes -1 on end-of-file
real ln_Z, ln_Zliq, ln_Zgas		! partition function, also for gas and liquid
real, allocatable :: dens(:) ! dens stores the probability of a certain number
real avgnum	! mean number
real ene_liq,ene_gas	! energy of liquid and gas
real nmid	! mid-point number of particles on coexistence
real nliq,ngas,nbelow	! liquid and gas densities on coexistence, accumulator of densities for calculation
								! of coexistence
integer iter
integer, parameter :: maxiter=100	! maximum iterations for phase equil.
real tmin,tmax,t_incr	! for calculation of phase diagram
real nbelow1	! for Newton's method
real muincr ! for Newton's derivative eval.
real, parameter :: maxexponent = 300	! for avoiding overflows
real mindens ! minimum accepted peak separation
character(2), allocatable :: suffix(:)	! for dealing with series of runs
integer nsuffix,isuffix	! number of different run series
logical do_phase	! indicates if a phase equil. calculation is to be done
real width1, xdim1, ydim1, zdim1	! these are only read, not used anywhere
real zfactor ! auxiliary variable in calculation of weight and Z
logical is_2D, read_s	! determines if calculation is against 2D distribution, if s is read from screen
character*60 directory	! location of ising.dat files
integer		nising ! = 168	(3D) 224 (2D) ! number of points in ising.dat file
real, allocatable :: xising(:),yising(:)	! from ising.dat
real specexp,tmp1,auxreal
specexp(auxreal,tmp1) = max(auxreal,tmp1) + log(1.+exp(-abs(auxreal-tmp1)))

integer ipoint	! for counting on ising data
real funk
external funk
integer max_index

real iPrecise 	! added by NAM
iPrecise = 1.0d15

open (1,file='entr2_par.dat')
read (1,*) directory 
read (1,*) muincr, is_2D, read_s 
read (1,*) tolerance,accfactor
read (1,*) mindens
read (1,*) smix,ncrit
close (1)

if (is_2D) then
	open (1,file=trim(directory)//'ising2D.dat')
	nising = 224
else
	open (1,file=trim(directory)//'ising.dat')
	nising = 168
endif
allocate (xising(nising),yising(nising))
do ipoint =1,nising
	read (1,*) xising(ipoint),yising(ipoint)
enddo
close(1)

open (1,file='input_hs2.dat')
read (1,*) nsuffix
allocate (suffix(nsuffix))
do isuffix = 1,nsuffix
	read (1,*) suffix(isuffix)
enddo
read (1,*) nfiles, npoints, maxpart
close (1)

allocate (entropy(-1:npoints,0:maxpart),omega(-1:npoints,0:maxpart))
allocate (counts(npoints),counts_file(nfiles,nsuffix),g_n(nfiles,nsuffix))
allocate (beta(nfiles,nsuffix),mu(nfiles,nsuffix))
allocate (filename(nfiles),weight(nfiles,nsuffix),oldweight(nfiles,nsuffix))
allocate (dens(0:maxpart))

write (*,'(a)') '  -----------------------------------------------------------------'
if (is_2D) then
	write (*,'(a)') '  2D Ferrenberg-Swedsen patching - version 2.39 (NAM, April 2014)'
else
	write (*,'(a)') '  Ferrenberg-Swedsen patching - version 2.39 (NAM, April 2014)'
endif
write (*,'(a)') '  -----------------------------------------------------------------'

success = .false.
do while (.not. success) 
	!entropy = min_entr
	counts = 0
	entropy = -iPrecise	! program looks at counts, so a zero value is not taken seriously
					! unless backed up by data -iPrecise corresponds to a small count
	entropy(0,0:maxpart) = min_ene
	entropy(-1,0:maxpart) = min_ene

	open (1,file='input_hs2.dat')
	read (1,*) 
	do isuffix = 1,nsuffix
		read (1,*) 
	enddo
	read (1,*) 
	do ifile = 1,nfiles
		do isuffix = 1,nsuffix
			read (1,*) filename(ifile), counts_file(ifile,isuffix), &
					   weight(ifile,isuffix),beta(ifile,isuffix),mu(ifile,isuffix),width,g_n(ifile,isuffix)
			beta(ifile,isuffix) = 1./beta(ifile,isuffix)
			if (.not.(ifile==1 .and. isuffix==1)) weight(ifile,isuffix) = weight(ifile,isuffix) - weight(1,1)
		enddo
	enddo
	weight(1,1) = 0.
	close (1)

	min_part = maxpart+1	! so that these are later set to correct values
	max_part = -1	! so that these are later set to correct values

	do ifile = 1,nfiles

		do isuffix = 1,nsuffix	! no ident for this do loop
		open (23,file='his'//trim(filename(ifile))//trim(suffix(isuffix))//'.dat')
		read (23,*)
		read (23,*) beta1,mu1,width1,xdim1,ydim1,zdim1	! although not actually used, these need to 
														! be read so that the correct number of lines are
														! skipped in the input file.
		min_part1 = maxpart + 1
		max_part1 = - 1
		iostat = 0
		omega = -iPrecise
		omega(0,0:maxpart) = min_ene
		omega(-1,0:maxpart) = min_ene
		counts = 0
		dens = 0.
		do while (iostat/=-1)
			read (23,*,iostat=iostat) ipart,ncount,energy_start
			if (iostat/=-1) then
				if (ipart>maxpart.or.ncount>npoints) then
					write (*,*) 'Maxpart or Npoints exceeded, file, ipart, ncount = ',trim(filename(ifile))//trim(suffix(isuffix)),ipart,ncount
					stop
				endif
				omega(0,ipart) = energy_start
				omega(-1,ipart) = energy_start + ncount*width ! ncout-1 ?? but turns out last (superfluous point) isnt used
				read (23,*) (counts(i),i=1,ncount)
				do i = 1,ncount
					if (counts(i)>0) then
						ln_num_P = -iPrecise
						do jfile = 1,nfiles
							do jsuffix = 1,nsuffix
								! to keep exponentials within reason, all are referenced to beta(1,1) and
								! mu(1,1)
								zfactor =  + beta(1,1)*(omega(0,ipart)+(i-1)*width) - &
								   beta(1,1)*mu(1,1)*ipart &
								   - beta(jfile,jsuffix)*(omega(0,ipart)+(i-1)*width) + &
								   beta(jfile,jsuffix)*mu(jfile,jsuffix)*ipart - weight(jfile,jsuffix)
								if (counts_file(jfile,jsuffix)>0) &
									ln_num_P  = specexp(ln_num_P,log(counts_file(jfile,jsuffix)/g_n(jfile,jsuffix))+zfactor)
							enddo
						enddo
						omega(i,ipart) = log(counts(i))-log(g_n(ifile,isuffix))-ln_num_P
					endif
				enddo
				min_part = min(min_part,ipart)
				max_part = max(max_part,ipart)
				min_part1 = min(min_part1,ipart)
				max_part1 = max(max_part1,ipart)
			endif
		enddo
		close(23)
		
		if (ifile==1.and.isuffix==1) then
			entropy = omega
		else ! patch current data with previous data 
			do ipart = min_part,max_part
				idiff = 0
				if (omega(0,ipart)/=min_ene.and.entropy(0,ipart)/=min_ene) idiff = nint (( entropy(0,ipart) - omega(0,ipart)  ) / width)
				if ( idiff > 0) then	! shift entropy
					if ((max(entropy(-1,ipart),omega(-1,ipart))-omega(0,ipart))/width <= npoints) then
						entropy(0,ipart) = omega(0,ipart)
						entropy(-1,ipart) = max(entropy(-1,ipart),omega(-1,ipart))
						do i = npoints,1+idiff,-1
							entropy(i,ipart) = entropy(i-idiff,ipart)
						enddo
						do i = 1,1+idiff-1
							entropy(i,ipart) = -iPrecise
						enddo
					else
						write (*,*) ' Width exceeded (low), ipart,file = ',ipart,trim(filename(ifile))//trim(suffix(isuffix))
						stop
					endif
				else if (idiff < 0) then
					if ((max(entropy(-1,ipart),omega(-1,ipart))-entropy(0,ipart))/width <= npoints) then
						entropy(-1,ipart) = max(entropy(-1,ipart),omega(-1,ipart))
						omega(0,ipart) = entropy(0,ipart)
						do i = npoints,1-idiff,-1
							omega(i,ipart) = omega(i+idiff,ipart)
						enddo
						do i = 1,1-idiff-1
							omega(i,ipart) = -iPrecise
						enddo
					else
						write (*,*) ' Width exceeded (high), ipart,file = ',ipart,trim(filename(ifile))//trim(suffix(isuffix))
						stop
					endif
				endif
			enddo	
		
			do ipart = min_part,max_part
				do i = -1,0
					if (entropy(i,ipart) == min_ene) then
						entropy(i,ipart) = omega(i,ipart)
					endif
				enddo
				do i = 1,npoints
					entropy(i,ipart) = specexp(entropy(i,ipart),omega(i,ipart))
				enddo
			enddo						

		endif

	enddo	! over isuffix (no ident)
	enddo	! over ifile

	auxmin = entropy(0,min_part)	! find true minimum energy
	do ipart = min_part+1,max_part
		auxmin = min(auxmin,entropy(0,ipart))
	enddo

	oldweight = weight	
	weight = -iPrecise
	success = .true.
	do ifile = 1,nfiles
		do isuffix = 1,nsuffix
			do ipart = min_part,max_part
				do i = 1,npoints
					zfactor =    + beta(1,1)*(entropy(0,ipart)+(i-1)*width) - &
								   beta(1,1)*mu(1,1)*ipart &
							-beta(ifile,isuffix)*(entropy(0,ipart)+(i-1)*width) + &
					beta(ifile,isuffix)*mu(ifile,isuffix)*ipart- oldweight(ifile,isuffix)
					if (entropy(i,ipart)>-iPrecise) &
						weight(ifile,isuffix)  = specexp(weight(ifile,isuffix),entropy(i,ipart)+ zfactor)
				enddo
			enddo
			if (abs(weight(ifile,isuffix)) > tolerance) then
				if (success) write (*,'(a)') '  ---------------------------'
				write (*,'(3a,f10.4)') '  Not Converged: ',trim(filename(ifile)),trim(suffix(isuffix)),weight(ifile,isuffix)
				success = .false.
			endif
			weight(ifile,isuffix) = accfactor*weight(ifile,isuffix) + oldweight(ifile,isuffix)
			if (.not.(ifile==1 .and. isuffix==1)) weight(ifile,isuffix) = weight(ifile,isuffix) - weight(1,1)
		enddo
	enddo
	weight(1,1) = 0.
	!write (*,'(a,8f8.3/<nsuffix*nfiles>f8.3)') ' New Weights=',((weight(ifile,isuffix),isuffix=1,nsuffix),ifile=1,nfiles)

	open (1,file='input_hs2.dat')
	write (1,*) nsuffix
	do isuffix = 1,nsuffix
		write (1,*) suffix(isuffix)
	enddo
	write (1,*) nfiles, npoints, maxpart
	do ifile = 1,nfiles
		do isuffix = 1,nsuffix
			write (1,'(A,g12.4,2f16.6,f13.6,f12.6,f7.1)') trim(filename(ifile)), counts_file(ifile,isuffix), &
				weight(ifile,isuffix),1./beta(ifile,isuffix),mu(ifile,isuffix),width, &
				g_n(ifile,isuffix)
		enddo
	enddo
	close (1)

enddo 

open (25,file='phase.dat')
write (25,'(A,20A5,3h *//3h /*,600A,3h *//3h /*,600A,3h */)') '/* Suffixes and Files = '&
	,(suffix(isuffix),isuffix=1,nsuffix),(trim(filename(ifile)),ifile=1,nfiles)
write (25,'(A)') '/*  T       mu       N_liq     N_gas     U_liq        U_gas     ln(Zliq)     ln(Zgas)'
open (28,file='pvt.dat')
write (28,'(A,20A5,3h *//3h /*,600A,3h *//3h /*,600A,3h */)') '/* Suffixes and Files = '&
	,(suffix(isuffix),isuffix=1,nsuffix),(trim(filename(ifile)),ifile=1,nfiles)
write (28,'(A)') '/*    T          mu           <N>           <U>          lnZ      ln(Zliq)   ln(Zgas)'
	
100 continue

if (read_s) then
	write (*,*) 'Enter T, Mu, Nmid, (lett. to stop, nmid<0 for phase eq.),smix, ncrit'
	read (*,*,err=200) beta1,mu1,nmid,smix,ncrit
else
	write (*,*) 'Enter T, Mu, Nmid, (lett. to stop, nmid<0 for phase eq.)'
	read (*,*,err=200) beta1,mu1,nmid
endif
if (ncrit>maxcrit) then
	write (*,*) '  Maxcrit < ncrit '
	stop
endif
beta1 = 1./beta1
do_phase = .false.
if (nmid<0) then
	do_phase = .true.
endif
nmid = abs(nmid)

write (*,'(A,f10.3)') ' Deviation = ', 100*funk(beta1,mu1,smix,ncrit,npoints,maxpart,nfiles,nsuffix, &
					min_part,max_part,entropy,width,nising,xising,yising,beta,mu)

ln_Z = -iPrecise
ln_Zliq = -iPrecise
ln_Zgas = -iPrecise	

do ipart = min_part,max_part
	do i = 1,npoints
		zfactor =    + beta(1,1)*(entropy(0,ipart)+(i-1)*width) - &
					   beta(1,1)*mu(1,1)*ipart &
					 - beta1*(entropy(0,ipart)+(i-1)*width)+beta1*mu1*ipart
		if (entropy(i,ipart)>-iPrecise) then
			ln_Z = specexp(ln_Z,entropy(i,ipart)+zfactor)
			if (ipart<= nmid) then
				ln_Zgas = specexp(ln_Zgas,entropy(i,ipart)+zfactor)
			else
				ln_Zliq = specexp(ln_Zliq,entropy(i,ipart)+zfactor)
			endif
		endif
	enddo
enddo

open (24,file='dens.dat')
write (24,'(A,f6.2,A,f8.2)') '/* T=',1/beta1,'; mu=',mu1
!write (24,'(A,f6.2,A,f6.3,A,20A5,3h *//3h /*,600A,3h *//3h /*,900A,3h */)') '/* T=',1/beta1,'; mu=',mu1,'; Suffixes and Files = '&
!	,(suffix(isuffix),isuffix=1,nsuffix),(trim(filename(ifile)),ifile=1,nfiles)
iter = 0
dens = 0.
avgnum = 0.
nbelow = 0.
ngas = 0.
nliq = 0.
ene_liq = 0.
do ipart = min_part,max_part
	do i = 1,npoints
		zfactor = + beta(1,1)*(entropy(0,ipart)+(i-1)*width) - &
					   beta(1,1)*mu(1,1)*ipart &
			 	  - beta1*(entropy(0,ipart)+(i-1)*width)+beta1*mu1*ipart
		ln_probab = -iPrecise
		if (entropy(i,ipart)>-iPrecise) ln_probab = specexp(ln_probab,entropy(i,ipart)+zfactor)
		dens(ipart) = dens(ipart) + exp(ln_probab-ln_Z)
		ene_liq = ene_liq + (entropy(0,ipart)+(i-1)*width)*exp(ln_probab-ln_Z)
	enddo
	avgnum = avgnum + ipart * dens(ipart)
	if (ipart<=nmid) then
		nbelow = nbelow + dens(ipart)
		ngas = ngas + ipart * dens(ipart)
	else
		nliq = nliq + ipart * dens(ipart)
	endif
	if (dens(ipart)>1e-30) write (24,'(i5,g12.3)') ipart,dens(ipart)
enddo

write (*,'(A,f7.3,A,f9.2,A,f4.0,A,f5.2,A,f8.2,A,f8.3,A,f8.3)') '  <N> =',avgnum,'  <E>/<N> =',ene_liq/avgnum,'; Frac. < ',nmid,'=',nbelow,'; nliq =',nliq/(1.-nbelow),'; ngas=',ngas/nbelow,' lnZ=',ln_Z
write (24,'(A,f7.3,A,f8.3,A,f4.0,A,f5.2,A,f9.3,A,f9.4,A,f9.4,3h */)') '/* <N> =',avgnum,'  <E>/<N> =',ene_liq/avgnum,'; Frac. < ',nmid,'=',nbelow,'; nliq =',nliq/(1.-nbelow),'; ngas=',ngas/nbelow,' lnZ=',ln_Z
close(24)
!write (29,'(A,f7.3,A,f8.3,A,f4.0,A,f5.3,A,f9.3,A,f9.4,A,f9.4,3h */)') '/* <N> =',avgnum,'  <E>/<N> =',ene_liq/avgnum,'; Frac. < ',nmid,'=',nbelow,'; nliq =',nliq/(1.-nbelow),'; ngas=',ngas/nbelow,' lnZ=',ln_Z
!close(29)

write (28,'(9g13.5)') 1/beta1,mu1,avgnum,ene_liq,ln_Z,ln_Zliq,ln_Zgas

open (24,file='crit.dat')
write (24,'(A,f8.4,A,f10.4)') '/* T=',1/beta1,'; mu=',mu1
write (24,'(A,f8.4,A,f8.4,A,20A5,3h *//3h /*,600A,3h *//3h /*,900A,3h */)') '/* T=',1/beta1,'; mu=',mu1,'; Suffixes and Files = '&
	,(suffix(isuffix),isuffix=1,nsuffix),(trim(filename(ifile)),ifile=1,nfiles),' */'
p_l = 0.
avg_x = 0.
min_x = min_part - smix*entropy(-1,min_part)	! we are looking for extreme values
max_x = min_x
do ipart = min_part,max_part
	do i = 1,npoints
		min_x = MIN(min_x, ipart - smix*(entropy(0,ipart)+(i-1)*width))
		max_x = MAX(max_x, ipart - smix*(entropy(0,ipart)+(i-1)*width))
	enddo
enddo
!min_x = min(min_x, min_part - smix*entropy(0,min_part))
!max_x = max(max_x, min_part - smix*entropy(0,min_part))
!min_x = min(min_x, min_part - smix*(entropy(0,min_part)+npoints*width))
!max_x = max(max_x, min_part - smix*(entropy(0,min_part)+npoints*width))
!min_x = min(min_x, max_part - smix*entropy(0,max_part))
!max_x = max(max_x, max_part - smix*entropy(0,max_part))
!min_x = min(min_x, max_part - smix*(entropy(0,min_part)+npoints*width))
!max_x = max(max_x, max_part - smix*(entropy(0,min_part)+npoints*width))

max_index = 0
do ipart = min_part,max_part
	do i = 1,npoints
		x = ipart - smix*(entropy(0,ipart)+(i-1)*width)
		zfactor = + beta(1,1)*(entropy(0,ipart)+(i-1)*width) - &
					beta(1,1)*mu(1,1)*ipart &
 				  - beta1*(entropy(0,ipart)+(i-1)*width)+beta1*mu1*ipart
		ln_probab = -iPrecise
		if (entropy(i,ipart)>-iPrecise) ln_probab = specexp(ln_probab,entropy(i,ipart)+zfactor)
		index = nint(ncrit*(x-min_x)/(max_x-min_x))
		max_index = max(index,max_index)
		p_l(index) = p_l(index) + exp(ln_probab-ln_Z)
	enddo
enddo

do i = 0,max_index
	avg_x = avg_x + p_l(i)*i
enddo

variance = 0
do i = 0,max_index
	variance = variance + p_l(i)*(i-avg_x)**2
enddo

do i = 0,max_index
	write (24,'(2g12.3)') (i-avg_x)/sqrt(variance), &
	                       p_l(i)*sqrt(variance)
enddo
write (24,'(A,2g14.3,A)') '/* <x> =',avg_x,variance,' */'
close(24)

if (do_phase) then	! no ident for this if
write (*,*) 'Enter Tmin, Tmax, Tincr (Tmin < Tmax; Tincr can be + or -)'
read (*,*) tmin,tmax,t_incr
if (t_incr>0) then
	beta1 = 1/tmin
else
	beta1 = 1/tmax
endif
iter = 0
do while (1./beta1<tmax+.00001.and.1./beta1>tmin-.00001.and.iter<=maxiter)
	iter = 0
	nbelow = 0.
	nmid = ngas+nliq
	do while (abs(0.5-nbelow)>0.0005.and.iter<=maxiter)
		if (iter/=0) then
			mu1 = mu1 + muincr	! arbitrary formula
			ln_Z = -iPrecise
			ln_Zliq = -iPrecise
			ln_Zgas = -iPrecise	
			do ipart = min_part,max_part
				do i = 1,npoints
					zfactor = + beta(1,1)*(entropy(0,ipart)+(i-1)*width) - &
								beta(1,1)*mu(1,1)*ipart &
			 			      - beta1*(entropy(0,ipart)+(i-1)*width)+beta1*mu1*ipart
					ln_probab = -iPrecise
					if (entropy(i,ipart)>-iPrecise) then
						ln_Z = specexp(ln_Z,entropy(i,ipart)+zfactor)
						if (ipart<= nmid) then
							ln_Zgas = specexp(ln_Zgas,entropy(i,ipart)+zfactor)
						else
							ln_Zliq = specexp(ln_Zliq,entropy(i,ipart)+zfactor)
						endif
					endif
				enddo
			enddo
			dens = 0.
			nbelow1 = 0.
			do ipart = min_part,max_part
				do i = 1,npoints
					zfactor = + beta(1,1)*(entropy(0,ipart)+(i-1)*width) - &
								beta(1,1)*mu(1,1)*ipart &
			 				  - beta1*(entropy(0,ipart)+(i-1)*width)+beta1*mu1*ipart
					ln_probab = -iPrecise
					if (entropy(i,ipart)>-iPrecise) ln_probab =  specexp(ln_probab,entropy(i,ipart)+zfactor)
					dens(ipart) = dens(ipart) + exp(ln_probab-ln_Z)
				enddo
				if (ipart<=nmid) then
					nbelow1 = nbelow1 + dens(ipart)
				endif
			enddo	
			! the function being set to zero is ln(x) - ln(1-x), which
			! behaves much better in convergence
			mu1 = mu1 - (log(nbelow1)-log(1.-nbelow1))*muincr &
				/ (log(nbelow1) - log(nbelow) - log(1.-nbelow1) + log(1.-nbelow))	! Newton's formula

		endif
		ln_Z = -iPrecise
		ln_Zliq = -iPrecise
		ln_Zgas = -iPrecise	
		do ipart = min_part,max_part
			do i = 1,npoints
				zfactor = + beta(1,1)*(entropy(0,ipart)+(i-1)*width) - &
							beta(1,1)*mu(1,1)*ipart &
			 			  - beta1*(entropy(0,ipart)+(i-1)*width)+beta1*mu1*ipart
				ln_probab = -iPrecise
				if (entropy(i,ipart)>-iPrecise) then
					ln_Z = specexp(ln_Z,entropy(i,ipart)+zfactor)
					if (ipart<= nmid) then
						ln_Zgas = specexp(ln_Zgas,entropy(i,ipart)+zfactor)
					else
						ln_Zliq = specexp(ln_Zliq,entropy(i,ipart)+zfactor)
					endif
				endif
			enddo
		enddo
		dens = 0.
		ngas = 0.
		nliq = 0.
		nbelow = 0.
		iter = iter + 1
		do ipart = min_part,max_part
			do i = 1,npoints
				zfactor = + beta(1,1)*(entropy(0,ipart)+(i-1)*width) - &
							beta(1,1)*mu(1,1)*ipart &
			 			  - beta1*(entropy(0,ipart)+(i-1)*width)+beta1*mu1*ipart
				ln_probab = -iPrecise
				if (entropy(i,ipart)>-iPrecise)	ln_probab = specexp(ln_probab,entropy(i,ipart)+zfactor)
				dens(ipart) = dens(ipart) + exp(ln_probab-ln_Z)
			enddo
			if (ipart<=nmid) then
				nbelow = nbelow + dens(ipart)
				ngas = ngas + ipart * dens(ipart)
			else
				nliq = nliq + ipart * dens(ipart)
			endif
		enddo
	enddo
	if (iter>maxiter.or.nbelow==0) then
		write (*,'(A,f7.3,A,I3,A,f7.2)') ' No convergence at T=',1/beta1,' Iter =',iter,'; Mu=',mu1
	else
		if (dens(int(nmid))>mindens*dens(int(nliq/(1.-nbelow))).or.dens(int(nmid))>mindens*dens(int(ngas/nbelow))) then
			write (*,'(A,f7.3,A,i4/A,3(I4,f6.3))') ' T=',1/beta1,' too close to Tc: peaks overlap at nmid = ',int(nmid),&
				' Nliq, Nmid, Ngas and peak heights:',int(nliq/(1.-nbelow)),dens(int(nliq/(1.-nbelow))), &
				int(nmid),dens(int(nmid)),int(ngas/nbelow),dens(int(ngas/nbelow))
		!	iter = maxiter+1
		else
			! calculate average energy
			ene_liq = 0.
			ene_gas = 0.
			do ipart = min_part,max_part
				do i = 1,npoints
					zfactor = + beta(1,1)*(entropy(0,ipart)+(i-1)*width) - &
								beta(1,1)*mu(1,1)*ipart &
			 				  - beta1*(entropy(0,ipart)+(i-1)*width)+beta1*mu1*ipart
					ln_probab = -iPrecise
					if (entropy(i,ipart)>-iPrecise) ln_probab = specexp(ln_probab,entropy(i,ipart)+zfactor)
					if (ipart<=nmid) then
						ene_gas = ene_gas + (entropy(0,ipart)+(i-1)*width)*exp(ln_probab-ln_Z)
					else
						ene_liq = ene_liq + (entropy(0,ipart)+(i-1)*width)*exp(ln_probab-ln_Z)				
					endif
				enddo
			enddo
			if (1/beta1 < 50) then				
				write (*,'(A,f7.4,A,I3,A,f10.4,A,f10.3,A,f10.4)') '  T=',1/beta1,' Iter =',iter,'; Mu=',mu1,'; nliq =',nliq/(1.-nbelow),'; ngas=',ngas/nbelow
				write (25,'(f8.5,F11.5,f10.3,f10.4,2g12.4,2f11.4)') 1/beta1,mu1, nliq/(1.-nbelow),ngas/nbelow, ene_liq/nliq,ene_gas/ngas,ln_Zliq,ln_Zgas
			else
				write (*,'(A,f7.2,A,I3,A,f10.4,A,f10.3,A,f10.4)') '  T=',1/beta1,' Iter =',iter,'; Mu=',mu1,'; nliq =',nliq/(1.-nbelow),'; ngas=',ngas/nbelow
				write (25,'(f7.2,F9.2,f10.3,f10.4,f11.3,3f11.4)') 1/beta1,mu1, nliq/(1.-nbelow),ngas/nbelow, ene_liq/nliq,ene_gas/ngas,ln_Zliq,ln_Zgas
			endif
		endif
	endif
	beta1 = 1./(1./beta1 + t_incr)
enddo
endif	! end of no ident if loop

goto 100
200 close(25)

stop

end

real function funk(beta1,mu1,smix,ncrit,npoints,maxpart,nfiles,nsuffix, &
					min_part,max_part,entropy,width,nising,xising,yising,beta,mu)

integer maxpart,npoints,nfiles,nsuffix
real beta1,mu1,smix,ncrit
real ln_Z					! partition function
integer min_part,max_part,ipart,i	! local counters
real entropy(-1:npoints,0:maxpart)	! contains counts for a certain energy+number
									! entropy (0,maxpart) is the starting (lowest) value of the energy
									! entropy (-1,maxpart) is the highest active value of the energy
integer, parameter :: maxcrit	= 5000 ! discretization of field mixing parameter
real p_l(0:maxcrit)			! ordering operator distribution
real x						! ordering operator distribution x axis
real max_x, min_x, avg_x	! limits of x, average x
integer nising 	! number of points in ising.dat file
integer ipoint	! for counting on ising data
real xising(nising),yising(nising)	! from ising.dat
real xpoint, ypoint, yvalue
real beta(nfiles,nsuffix),mu(nfiles,nsuffix),width	! temperature, chemical potential, width of energy distributions in each
							! histogram file
!common / c1 / min_part,max_part,entropy,width,xising,yising,beta,mu
! Function calculates the difference between the ising and the real distribution
integer inonzero	! number of actual points in distribution
real zfactor
real specexp,tmp1,auxreal
specexp(auxreal,tmp1) = max(auxreal,tmp1) + log(1.+exp(-abs(auxreal-tmp1)))
real ln_probab,variance	! auxiliary variables
integer index,max_index	! counter

real iPrecise 	! added by NAM
iPrecise = 1.0d15

ln_Z = -iPrecise	

do ipart = min_part,max_part
	do i = 1,npoints
		zfactor = + beta(1,1)*(entropy(0,ipart)+(i-1)*width) - &
					beta(1,1)*mu(1,1)*ipart &
 				  - beta1*(entropy(0,ipart)+(i-1)*width)+beta1*mu1*ipart
		if (entropy(i,ipart) > -iPrecise) ln_Z = specexp(ln_Z,entropy(i,ipart)+zfactor)
	enddo
enddo

p_l = 0.
avg_x = 0.
min_x = min_part - smix*entropy(-1,min_part)	! we are looking for extreme values
max_x = min_x
do ipart = min_part,max_part
	do i = 1,npoints
		min_x = MIN(min_x, ipart - smix*(entropy(0,ipart)+(i-1)*width))
		max_x = MAX(max_x, ipart - smix*(entropy(0,ipart)+(i-1)*width))
	enddo
enddo
!min_x = min(min_x, min_part - smix*entropy(0,min_part))
!max_x = max(max_x, min_part - smix*entropy(0,min_part))
!min_x = min(min_x, min_part - smix*(entropy(0,min_part)+npoints*width))
!max_x = max(max_x, min_part - smix*(entropy(0,min_part)+npoints*width))
!min_x = min(min_x, max_part - smix*entropy(0,max_part))
!max_x = max(max_x, max_part - smix*entropy(0,max_part))
!min_x = min(min_x, max_part - smix*(entropy(0,min_part)+npoints*width))
!max_x = max(max_x, max_part - smix*(entropy(0,min_part)+npoints*width))

max_index = 0
do ipart = min_part,max_part
	do i = 1,npoints
		x = ipart - smix*(entropy(0,ipart)+(i-1)*width)
		zfactor = + beta(1,1)*(entropy(0,ipart)+(i-1)*width) - &
					beta(1,1)*mu(1,1)*ipart &
 				  - beta1*(entropy(0,ipart)+(i-1)*width)+beta1*mu1*ipart
		ln_probab = -iPrecise
		if (entropy(i,ipart)>-iPrecise) ln_probab = specexp(ln_probab,entropy(i,ipart)+zfactor)
		index = nint(ncrit*(x-min_x)/(max_x-min_x))
		if (index.gt.maxcrit.or.index.lt.0) then
			write(*,*) 'Invalid index in funk',index,maxcrit
			stop
		endif
		max_index = max(index,max_index)
		p_l(index) = p_l(index) + exp(ln_probab-ln_Z)
	enddo
enddo

do i = 0,max_index
	avg_x = avg_x + p_l(i)*i
enddo

variance = 0
do i = 0,max_index
	variance = variance + p_l(i)*(i-avg_x)**2
enddo

funk = 0.
ipoint = 1
inonzero = 0
do i = 0,max_index
	xpoint = (i-avg_x)/sqrt(variance)
	ypoint = p_l(i)*sqrt(variance)
	do while (ipoint <= nising-2 .and. xising(ipoint+1) < xpoint)
		ipoint = ipoint + 1
	enddo
	yvalue = yising(ipoint) + &
	(yising(ipoint+1) - yising(ipoint))/(xising(ipoint+1)-xising(ipoint)) * &
	(xpoint - xising(ipoint))
	funk = funk + sqrt(ypoint)*abs(ypoint - yvalue)
	if (ypoint > 1e-6 .and. yvalue > 1e-6) inonzero = inonzero + 1
enddo

funk = funk/inonzero

return 
end
