program order

! version 1.1 - July '96 (AZP)
! program takes as input observations of histograms of N and E
! and orders them to maximize overlaps for the calculation of entropy.


implicit none

integer, parameter :: maxpart = 5000	! maximum number of particles
integer, parameter :: npoints = 5000 ! number of points for entropy discretization
									! for every particle number 
integer, parameter :: maxfiles = 1000	! maximum number of histogram files
integer, parameter :: maxsuffix = 9	! maximum number of run series 
real entropy(-1:0,0:maxpart)	! contains counts for a certain energy+number
									! entropy (0,maxpart) is the starting (lowest) value of the energy
									! entropy (-1,maxpart) is the highest active value of the energy
real tot_counts(npoints,0:maxpart)	! contains counts (observations) for a certain energy+number
									! reference for energy is the same as for entropy
!real, parameter :: min_entr = -999.	! nominal minimum entropy value (to avoid ln(0))
real, parameter :: min_ene = -99999.	! minimum energy
!real, parameter :: min_obs = 0.5	! how many observations is minimum to take them "sriously"
real omega(-1:0,0:maxpart),energy_start	! temporary holding variable, reads contr. to entropy from
												! each file being read
real counts(npoints,0:maxpart)	! counts for every file being read
integer ncount	! counts how many observations there are for a certain number of particles
integer min_part,max_part,ipart,i,idiff	! local counters
integer min_part1,max_part1		! min. and max. number of particles in file
real, parameter :: incr_width = 1.	! width for energy/N in counts file
real*8 beta,mu,width,width1	! temperature, chemical potential, width of energy distributions in each
							! histogram file; width1 makes sure that all widths are the same
character(20) filename(maxfiles)		! histogram file name - if the form his?????.dat; only the ????? part is stored
integer nfiles,ifile		! how many files being read
integer iostat				! iostat becomes -1 on end-of-file
integer noverlap(maxfiles),noverlap1,maxoverlap	! counters of overlaps between histograms 
character(2) suffix(maxsuffix)	! for dealing with series of runs
integer nsuffix,isuffix	! number of different run series
real xdim,ydim,zdim
integer filestart,maxindex


counts = 0
tot_counts = 0
entropy = 0.	! program looks at counts, so a zero value is not taken seriously
				! unless backed up by data
entropy(0,0:maxpart) = min_ene
entropy(-1,0:maxpart) = min_ene

open (1,file='input_hs.dat')
read (1,*) nsuffix
if (nsuffix>maxsuffix) then
	write (*,*) 'Nsuffix is greater than maxsuffix'
	stop
endif
do isuffix = 1,nsuffix
	read (1,*) suffix(isuffix)
enddo

nfiles = 1
iostat = 0
do while (iostat/=-1)
	read (1,*,iostat=iostat) filename(nfiles)
	nfiles=nfiles+1
enddo
nfiles = nfiles-2

if (nfiles*nsuffix>maxfiles) then
	write (*,*) 'Nfiles*nsuffix is greater than maxfiles'!, nfiles*nsuffix
	stop
endif
close (1)

min_part = maxpart+1	! so that these are later set to correct values
max_part = -1	! so that these are later set to correct values

ifile = 1
isuffix = 1	! we only order the first suffix
open (23,file='his'//trim(filename(ifile))//trim(suffix(isuffix))//'.dat')
read (23,*)
read (23,*) beta,mu,width1,xdim,ydim,zdim
min_part1 = maxpart + 1
max_part1 = - 1
beta = 1./beta
iostat = 0
omega = 0.
omega(0,0:maxpart) = min_ene
omega(-1,0:maxpart) = min_ene
counts = 0
open(26,file='cov'//trim(filename(ifile))//trim(suffix(isuffix))//'.dat')
do while (iostat/=-1)
	read (23,*,iostat=iostat) ipart,ncount,energy_start
	if (iostat/=-1) then
		if (ipart>maxpart.or.ncount>npoints) then
			write (*,*) 'Maxpart or Npoints exceeded, file, ipart, ncount = ',trim(filename(ifile))//trim(suffix(isuffix)),ipart,ncount
			stop
		endif
		omega(0,ipart) = energy_start
		omega(-1,ipart) = energy_start + (ncount-1)*width1
		write (26,'(i7,2f10.0)') ipart,omega(0,ipart),omega(-1,ipart)
		read (23,*) (counts(i,ipart),i=1,ncount)
		min_part = min(min_part,ipart)
		max_part = max(max_part,ipart)
		min_part1 = min(min_part1,ipart)
		max_part1 = max(max_part1,ipart)
	endif
enddo
close(23)
close(26)

open (24,file="input_hs.dat")
write (24,'(i1)') nsuffix
do isuffix = 1,nsuffix
	write (24,'(A)') trim(suffix(isuffix))
enddo
write (24,'(A)') trim(filename(ifile))

entropy = omega
tot_counts = counts
width = width1
isuffix = 1
filestart = 2

do while (filestart<nfiles)

	do ifile = filestart,nfiles

		open (23,file='his'//trim(filename(ifile))//trim(suffix(isuffix))//'.dat')
		read (23,*)
		read (23,*) beta,mu,width1,xdim,ydim,zdim
		min_part1 = maxpart + 1
		max_part1 = - 1
		beta = 1./beta
		iostat = 0
		omega = 0.
		omega(0,0:maxpart) = min_ene
		omega(-1,0:maxpart) = min_ene
		counts = 0
		open(26,file='cov'//trim(filename(ifile))//trim(suffix(isuffix))//'.dat')
		do while (iostat/=-1)
			read (23,*,iostat=iostat) ipart,ncount,energy_start
			if (iostat/=-1) then
				if (ipart>maxpart.or.ncount>npoints) then
					write (*,*) 'Maxpart or Npoints exceeded, file, ipart, ncount = ',trim(filename(ifile))//trim(suffix(isuffix)),ipart,ncount
					stop
				endif
				omega(0,ipart) = energy_start
				omega(-1,ipart) = energy_start + (ncount-1)*width1
				write (26,'(i7,2f10.0)') ipart,omega(0,ipart),omega(-1,ipart)
				read (23,*) (counts(i,ipart),i=1,ncount)
				min_part = min(min_part,ipart)
				max_part = max(max_part,ipart)
				min_part1 = min(min_part1,ipart)
				max_part1 = max(max_part1,ipart)
			endif
		enddo
		close(23)
		close(26)

		if (width1 /= width) then
			write (*,*) 'Width not consistent with previous, file = ',trim(filename(ifile))//trim(suffix(isuffix))
			stop
		endif
		noverlap(ifile) = 0
		do ipart = min_part,max_part
			idiff = 0
			if (omega(0,ipart)/=min_ene.and.entropy(0,ipart)/=min_ene) idiff = int (( entropy(0,ipart) - omega(0,ipart)  ) / width)
			do i = 1,npoints
				if (i-idiff>=1 .and. i-idiff<=npoints) then
					if (counts(i,ipart)>0.and.tot_counts(i-idiff,ipart)>0) then
						noverlap1 = 2./(1./counts(i,ipart)+1./tot_counts(i-idiff,ipart))
						noverlap(ifile) = noverlap(ifile) + noverlap1
					endif
				endif
			enddo
			if ( idiff > 0) then	! shift entropy
				if ((max(entropy(-1,ipart),omega(-1,ipart))-omega(0,ipart))/width <= npoints) then
					entropy(0,ipart) = omega(0,ipart)
					entropy(-1,ipart) = max(entropy(-1,ipart),omega(-1,ipart))
					do i = npoints,1+idiff,-1
						tot_counts(i,ipart) = tot_counts(i-idiff,ipart)
					enddo
					do i = 1,1+idiff-1
						tot_counts(i,ipart) = 0.
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
						counts(i,ipart) = counts(i+idiff,ipart)
					enddo
					do i = 1,1-idiff-1
						counts(i,ipart) = 0
					enddo
				else
					write (*,*) ' Width exceeded (high), ipart,file = ',ipart,trim(filename(ifile))//trim(suffix(isuffix))
					stop
				endif
			endif
		enddo	
	
!		write (*,'(3A,f6.2,A,f6.2,A,2i5,A,f7.3)') ' File: ',trim(filename(ifile))//trim(suffix(isuffix)),'; T=',1/beta,'; mu=',mu,'; min/max part.:',min_part1,max_part1,'; Overl (M):',real(noverlap(ifile))/1000000.
	enddo	! over ifile

	! order files with respect to overlaps

	maxoverlap = 0

	do ifile = filestart,nfiles 
		if (noverlap(ifile)>maxoverlap) then
			maxoverlap = noverlap(ifile)
			maxindex = ifile
		endif
	enddo

	if (maxoverlap > 0) then
		write (24,'(A)') trim(filename(maxindex))
		open (23,file='his'//trim(filename(maxindex))//trim(suffix(isuffix))//'.dat')
		read (23,*)
		read (23,*) beta,mu,width1,xdim,ydim,zdim
		min_part1 = maxpart + 1
		max_part1 = - 1
		beta = 1./beta
		iostat = 0
		omega = 0.
		omega(0,0:maxpart) = min_ene
		omega(-1,0:maxpart) = min_ene
		counts = 0
		open(26,file='cov'//trim(filename(maxindex))//trim(suffix(isuffix))//'.dat')
		do while (iostat/=-1)
			read (23,*,iostat=iostat) ipart,ncount,energy_start
			if (iostat/=-1) then
				omega(0,ipart) = energy_start
				omega(-1,ipart) = energy_start + (ncount-1)*width1
				write (26,'(i7,2f10.0)') ipart,omega(0,ipart),omega(-1,ipart)
				read (23,*) (counts(i,ipart),i=1,ncount)
				min_part = min(min_part,ipart)
				max_part = max(max_part,ipart)
				min_part1 = min(min_part1,ipart)
				max_part1 = max(max_part1,ipart)
			endif
		enddo
		close(23)
		close(26)
		do ipart = min_part,max_part
			idiff = 0
			if (omega(0,ipart)/=min_ene.and.entropy(0,ipart)/=min_ene) idiff = int (( entropy(0,ipart) - omega(0,ipart)  ) / width)
			if ( idiff > 0) then	! shift entropy
				if ((max(entropy(-1,ipart),omega(-1,ipart))-omega(0,ipart))/width <= npoints) then
					entropy(0,ipart) = omega(0,ipart)
					entropy(-1,ipart) = max(entropy(-1,ipart),omega(-1,ipart))
					do i = npoints,1+idiff,-1
						tot_counts(i,ipart) = tot_counts(i-idiff,ipart)
					enddo
					do i = 1,1+idiff-1
						tot_counts(i,ipart) = 0.
					enddo
				endif
			else if (idiff < 0) then
				if ((max(entropy(-1,ipart),omega(-1,ipart))-entropy(0,ipart))/width <= npoints) then
					entropy(-1,ipart) = max(entropy(-1,ipart),omega(-1,ipart))
					omega(0,ipart) = entropy(0,ipart)
					do i = npoints,1-idiff,-1
						counts(i,ipart) = counts(i+idiff,ipart)
					enddo
					do i = 1,1-idiff-1
						counts(i,ipart) = 0
					enddo
				endif
			endif
		enddo	
	
		write (*,'(3A,f6.2,A,f6.2,A,2i5,A,f7.3)') ' MAX: ',trim(filename(maxindex))//trim(suffix(isuffix)),'; T=',1/beta,'; mu=',mu,'; min/max part.:',min_part1,max_part1,'; Overl (M):',real(noverlap(maxindex))/1000000.
		do ipart = min_part,max_part
			do i = -1,0
				if (entropy(i,ipart) == min_ene) then
					entropy(i,ipart) = omega(i,ipart)
				endif
			enddo
			do i = 1,npoints
				if (tot_counts(i,ipart) == 0) then
					tot_counts(i,ipart) = counts(i,ipart)
				else 
					if (counts(i,ipart)>0) then
						tot_counts(i,ipart) = tot_counts(i,ipart) + counts(i,ipart)
					endif
				endif
			enddo
		enddo						
		filename(maxindex) = filename(filestart)
		filestart = filestart + 1

	else
		write (*,*) ' No overlap with any file ', filename(filestart)
		stop
	endif

enddo
write (24,'(A)') trim(filename(nfiles))

stop
end
