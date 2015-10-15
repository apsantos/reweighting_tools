Module data

integer, public :: maxpart,npoints,min_part,max_part,nfiles, nsuffix
real, public, allocatable :: entropy(:,:)   ! contains logs of counts for a certain energy+number
									! entropy (0,maxpart) is the starting (lowest) value of the energy
									! entropy (-1,maxpart) is the highest active value of the energy
character(20), public, allocatable :: filename(:)	! histogram file name - if the form his?????.dat; only the ????? part is stored
character(2), public, allocatable :: suffix(:)	! for dealing with series of runs
real, public :: width				! width of energy distributions in each histogram file
real smix					! s field mixing parameter
real, public, allocatable :: beta(:,:),mu(:,:) ! temperature, chemical potential
real, public :: dmu_dT	! Change of mu vs. T at rho = rho_c

end module


program entropy4

! Version 3.13 - April '14: Fixed +1 indexing error (NAM)
! Version 3.12 - Dec. '04: Formats
! Version 3.10 - Nov. '04: based on v. 2.37 of entropy3; 3D only
! Version 3.00 - Dec. '04: based on v. 2.37 of entropy3; 3D only
! program takes as input observations of histograms of N and E
! and optimizes the critical temperature, chemical potential and field mixing parameters
! works only for 3-D systems, using the Tsypin and Bloete [PRE 62:73, 2000] function
!
! Program takes as input observations of histograms of N and E
! and weights calculated by the
! Ferrenberg and Swedsen algorithm (PRL, 1989) (from program entropy2)
! and optimizes the Critical Temperature, Critical Chem. Pot., 
! field mixing parameter S and the number of bins.

use data

implicit none
real, parameter :: min_ene = -999999.	! minimum energy
real, allocatable :: omega(:,:)   ! temporary holding variable, reads contr. to entropy from
								! each file being read
real energy_start

real, allocatable :: counts(:)		! counts for every file being read per (E,N) point
real, allocatable :: counts_file(:,:)	! total counts for each file
real, allocatable :: g_n(:,:)		! related to "correlation time" for each file - g_n in S&W paper
								! if g_n=1 every configuration is considered independent, 
								! greater values mean correlated configurations.
integer ncount	! counts how many observations there are for a certain number of particles
integer min_part1,max_part1		! min. and max. number of particles in file
real, parameter :: incr_width = 0.1	! width for energy/N in counts file
real beta1,mu1			! temperature, chem. pot., field mixing
integer, parameter :: maxcrit	= 150000 ! discretization of field mixing parameter
real, allocatable :: weight(:,:)		! f_n in Ferr. + Swed. paper
real ln_num_P				! Numerator for P(S,K) in eq. (3) of Ferr. + Swed. paper
integer ifile,jfile,jsuffix	! how many files being read
integer iostat				! iostat becomes -1 on end-of-file
real specexp,tmp1,auxreal
specexp(auxreal,tmp1) = max(auxreal,tmp1) + log(1.+exp(-abs(auxreal-tmp1))) ! special exponential in-line function
real, allocatable :: dens(:) ! dens stores the probability of a certain number
integer iter,i,ipart,idiff
integer isuffix	! number of different run series
real width1, xdim1, ydim1, zdim1	! these are only read, not used anywhere
integer, parameter :: nopt_dim = 3	! number of optimizing directions:
									! T_c, mu0, smix
real ptry(nopt_dim)
real p(nopt_dim+1,nopt_dim)		! input to optimization function
real y(nopt_dim+1)				! input to optimization function
real ftol ! tolerance for convergence of optimization
real trial_vec(nopt_dim)		! factors for initial search in Tcrit
real funk
external funk
real zfactor	! auxiliary variable in calculation of weight and Z
character(len=32) :: nfile_fmt

open (1,file='entr4_par.dat')
read (1,*) (trial_vec(i),i=1,nopt_dim)
read (1,*) ftol, dmu_dT, smix
close (1)

open (1,file='input_hs2.dat')
read (1,*) nsuffix
allocate (suffix(nsuffix))
do isuffix = 1,nsuffix
	read (1,*) suffix(isuffix)
enddo
read (1,*) nfiles, npoints, maxpart
close(1)

allocate (entropy(-1:npoints,0:maxpart),omega(-1:npoints,0:maxpart))
allocate (counts(npoints),counts_file(nfiles,nsuffix),g_n(nfiles,nsuffix))
allocate (beta(nfiles,nsuffix),mu(nfiles,nsuffix))
allocate (filename(nfiles),weight(nfiles,nsuffix))
allocate (dens(0:maxpart))

write (*,'(a)') '  -----------------------------------------------------------------'
write (*,'(a)') '  3D critical parameters - version 3.13 (NAM, April 2014)'
write (*,'(a)') '  -----------------------------------------------------------------'

!entropy = min_entr
counts = 0
entropy = -1.d6	! program looks at counts, so a zero value is not taken seriously
				! unless backed up by data -1.d6 corresponds to a small count
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
	omega = -1.d6
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
			omega(-1,ipart) = energy_start + ncount*width
			read (23,*) (counts(i),i=1,ncount)
			do i = 1,ncount
				if (counts(i)>0) then
					ln_num_P = -1.d6
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
						entropy(i,ipart) = -1.d6
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
						omega(i,ipart) = -1.d6
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


write(nfile_fmt, '(a, i0, a, i0, a)') '(A,', nsuffix+nfiles, 'A,A,8f8.3/',nsuffix*nfiles,'f8.3)'
!write (*,'(a,<nsuffix+nfiles>a,a,8f8.3/<nsuffix*nfiles>f8.3)') '  Weights (optimized in entropy2!): ', &
write (*,nfile_fmt) '  Weights (optimized in entropy2!): ', &
	( trim(filename(ifile)),ifile=1,nfiles ), ( trim(suffix(isuffix)),isuffix=1,nsuffix),&
	' =',((weight(ifile,isuffix),isuffix=1,nsuffix),ifile=1,nfiles)

open (1,file='input_hs2.dat')
write (1,*) nsuffix
do isuffix = 1,nsuffix
	write (1,*) suffix(isuffix)
enddo
write (1,*) nfiles, npoints, maxpart
do ifile = 1,nfiles
	do isuffix = 1,nsuffix
		write (1,'(A,g12.4,2f12.6,f13.6,f12.6,f7.1)') trim(filename(ifile)), counts_file(ifile,isuffix), &
			weight(ifile,isuffix),1./beta(ifile,isuffix),mu(ifile,isuffix),width, &
			g_n(ifile,isuffix)
	enddo
enddo
close (1)

write (*,*) '--------------------------------------------------------------------------'
write (*,*) ' Iter       T          Mu        <N>        Smix          Dev.'
write (*,*) '--------------------------------------------------------------------------'

p(1,1) = beta(1,1)
p(1,2) = mu(1,1)
p(1,3) = smix
ptry(1:nopt_dim) = p(1,1:nopt_dim)
y(1) = funk(ptry)
p(2,1) = beta(1,1) - trial_vec(1)*beta(1,1)**2
p(2,2) = mu(1,1)
p(2,3) = smix
ptry(1:nopt_dim) = p(2,1:nopt_dim)
y(2) = funk(ptry)
p(3,1) = beta(1,1)
p(3,2) = mu(1,1)  + trial_vec(2)
p(3,3) = smix
ptry(1:nopt_dim) = p(3,1:nopt_dim)
y(3) = funk(ptry)
p(4,1) = beta(1,1)
p(4,2) = mu(1,1)  
p(4,3) = smix + trial_vec(3)
ptry(1:nopt_dim) = p(4,1:nopt_dim)
y(4) = funk(ptry)

call amoeba(p,y,nopt_dim+1,nopt_dim,nopt_dim,ftol,funk,iter)

write (*,'(A,f12.5,f12.5,9x,f13.6)') '   +-', 1./minval(p(1:nopt_dim+1,1)) - 1./maxval(p(1:nopt_dim+1,1)),&
	maxval(p(1:nopt_dim+1,2))-minval(p(1:nopt_dim+1,2)),maxval(p(1:nopt_dim+1,3))-minval(p(1:nopt_dim+1,3))

stop

end

real function funk(ptry)
! Function calculates the difference between the ising and the real distribution
! Ising distribution from Tsypin and Bloete [PRE 62:73, 2000] 
! Objective function is the sum of sqrt(ypoint)*abs(ypoint - yvalue) divided by the number of points.
! The value reported from entropy4 is 1000x the value of the objective function.

use data
real beta1,mu1
real ln_Z					! partition function
integer ipart,i	! local counters
real, allocatable:: p_l(:)	! ordering operator distribution
real x						! ordering operator distribution x axis
real max_x, min_x, avg_x	! limits of x, average x
real xpoint, ypoint, yvalue
real zfactor
real specexp,tmp1,auxreal
specexp(auxreal,tmp1) = max(auxreal,tmp1) + log(1.+exp(-abs(auxreal-tmp1)))
real ln_probab,variance	! auxiliary variables
integer index,isuffix,ifile	! counters
integer ncrit	! set to max_part - min_part + 1
real dx			! increment in  order parameter x
real n_pos		! counts how many real points we had
real fvalue		! value of function
integer, save:: count
real ptry(3)
real avg_N

ncrit = max_part - min_part + 1
allocate (p_l(0:ncrit))

ln_Z = -1.d6

beta1 = ptry(1)
mu1 = ptry(2) + dmu_dT*(1./beta1 - 1./beta(1,1))
smix = ptry(3)
	
do ipart = min_part,max_part
	do i = 1,npoints
		zfactor = + beta(1,1)*(entropy(0,ipart)+(i-1)*width) - &
					beta(1,1)*mu(1,1)*ipart &
 				  - beta1*(entropy(0,ipart)+(i-1)*width)+beta1*mu1*ipart
		if (entropy(i,ipart) > -1.d6) ln_Z = specexp(ln_Z,entropy(i,ipart)+zfactor)
	enddo
enddo

open (24,file='crit.dat')
write (24,'(A,f8.4,A,g14.8,A,g10.4,A,20A5,3h *//3h /*,50A,3h *//3h /*,50A,3h */)') '/* T=' &
             ,1/beta1,'; mu=',mu1,'; smix=',smix,'; Suffixes and Files = '&
             ,(suffix(isuffix),isuffix=1,nsuffix),(trim(filename(ifile)),ifile=1,nfiles),' */'

p_l = 0.
avg_x = 0.
avg_N = 0.
min_x = min_part - smix*entropy(-1,min_part)	! we are looking for extreme values
max_x = min_x
min_x = min(min_x, min_part - smix*entropy(0,min_part))
max_x = max(max_x, min_part - smix*entropy(0,min_part))
min_x = min(min_x, min_part - smix*(entropy(0,min_part)+npoints*width))
max_x = max(max_x, min_part - smix*(entropy(0,min_part)+npoints*width))
min_x = min(min_x, max_part - smix*entropy(0,max_part))
max_x = max(max_x, max_part - smix*entropy(0,max_part))
min_x = min(min_x, max_part - smix*(entropy(0,max_part)+npoints*width))
max_x = max(max_x, max_part - smix*(entropy(0,max_part)+npoints*width))
dx = (max_x - min_x) / ncrit
do ipart = min_part,max_part
	do i = 1,npoints
		x = ipart - smix*(entropy(0,ipart)+(i-1)*width)
		zfactor = + beta(1,1)*(entropy(0,ipart)+(i-1)*width) - &
					beta(1,1)*mu(1,1)*ipart &
 				  - beta1*(entropy(0,ipart)+(i-1)*width)+beta1*mu1*ipart
		ln_probab = -1.d6
		if (entropy(i,ipart)>-1.d6) ln_probab = specexp(ln_probab,entropy(i,ipart)+zfactor)
		index = nint((x-min_x)/dx)
		if (index<=ncrit .and. index >= 0) then
			p_l(index) = p_l(index) + exp(ln_probab-ln_Z)
			avg_N = avg_N + ipart*exp(ln_probab-ln_Z)
		endif
	enddo
enddo

do i = 0,ncrit
	avg_x = avg_x + p_l(i)*i
enddo

variance = 0
do i = 0,ncrit
	variance = variance + p_l(i)*(i-avg_x)**2
enddo

fvalue = 0.
n_pos = 1e-6
do i = 0,ncrit
	xpoint = (i-avg_x)/sqrt(variance)
	ypoint = p_l(i)*sqrt(variance)
	yvalue = exp(-(0.7774*xpoint*xpoint-1)**2*(0.158*0.7774*xpoint*xpoint+0.776))/2.3306
	if (yvalue > 1e-4) n_pos = n_pos + 1.
	fvalue = fvalue + sqrt(ypoint)*abs(ypoint - yvalue)
	write(24,'(g9.3,2g12.3)') xpoint,ypoint,yvalue
	if (ypoint > 1e-6 .and. yvalue > 1e-6) then
		n_pos = n_pos + 1.
	endif
enddo
fvalue = fvalue/n_pos
write (24,'(A,3g14.3,A)') '/* n_pos , var, deviation =',n_pos,variance,fvalue,' */'
close(24)

funk = fvalue
count = count + 1
write (*,'(i5,2f12.5,f10.3,f13.7,1x,f10.3)') count,1./beta1,mu1,avg_N,smix,100.*fvalue

return 
end


      FUNCTION amotry(p,y,psum,mp,np,ndim,funk,ihi,fac)
      INTEGER ihi,mp,ndim,np,NMAX
      REAL amotry,fac,p(mp,np),psum(np),y(mp),funk
      PARAMETER (NMAX=3)
      EXTERNAL funk
!    USES funk
      INTEGER j
      REAL fac1,fac2,ytry,ptry(NMAX)
	  integer npoints, maxpart,nfiles,nsuffix,min_part,max_part
	  real entropy,width,beta,mu

      fac1=(1.-fac)/ndim
      fac2=fac1-fac
      do 11 j=1,ndim
        ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
11    continue
      ytry=funk(ptry,npoints,maxpart,nfiles,nsuffix,&
					min_part,max_part,entropy,width,beta,mu)
      if (ytry.lt.y(ihi)) then
        y(ihi)=ytry
        do 12 j=1,ndim
          psum(j)=psum(j)-p(ihi,j)+ptry(j)
          p(ihi,j)=ptry(j)

12      continue
      endif
      amotry=ytry
!	  write(*,'(2g13.5,g13.6,f10.2)') ytry*100,ptry(1),ptry(2),ptry(3)
      return
      END

SUBROUTINE amoeba(p,y,mp,np,ndim,ftol,funk,iter)
INTEGER iter,mp,ndim,np,NMAX,ITMAX
REAL ftol,p(mp,np),y(mp),funk
PARAMETER (NMAX=3,ITMAX=5000)
EXTERNAL funk

!     USES amotry,funk
      INTEGER i,ihi,ilo,inhi,j,m,n
      REAL rtol,sum,swap,ysave,ytry,psum(NMAX),amotry
      iter=0
1     do 12 n=1,ndim
        sum=0.
        do 11 m=1,ndim+1
          sum=sum+p(m,n)
11      continue
        psum(n)=sum
12    continue
2     ilo=1
      if (y(1).gt.y(2)) then
        ihi=1
        inhi=2
      else
        ihi=2

        inhi=1
      endif
      do 13 i=1,ndim+1
        if(y(i).le.y(ilo)) ilo=i
        if(y(i).gt.y(ihi)) then
          inhi=ihi
          ihi=i
        else if(y(i).gt.y(inhi)) then
          if(i.ne.ihi) inhi=i
        endif
13    continue
      rtol=2.*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))
      if (rtol.lt.ftol) then
        swap=y(1)
        y(1)=y(ilo)
        y(ilo)=swap
        do 14 n=1,ndim
          swap=p(1,n)
          p(1,n)=p(ilo,n)
          p(ilo,n)=swap
14      continue

        return
      endif
      if (iter.ge.ITMAX) then
        write(*,*) 'ITMAX exceeded in amoeba, type [enter] to continue'
        read (*,*) 
      endif
      iter=iter+2
      ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,-1.0)
      if (ytry.le.y(ilo)) then
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,2.0)
      else if (ytry.ge.y(inhi)) then
        ysave=y(ihi)
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,0.5)
        if (ytry.ge.ysave) then
          do 16 i=1,ndim+1
            if(i.ne.ilo)then
              do 15 j=1,ndim
                psum(j)=0.5*(p(i,j)+p(ilo,j))
                p(i,j)=psum(j)
15            continue
              y(i)=funk(psum)
            endif
16        continue
          iter=iter+ndim
          goto 1

        endif
      else
        iter=iter-1
      endif
      goto 2
      END
