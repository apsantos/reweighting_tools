
program fspatch

! v. 2.11  6/24/08 - formats (from Antti-Pekka's version)
! v. 2.1  7/12/99 - weight = log(weight) to avoid overflows
! v. 2.0  7/5/99 - patching of new histograms done incrementally
! v. 1.4  6/26/99 - all arrays allocatable, different names/format of I/O files (AZP)

! Initial code written by Jeff Potoff
! version 1.0 - 10/31/97  This program iterates to find the weights for
! patching the individual histograms together.  The histogram data is read
! in the form of N1, N2, E for each configuration.
! 2/11/98 - Histogram patching rewritten to avoid overflows 
! 4/29/98 - New data type used to cut memory usage


!Notes: The low density reference state must be the first list read
! When attemping to patch new data into previously patch runs, the 
! new data must go at the *end* of the list of lists in 'inp_fs.dat'


!USE PORTLIB    
!use IFPORT ! for hostnam

implicit none

integer mentry    ! maximum number of observations
real tolerance    ! tolerance for convergence of weights
real ksee        ! acceleration factor
 
integer i    ! local counter
integer, allocatable :: min_part1(:),max_part1(:)    ! min. and max. number of particles in file
integer, allocatable :: min_part2(:),max_part2(:)    ! min. and max. number of particles in file
integer nfiles,ifile,k        ! how many files being read
integer iostat            ! iostat becomes -1 on end-of-file
                                
real*8 xdim,ydim,zdim    ! linear dimensions of box on which data were taken
real*8 xdim1,ydim1,zdim1    ! linear dimensions of box for current file

! variables for new histogram format
integer jfile,kfile,nfiles_old,iter,iter_old
integer, allocatable, target :: n1_temp(:)
integer, allocatable, target :: n2_temp(:)
real*8, allocatable :: e_temp(:)
real*8 prob,maxd,maxd_old,val,y,ylog,tmp1
character*20, dimension(:), allocatable :: filename
real*8, dimension(:), allocatable :: beta
real*8, dimension(:), allocatable :: mu1
real*8, dimension(:), allocatable :: mu2
real*8, dimension(:), allocatable :: nentry
real*8, dimension(:), allocatable :: weight
real*8, dimension(:), allocatable :: oldweight

real trial_vec(3)        ! trial_vec for fscrm

!define data type - this is to required to reduce memory usage
type column1
  integer, pointer ::flength(:) ! Number of entries in the file
end type

type column2
  real*8, pointer :: flength2(:) ! Number of entries in the file
end type

type(column1), dimension(:), allocatable :: n1,n2
type(column2), dimension(:), allocatable :: e

logical oldrun

integer                    istat                ! aux variable
!character*20               hostname
character*20               auxfile
real(4)                    itim, ttt, tarr(2)    ! timing stuff
real*8 beta_tmp,mu1_tmp,mu2_tmp ! temporary variables
integer side1,side2,oldfile,overlap    ! counters
integer maxoverlap    ! file with which a new file has the maximum 

!-------------Executable part----------------------------

istat = 0
iostat = 0
nfiles = 0
nfiles_old = 0
ifile = 0
jfile = 0
kfile = 0
mentry = 0
i = 0
k = 0
side1 = 0
side2 = 0
oldfile = 0
overlap = 0
maxoverlap = 0
tolerance = 0.0
ksee = 0.0
xdim = 0.0
xdim1 = 0.0
ydim = 0.0
ydim1 = 0.0
zdim = 0.0
zdim1 = 0.0
prob = 0.0
maxd = 0.0
val = 0.0
y = 0.0
ylog = 0.0
tmp1 = 0.0
itim = 0.0
ttt = 0.0
beta_tmp = 0.0
mu1_tmp = 0.0
mu2_tmp = 0.0

!ttt = timef()  ! total elapsed time
!itim = dtime(tarr)  ! CPU time

write (*,*) '--------------------------------------------------------------------------'
write (*,*) ' Ferrenberg-Swedsen histogram patching, v. 2.11, July 08 (AZP)'
write (*,*) '--------------------------------------------------------------------------'

open (1,file='inp_fs.dat')
read (1,*) 
do while (iostat/=-1)
    read (1,*,iostat=iostat) auxfile
    nfiles=nfiles+1
enddo
nfiles = nfiles-1
close(1)

write(*,*) nfiles

oldrun = .false.
open(2,file='inp_fs2.dat',status='old',iostat=iostat)
if (iostat==0) then
    oldrun = .true.
    read(2,*) nfiles_old
    nfiles = nfiles_old + nfiles
endif
close(2)

open (1,file='inp_fs.dat')
read (1,*) ksee, tolerance, mentry, trial_vec(1:3)
allocate(n1_temp(mentry),n2_temp(mentry),e_temp(mentry))
allocate(min_part1(nfiles),max_part1(nfiles),min_part2(nfiles),max_part2(nfiles))
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
n1_temp(:) = 0
n2_temp(:) = 0
e_temp(:) = 0.0
beta(:) = 0.0
mu1(:) = 0.0
mu2(:) = 0.0
nentry(:) = 0.0
weight(:) = 0.0
oldweight(:) = 0.0

filename(:) = ""
do ifile = nfiles_old+1,nfiles
    read (1,*) filename(ifile)
enddo
close (1)

iter = 0
iter_old = 0
maxd_old = 100
open(2,file='inp_fs2.dat',status='old',iostat=iostat)
if (iostat==0) then
    read(2,*) 
    do ifile = 1,nfiles_old
        read(2,*) filename(ifile),nentry(ifile),oldweight(ifile),beta(ifile),&
                mu1(ifile),mu2(ifile)
        beta(ifile) = 1./beta(ifile)
    enddo
    read(2,*,iostat=iostat) iter_old,maxd_old
endif
close(2)

do ifile = 1,nfiles_old
    do jfile = nfiles_old+1,nfiles
        if (trim(filename(ifile))==trim(filename(jfile))) then
            write(*,'(3A)') ' File ',trim(filename(jfile)),' already in inp_fs2'
            stop
        endif
    enddo
enddo

write(*,'(A)') '         T          mu1        mu2    minp1   minp2  maxp1  maxp2  npoints'
                                        
do kfile= 1,nfiles
    min_part1(kfile) = 1000000    ! so that these are later set to correct values
    max_part1(kfile) = -1    
    min_part2(kfile) = 1000000
    max_part2(kfile) = -1
        
    open (23,file='his'//trim(filename(kfile))//'.dat')
    read (23,*)
    read (23,*) beta_tmp,mu1_tmp,mu2_tmp,xdim1,ydim1,zdim1
    beta_tmp = 1./beta_tmp
    if (kfile<=nfiles_old) then
        if (abs(beta_tmp-beta(kfile))+abs(mu1_tmp-mu1(kfile))+abs(mu2_tmp-mu2(kfile)) > 1.e-5) then
            write(*,*) 'Histogram file ',trim(filename(kfile)),' inconsistent with inp_fs2.dat'
            stop
        endif
    else
        mu1(kfile)=mu1_tmp
        mu2(kfile)=mu2_tmp
        beta(kfile) = beta_tmp
    endif
    if (kfile==1) then
        xdim = xdim1
        ydim = ydim1
        zdim = zdim1
    else
        if (xdim1/=xdim.or.ydim1/=ydim.or.zdim1/=zdim) then
            write (*,*) 'System size in file ',trim(filename(kfile)),' inconsistent with first file:'
            write (*,'(A,3F8.2,A,3F8.2,A)') '(',xdim, ydim, zdim,') /= (',xdim1, ydim1, zdim1,')'
            stop
        endif
    endif
    iostat = 0 
    i = 0
    do while (iostat/=-1) 
        i = i + 1    
        if (i>mentry) then
            write(*,*) '  Maximum number of entries exceeded in file his'//trim(filename(kfile))//'.dat'
            stop
        endif
        read (23,*,iostat=iostat) n1_temp(i),n2_temp(i),e_temp(i)  
        if(iostat.ne.-1) then
            min_part1(kfile) = min(n1_temp(i),min_part1(kfile))
            max_part1(kfile) = max(n1_temp(i),max_part1(kfile))
            min_part2(kfile) = min(n2_temp(i),min_part2(kfile)) 
            max_part2(kfile) = max(n2_temp(i),max_part2(kfile))
        endif
    enddo
    nentry(kfile) = i - 1

    allocate(n1(kfile)%flength(1:i-1))
    allocate(n2(kfile)%flength(1:i-1))
    allocate(e(kfile)%flength2(1:i-1))
    
    n1(kfile)%flength(1:i-1) = n1_temp(1:nentry(kfile))
    n2(kfile)%flength(1:i-1) = n2_temp(1:nentry(kfile))
    e(kfile)%flength2(1:i-1) = e_temp(1:nentry(kfile))
    close(23)        
    write(*,'(1x,A,1x,f9.5,2x,f10.5,2x,f10.5,1x,i5,2x,i5,2x,i5,2x,i5,1x,f10.0)') &
        trim(filename(kfile)), &
        1./beta(kfile),mu1(kfile),mu2(kfile),min_part1(kfile),min_part2(kfile),&
        max_part1(kfile),max_part2(kfile),nentry(kfile)
    weight(1) = 0.
    oldweight(1) = 0.
    if((oldrun==.false. .or. kfile>nfiles_old) .and. kfile > 1) then
        ! find decent quess for weight of new file by patching it with previous file
        ! of maximum overlap
        overlap = 0
        do oldfile = 1,kfile-1
            side1 = min(max_part1(kfile),max_part1(oldfile)) - max(min_part1(kfile),min_part1(oldfile))
            side2 = min(max_part2(kfile),max_part2(oldfile)) - max(min_part2(kfile),min_part2(oldfile))
            if (max_part2(kfile)==0.and.min_part2(kfile)==0) side2=1
            if (side1<=0.or.side2<=0) side2 = 0
            if (side1*side2 > overlap) then
                overlap = side1*side2
                maxoverlap = oldfile
            endif
        enddo
        weight(kfile) = 0.
        if (overlap==0)    then
            write(*,'(3A,I4,A,I4,A)') 'Warning: no overlap for file ',trim(filename(kfile)), ' side1=',side1,' side2=',side2
            oldweight(kfile) = 0.
        else
            oldweight(kfile) = 0.
            iter = 0
            maxd = 100
            do while(maxd.gt.tolerance.and.iter.lt.1e4)
                iter = iter +1
                maxd = 0.0  
                weight(kfile) = -1e9
                do ifile = maxoverlap,kfile,kfile-maxoverlap    ! just does it twice      
                    do i=1,nentry(ifile)
                        y = 0.0
                        ylog = -1e9  ! exp(ylog) = y = 0.0
                        do jfile = maxoverlap,kfile,kfile-maxoverlap     
                            prob = -(beta(jfile)-beta(kfile))*e(ifile)%flength2(i) +&
                                (beta(jfile)*mu1(jfile)-beta(kfile)*mu1(kfile))*n1(ifile)%flength(i)+&
                                (beta(jfile)*mu2(jfile)-beta(kfile)*mu2(kfile))*n2(ifile)%flength(i)
                            tmp1 = prob + log(nentry(jfile)) - oldweight(jfile)      
                            ylog = max(ylog,tmp1) + log(1.+exp(-abs(ylog-tmp1)))
                        enddo
                        weight(kfile) = max(weight(kfile),-ylog) + log(1.+exp(-abs(weight(kfile)+ylog)))
                    enddo
                enddo
                maxd = abs(weight(kfile)-oldweight(kfile))
                if (iter>1) weight(kfile) = ksee*weight(kfile)+(1-ksee)*oldweight(kfile)
                oldweight(kfile)=weight(kfile)
                if (iter < 5 .or. mod(iter,10)==0) write(*,*) iter,weight(kfile),maxd
            enddo
            write(*,'(5a,i6,a,g11.3)') ' Patched ',trim(filename(kfile)),' with ',trim(filename(maxoverlap)),&
                        ', overlap =',overlap,', weight=',weight(kfile)

        endif    
    endif
enddo
! iterate to find the weights
if(nfiles==nfiles_old) then 
    maxd = maxd_old
    iter = iter_old
else
    maxd = 100
    iter = 0
endif

do while(maxd.gt.tolerance.and.iter.lt.1e4)
    if (mod(iter,1).eq.0) then
!        write(*,'(A,i5,5x,A,f12.8)') ' iteration = ',iter,'deviation = ',maxd
        open(24,file='inp_fs2.dat')
        write(24,*) nfiles, trial_vec(1:3)
        do ifile=1,nfiles
            write(24,'(A12,f10.0,g15.6,3f14.5)') filename(ifile),nentry(ifile),&
                oldweight(ifile),1./beta(ifile),mu1(ifile),mu2(ifile)
        enddo
        write(24,'(i5,f9.4,A)') iter,maxd,'  ! total iterations, deviation '
        close(24)
    end if    
    iter = iter +1
    maxd = 0.0  
    do kfile = 1,nfiles
        weight(kfile) = -1e9
        do ifile = 1,nfiles      
            do i=1,nentry(ifile)
                y = 0.0
                ylog = -1e9  ! exp(ylog) = y = 0.0
                do jfile = 1,nfiles      
                    prob = -(beta(jfile)-beta(kfile))*e(ifile)%flength2(i) +&
                        (beta(jfile)*mu1(jfile)-beta(kfile)*mu1(kfile))*n1(ifile)%flength(i)+&
                        (beta(jfile)*mu2(jfile)-beta(kfile)*mu2(kfile))*n2(ifile)%flength(i)
                    ! y = y+nentry(jfile)*exp(prob-log(oldweight(jfile)))
                    ! The following is a code rearrangement to avoid overflow problems
                    ! This trick was compliments of Peter-Lawrence Montgomery 2/12/98
                    ! Replace y with log(y)
                    ! Define log_expsum(y,z) = log(exp(y)+exp(z))
                    ! Evaluate log_expsum(y,z) = max(y,z) + log(1+exp(-abs(y-z)))
                    tmp1 = prob + log(nentry(jfile)) - oldweight(jfile)
                    ylog = max(ylog,tmp1) + log(1.+exp(-abs(ylog-tmp1)))
                enddo
                weight(kfile) = max(weight(kfile),-ylog) + log(1.+exp(-abs(weight(kfile)+ylog)))
            enddo
        enddo
        val = abs(weight(kfile)-oldweight(kfile))
        if (val.gt.maxd) maxd=val
    enddo
    do i = 2,nfiles
        if (iter==1) then
            weight(i) = weight(i)-weight(1)
        else
            weight(i) = ksee*(weight(i)-weight(1))+(1-ksee)*oldweight(i)
        endif
    enddo
    weight(1) = 0.
    oldweight(1:nfiles)=weight(1:nfiles)
enddo

write(*,'(A,<nfiles>g11.3)') ' Final weights: ', weight(1:nfiles)
open(unit=24, file='inp_fs2.dat',status='unknown')
write(24,*) nfiles, trial_vec(1:3)
do ifile=1,nfiles
    write(24,'(A5,f10.0,g15.6,3f14.5)') filename(ifile),nentry(ifile),weight(ifile),&
    1./beta(ifile),mu1(ifile),mu2(ifile)
enddo
write(24,'(i5,f8.4,A)') iter,maxd,'  ! total iterations, deviation '
close(unit = 24)

! Calculates average properties such as 
!do kfile = 1,nfiles
!    do i=1,nentry(ifile)
!        y = 0.0
!        ylog = -1e9  ! exp(ylog) = y = 0.0
!        do jfile = 1,nfiles      
!            prob = -(beta(jfile)-beta(kfile))*e(ifile)%flength2(i) +&
!                (beta(jfile)*mu1(jfile)-beta(kfile)*mu1(kfile))*n1(ifile)%flength(i)+&
!                (beta(jfile)*mu2(jfile)-beta(kfile)*mu2(kfile))*n2(ifile)%flength(i)
!            tmp1 = prob + log(nentry(jfile)) - oldweight(jfile)
!            ylog = max(ylog,tmp1) + log(1.+exp(-abs(ylog-tmp1)))
!        enddo
!        weight(kfile) = max(weight(kfile),-ylog) + log(1.+exp(-abs(weight(kfile)+ylog)))
!    enddo
!    val = abs(weight(kfile)-oldweight(kfile))
!    if (val.gt.maxd) maxd=val
!enddo


!ttt = timef()
!itim = dtime(tarr)
!istat = hostnam(hostname)
!
!if (itim<300) then
!    write(*,'(A,f6.1,A,f5.1,2A)') ' Program run in', itim, ' CPU s, using ',itim/ttt*100,'% of one proc. on ',trim(hostname)
!    write(24,'(A,f6.1,A,f5.1,2A)') ' Program run in', itim, ' CPU s, using ',itim/ttt*100,'% of one proc. on ',trim(hostname)
!else if (itim<3600) then
!    write(*,'(A,f6.1,A,f5.1,2A)') ' Program run in', itim/60, ' CPU min, using ',itim/ttt*100,'% of one proc. on ',trim(hostname)
!    write(24,'(A,f6.1,A,f5.1,2A)') ' Program run in', itim/60, ' CPU min, using ',itim/ttt*100,'% of one proc. on ',trim(hostname)
!else
!    write(*,'(A,f6.2,A,f5.1,2A)') ' Program run in', itim/3600, ' CPU hr, using ',itim/ttt*100,'% of one proc. on ',trim(hostname)
!    write(24,'(A,f6.2,A,f5.1,2A)') ' Program run in', itim/3600, ' CPU hr, using ',itim/ttt*100,'% of one proc. on ',trim(hostname)
!endif


end

   

