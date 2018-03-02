program fspvt


! v. 2.11  6/24/08 - formats (from Antti-Pekka's version)
! 7/13/99 - weights are logarithms
! 6/29/99 - eliminated patching part; now requires patched data in 
! 6/26/99 - all arrays allocatable, different names/format of I/O files (AZP)

! Initial code written by Jeff Potoff
! version 1.0 - 10/31/97  This program iterates to find the weights for
! patching the individual histograms together.  The histogram data is read
! in the form of N1, N2, E for each configuration.

! 2/11/98 - Histogram patching rewritten to avoid overflows 
! 4/29/98 - New data type used to cut memory usage


!USE PORTLIB    
!use IFPORT ! for hostnam

implicit none

integer :: mentry                                       ! maximum number of observations
 
integer :: i                                            ! local counter
integer, allocatable :: min_part1(:),max_part1(:)       ! min. and max. number of particles in file
integer, allocatable :: min_part2(:),max_part2(:)       ! min. and max. number of particles in file
integer :: glob_min1, glob_max1, glob_min2, glob_max2   ! global minimum and maximum particle numbers
integer :: nfiles, ifile, k                             ! how many files being read
integer :: iostat                                       ! iostat becomes -1 on end-of-file
                                
real*8 :: xdim, ydim, zdim                              ! linear dimensions of box on which data were taken
real*8 :: xdim1, ydim1, zdim1                           ! linear dimensions of box for current file

! variables for new histogram format
integer :: jfile, kfile, iter
integer, allocatable, target :: n1_temp(:)
integer, allocatable, target :: n2_temp(:)
real*8, allocatable :: e_temp(:)
real*8 :: prob, maxd, val, y, ylog, tmp1
character*20, dimension(:), allocatable :: filename
real*8, dimension(:), allocatable :: beta
real*8, dimension(:), allocatable :: mu1
real*8, dimension(:), allocatable :: mu2
real*8, dimension(:), allocatable :: nentry
real*8, dimension(:), allocatable :: weight

!define data type - this is to required to reduce memory usage
type column1
  integer, pointer :: flength(:)                        ! Number of entries in the file
end type

type column2
  real*8, pointer :: flength2(:)                        ! Number of entries in the file
end type

type(column1), dimension(:), allocatable :: n1,n2
type(column2), dimension(:), allocatable :: e

character*20 ::auxfile
!character*20 :: hostname
!integer :: istat                                        ! aux variable
!real(4) :: itim, ttt, tarr(2)                           ! timing stuff


real :: betap, mu1p, mu2p                               ! points for evaluation of pvt
real*8 :: logZ                                          ! log of partition function
real*8, allocatable :: dens(:,:)                        ! [glob_min1:glob_max1,glob_min2:glob_max2] (n1,n2) log of density distr.
real*8 :: avgnum1, avgnum2, avgene                      ! mean numbers and energy
integer :: ipart, jpart                                 ! counters
real :: beta_tmp, mu1_tmp, mu2_tmp                      ! temp. variables

!-----------------executable part-----------------------
write (*,*) '--------------------------------------------------------------------------'
write (*,*) ' FS PVT properties, v. 2.11, July 08 (AZP) [run fspatch first]'
write (*,*) '--------------------------------------------------------------------------'

!ttt = timef()  ! total elapsed time
!itim = dtime(tarr)  ! CPU time

open(1,file='inp_fs2.dat',status='old')
read(1,*) nfiles
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

min_part1(:) = 0
max_part1(:) = 0
min_part2(:) = 0
max_part2(:) = 0
filename(:) = ""
beta(:) = 0.0
mu1(:) = 0.0
mu2(:) = 0.0
nentry(:) = 0.0
weight(:) = 0.0

do ifile = 1,nfiles
    read(1,*) filename(ifile),nentry(ifile),weight(ifile),beta(ifile),&
            mu1(ifile),mu2(ifile)
    beta(ifile) = 1./beta(ifile)
enddo
close(1)

mentry=maxval(nentry)
allocate(n1_temp(mentry),n2_temp(mentry),e_temp(mentry))

n1_temp(:) = 0
n2_temp(:) = 0
e_temp(:) = 0.0

write(*,'(A)') '         T        mu1       mu2    minp1   minp2  maxp1  maxp2  npoints'
                                        
do ifile= 1,nfiles
    min_part1(ifile) = 1000000    ! so that these are later set to correct values
    max_part1(ifile) = -1    
    min_part2(ifile) = 1000000
    max_part2(ifile) = -1
        
    open (23,file='his'//trim(filename(ifile))//'.dat')
    read (23,*)
    read (23,*) beta_tmp,mu1_tmp,mu2_tmp,xdim1,ydim1,zdim1
    beta_tmp = 1./beta_tmp
    if (abs((beta_tmp-beta(ifile)))+abs(mu1_tmp-mu1(ifile))+abs(mu2_tmp-mu2(ifile)) > 1.e-3) then
        write(*,*) 'Histogram file ',trim(filename(ifile)),' inconsistent with inp_fs2.dat'
        write (*,*) beta_tmp,beta(ifile),mu1_tmp,mu1(ifile),mu2_tmp,mu2(ifile)
        stop
    endif
    if (ifile==1) then
        xdim = xdim1
        ydim = ydim1
        zdim = zdim1
    else
        if (xdim1/=xdim.or.ydim1/=ydim.or.zdim1/=zdim) then
            write (*,*) 'System size in file ',trim(filename(ifile)),' inconsistent with first file'
            stop
        endif
    endif
    do i=1,nentry(ifile)
        read (23,*,iostat=iostat) n1_temp(i),n2_temp(i),e_temp(i)  
        min_part1(ifile) = min(n1_temp(i),min_part1(ifile))
        max_part1(ifile) = max(n1_temp(i),max_part1(ifile))
        min_part2(ifile) = min(n2_temp(i),min_part2(ifile)) 
        max_part2(ifile) = max(n2_temp(i),max_part2(ifile))
    enddo

    allocate(n1(ifile)%flength(1:i-1))
    allocate(n2(ifile)%flength(1:i-1))
    allocate(e(ifile)%flength2(1:i-1))

    n1(ifile)%flength(1:i-1) = n1_temp(1:nentry(ifile))
    n2(ifile)%flength(1:i-1) = n2_temp(1:nentry(ifile))
    e(ifile)%flength2(1:i-1) = e_temp(1:nentry(ifile))
    
    write(*,'(2x,A,f7.2,2x,f8.2,2x,f8.2,2x,i5,2x,i5,2x,i5,2x,i5,2x,f10.0)') &
      trim(filename(ifile)), &
        1./beta(ifile),mu1(ifile),mu2(ifile),min_part1(ifile),min_part2(ifile),&
        max_part1(ifile),max_part2(ifile),nentry(ifile)

enddo
close(23)
glob_max1 = maxval(max_part1)
glob_max2 = maxval(max_part2)
glob_min1 = minval(min_part1)
glob_min2 = minval(min_part2)
allocate (dens(glob_min1:glob_max1,glob_min2:glob_max2))
dens(:,:) = 0.0

open (28,file='pvt.dat')
write (28,'(A,20A5/3h /*,50A/3h /*,50A)') '/* Files = '&
  ,(trim(filename(ifile)),ifile=1,nfiles)
write (28,'(A)') '/*  T       mu1       mu2         <E>        <N1>        <N2>         lnZ'        
write (*,'(A)')  '    T       mu1       mu2         <E>        <N1>        <N2>         lnZ'        
    
open(unit=1,file='inp_pvt.dat',status='old')
iostat = 0
do while(iostat.ne.-1)
    read(1,*,iostat=iostat) betap,mu1p,mu2p
    if(iostat.ne.-1) then
        betap = 1./betap
        logZ = -1.e9
        dens = -1.e9
        avgnum1 = 0.
        avgnum2 = 0.
        avgene = 0.
        ! calculate the density distributions
        do ifile = 1,nfiles      
            do i=1,nentry(ifile)
                ylog = -1e9  ! exp(ylog) = y = 0.0
                do jfile = 1,nfiles      
                    prob = -(beta(jfile)-betap)*e(ifile)%flength2(i) +&
                        (beta(jfile)*mu1(jfile)-betap*mu1p)*n1(ifile)%flength(i)+&
                        (beta(jfile)*mu2(jfile)-betap*mu2p)*n2(ifile)%flength(i)
                    tmp1 = prob + log(nentry(jfile)) - weight(jfile)
                    ylog = max(ylog,tmp1) + log(1.+exp(-abs(ylog-tmp1)))
                enddo
                logZ = max(logZ,-ylog) + log(1.+exp(-abs(logZ+ylog)))
                dens(n1(ifile)%flength(i),n2(ifile)%flength(i)) = &
                    max(dens(n1(ifile)%flength(i),n2(ifile)%flength(i)),-ylog) + &
                    log(1.+exp(-abs(dens(n1(ifile)%flength(i),n2(ifile)%flength(i))+ylog)))
                avgene = avgene + e(ifile)%flength2(i)*exp(-ylog)
            enddo
        enddo
        avgene  = avgene / exp(logZ)
        do ipart = glob_min1,glob_max1
            do jpart = glob_min2,glob_max2
                avgnum1 = avgnum1 + ipart * exp(dens(ipart,jpart) - logZ)
                avgnum2 = avgnum2 + jpart * exp(dens(ipart,jpart) - logZ)
            enddo
        enddo
        write (28,'(f7.3,2f10.3,3x,4g12.4)') 1/betap,mu1p,mu2p,avgene,avgnum1,avgnum2,logZ
        write (*,'(f7.3,2f10.3,3x,4g12.4)') 1/betap,mu1p,mu2p,avgene,avgnum1,avgnum2,logZ
    endif
enddo

!ttt = timef()
!itim = dtime(tarr)
!istat = hostnam(hostname)
!
!if (itim<300) then
!    write(*,'(A,f6.1,A,f5.1,2A)') '  Program run in', itim, ' CPU s, using ',itim/ttt*100,'% of one proc. on ',trim(hostname)
!else if (itim<3600) then
!    write(*,'(A,f6.1,A,f5.1,2A)') '  Program run in', itim/60, ' CPU min, using ',itim/ttt*100,'% of one proc. on ',trim(hostname)
!else
!    write(*,'(A,f6.2,A,f5.1,2A)') '  Program run in', itim/3600, ' CPU hr, using ',itim/ttt*100,'% of one proc. on ',trim(hostname)
!endif

end

