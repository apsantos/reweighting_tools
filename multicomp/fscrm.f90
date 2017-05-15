
MODULE data
    integer nfiles          ! how many files being read
    integer, parameter :: maxcrit = 3000    ! discretization of the field mixing parameter
    integer, parameter :: nising = 168  ! number of points in ising.dat file
    real xising(nising),yising(nising)  ! from ising.dat
    real*8 betap,mu1p,mu2p      ! trial beta,mu1,mu2
    real*8, dimension(:), allocatable :: beta
    real*8, dimension(:), allocatable :: mu1
    real*8, dimension(:), allocatable :: mu2
    real*8, dimension(:), allocatable :: nentry
    real*8, dimension(:), allocatable :: weight
  
    !define data type - this is to required to reduce memory usage
    type column1
        integer, pointer ::flength(:) ! Number of entries in the file
    end type

    type column2
        real*8, pointer :: flength2(:) ! Number of entries in the file
    end type

    type(column1), dimension(:), allocatable :: n1,n2
    type(column2), dimension(:), allocatable :: e

! *** Moved from main program 16-2-2004
    integer, parameter :: ebin = 1001   ! number of energy bins
    real*8 eng(ebin),eng_ave,etemp  ! energy distribution
    real*8 eng_min,eng_max      !min and max of energies stored in the lists
    real*8, dimension(:), allocatable :: dens1 ! n1 distribution
    real*8, dimension(:), allocatable :: dens2 ! n2 distribution
    integer h
    real*8 deltau
    real*8 Z            ! partition function
    real*8 Z2
    real*8 smix1,smix2      ! field mixing parameters for e and n2
    real*8 x            ! Order parameter distribution x axis
! critical point determination
    integer index
    real*8 ncrit            ! number of points for p_order
    real*8 min_x,max_x,avg_x    ! Min, max and average of x
    real*8 p_order(0:maxcrit)   ! Order parameter distribution
    real*8 variance
! ***
END MODULE data

program crtm

! v. 2.12  6/24/08 - formats (from Antti-Pekka's version)
! v. 2.11 2/12/03 : unformatted data read added (activated with unform=.true.) (APH)
! v. 2.1  8/01 - weight = log(weight) to avoid overflows (AZP)
! Version 2.0 - 4/30/98
! This program patches the histograms together and matches the
! critical density distribution to the critical ordering operator
! distribution.  "fspatch" should be used to patch find the weights
! before using this program.

use data
implicit none

integer, parameter :: mentry = 6e6  ! maximum number of observations 

integer ifile,jfile,kfile,i,ipart   !local counters
integer iter,k,j        ! # of iter used in finding the patching constants
integer min_part1f,max_part1f   ! min. and max. number of particles in file
integer min_part2f,max_part2f   ! min. and max. number of particles in file
integer min_part1,max_part1 ! min. and max. number of particles - global
integer min_part2,max_part2 ! min. and max. number of particles - global

integer iostat          ! iostat becomes -1 on end-of-file
real*8 nliq1,ngas1,nliq2,ngas2  ! liquid and gas densities on coexistence, accumulator of densities for calculation
                ! of coexistence
!real*8 tmin,tmax,tincr     ! for calculation of phase diagram
real*8 xdim,ydim,zdim       ! linear dimensions of box on which data were taken
real*8 xdim1,ydim1,zdim1    ! linear dimensions of box for current file

! variables for list data format
integer, target :: n1_temp(mentry)  ! temporary storage for the list
integer, target :: n2_temp(mentry)  ! temporary storage for the list
real*8, target :: e_temp(mentry)    ! temporary storage for the list
real*8 prob,yz,ylog,tmp1    ! temporary variables
real*8 avgnum1,avgnum2      ! average number of particles 1 and 2

character*20, dimension(:), allocatable :: filename

! amoeba search
integer, parameter :: n_dim = 3    ! number of optimizing directions: mu1, smix1, smix2
integer ipoint              ! for counting ising data
real ptry(n_dim)            ! vector of trial parameters
real p(n_dim+1,n_dim)       ! input to optimization function
real y(n_dim+1)         ! input to optimization function
real, parameter :: ftol = 4e-3      ! tolerance for convergence of optimization
real trial_vec(n_dim)       ! factors for initial search in Tcrit
!data trial_vec / 0.0005, -0.05, -0.05/ ! initial values 
real Tstart,Tend,Tincr,dmu_dT   ! Temperatures to be studied, change of mu vs. T at rho = rho_c
integer Ntemp, itemp
real    devmin, betamin, mu1min, mu2min, smix1min, smix2min, ncritmin   ! minimum
integer itermin

real funk
external funk

logical unform

!---------------------------------------------------------------------------------------------
write (*,*) '--------------------------------------------------------------------------'
write (*,*) ' FS Critical Point Search, v. 2.12, July 08 (AZP) [run fspatch first]'
write (*,*) '--------------------------------------------------------------------------'

unform = .false.

if (unform) then
    write (*,*) '  Reading in unformatted data (unform=.true.)'
endif

!open (1,file='e:\pctn\ising.dat',status='old')
open (1,file='ising.dat',status='old')
do ipoint =1,nising
    read (1,*) xising(ipoint),yising(ipoint)
enddo
close(1)

open(unit=2,file='inp_fs2.dat',status='old')
read(2,*) nfiles, trial_vec(1:n_dim)

allocate(n1(nfiles))
allocate(n2(nfiles))
allocate(e(nfiles))
allocate(filename(nfiles))
allocate(beta(nfiles))
allocate(mu1(nfiles))
allocate(mu2(nfiles))
allocate(nentry(nfiles))
allocate(weight(nfiles))

do ifile = 1,nfiles
    read(2,*) filename(ifile),nentry(ifile),weight(ifile),beta(ifile),&
                mu1(ifile),mu2(ifile)
    beta(ifile) = 1./beta(ifile)
enddo
close(unit=2)

min_part1 = 1e5 ! so that these are later set to correct values
max_part1 = -1  
min_part2 = 1e5 
max_part2 = -1  
write(*,'(10x,A)') ' T         mu1        mu2     minp1 maxp1 minp2 maxp2 npoints'
                                        
do ifile= 1,nfiles
  
    min_part1f = 1e5 ! so that these are later set to correct values
    max_part1f = -1 
    min_part2f = 1e5    
    max_part2f = -1
    
    if (unform) then    
        open (23,file='his'//trim(filename(ifile))//'.dat',form='unformatted')
        read (23) beta(ifile),mu1(ifile),mu2(ifile),xdim1,ydim1,zdim1
    else
        open (23,file='his'//trim(filename(ifile))//'.dat')
        read (23,*)
        read (23,*) beta(ifile),mu1(ifile),mu2(ifile),xdim1,ydim1,zdim1
    endif
    beta(ifile) = 1./beta(ifile)
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
               
        if(i.gt.mentry) then
            write(*,'(A,f12.0,A,i12)') 'maximum number of entries exceeded',&
                nentry(ifile),'>',mentry
            stop
        endif

        if (unform) then
            read (23,iostat=iostat) n1_temp(i),n2_temp(i),e_temp(i)
        else
            read (23,*,iostat=iostat) n1_temp(i),n2_temp(i),e_temp(i)
        endif

        ! do not read zero entries (APH 25th March 2004)
!        if (n1_temp(i)==0.and.n2_temp(i)==0) then
!            i = i - 1
!        endif

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
      
    allocate(n1(ifile)%flength(1:i-1))
    allocate(n2(ifile)%flength(1:i-1))
    allocate(e(ifile)%flength2(1:i-1))
    
    n1(ifile)%flength(1:i-1) = n1_temp(1:nentry(ifile))
    n2(ifile)%flength(1:i-1) = n2_temp(1:nentry(ifile))
    e(ifile)%flength2(1:i-1) = e_temp(1:nentry(ifile))
  
    write(*,'(2x,A,2x,f8.5,2x,f9.5,2x,f9.5,2x,i5,1x,i5,1x,i5,1x,i5,2x,f8.0)') &
        trim(filename(ifile)),1./beta(ifile),mu1(ifile),mu2(ifile),min_part1f,max_part1f,&
        min_part2f,max_part2f,nentry(ifile)

enddo
close(unit=23)

deltau = (eng_min-eng_max)/dble(ebin-1)
allocate(dens1(min_part1:max_part1))
allocate(dens2(min_part2:max_part2))

open(unit=25, file='crit2.dat', status = 'unknown')
write (25,'(A,20A5/3h /*,50A/3h /*,50A)') '/* Files = '&
    ,(trim(filename(ifile)),ifile=1,nfiles)
    
write (25,'(A,4x,A)') '/* T        mu1        mu2         N_1          N_2    ln(Z)',&
       'Smix1     Smix2     <E>        DEV'     
  
k = 0
! start ======
do while (k.eq.0) 
            
    write(*,*) ' Enter T, mu1, mu2,smix1,smix2,ncrit (T<0 for optimization)'
    read(*,*) betap,mu1p,mu2p,smix1,smix2,ncrit
    betap = 1./betap
    if (betap < 0.) goto 100
    ! set initial conditions for the amoeba
    p(1,1) = mu1p
    p(1,2) = smix1
    p(1,3) = smix2
    ptry(1:n_dim) = p(1,1:n_dim)
    y(1) = funk(ptry)
     
    open (24,file='crit.dat',status = 'unknown')
    write (24,'(A,f8.4,A,f10.4,A,f10.4,A,80A5/3h /*,50A/3h /*,50A)') &
            '/* T=',1/betap,';mu1=',mu1p,';mu2=',mu2p,'; Files = ', &
            (trim(filename(ifile)),ifile=1,nfiles) 
    do i = 0,int(ncrit)
        write(24,'(2g12.3)') &
                (min_x+i*(max_x-min_x)/ncrit-avg_x)/sqrt(variance), &
                  p_order(i)*(ncrit/(max_x-min_x))*sqrt(variance)
    enddo
    write(24,'(A,2g14.3)') '/* <x> =',avg_x,variance
    close(24)

    ! Calculate density averages and write them to dens1.dat and dens2.dat
    avgnum1 = 0.
    ngas1 = 0.
    open (24,file='dens1.dat')
    write (24,'(A,f8.4,A,f8.4,A,f8.4,A,f9.5,A,f9.5,A)') &
        '/* T=',1/betap,'; mu1=',mu1p,'; mu2=',mu2p,'; Smix1=',smix1,'; Smix2=',smix2,' */'
    do ipart = min_part1,max_part1
        avgnum1 = avgnum1 + ipart * dens1(ipart)/Z
        ngas1 = ngas1 + ipart * dens1(ipart)/Z
        if (dens1(ipart)>1e-99) write (24,'(i5,g12.3)') ipart,dens1(ipart)
    enddo
    close(24)

    avgnum2 = 0.
    ngas2 = 0.
    open (24,file='dens2.dat')
    write (24,'(A,f8.4,A,f8.4,A,f8.4,A,f9.5,A,f9.5,A)') &
        '/* T=',1/betap,'; mu1=',mu1p,'; mu2=',mu2p,'; Smix1=',smix1,'; Smix2=',smix2,' */'
    do ipart = min_part2,max_part2
        avgnum2 = avgnum2 + ipart * dens2(ipart)/Z 
        ngas2 = ngas2 + ipart * dens2(ipart)/Z
        if (dens2(ipart)>1e-99) write (24,'(i5,g12.3)') ipart,dens2(ipart)
    enddo
    close(24)
  
    ! calculate the average energy 
    eng_ave = 0.0
    open(unit=27,file='energy.dat',status='unknown')
    do i=1,ebin-1 
        etemp = (i-1)*deltau+eng_max
        write(27,*) etemp, eng(i)/Z
        eng_ave = eng_ave + etemp*eng(i)/Z      
    enddo 
    close(27)
    write(*,'(A,f10.6,A,f8.3,A,f8.3,A,F14.8)') 'Deviation = ',real(y(1)*100),' <N1> = ',avgnum1,' <N2> = ',avgnum2,' ln(Z) = ',log(Z2)
    write (25,'(f8.4,2f12.4,3f10.3,2f10.6,f12.3,f8.4)') &
              1./betap,mu1p,mu2p,ngas1,ngas2,log(Z),smix1,smix2,eng_ave,y(1)*100
       
enddo 

100 continue

! *** Perform Amoeba optimization 16-2-2004 (APH)
write (*,*) ' Enter Tstart, Tend, Tincr, Mu1, Mu2, Smix1, Smix2, ncrit, dmu_dT'
read (*,*) Tstart,Tend,Tincr,mu1p,mu2p,smix1,smix2,ncrit,dmu_dT
write (*,*) ' Perturb.:      Mu1        Smix1        Smix2      Dev'

betap = 1./Tstart

p(1,1) = mu1p
p(1,2) = smix1
p(1,3) = smix2
ptry(1:n_dim) = p(1,1:n_dim)
y(1) = funk(ptry)
write(*,'(10x,3f12.5,f10.3)') ptry(1),ptry(2),ptry(3),y(1)*100
p(2,1) = mu1p + trial_vec(1)
p(2,2) = smix1
p(2,3) = smix2
ptry(1:n_dim) = p(2,1:n_dim)
y(2) = funk(ptry)
write(*,'(10x,3f12.5,f10.3)') ptry(1),ptry(2),ptry(3),y(2)*100
p(3,1) = mu1p
p(3,2) = smix1 + trial_vec(2)
p(3,3) = smix2
ptry(1:n_dim) = p(3,1:n_dim)
y(3) = funk(ptry)
write(*,'(10x,3f12.5,f10.3)') ptry(1),ptry(2),ptry(3),y(3)*100
p(4,1) = mu1p
p(4,2) = smix1
p(4,3) = smix2 + trial_vec(3)
ptry(1:n_dim) = p(4,1:n_dim)
y(4) = funk(ptry)
write(*,'(10x,3f12.5,f10.3)') ptry(1),ptry(2),ptry(3),y(4)*100
write (*,*) '--------------------------------------------------------------------------'
write (*,*) '    T         Mu1       Mu2      Smix1     Smix2  Ncrit  iter    Dev.  '
write (*,*) '--------------------------------------------------------------------------'

Ntemp = int((Tend - Tstart)/Tincr) + 1


devmin = 100000.0

do itemp = 1,Ntemp
	open (26,file='crit.dat')
    if (itemp /= 1) then    ! Set initial amoeba for new temperature
        betap = 1./(Tstart + (itemp-1)*Tincr)
        p(1:4,1) = p(1:4,1) + dmu_dT*Tincr  ! set chemical potential to a more reasonable value
        ptry(1:n_dim) = p(1,1:n_dim)
        y(1) = funk(ptry)
        p(2,1) = p(2,1) + trial_vec(1)/3. ! divide by 3 since already near optimum
        ptry(1:n_dim) = p(2,1:n_dim)
        y(2) = funk(ptry)
        p(3,2) = p(3,2) + trial_vec(2)/3.
        ptry(1:n_dim) = p(3,1:n_dim)
        y(3) = funk(ptry)
        p(4,3) = p(4,3) + trial_vec(3)/3.
        ptry(1:n_dim) = p(4,1:n_dim)
        y(4) = funk(ptry)
    endif

!   write (*,'(4f12.5)') 1./betap,mu1p,smix1,smix2

    call amoeba(p,y,n_dim+1,n_dim,n_dim,ftol,funk,iter)

    mu1p = p(1,1)
    smix1 = p(1,2)
    smix2 = p(1,3)
    ptry(1:n_dim) = p(1,1:n_dim)
    y(1) = funk(ptry)

!   write (*,*) ptry(1:n_dim)
    write (*,'(5f10.5,f6.1,i5,f10.3)') 1./betap,mu1p,mu2p,smix1,smix2,ncrit,iter,y(1)*100

    if (y(1) < devmin) then
        devmin = y(1)
        betamin = betap
        mu1min = mu1p
        mu2min = mu2p
        smix1min = smix1
        smix2min = smix2
        ncritmin = ncrit
        itermin = iter
    endif

    ! Write crit.dat
    write (26,'(A,f8.4,A,f10.4,A,f10.4,A,20A5/3h /*,50A/3h /*,50A)') &
            '/* T=',1/betap,';mu1=',mu1p,';mu2=',mu2p,'; Files = ', &
            (trim(filename(ifile)),ifile=1,nfiles)

    do i = 0,int(ncrit)
        write(26,'(2g12.3)') &
                (min_x+i*(max_x-min_x)/ncrit-avg_x)/sqrt(variance), &
                  p_order(i)*(ncrit/(max_x-min_x))*sqrt(variance)
    enddo
    write(26,'(A,2g14.3)') '/* <x> =',avg_x,variance
    close(26)

enddo
! Set system back to minimum
betap = betamin
mu1p = mu1min
mu2p = mu2min
smix1 = smix1min
smix2 = smix2min
ncrit = ncritmin
iter = itermin
write (*,*) ' minimum:'
write (*,'(5f10.5,f6.1,i5,f10.3)') 1./betap,mu1p,mu2p,smix1,smix2,ncrit,iter,devmin*100

close(25)
end

! ===========================================================================  

real function funk(ptry)
use data

integer ipoint              ! for counting on ising data
integer i,ifile,jfile       ! local counters
integer inonzero        ! number of actual points in distribution
real xpoint, ypoint, yvalue
real ptry(5)            ! vector of trial parameters
real*8 prob,yz,ylog,tmp1    ! temporary variables

! *** Commented out 16-2-2004
!integer index
!real*8 ncrit           ! number of points in the ordering op. distribution
!real*8 p_order(0:maxcrit)  ! ordering operator distribution
!real*8 max_x, min_x, avg_x ! limits of x, average x
!real*8 variance
!real*8 Z           ! partition function
!real*8 smix1,smix2     ! field mixing parameters
!real*8 x               ! ordering operator distribution x axis
! ***

! Function calculates the difference between the ising and the real distribution

! mu2p and ncrit are assumed to be set to correct values already
mu1p = ptry(1)
smix1 = ptry(2)
smix2 = ptry(3)

Z = 0.
Z2 = 0.
eng = 0.0
avgnum1 = 0.
avgnum2 = 0.
dens1 = 0.0
dens2 = 0.0
ngas1 = 0.
nliq1 = 0.
ngas2 = 0.
nliq2 = 0.
p_order = 0.0
avg_x = 0.0
! Find min and max of critical point distribution
min_x = 1e6
max_x = -1e6
do ifile = 1,nfiles
    do j=1,nentry(ifile)
        x = n1(ifile)%flength(j) - smix1*e(ifile)%flength2(j) - smix2*n2(ifile)%flength(j)
        min_x = min(min_x,x)
        max_x = max(max_x,x)
    enddo
enddo

!write (*,*) 'min_x = ',min_x,' max_x = ',max_x

! calculate Z    
do ifile = 1,nfiles
    do i=1,nentry(ifile)
!    yz = 0.0
        ylog = -1e9  ! exp(ylog) = y = 0.0
        do jfile=1,nfiles
            prob =(-(beta(jfile)-betap)*e(ifile)%flength2(i) + &
                    (beta(jfile)*mu1(jfile)-betap*mu1p)*n1(ifile)%flength(i) + &
                    (beta(jfile)*mu2(jfile)-betap*mu2p)*n2(ifile)%flength(i))
            !     yz = yz + nentry(jfile)*exp(prob-log(weight(jfile)))
            ! The following is a code rearrangement to avoid overflow problems
            ! This trick was compliments of Peter-Lawrence Montgomery 2/12/98
            ! Replace y with log(y)
            ! Define log_expsum(y,z) = log(exp(y)+exp(z)
            ! Evaluate log_expsum(y,z) = max(y,z) + log(1+exp(-abs(y-z)))
            tmp1 = prob + log(nentry(jfile)) - weight(jfile)
            ylog = max(ylog,tmp1) + log(1.+exp(-abs(ylog-tmp1)))
        enddo
        Z2 = Z2 + exp(-ylog)
        if (n1(ifile)%flength(i) /= 0.or.n2(ifile)%flength(i) /= 0) then
        Z = Z + exp(-ylog)
        x = n1(ifile)%flength(i) - smix1*e(ifile)%flength2(i) - smix2*n2(ifile)%flength(i) ! calc pos on x axis
        index = int(ncrit*(x-min_x)/(max_x-min_x))
        p_order(index) = p_order(index) + exp(-ylog)
        avg_x = avg_x + x*exp(-ylog)
! *** Added 16-2-2004
        h = int((e(ifile)%flength2(i)-eng_max)/deltau) + 1
        eng(h) = eng(h) + exp(-ylog)
        dens1(n1(ifile)%flength(i)) = dens1(n1(ifile)%flength(i)) + exp(-ylog)
        dens2(n2(ifile)%flength(i)) = dens2(n2(ifile)%flength(i)) + exp(-ylog) 
! ***
        endif
    enddo  
enddo
 
! calculate the variance
avg_x = avg_x/Z
p_order = p_order/Z
variance = 0.0
do i = 0,maxcrit
    variance = variance + p_order(i)*(min_x+i*(max_x-min_x)/ncrit-avg_x)**2.    
enddo

!write (*,*) 'avg_x = ',avg_x

funk = 0.
ipoint = 1
inonzero = 0
do i = 0,int(ncrit)+1
    xpoint = (min_x+i*(max_x-min_x)/ncrit-avg_x)/sqrt(variance)
    ypoint = p_order(i)*sqrt(variance)/(max_x-min_x)*ncrit
    do while (ipoint <= nising-2 .and. xising(ipoint+1) < xpoint)
        ipoint = ipoint + 1
    enddo
    yvalue = yising(ipoint) + (yising(ipoint+1) - yising(ipoint))/(xising(ipoint+1)-xising(ipoint))*(xpoint - xising(ipoint))
    funk = funk + abs(ypoint - yvalue)
    if (ypoint > 1e-6 .and. yvalue > 1e-6) inonzero = inonzero + 1
enddo

funk = funk/inonzero

return 
end

! Routines from numerical recipies below.
!================================================================================
      FUNCTION amotry(p,y,psum,mp,np,ndim,funk,ihi,fac)
      INTEGER ihi,mp,ndim,np,NMAX
      REAL amotry,fac,p(mp,np),psum(np),y(mp),funk
      PARAMETER (NMAX=5)
      EXTERNAL funk
!    USES funk
      INTEGER j
      REAL fac1,fac2,ytry,ptry(NMAX)
      fac1=(1.-fac)/ndim
      fac2=fac1-fac
      do 11 j=1,ndim
        ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
11    continue
      ytry=funk(ptry)
      if (ytry.lt.y(ihi)) then
        y(ihi)=ytry
        do 12 j=1,ndim
          psum(j)=psum(j)-p(ihi,j)+ptry(j)
          p(ihi,j)=ptry(j)

12      continue
      endif
      amotry=ytry
!      write(*,'(3f12.5,2f15.6,f12.5)') &
!               ytry*100,ptry(1),ptry(2),ptry(3)
      return
      END

SUBROUTINE amoeba(p,y,mp,np,ndim,ftol,funk,iter)
INTEGER iter,mp,ndim,np,NMAX,ITMAX
REAL ftol,p(mp,np),y(mp),funk
PARAMETER (NMAX=5,ITMAX=5000)
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
      if (iter.ge.ITMAX) pause 'ITMAX exceeded in amoeba'
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
