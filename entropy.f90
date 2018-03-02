program calc_entropy

! version 2.02 - Oct. '05 (AZP) - specexp, moved to lnZ rather than Z
! version 2.01 - July '03 (AZP) - overlap problem fixed
! version 2.00 - July '02 (AZP) - complete rewrite of I/O (no shifting of arrays etc)
! version 1.31 - July '02 (AZP) Nmid=0 for scan; energy histogram
! version 1.30 - May '00 (AZP)
! calculation of cmc's
! version 1.25 - June '99 (AZP)
! calc. of Cv
! version 1.24 - Apr '99 (AZP)
! formats for mu
! version 1.23 - August '98 (AZP)
! fixed problem with ene_liq
! version 1.22 - August '98 (AZP)
! fixed problem with overestimation of npoints; moved shift_Z=0 statement
! version 1.21 - July '98 (AZP)
! differences from version 1.1:  (a) better protection against over/under flow 
! (b) ability to do a "scan" over chemical potentials
! (c) dynamic allocation of arrays
! (d) individual suffix coverage and density files are not produced 
!  (commented out) - units 26 and 27
! version 1.1 - July '96 (AZP)
! program takes as input observations of histograms of N and E
! and calculates the entropy S(N,E) = ln(Omega) by "patching"


implicit none

integer             minpart,maxpart ! global min / max number of particles 
real                min_ene,max_ene ! global min / max energy observed 
integer             minpart1,maxpart1   ! min / max number of particles in a file
real                min_ene1,max_ene1   ! min / max energy observed in a file
real                dum_ene,pt_ene  ! start energy on read, point energy
integer             npoints         ! points along energy dimension
real, allocatable:: entropy(:,:)    ! [npoints,minpart:maxpart]  contains counts for a certain energy+number
integer, allocatable:: tot_counts(:,:)  ! [npoints,minpart:maxpart] contains counts (observations) for a certain energy+number
!integer, allocatable:: counts(:,:)     ! counts for every file being read, per (E,N) point
real*8, allocatable:: counts(:,:)       ! counts for every file being read, per (E,N) point
real, allocatable:: omega(:,:)      ! temporary holding variable, reads contr. to entropy from
real*8, allocatable:: counts_file(:)    ! total counts for each file being read
real, allocatable:: sum_counts(:)   ! 
real*8              beta,mu         ! temperature, chemical potential
real                width           ! increment of energy distributions in each file
real                width1          ! width1 makes sure that all widths are the same
character(30),allocatable:: filename(:)     ! histogram file name - if the form his?????.dat; only the ????? part is stored
character(30) file_aux  ! auxiliary filename
integer nfiles,ifile        ! how many files being read
integer iostat              ! iostat becomes -1 on end-of-file

!real energy_start
!integer, parameter :: maxglobal = 5000 ! "global" max number of particles
!real energy1,energy2(0:maxglobal),energy3(0:maxglobal) ! aux. variable to determine npoints

!integer ncount ! counts how many observations there are for a certain number of particles
!integer min_part,max_part,ipart,i,idiff,iend   ! local counters
!integer min_part1,max_part1        ! min. and max. number of particles in file
!integer j  ! local counter
!real, parameter :: incr_width = 150    ! width for energy/N in counts file

real*8              scale           ! scale determines how different files are "patched" together
real*8              lnZ,exponent        ! ln of partition function, exponent of calculation
integer*8               noverlap,noverlap1  ! counters of overlaps between histograms 
real*8, allocatable:: dens(:)       ! stores the probability of a certain number
real*8, allocatable:: edens(:)      ! stores the probability of a certain energy
real*8              avgnum          ! is the mean number of particles
real*8              ene_liq,ene_gas ! energy of liquid and gas
real*8              nmid            ! mid-point number of particles on coexistence
real*8              nliq,ngas,nbelow    ! liquid and gas densities on coexistence, accumulator of densities for calculation
                                        ! of coexistence
integer             iter            ! iteration number
integer, parameter :: maxiter=100   ! maximum iterations for phase equil.
real*8              tmin,tmax,t_incr    ! for calculation of phase diagram
real                xdim,ydim,zdim  ! linear dimensions of box on which data were taken
real                xdim1,ydim1,zdim1   ! linear dimensions of box for current file
real*8              nbelow1         ! for Newton's method
real*8, parameter :: muincr = 0.001 ! for Newton's derivative eval.
real, parameter ::  maxexponent = 300   ! for avoiding overflows
real*8, parameter :: mindens = .45  ! minimum accepted peak separation
character(2),allocatable:: suffix(:)    ! for dealing with series of runs
integer             nsuffix,isuffix ! number of different run series
logical             do_phase    ! indicates if a phase equil. calculation is to be done
logical             do_scan ! indicates if a scan over pvt values is to be done
real*8              end_mu,incr_mu ! for scanning mu
!real*8             shift_Z ! amount by which to shift exponentials
real*8              idum    ! to read in dummy values in initial pass
real*8              ene2,num2   ! square of energy, number of particles
real*8              ene_num ! energy*number
real*8              cv_part ! heat capacity per particle
integer             nfiles_suf
real*8              min_slope, slope, slope0, slope1, slope2, slope2_min    ! dZ/dN 
real*8              mu_h0,mu_h1,mu_h2,num_h0,num_h1,num_h2,Z_h0,Z_h1,Z_h2   ! for derivatives
real*8              mu_cmc, num_cmc ! cmc parameters
integer             count_pass, index, ipart, ncount, i ! do loop indices
real*8 specexp,tmp1,auxreal
character(len=32) :: nfile_fmt
specexp(auxreal,tmp1) = max(auxreal,tmp1) + log(1.+dexp(-abs(auxreal-tmp1)))


open (1,file='input_hs.dat')
read (1,*) nsuffix
allocate (suffix(nsuffix))

do isuffix = 1,nsuffix
    read (1,*) suffix(isuffix)
enddo

nfiles = 0
do while (iostat/=-1)
    read (1,*,iostat=iostat) file_aux
    nfiles=nfiles+1
enddo
nfiles = nfiles-1

nfiles_suf = nfiles+nsuffix
allocate (counts_file(nfiles*nsuffix),filename(nfiles*nsuffix))
close (1)

open (1,file='input_hs.dat')
read (1,*)
do isuffix = 1,nsuffix
    read (1,*) 
enddo
do ifile = 1,nfiles
    read (1,*) filename(ifile)
enddo
close (1)

open (1,file='input_hs2.dat') ! prepare file for input to FS program
write (1,*) nsuffix
do isuffix = 1,nsuffix
    write (1,*) trim(suffix(isuffix))
enddo

write (*,'(a)') '  ------------------------------------------------------------------------'
write (*,'(a)') '  Approximate histogram patch program - version 2.02 (AZP, Oct. 05)'
write (*,'(a)') '  ------------------------------------------------------------------------'

maxpart = 0
minpart = huge(minpart)
max_ene = -huge(max_ene)
min_ene = huge(min_ene)
npoints = 0
do ifile = 1,nfiles  ! initial read to check and get maxpart and npoints
    do isuffix = 1,nsuffix  
        iostat = 0
        open (23,file='his'//trim(filename(ifile))//trim(suffix(isuffix))//'.dat',status="old")
        read (23,*)
        read (23,*) beta,mu,width1,xdim1,ydim1,zdim1
        if (ifile==1.and.isuffix==1) then
            xdim = xdim1
            ydim = ydim1
            zdim = zdim1
        else if (xdim1/=xdim.or.ydim1/=ydim.or.zdim1/=zdim) then
            write (*,*) 'System size in file ',trim(filename(ifile))//trim(suffix(isuffix)),' inconsistent with first file'
            stop
        endif
        do while (iostat/=-1)
            read (23,*,iostat=iostat) ipart,ncount,dum_ene
            if (iostat/=-1) then
                min_ene = min(min_ene,dum_ene)
                max_ene = max(max_ene,dum_ene+(ncount-1)*width1)
                minpart = min(minpart,ipart) ; maxpart = max(maxpart,ipart)
                read (23,*) (idum,i=1,ncount)
            endif
        enddo
        close(23)
    enddo
enddo
npoints = (max_ene-min_ene)/width1 + 2

! allocate arrays
allocate (entropy(npoints,minpart:maxpart),tot_counts(npoints,minpart:maxpart))
allocate (omega(npoints,minpart:maxpart),counts(npoints,minpart:maxpart),sum_counts(npoints))
allocate (dens(minpart:maxpart),edens(npoints))

write (*,'(a,3i7,A,f7.2)') '  Min/Maxpart, Npoints = ',minpart,maxpart,npoints,' Mem usage (MB)=', &
                ( kind(entropy)*size(entropy) + kind(tot_counts)*size(tot_counts) + &
                  kind(omega)*size(omega) + kind(counts)*size(counts) + &
                  kind(sum_counts)*size(sum_counts) + 2*kind(dens)*size(dens) )/(2.**20)
write (1,*) nfiles, npoints, maxpart

counts = 0
tot_counts = 0
entropy = 0.    ! program looks at counts, so a zero value is not taken seriously
                ! unless backed up by data

do ifile = 1,nfiles
    do isuffix = 1,nsuffix  ! no ident for this do loop
    open (23,file='his'//trim(filename(ifile))//trim(suffix(isuffix))//'.dat',status="old")
    read (23,*)
    read (23,*) beta,mu,width1,xdim1,ydim1,zdim1
    beta = 1./beta
    iostat = 0
    omega = 0.
    counts = 0
    counts_file(ifile) = 0.
    dens = 0.
    maxpart1 = 0 ;  minpart1 = huge(minpart1) ; max_ene1 = -huge(max_ene1) ; min_ene1 = huge(min_ene1)
    do while (iostat/=-1)
        read (23,*,iostat=iostat) ipart,ncount,dum_ene
        if (iostat/=-1) then
            maxpart1 = max(maxpart1,ipart) ; minpart1 = min(minpart1,ipart)
            index = (dum_ene - min_ene) / width1
            read (23,*) (counts(i,ipart),i=index+1,index+ncount)
            do i = index+1,index+ncount
                if (counts(i,ipart)>0) then
                    pt_ene = dum_ene+(i-index-1)*width1
                    omega(i,ipart) = log(real(counts(i,ipart))) + beta*pt_ene - beta*mu*ipart
                    dens(ipart) = dens(ipart) + counts(i,ipart)
                    counts_file(ifile) = counts_file(ifile)+counts(i,ipart)
                    max_ene1 = max(max_ene1,pt_ene) ; min_ene1 = min(min_ene1,pt_ene)
                endif
            enddo
        endif
    enddo
    close(23)

    if (ifile==1.and.isuffix==1) then
        entropy = omega
        tot_counts = counts
        width = width1
        write (*,'(A)') '  File     T         mu    min/max: N             E        Overlap (10^6)'
        write (*,'(2A,X,f9.4,X,f12.4,X,2i8,X,2f8.0)') '  ',trim(filename(ifile))//trim(suffix(isuffix)),1/beta,mu,minpart1,maxpart1,&
                                                   min_ene1,max_ene1
        write (1,'(A,g12.4,X,2f15.6,X,f13.6,X,f12.6,X,f7.1)') trim(filename(ifile)), counts_file(ifile), 0.,1./beta,mu,width,1.
    else ! patch current data with previous data 
        if (width1 /= width) then
            write (*,*) 'Width not consistent with previous, file = ',trim(filename(ifile))//trim(suffix(isuffix))
            stop
        endif
        scale = 0.
        noverlap = 0
        do ipart = minpart,maxpart
            do i = 1,npoints
                if (counts(i,ipart)>0.and.tot_counts(i,ipart)>0) then
                    ! scale is weighted by the mean of the inverse of the number of observations in the
                    ! existing and new files
                    noverlap1 = 0
                    noverlap1 = 2./(1./counts(i,ipart)+1./tot_counts(i,ipart))
                    scale = scale + (omega(i,ipart)-entropy(i,ipart)) * noverlap1
!                       write(30,'(2i4,i8,f10.2)') ipart,i,noverlap1,(omega(i,ipart)-entropy(i-idiff,ipart))
                    noverlap = noverlap + noverlap1
                endif
            enddo
        enddo   
        write (*,'(2A,X,f9.4,X,f16.4,X,2i8,X,2f8.0,X,f12.3)') '  ',trim(filename(ifile))//trim(suffix(isuffix)),1/beta,mu,minpart1,maxpart1&
                                                        ,min_ene1,max_ene1,real(noverlap)/1000000.
    
        if (noverlap == 0) then
            write (*,*) ' No overlap, file ',trim(filename(ifile))//trim(suffix(isuffix))
            stop
        endif

        scale = scale / noverlap
        do ipart = minpart,maxpart
            do i = 1,npoints
                if (counts(i,ipart)>0) omega(i,ipart) = omega(i,ipart) - scale
            enddo
        enddo

        write (1,'(A,X,g12.4,X,2f20.6,X,f13.6,X,f12.6,X,f7.1)') trim(filename(ifile)), counts_file(ifile), -scale, 1./beta, mu, width, 1.

        do ipart = minpart,maxpart
            do i = 1,npoints
                if (tot_counts(i,ipart) == 0) then
                    entropy(i,ipart) = omega(i,ipart)
                    tot_counts(i,ipart) = counts(i,ipart)
                else 
                    if (counts(i,ipart)>0) then
                        ! new entropy is weighted mean
                        entropy(i,ipart) = ( tot_counts(i,ipart)*entropy(i,ipart) &
                            + counts(i,ipart)*omega(i,ipart) ) / (tot_counts(i,ipart) + counts(i,ipart))
                        tot_counts(i,ipart) = tot_counts(i,ipart) + counts(i,ipart)
                        if (ipart==0.and.i==3264) write (*,*) tot_counts(i,ipart)
                    endif
                endif
            enddo
        enddo                       

    endif

enddo   ! over isuffix (no ident)
enddo   ! over ifile

close(1)

write(nfile_fmt, '(a, i0, a)') '(A,', nfiles_suf, 'A4,A)'

open (25,file='phase.dat')
write (25, nfile_fmt) '/* Suffixes and Files = '&
    ,(suffix(isuffix),isuffix=1,nsuffix),(trim(filename(ifile)),ifile=1,nfiles),' */'
write (25,'(A)') '/*  T       mu       N_liq      N_gas       U_liq        U_gas    lnZ  */'
open (28,file='pvt.dat')
write (28, nfile_fmt) '/* Suffixes and Files = '&
    ,(suffix(isuffix),isuffix=1,nsuffix),(trim(filename(ifile)),ifile=1,nfiles),' */'
write (28,'(A)') '/*     T         mu          <N>           <U>          lnZ           Cv        dlnZ/dN     d2lnZ/dN2 */'
open (31,file='cmc.dat')
write (31, nfile_fmt) '/* Suffixes and Files = '&
    ,(suffix(isuffix),isuffix=1,nsuffix),(trim(filename(ifile)),ifile=1,nfiles),' */'
write (31,'(A)') '/*    T      min_slope    mu_cmc        N_cmc  */'

!shift_Z = 0.

100 continue

write (*,*) 'Enter T, Mu, Nmid (letters to stop, 0/neg Nmid for scan/phase equil.)'
read (*,*,err=200) beta,mu,nmid
beta = 1./beta
do_phase = .false.
do_scan = .false.
if (nmid<0) do_phase = .true.
nmid = abs(nmid)
if (nmid==0) do_scan = .true.
!beta = abs(beta)
end_mu = mu
count_pass = 1
incr_mu = 1.e-6
if (do_scan) then
    write (*,*) 'Enter ending Mu, increment'
    read (*,*,err=200) end_mu,incr_mu
    write (*,'(A)') '      T       mu     <N>   <E>/<N>   lnZ     Cv   dlnZ/dN  d2lnZ/dN2'
endif

do while ((mu<=end_mu.and.incr_mu>0).or.(mu>=end_mu.and.incr_mu<0)) 

    lnZ = -huge(lnZ)    

    do ipart = minpart,maxpart
        do i = 1,npoints
            if (tot_counts(i,ipart)>0) then
                lnZ = specexp(lnZ,entropy(i,ipart)-beta*(min_ene+(i-1)*width)+beta*mu*ipart)
!               exponent = entropy(i,ipart)-beta*(min_ene+(i-1)*width)+beta*mu*ipart + shift_z
!               the idea here is to shift the exponent by fixed amount, to preserve accuracy
!               if (exponent>maxexponent) then
!                   write (*,*) 'Overflow 0 - i,ipart,exponent=',i,ipart,exponent, tot_counts(i,ipart)
!                   stop
!               endif
!               if (exponent<-maxexponent) then
!                   exponent = -maxexponent
!               endif
!               Z = Z + dexp(exponent)
            endif
        enddo
!       write (*,*) Z,exponent,entropy(1,ipart),tot_counts(0,ipart)
    enddo

    open (24,file='dens.dat')
    write (24,'(A,f8.4,A,f12.4,A,20A5,3h *//3h /*,600A,3h *//3h /*,600A,3h */)') '/* T=',1/beta &
        ,'; mu=',mu,'; Suffixes and Files = '&
        ,(suffix(isuffix),isuffix=1,nsuffix),(trim(filename(ifile)),ifile=1,nfiles),' */'
    open (29,file='edens.dat')
    write (29,'(A,f8.4,A,f12.4,A,20A5,3h *//3h /*,600A,3h *//3h /*,600A,3h */)') '/* T=',1/beta &
        ,'; mu=',mu,'; Suffixes and Files = '&
        ,(suffix(isuffix),isuffix=1,nsuffix),(trim(filename(ifile)),ifile=1,nfiles), ' */'
    iter = 0
    dens = 0.
    edens = 0.
    avgnum = 0.
    nbelow = 0.
    ngas = 0.
    nliq = 0.
    ene_liq = 0.
    ene2 = 0.
    num2 = 0.
    ene_num = 0.
    do ipart = minpart,maxpart
        do i = 1,npoints
            if (tot_counts(i,ipart)>0) then
                exponent = entropy(i,ipart)-beta*(min_ene+(i-1)*width)+beta*mu*ipart
                !exponent = entropy(i,ipart)-beta*(min_ene+(i-1)*width)+beta*mu*ipart + shift_z
                !if (exponent > maxexponent) then
                !   write (*,*) 'Overflow 1 - i,ipart,exponent=',i,ipart,exponent
                !   stop
                !endif
                dens(ipart) = dens(ipart) + dexp(exponent-lnZ)
                edens(i) = edens(i) + dexp(exponent-lnZ)
                ene_liq = ene_liq + (min_ene+(i-1)*width)*dexp(exponent-lnZ)
                ene2 = ene2 + (min_ene+(i-1)*width)**2*dexp(exponent-lnZ)
                ene_num = ene_num + (min_ene+(i-1)*width)*ipart*dexp(exponent-lnZ)
                num2 = num2 + ipart*ipart*dexp(exponent-lnZ)
            endif
        enddo
        avgnum = avgnum + ipart * dens(ipart)
        if (ipart<=nmid) then
            nbelow = nbelow + dens(ipart)
            ngas = ngas + ipart * dens(ipart)
        else
            nliq = nliq + ipart * dens(ipart)
        endif
        if (dens(ipart)>1.D-99) write (24,'(i5,g12.3)') ipart,dens(ipart)
    enddo

    do i=1,npoints
        if (edens(i)>1.D-99) write (29,'(g12.3,g12.3)') min_ene+(i-1)*width,edens(i)
    enddo

    if (count_pass > 2) then  ! third or more pass
        mu_h0 = mu_h1 ; num_h0 = num_h1 ; Z_h0 = Z_h1
        mu_h1 = mu_h2 ; num_h1 = num_h2 ; Z_h1 = Z_h2
        mu_h2 = mu    ; num_h2 = avgnum ; Z_h2 = lnZ
        slope1 = slope
        slope = (Z_h2 - Z_h1) / (num_h2 - num_h1)
        slope2 = (slope - slope1) / (num_h2 - num_h0) * 2.
        if (slope < min_slope) min_slope = slope
        if (slope2 < slope2_min) then
            mu_cmc = mu_h1  ! because of lag of calculation of derivative
            num_cmc = num_h1
            slope2_min = slope2
        else if (incr_mu<0) then
            if (slope2_min<0) end_mu = mu
        endif
    else if (count_pass > 1) then  ! second pass
        mu_h2 = mu
        num_h2 = avgnum
        Z_h2 = lnZ
        !Z_h2 = log(Z)- shift_Z
        slope0 = (Z_h2 - Z_h1) / (num_h2 - num_h1)
        slope = slope0 ; slope2 = 0.
        min_slope = slope0
        slope2_min = 9999.
    else ! first pass
        mu_h1 = mu
        num_h1 = avgnum
        Z_h1 = lnZ
!       Z_h1 = log(Z)- shift_Z
        slope = 0. ; slope2 = 0.
    endif

     cv_part = -(((ene_num - ene_liq*avgnum)**2)/(num2-avgnum*avgnum) - (ene2 - ene_liq*ene_liq))/avgnum*beta*beta + 1.5
    !cv_part = (num2-avgnum*avgnum)/avgnum ! valid for athermal only
    ! (-  ene2 + ene_liq*ene_liq + (ene_num - ene_liq*avgnum)**2/(num2 - avgnum*avgnum))/avgnum
    if (do_scan) then
        write (*,'(2f9.4,f8.2,3f8.3,2f8.3)') 1/beta,mu,avgnum,ene_liq/avgnum,lnZ, cv_part,slope,slope2*10
!       write (*,'(2f9.4,f8.2,3f8.3,2f8.3)') 1/beta,mu,avgnum,ene_liq/avgnum,log(Z)-shift_Z, cv_part,slope,slope2*10
    else
        if (avgnum<1000.) then
            write (*,'(A,f6.2,A,f9.2,A,f4.0,A,f6.2,A,f8.2,A,f8.3,A,f10.3)') ' <N>=',avgnum,' <E>/<N>=',ene_liq/avgnum &
                             ,'; Frac.<',nmid,'=',nbelow,'; nliq=',nliq/(1.-nbelow),'; ngas=',ngas/nbelow,' lnZ=',lnZ
        else
            write (*,'(A,f7.1,A,f9.2,A,f5.0,A,f6.2,A,f8.1,A,f8.2,A,f10.3)') ' <N>=',avgnum,' <E>/<N>=',ene_liq/avgnum &
                             ,'; Frac.<',nmid,'=',nbelow,'; nliq=',nliq/(1.-nbelow),'; ngas=',ngas/nbelow,' lnZ=',lnZ
        endif
    endif
    write (28,'(8g13.5)') 1/beta,mu,avgnum,ene_liq,lnZ,cv_part,slope,slope2*10
    write (24,'(A,f7.3,A,f8.3,A,f4.0,A,f7.3,A,f9.3,A,f9.4,A,f8.3)') '/* <N> =',avgnum,'  <E>/<N> =',ene_liq/avgnum &
                       ,'; Frac. < ',nmid,'=',nbelow,'; nliq =',nliq/(1.-nbelow),'; ngas=',ngas/nbelow,' lnZ=',lnZ
    close(24)
    write (29,'(A,f7.3,A,f8.3,A,f4.0,A,f7.3,A,f9.3,A,f9.4,A,f8.3)') '/* <N> =',avgnum,'  <E>/<N> =',ene_liq/avgnum &
                       ,'; Frac. < ',nmid,'=',nbelow,'; nliq =',nliq/(1.-nbelow),'; ngas=',ngas/nbelow,' lnZ=',lnZ
    close(29)

!   if (log(Z)>maxexponent*2/3) &
!       write (*,*) 'Possible loss of accuracy; consider decreasing Mu increment'
!   shift_Z = shift_Z - log(Z)   ! a good choice for shift_Z
    mu = mu + incr_mu
    count_pass = count_pass + 1
enddo

if (do_scan) then
    write (*,'(a,3f12.4)') '  Min_slope, Mu_cmc, Num_cmc =',min_slope,mu_cmc,num_cmc
    write (31,'(f9.4,3f12.4)') 1/beta,min_slope,mu_cmc,num_cmc
endif

if (do_phase) then  
    write (*,*) 'Enter Tmin, Tmax, Tincr (Tmin < Tmax; Tincr can be + or -)'
    read (*,*) tmin,tmax,t_incr
    if (t_incr>0) then
        beta = 1/tmin
    else
        beta = 1/tmax
    endif
    iter = 0
    write (*,'(A)') '    T    Iter    Mu      nliq      ngas      lnZ'
    do while (1./beta<tmax+.00001.and.1./beta>tmin-.00001.and.iter<=maxiter)
        iter = 0
        nbelow = 0.
        nmid = ngas+nliq
        do while (abs(0.5-nbelow)>0.0005.and.iter<=maxiter)
            if (iter/=0) then
                mu = mu + muincr    ! arbitrary formula
                lnZ = -huge(lnZ)
                do ipart = minpart,maxpart
                    do i = 1,npoints
                        if (tot_counts(i,ipart)>0) then
                            exponent = entropy(i,ipart)-beta*(min_ene+(i-1)*width)+beta*mu*ipart
                            lnZ = specexp(lnZ,exponent)
                        endif
                    enddo
                enddo
                dens = 0.
                nbelow1 = 0.
                do ipart = minpart,maxpart
                    do i = 1,npoints
                        if (tot_counts(i,ipart)>0) then
                            exponent = entropy(i,ipart)-beta*(min_ene+(i-1)*width)+beta*mu*ipart 
                            dens(ipart) = dens(ipart) + dexp(exponent-lnZ)
                        endif
                    enddo
                    if (ipart<=nmid) then
                        nbelow1 = nbelow1 + dens(ipart)
                    endif
                enddo   
                ! the function being set to zero is ln(x) - ln(1-x), which
                ! behaves much better in convergence
                mu = mu - (log(nbelow1)-log(1.-nbelow1))*muincr &
                    / (log(nbelow1) - log(nbelow) - log(1.-nbelow1) + log(1.-nbelow))   ! Newton's formula

            endif
            lnZ = -huge(lnZ)
            do ipart = minpart,maxpart
                do i = 1,npoints
                    if (tot_counts(i,ipart)>0) then
                        exponent =  entropy(i,ipart)-beta*(min_ene+(i-1)*width)+beta*mu*ipart
                        lnZ = specexp(lnZ,exponent)
                    endif
                enddo
            enddo
            dens = 0.
            ngas = 0.
            nliq = 0.
            nbelow = 0.
            iter = iter + 1
            do ipart = minpart,maxpart
                do i = 1,npoints
                    if (tot_counts(i,ipart)>0) then
                        exponent = entropy(i,ipart)-beta*(min_ene+(i-1)*width)+beta*mu*ipart 
                        dens(ipart) = dens(ipart) + dexp(exponent-lnZ)
                    endif
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
            write (*,'(A,f7.3,A,I3,A,f7.2)') ' No convergence at T=',1/beta,' Iter =',iter,'; Mu=',mu
        else
            if (dens(int(nmid))>mindens*dens(int(nliq/(1.-nbelow))).or.dens(int(nmid))>mindens*dens(int(ngas/nbelow))) then
                write (*,'(A,i4/A,3(I4,f6.3))') ' Too close to Tc (or no convergence): peaks overlapping at nmid = ',int(nmid),&
                    ' Nliq, Nmid, Ngas and peak heights:',int(nliq/(1.-nbelow)),dens(int(nliq/(1.-nbelow))), &
                    int(nmid),dens(int(nmid)),int(ngas/nbelow),dens(int(ngas/nbelow))
            !   iter = maxiter+1
            else
                ! calculate average energy
                ene_liq = 0.
                ene_gas = 0.
                do ipart = minpart,maxpart
                    do i = 1,npoints
                        if (tot_counts(i,ipart)>0) then
                            exponent = entropy(i,ipart)-beta*(min_ene+(i-1)*width)+beta*mu*ipart 
                            if (ipart<=nmid) then
                                ene_gas = ene_gas + (min_ene+(i-1)*width)*dexp(exponent-lnZ)
                            else
                                ene_liq = ene_liq + (min_ene+(i-1)*width)*dexp(exponent-lnZ)                            
                            endif
                        endif
                    enddo
                enddo
                if (1/beta < 50) then               
                    write (*,'(f7.3,I3,f9.2,f10.3,f10.4,f8.3)') 1/beta,iter,mu,nliq/(1.-nbelow),ngas/nbelow,lnZ
                    write (25,'(f8.5,F11.5,f10.3,g13.5,2g12.4,f8.3)') 1/beta,mu, nliq/(1.-nbelow),ngas/nbelow &
                                                                    , ene_liq/nliq,ene_gas/ngas,lnZ
                else
                    write (*,'(f7.1,I3,f9.1,f10.3,f10.4,f8.3)') 1/beta,iter,mu,nliq/(1.-nbelow),ngas/nbelow,lnZ
                    write (25,'(f7.1,F10.1,f10.3,g12.4,f11.3,f11.4,f8.3)') 1/beta,mu, nliq/(1.-nbelow),ngas/nbelow &
                                                                         , ene_liq/nliq,ene_gas/ngas,lnZ
                endif
                    
            endif
        endif
        beta = 1./(1./beta + t_incr)
    enddo
endif   

goto 100
200 close(25)

end
