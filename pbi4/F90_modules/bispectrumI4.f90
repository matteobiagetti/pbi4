subroutine measure_bk(nbins,grid,kbin,kcenter,iopen,kf,npool,itot,npart,counts,dcl,bout)

    include 'fftw3.f'

    integer, intent(in) :: nbins,kbin,kcenter,npool,grid
    integer(kind=8), intent(in) :: npart,itot
    logical, intent(in) :: iopen
    real, intent(in) :: kf, counts(:)
    complex(kind=8), intent(in) :: dcl(:,:,:)
    real(kind=8), intent(out) :: bout(9,itot)

    integer(kind=8) :: planf, gridD
    integer ::  status
    real(kind=8) :: rmax
    integer :: i, j, l, m, n, id, jd, ld, nmodes,idist
    integer, allocatable :: nk(:), nbk(:), indx(:)
    real(kind=8), allocatable :: m1(:,:),map1(:), map2(:), pow(:), tk(:)
    real(kind=8), allocatable, dimension(:,:) :: mapk
    real(kind=8) :: di, dj, dl, avg, kf3
    real(kind=8) ::  sum, powSN, powSN2, bispSN

    write(*,*) '---FORTRAN SUBROUTINE---'
    call flush()
    write(*,*)
    call flush()
    write(*,*) 'internal grid', grid
    call flush()
    write(*,*) 'number of bins', nbins
    call flush()
    
    rmax=kcenter+dble(nbins-1)*kbin ! largest k (in units of kf)

    ! k central values in units of the fundamental frequency
    allocate( tk(nbins) )
    do i=1,nbins
        tk(i)=kcenter+dble(i-1)*kbin
    enddo
    if (iopen) then
       ishift=int(kcenter/kbin)+1
    else
       ishift=int(kcenter/kbin)
    endif

    gridD=grid**3

    kf3=kf**3
    powSN=1/dble(npart)/kf3
    powSN2=powSN*powSN

    print *, 'Allocate all arrays'
    call flush()
    allocate( nk(0:2*grid), nbk(0:2*grid+1) )
    allocate( indx(grid**3), map1(grid**3), map2(grid**3)  )
    allocate( m1(2,grid**3))   
    allocate( pow(nbins))
    allocate( mapk(gridD,nbins),STAT=status )
    if (status .ne. 0) then
        print *,'allocation failed!'
        stop
    endif

    write(*,*) 'Find modes of amplitude |k|'
    call flush()

    nk=0
    do i = 0, grid -1
        do j = 0, grid  -1
            do l = 0, grid  -1
                di = dble(min(i,grid-i))
                dj = dble(min(j,grid-j))
                dl = dble(min(l,grid-l))
                !idist = int(dsqrt(di*di+dj*dj+dl*dl)/kbin+0.5d0)
                idist=nint(dsqrt(di*di+dj*dj+dl*dl)/kbin-kcenter/kbin-tiny(kf))+1
                nk(idist) = nk(idist) + 1
            enddo 
        enddo 
    enddo
      
    nbk(0) = 0 
    do i = 0, 2*grid
        nbk(i+1) = nbk(i) + nk(i)
        nk(i) = 0 
    enddo 

      
    write(*,*) 'Save coordinates'
    call flush()  
    m = 0
    do i = 0, grid -1 
        do j = 0, grid -1
            do l = 0, grid -1
                di = dble(min(i,grid-i))
                dj = dble(min(j,grid-j))
                dl = dble(min(l,grid-l))
                !idist = int(dsqrt(di*di + dj*dj +dl*dl)/kbin+0.5d0)
                idist=nint(dsqrt(di*di+dj*dj+dl*dl)/kbin-kcenter/kbin-tiny(kf))+1
                nk(idist) = nk(idist) + 1  
                n = nbk(idist) + nk(idist)               
                m = m+ 1
                indx(n) = m 
            enddo 
        enddo 
    enddo 


    write(*,*) 'Make FFT plans'
    call flush()

    call dfftw_init_threads(iret)
    call dfftw_plan_with_nthreads(npool)
    call dfftw_plan_dft_3d(planf,grid,grid,grid,m1,m1,FFTW_FORWARD,FFTW_ESTIMATE)
      
    write(*,*) 'FFT Threads: ', npool
    call flush()

    n = 0 
    do i = 1, grid 
        do j = 1, grid 
            do l = 1, grid 
                n = n+1
                if (i .le. grid/2+1) then 
                    map1(n) = dble(dcl(i,j,l))
                    map2(n) = dimag(dcl(i,j,l))
                else 
                    id = mod(grid-i+1,grid)+1
                    jd = mod(grid-j+1,grid)+1
                    ld = mod(grid-l+1,grid)+1
                    map1(n) = dble(dcl(id,jd,ld))
                    map2(n) = -dimag(dcl(id,jd,ld))
                endif
                if (mod(i-1,grid/2)+mod(j-1,grid/2)+mod(l-1,grid/2).eq.0) then
                    map2(n) = 0.d0
                endif  
            enddo
        enddo 
    enddo
!      map2(1)=0.
!    deallocate(dcl)

    write(*,*) 'Calculate k maps'

    do i = 1, nbins ! nmin,nmax
        m1=0.d0
        nmodes = 0
        do n = nbk(i)+1, nbk(i+1)
            m1(1,indx(n)) = map1(indx(n))
            m1(2,indx(n)) = map2(indx(n))
            nmodes = nmodes +1 
        enddo 
        call dfftw_execute(planf)
        avg = 0.d0 
        do n = 1, gridD
            avg = avg + m1(1,n)*m1(1,n)
            mapk(n,i) = m1(1,n)
        enddo
        pow(i)=kf3*avg/dble(gridd)/dble(nmodes)-powSN
    enddo


    write(*,*) 'Last sum'
    
    write(*,'(a)', advance='no') 'Progress: |'
    call flush()
    m=0
    do l = 1, nbins
        do j = 1, l
            do i = max(1,l-j+1-ishift),j
                m=m+1
                sum = 0.d0 
                do n = 1, gridD
                   sum = sum + mapk(n,i)*mapk(n,j)*mapk(n,l)
                enddo
                bispSN=(pow(i)+pow(j)+pow(l))*powSN+powSN2
                bout(1,m) = tk(l)
                bout(2,m) = tk(j)
                bout(3,m) = tk(i)
                bout(4,m) = pow(l)
                bout(5,m) = pow(j)
                bout(6,m) = pow(i)
                bout(7,m) = kf3*sum/counts(m)
                bout(8,m) = bispSN
                bout(9,m) = real(counts(m)/dble(gridD))
                call flush()
                if (itot.lt.100) then
                    write(*,'(i3,a1)', advance='no')  int(dble(m)/dble(itot)*100.), '% '
                else if (mod(m,int(itot/100)).eq.0) then
                    write(*,'(i3,a1)', advance='no')   int(dble(m)/dble(itot)*100.), '% '
                    call flush()
                endif
!write(*,*) 'Triangle:', int(tk(l)),int(tk(j)),int(tk(i))
            enddo
        enddo
    enddo
    write(*,'(a)', advance='no') '|'
    call flush()
    deallocate(m1)

end subroutine

subroutine measure_bk_multipoles(nbins,grid,kbin,kcenter,iopen,kf,npool,itot,npart,counts,dcl,dcl2,dcl4,bout)

    include 'fftw3.f'

    integer, intent(in) :: nbins,kbin,kcenter,npool,grid
    integer(kind=8), intent(in) :: npart,itot
    logical, intent(in) :: iopen
    real, intent(in) :: kf,counts(:)
    complex(kind=8), intent(in) :: dcl(:,:,:),dcl2(:,:,:),dcl4(:,:,:)
    real(kind=8), intent(out) :: bout(11,itot)

    integer(kind=8) :: planf, gridD
    integer ::  status
    real(kind=8) :: rmax
    integer :: i, j, l, m, n, id, jd, ld, nmodes,idist
    integer, allocatable :: nk(:), nbk(:), indx(:)
    !complex(kind=8), allocatable :: dcl(:,:,:)
    real(kind=8), allocatable :: m1(:,:),map1(:), map2(:), pow(:), tk(:)
    real(kind=8), allocatable, dimension(:,:) :: mapk, mapk2, mapk4
    real(kind=8) :: di, dj, dl, avg, kf3
    real(kind=8) ::  sum, sum2, sum4, powSN, powSN2, bispSN

    write(*,*) '---FORTRAN SUBROUTINE---'
    call flush()
    write(*,*)
    call flush()
    write(*,*) 'internal grid', grid
    call flush()
    write(*,*) 'number of bins', nbins
    call flush()
    
    rmax=kcenter+dble(nbins-1)*kbin ! largest k (in units of kf)

    ! k central values in units of the fundamental frequency
    allocate( tk(nbins) )
    do i=1,nbins
        tk(i)=kcenter+dble(i-1)*kbin
    enddo
    if (iopen) then
       ishift=int(kcenter/kbin)+1
    else
       ishift=int(kcenter/kbin)
    endif

    gridD=grid**3

    kf3=kf**3
    powSN=1/dble(npart)/kf3
    powSN2=powSN*powSN

    print *, 'Allocate all arrays'
    call flush()
    allocate( nk(0:2*grid), nbk(0:2*grid+1) )
    allocate( indx(grid**3), map1(grid**3), map2(grid**3)  )
    allocate( m1(2,grid**3))   
    allocate( pow(nbins) )
    allocate( mapk(gridD,nbins),STAT=status )
    if (status .ne. 0) then
        print *,'allocation failed!'
        stop
    endif

    write(*,*) 'Find modes of amplitude |k|'
    call flush()

    nk=0
    do i = 0, grid -1
        do j = 0, grid  -1
            do l = 0, grid  -1
                di = dble(min(i,grid-i))
                dj = dble(min(j,grid-j))
                dl = dble(min(l,grid-l))
                !idist = int(dsqrt(di*di+dj*dj+dl*dl)/kbin+0.5d0)
                idist=nint(dsqrt(di*di+dj*dj+dl*dl)/kbin-kcenter/kbin-tiny(kf))+1
                nk(idist) = nk(idist) + 1
            enddo 
        enddo 
    enddo
      
    nbk(0) = 0 
    do i = 0, 2*grid
        nbk(i+1) = nbk(i) + nk(i)
        nk(i) = 0 
    enddo 

      
    write(*,*) 'Save coordinates'
    call flush()  
    m = 0
    do i = 0, grid -1 
        do j = 0, grid -1
            do l = 0, grid -1
                di = dble(min(i,grid-i))
                dj = dble(min(j,grid-j))
                dl = dble(min(l,grid-l))
                !idist = int(dsqrt(di*di + dj*dj +dl*dl)/kbin+0.5d0)
                idist=nint(dsqrt(di*di+dj*dj+dl*dl)/kbin-kcenter/kbin-tiny(kf))+1
                nk(idist) = nk(idist) + 1  
                n = nbk(idist) + nk(idist)               
                m = m+ 1
                indx(n) = m 
            enddo 
        enddo 
    enddo 


    write(*,*) 'Make FFT plans'
    call flush()

    call dfftw_init_threads(iret)
    call dfftw_plan_with_nthreads(npool)
    call dfftw_plan_dft_3d(planf,grid,grid,grid,m1,m1,FFTW_FORWARD,FFTW_ESTIMATE)
      
    write(*,*) 'FFT Threads: ', npool
    call flush()

    n = 0 
    do i = 1, grid 
        do j = 1, grid 
            do l = 1, grid 
                n = n+1
                if (i .le. grid/2+1) then 
                    map1(n) = dble(dcl(i,j,l))
                    map2(n) = dimag(dcl(i,j,l))
                else 
                    id = mod(grid-i+1,grid)+1
                    jd = mod(grid-j+1,grid)+1
                    ld = mod(grid-l+1,grid)+1
                    map1(n) = dble(dcl(id,jd,ld))
                    map2(n) = -dimag(dcl(id,jd,ld))
                endif
                if (mod(i-1,grid/2)+mod(j-1,grid/2)+mod(l-1,grid/2).eq.0) then
                    map2(n) = 0.d0
                endif  
            enddo
        enddo 
    enddo
!      map2(1)=0.
!    deallocate(dcl)

    write(*,*) 'Calculate k maps'

    do i = 1, nbins ! nmin,nmax
        m1=0.d0
        nmodes = 0
        do n = nbk(i)+1, nbk(i+1)
            m1(1,indx(n)) = map1(indx(n))
            m1(2,indx(n)) = map2(indx(n))
            nmodes = nmodes +1 
        enddo 
        call dfftw_execute(planf)
        avg = 0.d0 
        do n = 1, gridD
            avg = avg + m1(1,n)*m1(1,n)
            mapk(n,i) = m1(1,n)
        enddo
        pow(i)=kf3*avg/dble(gridd)/dble(nmodes)-powSN
    enddo

    write(*,*) 'Compute maps for delta_2'

    deallocate(map1,map2,m1)
    call dfftw_destroy_plan(planf)
    call dfftw_cleanup()

    allocate( map1(grid**3), map2(grid**3)  )
    allocate( m1(2,grid**3))  
    allocate( mapk2(gridD,nbins),STAT=status )
    if (status .ne. 0) then
        print *,'allocation failed!'
        stop
    endif

    call dfftw_init_threads(iret)
    call dfftw_plan_with_nthreads(npool)
    call dfftw_plan_dft_3d(planf,grid,grid,grid,m1,m1,FFTW_FORWARD,FFTW_ESTIMATE)

    n = 0 
    do i = 1, grid 
        do j = 1, grid 
            do l = 1, grid 
                n = n+1
                if (i .le. grid/2+1) then 
                    map1(n) = dble(dcl2(i,j,l))
                    map2(n) = dimag(dcl2(i,j,l))
                else 
                    id = mod(grid-i+1,grid)+1
                    jd = mod(grid-j+1,grid)+1
                    ld = mod(grid-l+1,grid)+1
                    map1(n) = dble(dcl2(id,jd,ld))
                    map2(n) = -dimag(dcl2(id,jd,ld))
                endif
                if (mod(i-1,grid/2)+mod(j-1,grid/2)+mod(l-1,grid/2).eq.0) then
                    map2(n) = 0.d0
                endif  
            enddo
        enddo 
    enddo
!      map2(1)=0.

    write(*,*) 'Calculate k maps'

    do i = 1, nbins ! nmin,nmax
        m1=0.d0
        do n = nbk(i)+1, nbk(i+1)
            m1(1,indx(n)) = map1(indx(n))
            m1(2,indx(n)) = map2(indx(n))
        enddo 
        call dfftw_execute(planf)
        do n = 1, gridD
            mapk2(n,i) = m1(1,n)
        enddo
    enddo

    write(*,*) 'Compute maps for delta_4'

    deallocate(map1,map2,m1)

    call dfftw_destroy_plan(planf)
    call dfftw_cleanup()

    allocate( map1(grid**3), map2(grid**3)  )
    allocate( m1(2,grid**3))  
    allocate( mapk4(gridD,nbins),STAT=status )
    if (status .ne. 0) then
        print *,'allocation failed!'
        stop
    endif

    call dfftw_init_threads(iret)
    call dfftw_plan_with_nthreads(npool)
    call dfftw_plan_dft_3d(planf,grid,grid,grid,m1,m1,FFTW_FORWARD,FFTW_ESTIMATE)


    n = 0 
    do i = 1, grid 
        do j = 1, grid 
            do l = 1, grid 
                n = n+1
                if (i .le. grid/2+1) then 
                    map1(n) = dble(dcl4(i,j,l))
                    map2(n) = dimag(dcl4(i,j,l))
                else 
                    id = mod(grid-i+1,grid)+1
                    jd = mod(grid-j+1,grid)+1
                    ld = mod(grid-l+1,grid)+1
                    map1(n) = dble(dcl4(id,jd,ld))
                    map2(n) = -dimag(dcl4(id,jd,ld))
                endif
                if (mod(i-1,grid/2)+mod(j-1,grid/2)+mod(l-1,grid/2).eq.0) then
                    map2(n) = 0.d0
                endif  
            enddo
        enddo 
    enddo
!      map2(1)=0.
!    deallocate(dcl)

    write(*,*) 'Calculate k maps'

    do i = 1, nbins ! nmin,nmax
        m1=0.d0
        do n = nbk(i)+1, nbk(i+1)
            m1(1,indx(n)) = map1(indx(n))
            m1(2,indx(n)) = map2(indx(n))
        enddo 
        call dfftw_execute(planf)
        do n = 1, gridD
            mapk4(n,i) = m1(1,n)
        enddo
    enddo

    write(*,*) 'Last sum'
    
    write(*,'(a)', advance='no') 'Progress: |'
    call flush()
    m=0
    do l = 1, nbins
        do j = 1, l
            do i = max(1,l-j+1-ishift),j
                m=m+1
                sum = 0.d0 
                sum2 = 0.d0
                sum4 = 0.d0
                do n = 1, gridD
                    sum = sum + mapk(n,i)*mapk(n,j)*mapk(n,l)
                    sum2 = sum2 + mapk(n,i)*mapk(n,j)*mapk2(n,l)      ! change to  mapk2(n,i)*mapk(n,j)*mapk(n,l) if you want the angle wrt shorter k
                    sum4 = sum4 + mapk(n,i)*mapk(n,j)*mapk4(n,l)      ! change to  mapk4(n,i)*mapk(n,j)*mapk(n,l) if you want the angle wrt shorter k
                enddo
                bout(1,m) = tk(l)
                bout(2,m) = tk(j)
                bout(3,m) = tk(i)
                bout(4,m) = pow(l)
                bout(5,m) = pow(j)
                bout(6,m) = pow(i)
                bout(7,m) = kf3*sum/counts(m)
                bout(8,m) = bispSN
                bout(9,m) = real(counts(m)/dble(gridD))
                bout(10,m) = 5*kf3*sum2/counts(m)
                bout(11,m) = 9*kf3*sum4/counts(m)
                call flush()
                if (itot.lt.100) then
                    write(*,'(i3,a1)', advance='no')  int(dble(m)/dble(itot)*100.), '% '
                else if (mod(m,int(itot/100)).eq.0) then
                    write(*,'(i3,a1)', advance='no')   int(dble(m)/dble(itot)*100.), '% '
                    call flush()
                endif
!write(*,*) 'Triangle:', int(tk(l)),int(tk(j)),int(tk(i))
            enddo
        enddo
    enddo
    write(*,'(a)', advance='no') '|'
    call flush()
    deallocate(m1)

end subroutine

subroutine measure_bk_cross(nbins,grid,kbin,kcenter,iopen,kf,npool,itot,npart,counts,dcl1,dcl2,bout)

    include 'fftw3.f'

    integer, intent(in) :: nbins,kbin,kcenter,npool,grid
    integer(kind=8), intent(in) :: npart,itot
    logical, intent(in) :: iopen
    real, intent(in) :: kf,counts(:)
    complex(kind=8), intent(in) :: dcl1(:,:,:),dcl2(:,:,:)
    real(kind=8), intent(out) :: bout(5,itot)

    integer(kind=8) :: planf, gridD
    integer ::  status,it(3)
    real(kind=8) :: rmax
    integer :: i, j, l, m, n, id, jd, ld, nmodes,idist
    integer, allocatable :: nk(:), nbk(:), indx(:)
    !complex(kind=8), allocatable :: dcl(:,:,:)
    real(kind=8), allocatable :: m1(:,:),map1(:), map2(:)
    real(kind=8), allocatable :: tk(:)
    real(kind=8), allocatable, dimension(:,:) :: mapk, mapkB
    real(kind=8) :: di, dj, dl, avg, avgC, kf3
    real(kind=8) ::  sum

    write(*,*) '---FORTRAN SUBROUTINE---'
    call flush()
    write(*,*)
    call flush()
    write(*,*) 'internal grid', grid
    call flush()
    write(*,*) 'number of bins', nbins
    call flush()
    
    rmax=kcenter+dble(nbins-1)*kbin ! largest k (in units of kf)

    ! k central values in units of the fundamental frequency
    allocate( tk(nbins) )
    do i=1,nbins
        tk(i)=kcenter+dble(i-1)*kbin
    enddo
    if (iopen) then
       ishift=int(kcenter/kbin)+1
    else
       ishift=int(kcenter/kbin)
    endif

    gridD=grid**3

    kf3=kf**3

    print *, 'Allocate all arrays'
    call flush()
    allocate( nk(0:2*grid), nbk(0:2*grid+1) )
    allocate( indx(grid**3), map1(grid**3), map2(grid**3)  )
    allocate( m1(2,grid**3))   
    allocate( mapk(gridD,nbins),STAT=status )
    if (status .ne. 0) then
        print *,'allocation failed!'
        stop
    endif

    write(*,*) 'Find modes of amplitude |k|'
    call flush()

    nk=0
    do i = 0, grid -1
        do j = 0, grid  -1
            do l = 0, grid  -1
                di = dble(min(i,grid-i))
                dj = dble(min(j,grid-j))
                dl = dble(min(l,grid-l))
                !idist = int(dsqrt(di*di+dj*dj+dl*dl)/kbin+0.5d0)
                idist=nint(dsqrt(di*di+dj*dj+dl*dl)/kbin-kcenter/kbin-tiny(kf))+1
                nk(idist) = nk(idist) + 1
            enddo 
        enddo 
    enddo
      
    nbk(0) = 0 
    do i = 0, 2*grid
        nbk(i+1) = nbk(i) + nk(i)
        nk(i) = 0 
    enddo 

      
    write(*,*) 'Save coordinates'
    call flush()  
    m = 0
    do i = 0, grid -1 
        do j = 0, grid -1
            do l = 0, grid -1
                di = dble(min(i,grid-i))
                dj = dble(min(j,grid-j))
                dl = dble(min(l,grid-l))
                !idist = int(dsqrt(di*di + dj*dj +dl*dl)/kbin+0.5d0)
                idist=nint(dsqrt(di*di+dj*dj+dl*dl)/kbin-kcenter/kbin-tiny(kf))+1
                nk(idist) = nk(idist) + 1  
                n = nbk(idist) + nk(idist)               
                m = m+ 1
                indx(n) = m 
            enddo 
        enddo 
    enddo 


    write(*,*) 'Make FFT plans'
    call flush()

    call dfftw_init_threads(iret)
    call dfftw_plan_with_nthreads(npool)
    call dfftw_plan_dft_3d(planf,grid,grid,grid,m1,m1,FFTW_FORWARD,FFTW_ESTIMATE)
      
    write(*,*) 'FFT Threads: ', npool
    call flush()

    n = 0 
    do i = 1, grid 
        do j = 1, grid 
            do l = 1, grid 
                n = n+1
                if (i .le. grid/2+1) then 
                    map1(n) = dble(dcl1(i,j,l))
                    map2(n) = dimag(dcl1(i,j,l))
                else 
                    id = mod(grid-i+1,grid)+1
                    jd = mod(grid-j+1,grid)+1
                    ld = mod(grid-l+1,grid)+1
                    map1(n) = dble(dcl1(id,jd,ld))
                    map2(n) = -dimag(dcl1(id,jd,ld))
                endif
                if (mod(i-1,grid/2)+mod(j-1,grid/2)+mod(l-1,grid/2).eq.0) then
                    map2(n) = 0.d0
                endif  
            enddo
        enddo 
    enddo
!      map2(1)=0.
!    deallocate(dcl)

    write(*,*) 'Calculate k maps'

    do i = 1, nbins ! nmin,nmax
        m1=0.d0
        nmodes = 0
        do n = nbk(i)+1, nbk(i+1)
            m1(1,indx(n)) = map1(indx(n))
            m1(2,indx(n)) = map2(indx(n))
        enddo 
        call dfftw_execute(planf)
        avg = 0.d0 
        do n = 1, gridD
            avg = avg + m1(1,n)*m1(1,n)
            mapk(n,i) = m1(1,n)
        enddo
    enddo

    write(*,*) 'Compute maps for deltaB'

    deallocate(map1,map2,m1)
    call dfftw_destroy_plan(planf)
    call dfftw_cleanup()

    allocate( map1(grid**3), map2(grid**3)  )
    allocate( m1(2,grid**3))  
    allocate( mapkB(gridD,nbins),STAT=status )
    if (status .ne. 0) then
        print *,'allocation failed!'
        stop
    endif

    call dfftw_init_threads(iret)
    call dfftw_plan_with_nthreads(npool)
    call dfftw_plan_dft_3d(planf,grid,grid,grid,m1,m1,FFTW_FORWARD,FFTW_ESTIMATE)

    n = 0 
    do i = 1, grid 
        do j = 1, grid 
            do l = 1, grid 
                n = n+1
                if (i .le. grid/2+1) then 
                    map1(n) = dble(dcl2(i,j,l))
                    map2(n) = dimag(dcl2(i,j,l))
                else 
                    id = mod(grid-i+1,grid)+1
                    jd = mod(grid-j+1,grid)+1
                    ld = mod(grid-l+1,grid)+1
                    map1(n) = dble(dcl2(id,jd,ld))
                    map2(n) = -dimag(dcl2(id,jd,ld))
                endif
                if (mod(i-1,grid/2)+mod(j-1,grid/2)+mod(l-1,grid/2).eq.0) then
                    map2(n) = 0.d0
                endif  
            enddo
        enddo 
    enddo
!      map2(1)=0.

    write(*,*) 'Calculate k maps'

    do i = 1, nbins ! nmin,nmax
        m1=0.d0
        do n = nbk(i)+1, nbk(i+1)
            m1(1,indx(n)) = map1(indx(n))
            m1(2,indx(n)) = map2(indx(n))
        enddo 
        call dfftw_execute(planf)
        do n = 1, gridD
            mapkB(n,i) = m1(1,n)
        enddo
    enddo

    write(*,*) 'Last sum'
    
    write(*,'(a)', advance='no') 'Progress: |'
    call flush()
    m=0
    do l = 1, nbins
        do j = 1, l
            do i = max(1,l-j+1-ishift),min(nbins,j+l-1+ishift)
                m=m+1
                sum = 0.d0
                do n = 1, gridD
                    sum = sum + mapkB(n,i)*mapk(n,j)*mapk(n,l)
                enddo
                bout(1,m) = tk(l)
                bout(2,m) = tk(j)
                bout(3,m) = tk(i)
                bout(4,m) = kf3*sum/counts(m)
                bout(5,m) = real(counts(m)/dble(gridD))
                call flush()
                if (itot.lt.100) then
                    write(*,'(i3,a1)', advance='no')  int(dble(m)/dble(itot)*100.), '% '
                else if (mod(m,int(itot/100)).eq.0) then
                    write(*,'(i3,a1)', advance='no')   int(dble(m)/dble(itot)*100.), '% '
                    call flush()
                endif
!write(*,*) 'Triangle:', int(tk(l)),int(tk(j)),int(tk(i))
            enddo
        enddo
    enddo
    write(*,'(a)', advance='no') '|'
    call flush()
    deallocate(m1)

end subroutine

subroutine count_triangles(nbins,grid,kbin,kcenter,iopen,cross,npool,itot,counts)

      include 'fftw3.f'

      integer, intent(in) :: nbins,grid,kbin,kcenter,npool
      integer(kind=8), intent(in) :: itot
      logical, intent(in) :: iopen,cross
      real(kind=8), intent(out) :: counts(itot)

      integer(kind=8) :: planf,gridD
      integer :: status,ishift
      real(kind=8) :: rmax,di,dj,dl,sum
      integer :: i, j, l, m, n,idist, iret
      integer, allocatable :: nk(:), nbk(:), indx(:)
      real(kind=8), allocatable :: m1(:,:),mapk(:,:)
     
      rmax=kcenter+dble(nbins-1)*kbin ! largest k (in units of kf)

      if (iopen) then
         ishift=int(kcenter/kbin)+1
      else
         ishift=int(kcenter/kbin)
      endif

      gridD=grid**3

      write(*,*) '---FORTRAN SUBROUTINE---'
      write(*,*) 'FFT grid = ',grid
      call flush()

      print *, 'Allocate all arrays'
      call flush()
      allocate( nk(0:2*grid), nbk(0:2*grid+1) )
      allocate( indx(grid**3) )
      allocate( m1(2,grid**3))
      allocate( mapk(gridD,nbins),STAT=status )
      if (status .ne. 0) then
         print *,'allocation failed!'
         stop
      endif
      
      write(*,*) 'Find modes of amplitude |k|'
      call flush()

      nk=0
      do i = 0, grid -1
         do j = 0, grid  -1
            do l = 0, grid  -1
               di = dble(min(i,grid-i))
               dj = dble(min(j,grid-j))
               dl = dble(min(l,grid-l))
               !idist = int(dsqrt(di*di+dj*dj+dl*dl)/kbin+0.5d0)
               idist=nint(dsqrt(di*di+dj*dj+dl*dl)/kbin-kcenter/kbin-tiny(rmax))+1
               nk(idist) = nk(idist) + 1
            enddo 
         enddo 
      enddo
      
      nbk(0) = 0 
      do i = 0, 2*grid
         nbk(i+1) = nbk(i) + nk(i)
         nk(i) = 0 
      enddo 

      
      write(*,*) 'Save coordinates'
      call flush()  
      m = 0
      do i = 0, grid -1 
         do j = 0, grid -1
            do l = 0, grid -1
               di = dble(min(i,grid-i))
               dj = dble(min(j,grid-j))
               dl = dble(min(l,grid-l))
               !idist = int(dsqrt(di*di + dj*dj +dl*dl)/kbin+0.5d0)
               idist=nint(dsqrt(di*di+dj*dj+dl*dl)/kbin-kcenter/kbin-tiny(rmax))+1
               nk(idist) = nk(idist) + 1  
               n = nbk(idist) + nk(idist)               
               m = m+ 1
               indx(n) = m 
            enddo 
         enddo 
      enddo 
      deallocate(nk)

      write(*,*) 'Make FFT plans'
      call flush()

      call dfftw_init_threads(iret)
      call dfftw_plan_with_nthreads(npool)
      call dfftw_plan_dft_3d(planf,grid,grid,grid,m1,m1,FFTW_FORWARD,FFTW_ESTIMATE)
      
      write(*,*) 'FFT Threads: ', npool
      call flush()
        
      write(*,*) 'Calculate k maps'
      call flush()

      do i = 1,nbins ! nmin,nmax
         do n = 1, gridD 
            m1(1,n) = 0.
            m1(2,n) = 0.
         enddo
         do n = nbk(i)+1, nbk(i+1)
            m1(1,indx(n)) = 1.
         enddo
         call dfftw_execute(planf)
         do n = 1, gridD
             mapk(n,i) = m1(1,n)
         enddo
      enddo

      if (cross) then
        write(*,*) 'Compute counts for cross bispectrum'
        call flush()
      else
        write(*,*) 'Compute counts'
        call flush()
      endif
         
      m=0
      do l = 1, nbins
        do j = 1, l
            if (cross) then
                imax = min(nbins,j+l-1+ishift)
            else 
                imax = j
            endif
            do i = max(1,l-j+1-ishift),imax
               m=m+1
               sum=0.
               do n = 1, gridD
                  sum=sum+mapk(n,i)*mapk(n,j)*mapk(n,l)
               enddo
               counts(m)=sum
               call flush()
               if (itot.lt.100) then
                  write(*,'(i3,a1)', advance='no')  int(dble(m)/dble(itot)*100.), '% '
               else if (mod(m,int(itot/100)).eq.0) then
                  write(*,'(i3,a1)', advance='no')   int(dble(m)/dble(itot)*100.), '% '
                  call flush()
               endif
            enddo
         enddo
      enddo
      write(*,'(a)', advance='no') '|'
      call flush()

      write(*,*) 'Done'
      call flush()

end subroutine