! *******************************************************************
! *** Direct assignment *********************************************
! 
! 
! subroutine assign_direct(grid, r, dcl)
! 
!     integer, intent(in) :: grid
!     real(kind=8), intent(in) :: r(3)
!     complex(kind=8), intent(out) :: dcl(grid, grid, grid)
!     integer :: ix,iy,iz,icy,icz,ikx,iky,ikz
!     real(kind=8) :: rkx,rky,rkz,product,sp,cp,rtpi(3)
!     real, parameter :: twopi = 6.28319
!   
!     Nnyqx=grid/2+1
!     Nnyqy=grid/2+1
!     Nnyqz=grid/2+1
!     Nx=grid
!     Ny=grid
!     Nz=grid
!     rtpi=r*twopi
!    
!     do iz=1,Nnyqz
! 		ikz=mod(iz+Nz/2-2,Nz)-Nz/2+1
! 		rkz=dble(ikz)*rtpi(3)
! 		icz=mod(Nz-iz+1,Nz)+1
! 	  
! 		do iy=1,Nnyqy
! 			iky=mod(iy+Ny/2-2,Ny)-Ny/2+1
! 			rky=dble(iky)*rtpi(2)
! 			icy=mod(Ny-iy+1,Ny)+1
! 
! 			do ix=1,Nnyqx
! 				ikx=mod(ix+Nx/2-2,Nx)-Nx/2+1
! 				rkx=dble(ikx)*rtpi(1)
! 		   
! 				product=rkx+rky+rkz
! 				cp=dcos(product)
! 				sp=dsin(product)
! 				dcl(ix,iy ,iz ) = dcl(ix,iy ,iz ) + dcmplx(cp,sp)
! 
! 				if (iz.ne.Nnyqz .and. iz.ne.1) then                   
! 					product=rkx+rky-rkz
! 					cp=dcos(product)
! 					sp=dsin(product)
! 					dcl(ix,iy ,icz) = dcl(ix,iy ,icz) + dcmplx(cp,sp)
! 				endif
! 			   
! 				if (iy.ne.Nnyqy .and. iy.ne.1) then                   
! 					product=rkx-rky+rkz
! 					cp=dcos(product)
! 					sp=dsin(product)
! 					dcl(ix,icy,iz ) = dcl(ix,icy,iz ) + dcmplx(cp,sp)
! 				endif
! 			   
! 				if (iz.ne.Nnyqz .and. iy.ne.Nnyqy  .and. iz.ne.1 .and. iy.ne.1) then                   
! 					product=rkx-rky-rkz
! 					cp=dcos(product)
! 					sp=dsin(product)
! 					dcl(ix,icy,icz) = dcl(ix,icy,icz) + dcmplx(cp,sp)      
! 				endif
! 			
! 			enddo
! 		enddo
! 	enddo
! 	return
! end subroutine

subroutine assign_single_grid( grid,nintp, npart, pos, weight, dtl)

 integer, intent(in) ::  grid, nintp, npart
 real(kind=8), intent(in) :: pos(:,:),weight(:)
 real, parameter :: twopi = 6.28319
 real(kind=8), intent(out) :: dtl(grid, grid, grid)
 integer :: ix,iy,iz,ivr(3),iv(3,4),i,j,Nv(3), nintph
 real(kind=8) :: w(3,4),vr(3), h, h2
 
 if (nintp.ge.3) nintph=3
 if (nintp.lt.3) nintph=1
 
 Nv = (/grid,grid,grid/)
 w=0.d0
 
 do i=1,npart
 
    vr=dble(Nv)*pos(:,i)+1.d0

    if (nintp.eq.4) then ! PCS 
       ivr=int(vr)
       do j=1,nintp
          iv(:,j) = mod(ivr(:)-nintph+j+Nv(:),Nv(:))+1
       enddo
       do j=1,3
          h=vr(j)-int(vr(j))
          h2=h*h
          w(j,1)=(1.d0-h)**3/6.d0
          w(j,2)=4.d0/6.d0+(0.5d0*h-1.d0)*h2
          w(j,4)=h2*h/6.d0
          w(j,3)=1.d0-w(j,1)-w(j,2)-w(j,4)
       enddo
    endif   
    if (nintp.eq.3) then ! TSC
       ivr=nint(vr)
       do j=1,nintp
          iv(:,j) = mod(ivr(:)-nintph+j+Nv(:),Nv(:))+1
       enddo
       do j=1,3
          h=vr(j)-nint(vr(j))
          h2=h*h
          w(j,1)=0.5d0*(0.5d0-h)**2
          w(j,2)=0.75d0-h2
          w(j,3)=1.d0-w(j,1)-w(j,2)
       enddo
    endif   
    if (nintp.eq.2) then ! CIC
       ivr=int(vr)
       do j=1,nintp
          iv(:,j) = mod(ivr(:)-nintph+j+Nv(:),Nv(:))+1
       enddo
       do j=1,3
          h=vr(j)-nint(vr(j))
          w(j,1)=1.0d0-h
          w(j,2)=h
       enddo
    endif   
    if (nintp.eq.1) then ! CIC
       w=1.0d0
       ivr=nint(vr)
       iv(:,1) = mod(ivr(:)-1+Nv(:),Nv(:))+1
       
    endif
    
    do ix=1,nintp
       do iy=1,nintp
          do iz=1,nintp
             dtl(iv(1,ix),iv(2,iy),iv(3,iz))=dtl(iv(1,ix),iv(2,iy),iv(3,iz))+w(1,ix)*w(2,iy)*w(3,iz)*weight(i)
          enddo
       enddo
    enddo
 enddo
end subroutine

subroutine assign_double_grid( grid,nintp, npart, pos, weight, dcl)

 integer, intent(in) ::  grid, nintp, npart
 real(kind=8), intent(in) :: pos(:,:),weight(:)
 real, parameter :: twopi = 6.28319
 complex(kind=8), intent(out) :: dcl(grid, grid, grid)
 complex(kind=8) :: cpx
 integer :: ix,iy,iz,ivr(3), ivt(3), iv(3,4),jv(3,4), i,j,Nv(3), nintph
 real(kind=8) :: w1(3,4), w2(3,4), vr(3), vt(3), h1, h1sq, h2,h2sq
 
 if (nintp.ge.3) nintph=3
 if (nintp.lt.3) nintph=1
 
 Nv = (/grid,grid,grid/)
 w1=0.d0
 w2=0.d0
 
 do i=1,npart
 
    vr=dble(Nv)*pos(:,i)+1.d0
    vt=vr-0.5d0      !had to change from +0.5 to -0.5 from the fortran code, reason unknown

    if (nintp.eq.4) then ! PCS 
       ivr=int(vr)
       ivt=int(vt)
       do j=1,nintp
          iv(:,j) = mod(ivr(:)-nintph+j+Nv(:),Nv(:))+1
          jv(:,j) = mod(ivt(:)-nintph+j+Nv(:),Nv(:))+1
       enddo
       do j=1,3
          h1=vr(j)-int(vr(j))
          h2=vt(j)-int(vt(j))
          h1sq=h1*h1
          h2sq=h2*h2
          w1(j,1)=(1.d0-h1)**3/6.d0
          w1(j,2)=4.d0/6.d0+(0.5d0*h1-1.d0)*h1sq
          w1(j,4)=h1sq*h1/6.d0
          w1(j,3)=1.d0-w1(j,1)-w1(j,2)-w1(j,4)
          
          w2(j,1)=(1.d0-h2)**3/6.d0
          w2(j,2)=4.d0/6.d0+(0.5d0*h2-1.d0)*h2sq
          w2(j,4)=h2sq*h2/6.d0
          w2(j,3)=1.d0-w2(j,1)-w2(j,2)-w2(j,4)
       enddo
    endif   
    if (nintp.eq.3) then ! TSC
       ivr=nint(vr)
       ivt=nint(vt)
       do j=1,nintp
          iv(:,j) = mod(ivr(:)-nintph+j+Nv(:),Nv(:))+1
          jv(:,j) = mod(ivt(:)-nintph+j+Nv(:),Nv(:))+1
       enddo
       do j=1,3
          h1=vr(j)-nint(vr(j))
          h2=vt(j)-nint(vt(j))
          h1sq=h1*h1
          h2sq=h2*h2
          w1(j,1)=0.5d0*(0.5d0-h1)**2
          w1(j,2)=0.75d0-h1sq
          w1(j,3)=1.d0-w1(j,1)-w1(j,2)
          
          w2(j,1)=0.5d0*(0.5d0-h2)**2
          w2(j,2)=0.75d0-h2sq
          w2(j,3)=1.d0-w2(j,1)-w2(j,2)
       enddo
    endif   
    if (nintp.eq.2) then ! CIC
       ivr=int(vr)
       ivt=int(vt)
       do j=1,nintp
          iv(:,j) = mod(ivr(:)-nintph+j+Nv(:),Nv(:))+1
          jv(:,j) = mod(ivt(:)-nintph+j+Nv(:),Nv(:))+1
       enddo
       do j=1,3
          h1=vr(j)-int(vr(j))
          h2=vt(j)-int(vt(j))
          w1(j,1)=1.0d0-h1
          w1(j,2)=h1

          w2(j,1)=1.0d0-h2
          w2(j,2)=h2
       enddo
    endif   
    if (nintp.eq.1) then ! NGP
       w1=1.0d0
       w2=1.0d0
       ivr=nint(vr)
       ivt=nint(vt)
       iv(:,1) = mod(ivr(:)-1+Nv(:),Nv(:))+1
       jv(:,1) = mod(ivt(:)-1+Nv(:),Nv(:))+1
       
    endif
    
    do ix=1,nintp
       do iy=1,nintp
          do iz=1,nintp
          
             cpx=dcmplx(w1(1,ix)*w1(2,iy)*w1(3,iz)*weight(i),0.d0)
             dcl(iv(1,ix),iv(2,iy),iv(3,iz))=dcl(iv(1,ix),iv(2,iy),iv(3,iz))+cpx

             cpx=dcmplx(0.d0,w2(1,ix)*w2(2,iy)*w2(3,iz)*weight(i))
             dcl(jv(1,ix),jv(2,iy),jv(3,iz))=dcl(jv(1,ix),jv(2,iy),jv(3,iz))+cpx

          enddo
       enddo
    enddo
 enddo
end subroutine

subroutine fcomb(grid, nintp, interlacing, box, fdcl, dcl)

   integer, intent(in) :: grid, nintp
   real, intent(in) :: box
   logical, intent(in) :: interlacing
   complex(kind=8), intent(in) :: fdcl(:,:,:)
   complex(kind=8), intent(out) :: dcl(grid,grid,grid)
   integer :: Nnyqx,Nnyqy,Nnyqz,icx,icy,icz,ix,iy,iz,i,Nx,Ny,Nz
   real(kind=8) :: tpiNx,piNx,tpiNy,piNy,tpiNz,piNz,cf,rk,Wkx,Wky,Wkz,cfac,kF3
   complex(kind=8) :: recx,recy,recz,xrec,yrec,zrec
   complex(kind=8) :: c1,ci,c000,c001,c010,c011,cma,cmb,cmc,cmd
   real(kind=8) :: tWkx(grid/2+1),tWky(grid/2+1),tWkz(grid/2+1)
   real, parameter :: twopi = 6.28319

   Nnyqx=grid/2+1
   tpiNx=twopi/dble(grid)
   Nnyqy=grid/2+1
   tpiNy=twopi/dble(grid)
   Nnyqz=grid/2+1
   tpiNz=twopi/dble(grid)
   kF3=twopi*twopi*twopi/box/box/box
   Nx=grid
   Ny=grid
   Nz=grid
   
   dcl = fdcl

   tWkx(1)=1.d0
   do i=2,Nnyqx
      rk=tpiNx*dble(i-1)
      tWkx(i)=(dsin(rk/2.d0)/(rk/2.d0))**nintp
   enddo
   tWky(1)=1.d0
   do i=2,Nnyqy
      rk=tpiNy*dble(i-1)
      tWky(i)=(dsin(rk/2.d0)/(rk/2.d0))**nintp
   enddo
   tWkz(1)=1.d0
   do i=2,Nnyqz
      rk=tpiNz*dble(i-1)
      tWkz(i)=(dsin(rk/2.d0)/(rk/2.d0))**nintp
   enddo

   if (interlacing) then

      cf=1.d0/(4.d0*kF3)

      piNx=-tpiNx/2.d0
      recx=dcmplx(dcos(piNx),dsin(piNx))
      piNy=-tpiNy/2.d0
      recy=dcmplx(dcos(piNy),dsin(piNy))
      piNz=-tpiNz/2.d0
      recz=dcmplx(dcos(piNz),dsin(piNz))

      c1=dcmplx(1.d0,0.d0)
      ci=dcmplx(0.d0,1.d0)

      zrec=c1
      do iz=1,Nnyqz
         icz=mod(Nz-iz+1,Nz)+1
         Wkz=tWkz(iz)
         
         yrec=c1
         do iy=1,Nnyqy
            icy=mod(Ny-iy+1,Ny)+1
            Wky=tWky(iy)

            xrec=c1
            do ix=1,Nnyqx
               icx=mod(Nx-ix+1,Nx)+1
               Wkx=tWkx(ix)
               
               cfac=cf/(Wkx*Wky*Wkz)
               
               cma=ci*xrec*yrec*zrec
               cmb=ci*xrec*yrec*dconjg(zrec)
               cmc=ci*xrec*dconjg(yrec)*zrec
               cmd=ci*xrec*dconjg(yrec*zrec)
               
               c000=dcl(ix,iy ,iz )*(c1-cma)+dconjg(dcl(icx,icy,icz))*(c1+cma)
               c001=dcl(ix,iy ,icz)*(c1-cmb)+dconjg(dcl(icx,icy,iz ))*(c1+cmb)
               c010=dcl(ix,icy,iz )*(c1-cmc)+dconjg(dcl(icx,iy ,icz))*(c1+cmc)
               c011=dcl(ix,icy,icz)*(c1-cmd)+dconjg(dcl(icx,iy ,iz ))*(c1+cmd)
               
               dcl(ix,iy ,iz )=c000*cfac
               dcl(ix,iy ,icz)=c001*cfac
               dcl(ix,icy,iz )=c010*cfac
               dcl(ix,icy,icz)=c011*cfac
               dcl(icx,iy ,iz )=dconjg(dcl(ix,icy,icz))
               dcl(icx,iy ,icz)=dconjg(dcl(ix,icy,iz ))
               dcl(icx,icy,iz )=dconjg(dcl(ix,iy ,icz))
               dcl(icx,icy,icz)=dconjg(dcl(ix,iy ,iz ))
               
               xrec=xrec*recx
            enddo
            yrec=yrec*recy
         enddo
         zrec=zrec*recz
      enddo

   else ! no interlacing, only correcting for the window function and normalizing

   cf=1.d0/kF3

   do iz=1,Nnyqz
      icz=mod(Nz-iz+1,Nz)+1
      Wkz=tWkz(iz)
      
      do iy=1,Nnyqy
         icy=mod(Ny-iy+1,Ny)+1
         Wky=tWky(iy)
         
         do ix=1,Nnyqx
            Wkx=tWkx(ix)
            
            cfac=cf/(Wkx*Wky*Wkz)
            dcl(ix,iy ,iz )=dcl(ix,iy ,iz )*cfac
            if(iz.ne.icz) dcl(ix,iy,icz)=dcl(ix,iy,icz)*cfac
            if(iy.ne.icy) dcl(ix,icy,iz)=dcl(ix,icy,iz)*cfac
            if(iz.ne.icz .and. iy.ne.icy) dcl(ix,icy,icz)=dcl(ix,icy,icz)*cfac
            
         enddo
      enddo
   enddo

   endif

end subroutine

subroutine reduce_array(L,din,grid,dout)

    integer,intent(in) :: L,grid
    integer :: ix,iy,iz,icy,icz,jcy,jcz
    complex (kind=8), intent(in) :: din(:,:,:)
    complex (kind=8), intent(out) :: dout(grid,grid,grid)


    do iz=1,grid/2+1
        icz=mod(L-iz+1,L)+1
        jcz=mod(grid-iz+1,grid)+1

        do iy=1,grid/2+1
            icy=mod(L-iy+1,L)+1
            jcy=mod(grid-iy+1,grid)+1
            
            do ix=1,grid/2+1

               dout(ix,iy ,iz ) = din(ix,iy ,iz )
               dout(ix,jcy,iz ) = din(ix,icy,iz )
               dout(ix,iy ,jcz) = din(ix,iy ,icz)
               dout(ix,jcy,jcz) = din(ix,icy,icz)

            enddo
        enddo
    enddo

end subroutine reduce_array

subroutine compute_q2k(grid,box,deltak,i,j,Q2k)

    integer, intent(in) :: grid,i,j
    real, intent(in) :: box
    complex(kind=8), intent(in) :: deltak(:,:,:)
    complex(kind=8), intent(out) :: Q2k(grid,grid,grid)
    integer :: ix,iy,iz,ikx,iky,ikz,Nx, Ny, Nz
    real(kind=8) :: rk2,v(3),kFx, kFy, kFz, kNx, kNy, kNz, kF3
    real, parameter :: twopi = 6.28319

    Nx = grid
    Ny = grid
    Nz = grid

    kFx = twopi/box   ! fundamental frequencies
    kFy = twopi/box
    kFz = twopi/box
    kF3 = kFx*kFy*kFz
  
    kNx=kFx*dble(Nx)/2.d0    ! Nyquist frequencies
    kNy=kFy*dble(Ny)/2.d0 
    kNz=kFz*dble(Nz)/2.d0

    Q2k = deltak

    do iz=1,Nz
        ikz=mod(iz+Nz/2-2,Nz)-Nz/2+1
        v(3)=kFz*dble(ikz)

        do iy=1,Ny
            iky=mod(iy+Ny/2-2,Ny)-Ny/2+1
            v(2)=kFy*dble(iky)

            do ix=1,Nx/2+1
            ikx=mod(ix+Nx/2-2,Nx)-Nx/2+1
            !icx=mod(Nx-ix+1,Nx)+1
            v(1)=kFx*dble(ikx)
            
            if (ix+iy+iz.gt.3) then
                rk2=dot_product(v,v)
                Q2k(ix,iy,iz)=Q2k(ix,iy,iz)*v(i)*v(j)/rk2
            else
                Q2k(ix,iy,iz)=dcmplx(0.d0,0.d0)
            endif
            
            enddo
        enddo
    enddo

end subroutine

subroutine compute_q4k(grid,box,deltak,i,j,l,k,Q4k)

    integer, intent(in) :: grid,i,j,l,k
    real, intent(in) :: box
    complex(kind=8), intent(in) :: deltak(:,:,:)
    complex(kind=8), intent(out) :: Q4k(grid,grid,grid)
    integer :: ix,iy,iz,ikx,iky,ikz,Nx, Ny, Nz
    real(kind=8) :: rk2,v(3),kFx, kFy, kFz, kNx, kNy, kNz, kF3
    real, parameter :: twopi = 6.28319

    Nx = grid
    Ny = grid
    Nz = grid

    kFx = twopi/box   ! fundamental frequencies
    kFy = twopi/box
    kFz = twopi/box
    kF3 = kFx*kFy*kFz
  
    kNx=kFx*dble(Nx)/2.d0    ! Nyquist frequencies
    kNy=kFy*dble(Ny)/2.d0 
    kNz=kFz*dble(Nz)/2.d0

    Q4k = deltak

    do iz=1,Nz
        ikz=mod(iz+Nz/2-2,Nz)-Nz/2+1
        v(3)=kFz*dble(ikz)

        do iy=1,Ny
        iky=mod(iy+Ny/2-2,Ny)-Ny/2+1
        v(2)=kFy*dble(iky)

         do ix=1,Nx/2+1
               ikx=mod(ix+Nx/2-2,Nx)-Nx/2+1
               v(1)=kFx*dble(ikx)
               
               if (ix+iy+iz.gt.3) then
                  rk2=dot_product(v,v)**2
                  Q4k(ix,iy,iz)=Q4k(ix,iy,iz)*v(i)*v(j)*v(l)*v(k)/rk2
               else
                  Q4k(ix,iy,iz)=dcmplx(0.d0,0.d0)
               endif
            
         enddo
        enddo
    enddo

end subroutine

subroutine measure_pk(grid, box, kbin, kcenter, nbins, dcl1, dcl2, pout)

    integer, intent(in) :: grid, nbins
    real, intent(in) :: box, kbin, kcenter
    complex(kind=8), intent(in) :: dcl1(:,:,:), dcl2(:,:,:)
    real(kind=8), intent(out) :: pout(5,nbins)
    integer :: ix, iy, iz, ikx, iky, ikz, icx, icy, icz, imk, Nx, Ny, Nz
    real (kind=8), dimension(nbins) :: avgk, avgP, co
    real (kind=8) :: kFmin,rk,kNmax,dk,kFx, kFy, kFz, kNx, kNy, kNz, kF3
    real (kind=8) :: Le2, Le4, sp, sit1, cot1, rkx, rky, rkz, pk, coga, cp, cc
    complex (kind=8) :: ct1, ct2
    real, parameter :: twopi = 6.28319
  
    Nx = grid
    Ny = grid
    Nz = grid

    kFx = twopi/box   ! fundamental frequencies
    kFy = twopi/box
    kFz = twopi/box
    kF3 = kFx*kFy*kFz
  
    kNx=kFx*dble(Nx)/2.d0    ! Nyquist frequencies
    kNy=kFy*dble(Ny)/2.d0 
    kNz=kFz*dble(Nz)/2.d0
  
    kFmin=min(kFx,min(kFy,kFz))
    kNmax=max(kNx,max(kNy,kNz))

    ! bin size is determined based on the smallest fundamental frequency
    ! while the number of bins depends on the largest Nyquist frequency

    dk=kbin*kFmin   ! size of the wavenumbers bin
  
    avgk=0.d0
    avgP=0.d0
    co=0.d0

    do iz=1,Nz
        ikz=mod(iz+Nz/2-2,Nz)-Nz/2+1
        icz=mod(Nz+1-iz,Nz)+1
        rkz=kFz*float(ikz)

         do iy=1,Ny
            iky=mod(iy+Ny/2-2,Ny)-Ny/2+1
            icy=mod(Ny+1-iy,Ny)+1
            rky=kFy*float(iky)

            do ix=1,Nx/2+1
                ikx=mod(ix+Nx/2-2,Nx)-Nx/2+1
                icx=mod(Nx+1-ix,Nx)+1
                rkx=kFx*float(ikx)

	            if (ix.ne.1 .or. (ix.eq.1 .and. iky.gt.0) .or. (ix.eq.1 .and. iky.eq.0 .and. ikz.gt.0)) then

                    rk=dsqrt(rkx**2+rky**2+rkz**2)
                    imk=nint(rk/dk-kcenter/kbin-tiny(kFmin))+1
                    !imk=nint(rk/dk-kcenter/kbin-0.5-tiny(kFmin))+1  !if you want [kf,2kf) instead of [0.5kf,1,5kf)

	               if (imk.le.nbins .and. imk.gt.0) then

                        co(imk)=co(imk)+1.d0
                        avgk(imk)=avgk(imk)+rk
		                if (ix.le.Nx/2+1) then
                            ct1=dcl1(ix,iy,iz)
                            ct2=dcl2(ix,iy,iz)
		                else
                            ct1=dcl2(icx,icy,icz)
                            ct2=dcl1(icx,icy,icz)
		                endif
                        pk=conjg(ct1)*(ct2)
                        pk=real(pk)
                        avgP(imk)=avgP(imk)+pk
	               endif
	            endif
            enddo
         enddo
    enddo
  
   pout = 0.d0
   pout(1,:) = avgk(:)/co(:)
   pout(2,:) = avgP(:)*kF3/co(:)
   pout(3,:) = co(:)

end subroutine

subroutine measure_pk_multipoles(grid, box, kbin, kcenter, nbins, dcl, dcl2, dcl4, pout)

   integer, intent(in) :: grid, nbins
   real, intent(in) :: box, kbin, kcenter
   complex(kind=8), intent(in) :: dcl(:,:,:), dcl2(:,:,:), dcl4(:,:,:)
   real(kind=8), intent(out) :: pout(5,nbins)
   integer :: ix, iy, iz, icx,icy,icz,ikx,iky,ikz,imk, Nx, Ny, Nz
   real (kind=8), dimension(nbins) :: avgk, avgP, avgP2, avgP4, co
   real (kind=8) :: k,kFmin,rk,pSN,kNmax,dk, rkx, rky, rkz, pk
   real (kind=8) :: kFx, kFy, kFz, kNx, kNy, kNz, kF3
   complex (kind=8) :: ct, ct2, ct4
   real, parameter :: twopi = 6.28319

   Nx = grid
   Ny = grid
   Nz = grid

   kFx = twopi/box   ! fundamental frequencies
   kFy = twopi/box
   kFz = twopi/box
   kF3 = kFx*kFy*kFz

   kNx=kFx*dble(Nx)/2.d0    ! Nyquist frequencies
   kNy=kFy*dble(Ny)/2.d0 
   kNz=kFz*dble(Nz)/2.d0

   kFmin=min(kFx,min(kFy,kFz))
   kNmax=max(kNx,max(kNy,kNz))

   ! bin size is determined based on the smallest fundamental frequency
   ! while the number of bins depends on the largest Nyquist frequency

   dk=kbin*kFmin   ! size of the wavenumbers bin

   avgk=0.d0
   avgP=0.d0
   co=0.d0

   avgP2=0.d0
   avgP4=0.d0

   do iz=1,Nz
      ikz=mod(iz+Nz/2-2,Nz)-Nz/2+1
      icz=mod(Nz+1-iz,Nz)+1
      rkz=kFz*float(ikz)

      do iy=1,Ny
         iky=mod(iy+Ny/2-2,Ny)-Ny/2+1
         icy=mod(Ny+1-iy,Ny)+1
         rky=kFy*float(iky)

         do ix=1,Nx/2+1
            ikx=mod(ix+Nx/2-2,Nx)-Nx/2+1
            icx=mod(Nx+1-ix,Nx)+1
            rkx=kFx*float(ikx)

            if (ix.ne.1 .or. (ix.eq.1 .and. iky.gt.0) .or. (ix.eq.1 .and. iky.eq.0 .and. ikz.gt.0)) then

               rk=dsqrt(rkx**2+rky**2+rkz**2)
               imk=nint(rk/dk-kcenter/kbin-tiny(kFmin))+1

               if (imk.le.nbins .and. imk.gt.0) then

                  co(imk)=co(imk)+1.d0
                  avgk(imk)=avgk(imk)+rk
                  !if (ix.le.Nx/2+1) then
                  ct=dcl(ix,iy,iz)
                  ct2=dcl2(ix,iy,iz)
                  ct4=dcl4(ix,iy,iz)
                  !else
                  !   ct=dcl(icx,icy,icz)
                  !   ct2=dcl2(icx,icy,icz)
                  !endif
                  pk=(cdabs(ct))**2
                  avgP(imk)=avgP(imk)+pk
            
                  !pk=cdabs(conjg(ct)*ct2)
                  pk=dble(conjg(ct)*ct2)
                  avgP2(imk) = avgP2(imk) + pk
                  pk=dble(conjg(ct)*ct4)
                  !pk=(cdabs(ct2))**2
                  avgP4(imk) = avgP4(imk) + pk
               endif
            endif
         enddo
      enddo
   enddo

  
   pout = 0.d0
   pout(1,:) = avgk(:)/co(:)
   pout(2,:) = avgP(:)*kF3/co(:)
   pout(3,:) = 5.*avgP2(:)*kF3/co(:)
   pout(4,:) = 9.*avgP4(:)*kF3/co(:)
   !pout(4,:) = 17.5*avgP4(:)*kF3/co(:)-pout(3,:)-3.5*pout(2,:)    this is if we compute only delta2
   pout(5,:) = co(:)

end subroutine
