c
c     energy and forces
c
c
      subroutine edisp1
      use atoms
      use limits
      use deriv
      use ewald
      use energi
      use disp
      implicit none
      integer i
      real*8 pre
c
c     choose method of summing over pairwise interactions
c
      edis = 0.0d0
      do i = 1, n
         dedis(1,i) = 0.0d0
         dedis(2,i) = 0.0d0
         dedis(3,i) = 0.0d0
      end do
      if (use_dewald) then
         do i = 1, ndisp
            pre = adewald**6/12.0d0
            edis = edis + pre*csix(i)*csix(i)
         end do
         if (use_mlist) then
c            call edisp1d
         else
            call edisprecip1
c            call edisp1c
         end if
      else
         if (use_mlist) then
c            call edisp1b
         else
c            call edisp1a
         end if
      end if
      return
      end
c
c
c
c
c
c
c
c
      subroutine edisp1c
      use sizes
      use atoms
      use bound
      use boxes
      use couple
      use cell
      use disp
      use disppot
      use deriv
      use energi
      use ewald
      use group
      use shunt
      use usage
      use virial
      implicit none
      integer i,j,k
      integer in
      integer ii,kk,kkk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 e,de
      real*8 dedx,dedy,dedz
      real*8 ci,ck
      real*8 r2,r6,r8
      real*8 fgrp
      real*8 ralpha2,term,dterm,expterm
      real*8 damp,ddamp
      real*8 vxx,vyx,vzx
      real*8 vxy,vyy,vzy
      real*8 vxz,vyz,vzz
      real*8, allocatable :: disscale(:)
      logical proceed,usei,usek
      character*6 mode

c
c     zero out the dispersion interaction
c
      if (ndisp .eq. 0) return
c
c     perform dynamic allocation of some local arrays
c
      allocate (disscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         disscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      mode = 'DEWALD'
      call switch (mode)
c
c     compute the Ewald self-energy term
c




c
c     compute the reciprocal space part of the Ewald summation
c
      call edisprecip1
c
c     compute the real space portion of the Ewald summation
c
      do ii = 1, ndisp-1
         i = idisp(ii)
         ci = csix(ii)
         usei = use(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            disscale(i12(j,i)) = dis2scale
         end do
         do j = 1, n13(i)
            disscale(i13(j,i)) = dis3scale
         end do
         do j = 1, n14(i)
            disscale(i14(j,i)) = dis4scale
         end do
         do j = 1, n15(i)
            disscale(i15(j,i)) = dis5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii+1, ndisp
            k = idisp(kk)
            ck = csix(kk)
            usek = use(k)
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               call image (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  r6 = r2**3
                  r8 = r6*r2
                  ralpha2 = r2 * adewald**2
                  damp = 1.0d0
                  if (ralpha2 .lt. 50.0d0) then
                     expterm = exp(-ralpha2)
                     term = 1.0d0 + ralpha2 + 0.5d0*ralpha2**2
                     damp = term*expterm
                  end if
                  e = -ci*ck*(damp + (disscale(k) - 1.0d0))/r6
                  edis = edis + e
                  dterm = term + (ralpha2**3)/6.0d0
                  ddamp = dterm*expterm
                  de = ci*ck*6.0d0*(ddamp + (disscale(k) - 1.0d0))/r8
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
                  dedis(1,i) = dedis(1,i) + dedx
                  dedis(2,i) = dedis(2,i) + dedy
                  dedis(3,i) = dedis(3,i) + dedz
                  dedis(1,k) = dedis(1,k) - dedx
                  dedis(2,k) = dedis(2,k) - dedy
                  dedis(3,k) = dedis(3,k) - dedz
c
c     increment the internal virial tensor components
c
                     vxx = xr * dedx
                     vyx = yr * dedx
                     vzx = zr * dedx
                     vxy = xr * dedy
                     vyy = yr * dedy
                     vzy = zr * dedy
                     vxz = xr * dedz
                     vxy = yr * dedz
                     vzz = zr * dedz
                     vir(1,1) = vir(1,1) + vxx
                     vir(2,1) = vir(2,1) + vyx
                     vir(3,1) = vir(3,1) + vzx
                     vir(1,2) = vir(1,2) + vxy
                     vir(2,2) = vir(2,2) + vyy
                     vir(3,2) = vir(3,2) + vzy
                     vir(1,3) = vir(1,3) + vxz
                     vir(2,3) = vir(2,3) + vyz
                     vir(3,3) = vir(3,3) + vzz
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            disscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            disscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            disscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            disscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c
c
      print *,"you're using replica's dummy -barf"
      stop
      end
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine edisprecip1
      use sizes
      use boxes
      use bound
      use disp
      use deriv
      use energi
      use ewald
      use math
      use pme
      use virial
      implicit none
      integer i,j,k
      integer k1,k2,k3
      integer m1,m2,m3
      integer nf1,nf2,nf3
      integer nff,npoint
      real*8 e,f,denom
      real*8 term,expterm
      real*8 hsq,struc2
      real*8 h1,h2,h3
      real*8 r1,r2,r3
      real*8 h,hhh,b
      real*8 term1,vterm
      real*8 fac1,fac2,fac3,bfac
      real*8 denom0
      real*8 erfcterm
      real*8 de1,de2,de3
      real*8 dn1,dn2,dn3
      real*8 dt1,dt2,dt3
      real*8 fi
      integer i0,iatm,igrd0
      integer isite
      integer it1,it2,it3
      integer j0,jgrd0,k0,kgrd0
      real*8 t1,t2,t3
      real*8 pre
c
c
c     zero out the particle mesh Ewald grid values
c
      do k = 1, ndfft3
         do j = 1, ndfft2
            do i = 1, ndfft1
               qgrid(1,i,j,k) = 0.0d0
               qgrid(2,i,j,k) = 0.0d0
            end do
         end do
      end do
c
c     get B-spline coefficients and put charges onto grid
c
      call bspline_fill
      call table_fill
      call grid_csix
c
c     perform the 3-D FFT forward transformation
c
      call fftfront
c
c     use scalar sum to get the reciprocal space energy
c
      npoint = ndfft1 * ndfft2 * ndfft3
      nff = ndfft1 * ndfft2
      nf1 = (ndfft1+1) / 2
      nf2 = (ndfft2+1) / 2
      nf3 = (ndfft3+1) / 2
      bfac = pi / adewald
      fac1 = 2.0d0*pi**(3.5d0)
      fac2 = adewald**3
      fac3 = -2.0d0*adewald*pi**2
      denom0 = (6.0d0*volbox)/(pi**1.5d0)
      do i = 1, npoint-1
         k3 = i/nff + 1
         j = i - (k3-1)*nff
         k2 = j/ndfft1 + 1
         k1 = j - (k2-1)*ndfft1 + 1
         m1 = k1 - 1
         m2 = k2 - 1
         m3 = k3 - 1
         if (k1 .gt. nf1)  m1 = m1 - ndfft1
         if (k2 .gt. nf2)  m2 = m2 - ndfft2
         if (k3 .gt. nf3)  m3 = m3 - ndfft3
         r1 = dble(m1)
         r2 = dble(m2)
         r3 = dble(m3)
         h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
         h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
         h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
         hsq = h1*h1 + h2*h2 + h3*h3
         h = sqrt(hsq)
         b = h*bfac
         hhh = h*hsq
         term = -b*b
         expterm = 0.0d0
         erfcterm = erfc(b)
         denom = denom0*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
         if (term .gt. -50.0d0) then
            expterm = exp(term)
            erfcterm = erfc(b)
            if (.not. use_bounds) then
               expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
               erfcterm = erfcterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
            else if (octahedron) then
               if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
               if (mod(m1+m2+m3,2) .ne. 0)  erfcterm = 0.0d0
            end if
            term1 = fac1*erfcterm*hhh + expterm*(fac2 + fac3*hsq)
            struc2 = qgrid(1,k1,k2,k3)**2 + qgrid(2,k1,k2,k3)**2
            e = -(term1 / denom) * struc2
            edis = edis + e
            vterm = 3.0d0*(fac1*erfcterm*h + fac3*expterm)*struc2/denom
            vir(1,1) = vir(1,1) + h1*h1*vterm - e
            vir(2,1) = vir(2,1) + h1*h2*vterm
            vir(3,1) = vir(3,1) + h1*h3*vterm
            vir(1,2) = vir(1,2) + h2*h1*vterm
            vir(2,2) = vir(2,2) + h2*h2*vterm - e
            vir(3,2) = vir(3,2) + h2*h3*vterm
            vir(1,3) = vir(1,3) + h3*h1*vterm
            vir(2,3) = vir(2,3) + h3*h2*vterm
            vir(3,3) = vir(3,3) + h3*h3*vterm - e
         end if
         qgrid(1,k1,k2,k3) = -(term1/denom) * qgrid(1,k1,k2,k3) 
         qgrid(2,k1,k2,k3) = -(term1/denom) * qgrid(2,k1,k2,k3)
      end do
c
c     account for zero point term
c        
      pre = (adewald**3)/denom0
      do i = 1, ndisp
         do j = 1, ndisp
            edis = edis - csix(i)*csix(j)*pre
            vir(1,1) = vir(1,1) + csix(i)*csix(j)*pre
            vir(2,2) = vir(2,2) + csix(i)*csix(j)*pre
            vir(3,3) = vir(3,3) + csix(i)*csix(j)*pre
         end do
      end do
c
c     perform the 3-D FFT backward transformation
c
      qgrid(1,1,1,1) = 0.0d0
      qgrid(2,1,1,1) = 0.0d0
c
      call fftback
c
c     get first derivatives of the reciprocal space energy 
c
      dn1 = dble(ndfft1)
      dn2 = dble(ndfft2)
      dn3 = dble(ndfft3)
      do isite = 1, ndisp
         iatm = idisp(isite)
         igrd0 = igrid(1,iatm)
         jgrd0 = igrid(2,iatm)
         kgrd0 = igrid(3,iatm)
         fi = csix(isite)
         de1 = 0.0d0
         de2 = 0.0d0
         de3 = 0.0d0
         k0 = kgrd0
         do it3 = 1, bsorder
            k0 = k0 + 1
            k = k0 + 1 + (ndfft3-sign(ndfft3,k0))/2
            t3 = thetai3(1,it3,iatm)
            dt3 = dn3 * thetai3(2,it3,iatm)
            j0 = jgrd0
            do it2 = 1, bsorder
               j0 = j0 + 1
               j = j0 + 1 + (ndfft2-sign(ndfft2,j0))/2
               t2 = thetai2(1,it2,iatm)
               dt2 = dn2 * thetai2(2,it2,iatm)
               i0 = igrd0
               do it1 = 1, bsorder
                  i0 = i0 + 1
                  i = i0 + 1 + (ndfft1-sign(ndfft1,i0))/2
                  t1 = thetai1(1,it1,iatm)
                  dt1 = dn1 * thetai1(2,it1,iatm)
                  term = qgrid(1,i,j,k)
                  de1 = de1 + 2.0d0*term*dt1*t2*t3
                  de2 = de2 + 2.0d0*term*dt2*t1*t3
                  de3 = de3 + 2.0d0*term*dt3*t1*t2
               end do
            end do
         end do
         dedis(1,iatm) = dedis(1,iatm) + fi*(recip(1,1)*de1+recip(1,2)*
     &                                      de2+recip(1,3)*de3)
         dedis(2,iatm) = dedis(2,iatm) + fi*(recip(2,1)*de1+recip(2,2)*
     &                                      de2+recip(2,3)*de3)
         dedis(3,iatm) = dedis(3,iatm) + fi*(recip(3,1)*de1+recip(3,2)*
     &                                      de2+recip(3,3)*de3)
      end do
      return
      end

