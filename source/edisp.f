c
c
c
c
      subroutine edisp
      use limits
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
      if (use_dewald) then
         do i = 1, ndisp
            pre = adewald**6/12.0d0
            edis = edis + pre*csix(i)*csix(i)
         end do
         print *,"self",edis
         if (use_mlist) then
c            call edisp0d
         else
            call edisp0c
         end if
      else
         if (use_mlist) then
c            call edisp0b
         else
c            call edisp0a
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
      subroutine edisp0c
      use sizes
      use atoms
      use bound
      use boxes
      use couple
      use cell
      use disp
      use disppot
      use energi
      use ewald
      use group
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer in
      integer ii,kk,kkk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 e
      real*8 ci,ck
      real*8 r2,r6
      real*8 fgrp
      real*8 ralpha2,term,expterm
      real*8 damp
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
c      call edisprecip


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
                  ralpha2 = r2 * adewald**2
                  damp = 1.0d0
                  if (ralpha2 .lt. 50.0d0) then
                     expterm = exp(-ralpha2)
                     print *,"expterm",expterm
                     term = 1.0d0 + ralpha2 + 0.5d0*ralpha2**2
                     damp = term*expterm
                  end if
                  e = -ci*ck*(damp + (disscale(k) - 1.0d0))/r6
                  edis = edis + e
                  print *,"e",i,k,damp,e,edis
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
      subroutine edisprecip
      use sizes
      use bound
      use boxes
      use energi
      use ewald
      use math
      use pme
      implicit none
      integer i,j,k
      integer k1,k2,k3
      integer m1,m2,m3
      integer nf1,nf2,nf3
      integer nff,npoint
      real*8 e,f,denom
      real*8 term,expterm
      real*8 pterm,volterm
      real*8 hsq,struc2
      real*8 h1,h2,h3
      real*8 r1,r2,r3
c
c
c     zero out the particle mesh Ewald grid values
c
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
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
      call grid_pchg
c
c     perform the 3-D FFT forward transformation
c
      call fftfront
c
c     use scalar sum to get the reciprocal space energy
c
      f = 0.5d0 * electric / dielec
      npoint = nfft1 * nfft2 * nfft3
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nff = nfft1 * nfft2
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
      do i = 1, npoint-1
         k3 = i/nff + 1
         j = i - (k3-1)*nff
         k2 = j/nfft1 + 1
         k1 = j - (k2-1)*nfft1 + 1
         m1 = k1 - 1
         m2 = k2 - 1
         m3 = k3 - 1
         if (k1 .gt. nf1)  m1 = m1 - nfft1
         if (k2 .gt. nf2)  m2 = m2 - nfft2
         if (k3 .gt. nf3)  m3 = m3 - nfft3
         r1 = dble(m1)
         r2 = dble(m2)
         r3 = dble(m3)
         h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
         h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
         h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
         hsq = h1*h1 + h2*h2 + h3*h3
         term = -pterm * hsq
         expterm = 0.0d0
         if (term .gt. -50.0d0) then
            denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
            expterm = exp(term) / denom
            if (.not. use_bounds) then
               expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
            else if (octahedron) then
               if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
            end if
            struc2 = qgrid(1,k1,k2,k3)**2 + qgrid(2,k1,k2,k3)**2
            e = f * expterm * struc2
            ec = ec + e
         end if
      end do
c
c     account for zeroth grid point for nonperiodic system
c
      if (.not. use_bounds) then
         expterm = 0.5d0 * pi / xbox
         struc2 = qgrid(1,1,1,1)**2 + qgrid(2,1,1,1)**2
         e = f * expterm * struc2
         ec = ec + e
      end if
      return
      end
