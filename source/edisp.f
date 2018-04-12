c
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
c
c     compute self energy term
c
            edis = edis + pre*csix(i)*csix(i)
         end do
         if (use_mlist) then
c            call edisprecip
            call edisp0d
         else
            call edisp0c
c            call edisprecip
         end if
      else
         if (use_mlist) then
            call edisp0b
         else
            call edisp0a
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
      subroutine edisp0a
      use sizes
      use atoms
      use bound
      use boxes
      use couple
      use cell
      use disp
      use disppot
      use energi
c      use ewald
      use group
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer in
      integer ii,kk,kkk
      integer counter
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 e
      real*8 ci,ck
      real*8 r,r2,r6
      real*8 fgrp
      real*8 ralpha2,term,expterm
      real*8 damp
      real*8 alphai,alphak
      real*8 dampi,dampk
      real*8 expdampi,expdampk
      real*8 termi,termk
      real*8 damp3,damp5,damping
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
c
c     I AM USING THE MPOLE CUTOFF HERE TO AVOID CREATING A NEW CUTOFF
c
      mode = 'MPOLE'
      call switch (mode)
c
c     compute the real space portion of the Ewald summation
c
      counter = 0
      do ii = 1, ndisp-1
         i = idisp(ii)
c
c     read in dispersion parameters for atom i
c
         ci = csix(ii)
         alphai = adisp(ii)
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
            alphak = adisp(kk)
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
c               print *,"offset",r2,off2
               if (r2 .le. off2) then
                  counter = counter + 1
                  r6 = r2**3
                  r = sqrt(r2)
c                  ralpha2 = r2 * adewald**2
c                  damp = 1.0d0
c                  if (ralpha2 .lt. 50.0d0) then
c                     expterm = exp(-ralpha2)
c                     term = 1.0d0 + ralpha2 + 0.5d0*ralpha2**2
c                     damp = term*expterm
c                  end if
c                  print *,"disscale",disscale(k)
c                  e = -ci*ck*(damp + (disscale(k) - 1.0d0))/r6
                  e = -ci*ck*disscale(k)/r6
c
c     apply dispersion damping
c
                  if (dispdamp) then
c
c     gordon damping
c 
                     dampi = alphai*r
                     dampk = alphak*r
                     expdampi = exp(-dampi)
                     expdampk = exp(-dampk)
                     if (alphai.ne.alphak) then
                        termi = alphak**2/(alphak**2 - alphai**2)
                        termk = alphai**2/(alphai**2 - alphak**2)
                        damp3 = 1.0d0 - termi*(1.0d0 +dampi)*expdampi
     &                       - termk*(1.0d0 + dampk)*expdampk
                        damp5 = 1.0d0 - termi*(1.0d0 + dampi +
     &                       (1.0d0/3.0d0)*dampi**2)*expdampi -
     &                       termk*(1.0d0 + dampk +
     &                       (1.0d0/3.0d0)*dampk**2)*expdampk
                        damping = (3.0d0*damp5 - damp3)/2.0d0
                     else
                        damp3 = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2)
     &                       *expdampi
                        damp5 = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2
     &                       + (1.0d0/6.0d0)*dampi**3)*expdampi
                        damping = (3.0d0*damp5 - damp3)/2.0d0
                     end if
c                     e = -ci*ck*(damp +
c     &                    (disscale(k)*damping**2 - 1.0d0))/r6
                     e = -ci*ck*disscale(k)*(damping**2)/r6
c
c     insert gordon1 damping later
c
                  end if
c
c     accumulate dispersion energy
c
                  edis = edis + e
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
      print *,"brute force dispersion energy",edis,counter
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c
c
      print *,"you're using replicas dummy -barf"
      stop
      end
c
c
c
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
      subroutine edisp0b
      use sizes
      use atoms
      use bound
      use boxes
      use couple
      use cell
      use disp
      use disppot
      use energi
      use group
      use neigh
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer in
      integer ii,kk,kkk
      integer counter
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 e
      real*8 ci,ck
      real*8 r,r2,r6
      real*8 fgrp
      real*8 ralpha2,term,expterm
      real*8 damp
      real*8 alphai,alphak
      real*8 dampi,dampk
      real*8 expdampi,expdampk
      real*8 termi,termk
      real*8 damp3,damp5,damping
      real*8, allocatable :: disscale(:)
      logical proceed,usei,usek
      character*6 mode
c
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
      mode = 'MPOLE'
      call switch (mode)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(ndisp,idisp,counter,
!$OMP& csix,adisp,use,x,y,z,n12,n13,n14,n15,dispdamp,
!$OMP& i12,i13,i14,i15,nelst,elst,
!$OMP& dis2scale,dis3scale,dis4scale,dis5scale,
!$OMP& use_group,off2)
!$OMP& firstprivate(disscale)
!$OMP& shared(edis)
!$OMP DO reduction(+:edis) schedule(guided)
c
c     compute the real space portion of the Ewald summation
c
      do ii = 1, ndisp
         i = idisp(ii)
c
c     read in dispersion parameters for atom i
c
         ci = csix(ii)
         alphai = adisp(ii)
c
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
c
c     for now, use mpole neighborlist
c
         do kkk = 1, nelst(ii)
            kk = elst(kkk,ii)
            k = idisp(kk)
c
c     read in dispersion parameters for atom k
c
            ck = csix(kk)
            alphak = adisp(kk)
c
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
c                  r = sqrt(r2)
                  e = -ci*ck*disscale(k)/r6
c
c     apply dispersion damping
c
                  if (dispdamp) then
c
c     gordon damping
c
                     if (r2 .le. 16.0d0) then
                     r = sqrt(r2)
                     dampi = alphai*r
                     dampk = alphak*r
                     expdampi = exp(-dampi)
                     expdampk = exp(-dampk)
                     if (alphai.ne.alphak) then
                        termi = alphak**2/(alphak**2 - alphai**2)
                        termk = alphai**2/(alphai**2 - alphak**2)
                        damp3 = 1.0d0 - termi*(1.0d0 +dampi)*expdampi
     &                       - termk*(1.0d0 + dampk)*expdampk
                        damp5 = 1.0d0 - termi*(1.0d0 + dampi +
     &                       (1.0d0/3.0d0)*dampi**2)*expdampi -
     &                       termk*(1.0d0 + dampk +
     &                       (1.0d0/3.0d0)*dampk**2)*expdampk
                        damping = (3.0d0*damp5 - damp3)/2.0d0
                     else
                        damp3 = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2)
     &                       *expdampi
                        damp5 = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2
     &                       + (1.0d0/6.0d0)*dampi**3)*expdampi
                        damping = (3.0d0*damp5 - damp3)/2.0d0
                     end if
                     end if
                     e = -ci*ck*disscale(k)*(damping**2)/r6
c
c     insert gordon1 damping later
c
                  end if
c
c     accumulate dispersion energy
c
                  edis = edis + e
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
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (disscale)
c      print *,"neighborlist, no ewald dispersion energy",edis,counter
      return
      end
c
c
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
      real*8 r,r2,r6
      real*8 fgrp
      real*8 ralpha2,term,expterm
      real*8 damp
      real*8 alphai,alphak
      real*8 dampi,dampk
      real*8 expdampi,expdampk
      real*8 termi,termk
      real*8 damp3,damp5,damping
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
      mode = 'MPOLE'
      call switch (mode)
c      print *,"cutoff",off
c
c     compute the Ewald self-energy term
c


c     I BELIEVE THERE IS NO SELF-ENERGY TERM?

c
c     compute the recipocal space part of the Ewald summation
c
      call edisprecip
c
c     compute the real space portion of the Ewald summation
c
      do ii = 1, ndisp-1
         i = idisp(ii)
c
c     read in dispersion parameters for atom i
c
         ci = csix(ii)
         alphai = adisp(ii)
c
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
c
c     read in dispersion parameters for atom k
c
            ck = csix(kk)
            alphak = adisp(kk)
c
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
                  r = sqrt(r2)
                  ralpha2 = r2 * adewald**2
                  damp = 1.0d0
                  if (ralpha2 .lt. 50.0d0) then
                     expterm = exp(-ralpha2)
                     term = 1.0d0 + ralpha2 + 0.5d0*ralpha2**2
                     damp = term*expterm
                  end if
                  e = -ci*ck*(damp + (disscale(k) - 1.0d0))/r6
c
c     apply dispersion damping
c
                  if (dispdamp) then
c
c     gordon damping
c
                     dampi = alphai*r
                     dampk = alphak*r
                     expdampi = exp(-dampi)
                     expdampk = exp(-dampk)
                     if (alphai.ne.alphak) then
                        termi = alphak**2/(alphak**2 - alphai**2)
                        termk = alphai**2/(alphai**2 - alphak**2)
                        damp3 = 1.0d0 - termi*(1.0d0 +dampi)*expdampi
     &                       - termk*(1.0d0 + dampk)*expdampk
                        damp5 = 1.0d0 - termi*(1.0d0 + dampi +
     &                       (1.0d0/3.0d0)*dampi**2)*expdampi -
     &                       termk*(1.0d0 + dampk +
     &                       (1.0d0/3.0d0)*dampk**2)*expdampk
                        damping = (3.0d0*damp5 - damp3)/2.0d0
                     else
                        damp3 = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2)
     &                       *expdampi
                        damp5 = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2
     &                       + (1.0d0/6.0d0)*dampi**3)*expdampi
                        damping = (3.0d0*damp5 - damp3)/2.0d0
                     end if
c                     if ((r .gt. 4.0d0).and.(r .lt. 4.1d0)) 
c     &                    print *,"i,k",i,k,r,damping**2
                     e = -ci*ck*(damp +
     &                    (disscale(k)*damping**2 - 1.0d0))/r6
c
c     insert gordon1 damping later
c
                  end if
c
c     accumulate dispersion energy
c
                  edis = edis + e
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
c      print *,"non-list ewald dispersion energy",edis
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c
c
      print *,"you're using replicas dummy -barf"
      stop
      end
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
      subroutine edisp0d
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
      use neigh
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer in
      integer ii,kk,kkk
      integer counter
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 e
      real*8 ci,ck
      real*8 r,r2,r6
      real*8 fgrp
      real*8 ralpha2,term,expterm
      real*8 damp
      real*8 alphai,alphak
      real*8 dampi,dampk
      real*8 expdampi,expdampk
      real*8 termi,termk
      real*8 damp3,damp5,damping
      real*8, allocatable :: disscale(:)
      logical proceed,usei,usek
      character*6 mode
c
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
      mode = 'MPOLE'
      call switch (mode)
c
c     compute the Ewald self-energy term
c


c     I BELIEVE THERE IS NO SELF-ENERGY TERM? NOPE

c
c     compute the recipocal space part of the Ewald summation
c
      call edisprecip
c      counter = 0
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(ndisp,idisp,
!$OMP& csix,adisp,use,x,y,z,n12,n13,n14,n15,dispdamp,
!$OMP& i12,i13,i14,i15,nelst,elst,
!$OMP& dis2scale,dis3scale,dis4scale,dis5scale,
!$OMP& use_group,off2,adewald)
!$OMP& firstprivate(disscale)
!$OMP& shared(edis)
!$OMP DO reduction(+:edis) schedule(guided)

c
c     compute the real space portion of the Ewald summation
c
      do ii = 1, ndisp
         i = idisp(ii)
c
c     read in dispersion parameters for atom i
c
         ci = csix(ii)
         alphai = adisp(ii)
c
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
c
c     for now, use vdw neighborlist
c
         do kkk = 1, nelst(ii)
c            print *,"nelst",nelst(ii)
c            print *,"elst",elst(kkk,ii)
            kk = elst(kkk,ii)
            k = idisp(kk)
c
c     read in dispersion parameters for atom k
c
            ck = csix(kk)
            alphak = adisp(kk)
c
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
c                  counter = counter + 1
                  r6 = r2**3
c                  r = sqrt(r2)
                  ralpha2 = r2 * adewald**2
                  damp = 1.0d0
                  if (ralpha2 .lt. 50.0d0) then
                     expterm = exp(-ralpha2)
                     term = 1.0d0 + ralpha2 + 0.5d0*ralpha2**2
                     damp = term*expterm
                  end if
                  e = -ci*ck*(damp + (disscale(k) - 1.0d0))/r6
c
c     apply dispersion damping
c
                  if (dispdamp) then
c
c     gordon damping
c
                     damping = 1.0d0
                     if (r2 .lt. 16.0d0) then
                        r = sqrt(r2)
                        dampi = alphai*r
                        dampk = alphak*r
                        expdampi = exp(-dampi)
                        expdampk = exp(-dampk)
                        if (alphai.ne.alphak) then
                           termi = alphak**2/(alphak**2 - alphai**2)
                           termk = alphai**2/(alphai**2 - alphak**2)
                           damp3 = 1.0d0 - termi*(1.0d0 +dampi)*expdampi
     &                          - termk*(1.0d0 + dampk)*expdampk
                           damp5 = 1.0d0 - termi*(1.0d0 + dampi +
     &                          (1.0d0/3.0d0)*dampi**2)*expdampi -
     &                          termk*(1.0d0 + dampk +
     &                          (1.0d0/3.0d0)*dampk**2)*expdampk
                           damping = (3.0d0*damp5 - damp3)/2.0d0
                        else
                           damp3 = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2)
     &                          *expdampi
                           damp5 = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2
     &                          + (1.0d0/6.0d0)*dampi**3)*expdampi
                           damping = (3.0d0*damp5 - damp3)/2.0d0
                        end if
                     end if
                     e = -ci*ck*(damp +
     &                    (disscale(k)*damping**2 - 1.0d0))/r6
c
c     insert gordon1 damping later
c
                  end if
c
c     accumulate dispersion energy
c
                  edis = edis + e
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
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (disscale)
c      print *,"neighborlist AND ewald dispersion energy",edis
      return
      end

















c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine edisprecip
      use sizes
      use boxes
      use bound
      use disp
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
      real*8 hsq,struc2
      real*8 h1,h2,h3
      real*8 r1,r2,r3
      real*8 h,hhh,b
      real*8 term1
      real*8 fac1,fac2,fac3,bfac
      real*8 denom0
      real*8 erfcterm
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
         if (term .gt. -50.0d0) then
            denom = denom0*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
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
         end if
      end do
c
c     account for zero point term
c        
      do i = 1, ndisp
         do j = 1, ndisp
            edis = edis - csix(i)*csix(j)*(adewald**3)/denom0
         end do
      end do
c
      return
      end
