c
c
c     #############################################################
c     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine emoeba3  --  atomic multipole energy & analysis  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "emoeba3" calculates the electrostatic energy due to atomic
c     multipole interactions, and partitions the energy among atoms
c
c
      subroutine emoeba3
      use limits
      use mpole
      implicit none
c
c     save permanent electric field
c     
      savefield = .true.
c
c
c     choose the method for summing over multipole interactions
c
      if (use_ewald) then
         if (use_mlist) then
            call emoeba3d
         else
            call emoeba3c
         end if
      else
         if (use_mlist) then
            call emoeba3b
         else
            call emoeba3a
         end if
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine emoeba3a  --  double loop multipole analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "emoeba3a" calculates the atomic multipole interaction energy
c     using a double loop, and partitions the energy among atoms
c
c
      subroutine emoeba3a
      use sizes
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use cell
      use chgpot
      use couple
      use energi
      use group
      use inform
      use inter
      use iounit
      use math
      use molcul
      use mplpot
      use mpole
      use potent
      use shunt
      use usage
      use chgpen
      use disp
      use pauli
      implicit none
      integer i,j,k
      integer ii,kk
      integer ix,iy,iz
      integer kx,ky,kz
      integer rorder
      real*8 e,f,fgrp
      real*8 ecc,ecv,evc,evv
      real*8 e_ele,e_disp,e_pauli
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 alphai,alphak
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 corei,vali
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 corek,valk
      real*8 qrix,qriy,qriz
      real*8 qrkx,qrky,qrkz
      real*8 dri,drk,dik
      real*8 qrri,qrrk
      real*8 qrrik,qik
      real*8 diqrk,dkqri
      real*8 term1,term2,term3
      real*8 term4,term5
      real*8 term1ik,term2ik,term3ik
      real*8 term4ik,term5ik
      real*8 term1i,term2i,term3i
      real*8 term1k,term2k,term3k
      real*8 rr6
      real*8 c6i,c6k,c6ik
      real*8 displam
      real*8 pvali,pvalk
      real*8 overlapi,overlapk,oik
      real*8 apauli,apaulk
      real*8 fid(3),fkd(3)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: lambdai(:)
      real*8, allocatable :: lambdak(:)
      real*8, allocatable :: lambdaik(:)
      logical proceed
      logical header,huge
      logical usei,usek
      character*6 mode
c
c
c     zero out total atomic multipole energy and partitioning
c
      nem = 0
      nedis = 0
      nepr = 0
      em = 0.0d0
      edis = 0.0d0
      epr = 0.0d0
      do i = 1, n
         aem(i) = 0.0d0
         aedis(i) = 0.0d0
         aepr(i) = 0.0d0
         do j = 1, 3
            permfield(j,i) = 0.0d0
         end do
      end do
      if (npole .eq. 0)  return
c
c     set maximum power of 1/r for damping
c     (9 for quadrupole-quadrupole energy)
c
      rorder = 9
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. npole.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Atomic Multipole Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (lambdai(rorder))
      allocate (lambdak(rorder))
      allocate (lambdaik(rorder))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         mscale(i) = 1.0d0
      end do
c
c     initialize charge penetration scale factors
c
      do i = 1, rorder
         lambdai(i) = 1.0d0
         lambdak(i) = 1.0d0
         lambdaik(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     calculate the multipole interaction energy term
c
      do i = 1, npole-1
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         corei = monopole(1,i)
         vali = monopole(2,i)
         alphai = alphaele(i)
         c6i = csix(i)
         overlapi = overpauli(i)
         apauli = alphapauli(i)
         pvali = monopauli(i)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         usei = (use(ii) .or. use(iz) .or. use(ix) .or. use(iy))
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do k = i+1, npole
            kk = ipole(k)
            kz = zaxis(k)
            kx = xaxis(k)
            ky = yaxis(k)
            usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
            if (proceed) then
               xr = x(kk) - xi
               yr = y(kk) - yi
               zr = z(kk) - zi
               if (use_bounds)  call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  corek = monopole(1,k)
                  valk = monopole(2,k)
                  alphak = alphaele(k)
                  c6k = csix(k)
                  overlapk = overpauli(k)
                  apaulk = alphapauli(k)
                  pvalk = monopauli(k)
                  ck = rpole(1,k)
                  dkx = rpole(2,k)
                  dky = rpole(3,k)
                  dkz = rpole(4,k)
                  qkxx = rpole(5,k)
                  qkxy = rpole(6,k)
                  qkxz = rpole(7,k)
                  qkyy = rpole(9,k)
                  qkyz = rpole(10,k)
                  qkzz = rpole(13,k)
c
c     get reciprocal distance terms for this interaction
c
c                  rr1 = f * mscale(kk) / r
                  rr1 = 1.0d0 / r
                  rr3 = rr1 / r2
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr9 = 7.0d0 * rr7 / r2
c
c     intermediates involving moments and distance separation
c
                  dri = dix*xr + diy*yr + diz*zr
                  drk = dkx*xr + dky*yr + dkz*zr
                  dik = dix*dkx + diy*dky + diz*dkz
                  qrix = qixx*xr + qixy*yr + qixz*zr
                  qriy = qixy*xr + qiyy*yr + qiyz*zr
                  qriz = qixz*xr + qiyz*yr + qizz*zr
                  qrkx = qkxx*xr + qkxy*yr + qkxz*zr
                  qrky = qkxy*xr + qkyy*yr + qkyz*zr
                  qrkz = qkxz*xr + qkyz*yr + qkzz*zr
                  qrri = qrix*xr + qriy*yr + qriz*zr
                  qrrk = qrkx*xr + qrky*yr + qrkz*zr
                  qrrik = qrix*qrkx + qriy*qrky + qriz*qrkz
                  qik = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                     + qixx*qkxx + qiyy*qkyy + qizz*qkzz
                  diqrk = dix*qrkx + diy*qrky + diz*qrkz
                  dkqri = dkx*qrix + dky*qriy + dkz*qriz
c
c     calculate intermediate terms for multipole interaction
c
c
c     calculate valence - valence interaction intermediate terms
c
c                  term1 = ci*ck
                  term1ik = vali*valk
                  term2ik = valk*dri - vali*drk + dik
                  term3ik = vali*qrrk + valk*qrri - dri*drk
     &                       + 2.0d0*(dkqri-diqrk+qik)
                  term4ik = dri*qrrk - drk*qrri - 4.0d0*qrrik
                  term5ik = qrri*qrrk
c
c     calculate core - valence interaction intermediate terms
c
                  term1i = corek*vali
                  term2i = corek*dri
                  term3i = corek*qrri
c
c     calculate valence - core interaction intermediate terms
c
                  term1k = corei*valk
                  term2k = -corei*drk
                  term3k = corei*qrrk
c
c     calculate core - core interaction intermediate terms
c
                  term1 = corei*corek
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     electrostatic interaction energy!
c
                  call damphlike(r,rorder,alphai,alphak,
     &                 lambdai,lambdak,lambdaik)
c
c     compute the valence - valence energy contribution for this interaction
c
                  evv = term1ik*rr1*lambdaik(1) + 
     &                  term2ik*rr3*lambdaik(3) + 
     &                  term3ik*rr5*lambdaik(5) +
     &                  term4ik*rr7*lambdaik(7) + 
     &                  term5ik*rr9*lambdaik(9)
c
c     compute the core - valence energy contribution for this interaction
c
                  ecv = term1i*rr1*lambdai(1) +
     &                  term2i*rr3*lambdai(3) +
     &                  term3i*rr5*lambdai(5)
c
                  evc = term1k*rr1*lambdak(1) +
     &                  term2k*rr3*lambdak(3) +
     &                  term3k*rr5*lambdak(5)
c
c     compute the core - core energy contribution for this interaction
c
                  ecc = term1*rr1
c
c     add together energy components for total damped energy
c
                  e_ele = f * mscale(kk) * (evv + ecv + evc + ecc)
c
                  if (use_group)  e_ele = e_ele * fgrp
c
c     save permanent electric field for induced dipole calculation
c
               fid(1) = -xr*(rr3*corek + rr3*lambdak(3)*valk -
     &              rr5*lambdak(5)*drk + rr7*lambdak(7)*qrrk)
     &              - rr3*lambdak(3)*dkx + 2.0d0*rr5*lambdak(5)*qrkx
               fid(2) = -yr*(rr3*corek + rr3*lambdak(3)*valk -
     &              rr5*lambdak(5)*drk+rr7*lambdak(7)*qrrk)
     &              - rr3*lambdak(3)*dky + 2.0d0*rr5*lambdak(5)*qrky
               fid(3) = -zr*(rr3*corek + rr3*lambdak(3)*valk -
     &              rr5*lambdak(5)*drk+rr7*lambdak(7)*qrrk)
     &              - rr3*lambdak(3)*dkz + 2.0d0*rr5*lambdak(5)*qrkz
               fkd(1) = xr*(rr3*corei + rr3*lambdai(3)*vali +
     &              rr5*lambdai(5)*dri + rr7*lambdai(7)*qrri)
     &              - rr3*lambdai(3)*dix - 2.0d0*rr5*lambdai(5)*qrix
               fkd(2) = yr*(rr3*corei + rr3*lambdai(3)*vali +
     &              rr5*lambdai(5)*dri + rr7*lambdai(7)*qrri)
     &              - rr3*lambdai(3)*diy - 2.0d0*rr5*lambdai(5)*qriy
               fkd(3) = zr*(rr3*corei + rr3*lambdai(3)*vali +
     &              rr5*lambdai(5)*dri + rr7*lambdai(7)*qrri)
     &              - rr3*lambdai(3)*diz - 2.0d0*rr5*lambdai(5)*qriz
c
c     increment electric field on both sites
c
               do j = 1, 3
                  permfield(j,i) = permfield(j,i) + fid(j)*mscale(kk)
                  permfield(j,k) = permfield(j,k) + fkd(j)*mscale(kk)
               end do
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     dispersion energy!
c
                  rr6 = rr3**2
c
c     c6 multiplicative combining rule
c
                  c6ik = c6i*c6k
c
c     dispersion damping factor
c
                  displam = 0.5d0*(3.0d0*lambdaik(5) - lambdaik(3))
c
c     compute damped 1/r^6 energy term
c
                  e_disp = -c6ik * rr6 * mscale(kk) * displam**2
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     pauli repulsion energy!
c
                  call damppauli(r,r2,rr1,rr3,rr5,rr7,rr9,rr11,rorder,
     &                 apauli,apaulk,lambdaik)                  
c
c     recompute terms with number of pauli valence electrons
c
                  term1ik = pvali*pvalk
                  term2ik = pvalk*dri - pvali*drk + dik
                  term3ik = pvali*qrrk + pvalk*qrri - dri*drk
     &                       + 2.0d0*(dkqri-diqrk+qik)
c
c     compute valence - valence energy contribution for this interaction
c     (pauli repulsion has no terms involving the core)
c
                  evv = term1ik*lambdaik(1) +
     &                  term2ik*lambdaik(3) +
     &                  term3ik*lambdaik(5) +
     &                  term4ik*lambdaik(7) +
     &                  term5ik*lambdaik(9)
c
c     combining rule for pauli repulsion prefactor
c
                  oik = overlapi*overlapk
c
c     total pauli repulsion energy
c
                  e_pauli = oik * mscale(kk) * evv * rr1
c
                  if (use_group)  e_pauli = e_pauli * fgrp
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     increment the overall multipole energy components
c
                  if (e_ele .ne. 0.0d0) then
                     nem = nem + 1
                     em = em + e_ele
                     aem(ii) = aem(ii) + 0.5d0*e_ele
                     aem(kk) = aem(kk) + 0.5d0*e_ele
                     if (molcule(ii) .ne. molcule(kk))
     &                  einter = einter + e_ele
                  end if
c
c     increment the overall dispersion energy components
c
                  if (e_disp .ne. 0.0d0) then
                     nedis = nedis + 1
                     edis = edis + e_disp
                     aedis(ii) = aedis(ii) + 0.5d0*e_disp
                     aedis(kk) = aedis(kk) + 0.5d0*e_disp
                     if (molcule(ii) .ne. molcule(kk))
     &                  einter = einter + e_disp
                  end if
c
c     increment the overall pauli repulsion energy components
c
                  if (e_pauli .ne. 0.0d0) then
                     nepr = nepr + 1
                     epr = epr + e_pauli
                     aepr(ii) = aepr(ii) + 0.5d0*e_pauli
                     aepr(kk) = aepr(kk) + 0.5d0*e_pauli
                     if (molcule(ii) .ne. molcule(kk))
     &                  einter = einter + e_pauli
                  end if
c
c     print message if the energy of this interaction is large
c
c                  huge = (abs(e) .gt. 100.0d0)
c                  if ((debug.and.e.ne.0.0d0)
c     &                  .or. (verbose.and.huge)) then
c                     if (header) then
c                        header = .false.
c                        write (iout,20)
c   20                   format (/,' Individual Atomic Multipole',
c     &                             ' Interactions :',
c     &                          //,' Type',14x,'Atom Names',
c     &                             15x,'Distance',8x,'Energy',/)
c                     end if
c                     write (iout,30)  ii,name(ii),kk,name(kk),r,e
c   30                format (' M-Pole',4x,2(i7,'-',a3),9x,
c     &                          f10.4,2x,f12.4)
c                  end if
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (use_replica) then
c
c     calculate interaction energy with other unit cells
c
         do i = 1, npole
            ii = ipole(i)
            iz = zaxis(i)
            ix = xaxis(i)
            iy = yaxis(k)
            xi = x(ii)
            yi = y(ii)
            zi = z(ii)
            ci = rpole(1,i)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            qixx = rpole(5,i)
            qixy = rpole(6,i)
            qixz = rpole(7,i)
            qiyy = rpole(9,i)
            qiyz = rpole(10,i)
            qizz = rpole(13,i)
            usei = (use(ii) .or. use(iz) .or. use(ix) .or. use(iy))
            do j = 1, n12(ii)
               mscale(i12(j,ii)) = m2scale
            end do
            do j = 1, n13(ii)
               mscale(i13(j,ii)) = m3scale
            end do
            do j = 1, n14(ii)
               mscale(i14(j,ii)) = m4scale
            end do
            do j = 1, n15(ii)
               mscale(i15(j,ii)) = m5scale
            end do
c
c     evaluate all sites within the cutoff distance
c
            do k = i, npole
               kk = ipole(k)
               kz = zaxis(k)
               kx = xaxis(k)
               ky = yaxis(k)
               usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
               if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
               proceed = .true.
               if (proceed)  proceed = (usei .or. usek)
               if (proceed) then
                  do j = 1, ncell
                     xr = x(kk) - xi
                     yr = y(kk) - yi
                     zr = z(kk) - zi
                     call imager (xr,yr,zr,j)
                     r2 = xr*xr + yr* yr + zr*zr
                     if (.not. (use_polymer .and. r2.le.polycut2))
     &                  mscale(kk) = 1.0d0
                     if (r2 .le. off2) then
                        r = sqrt(r2)
                        ck = rpole(1,k)
                        dkx = rpole(2,k)
                        dky = rpole(3,k)
                        dkz = rpole(4,k)
                        qkxx = rpole(5,k)
                        qkxy = rpole(6,k)
                        qkxz = rpole(7,k)
                        qkyy = rpole(9,k)
                        qkyz = rpole(10,k)
                        qkzz = rpole(13,k)
c
c     get reciprocal distance terms for this interaction
c
                        rr1 = f / r
                        rr3 = rr1 / r2
                        rr5 = 3.0d0 * rr3 / r2
                        rr7 = 5.0d0 * rr5 / r2
                        rr9 = 7.0d0 * rr7 / r2
c
c     intermediates involving moments and distance separation
c
                        dri = dix*xr + diy*yr + diz*zr
                        drk = dkx*xr + dky*yr + dkz*zr
                        dik = dix*dkx + diy*dky + diz*dkz
                        qrix = qixx*xr + qixy*yr + qixz*zr
                        qriy = qixy*xr + qiyy*yr + qiyz*zr
                        qriz = qixz*xr + qiyz*yr + qizz*zr
                        qrkx = qkxx*xr + qkxy*yr + qkxz*zr
                        qrky = qkxy*xr + qkyy*yr + qkyz*zr
                        qrkz = qkxz*xr + qkyz*yr + qkzz*zr
                        qrri = qrix*xr + qriy*yr + qriz*zr
                        qrrk = qrkx*xr + qrky*yr + qrkz*zr
                        qrrik = qrix*qrkx + qriy*qrky + qriz*qrkz
                        qik = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                           + qixx*qkxx + qiyy*qkyy + qizz*qkzz
                        diqrk = dix*qrkx + diy*qrky + diz*qrkz
                        dkqri = dkx*qrix + dky*qriy + dkz*qriz
c
c     calculate intermediate terms for multipole interaction
c
                        term1 = ci*ck
                        term2 = ck*dri - ci*drk + dik
                        term3 = ci*qrrk + ck*qrri - dri*drk
     &                             + 2.0d0*(dkqri-diqrk+qik)
                        term4 = dri*qrrk - drk*qrri - 4.0d0*qrrik
                        term5 = qrri*qrrk
c
c     compute the energy contribution for this interaction
c
                        e = term1*rr1 + term2*rr3 + term3*rr5
     &                         + term4*rr7 + term5*rr9
                        e = e * mscale(kk)
                        if (use_group)  e = e * fgrp
                        if (ii .eq. kk)  e = 0.5d0 * e
c
c     increment the overall multipole energy components
c
                        if (e .ne. 0.0d0) then
                           nem = nem + 1
                           em = em + e
                           aem(ii) = aem(ii) + 0.5d0*e
                           aem(kk) = aem(kk) + 0.5d0*e
                           einter = einter + e
                        end if
c
c     print message if the energy of this interaction is large
c
                        huge = (abs(e) .gt. 100.0d0)
                        if ((debug.and.e.ne.0.0d0)
     &                        .or. (verbose.and.huge)) then
                           if (header) then
                              header = .false.
                              write (iout,40)
   40                         format (/,' Individual Atomic Multipole',
     &                                   ' Interactions :',
     &                                //,' Type',14x,'Atom Names',
     &                                   15x,'Distance',8x,'Energy',/)
                           end if
                           write (iout,50)  ii,name(ii),kk,name(kk),r,e
   50                      format (' M-Pole',4x,2(i7,'-',a3),1x,
     &                                '(XTAL)',2x,f10.4,2x,f12.4)
                        end if
                     end if
                  end do
               end if
            end do
c
c     reset exclusion coefficients for connected atoms
c
            do j = 1, n12(ii)
               mscale(i12(j,ii)) = 1.0d0
            end do
            do j = 1, n13(ii)
               mscale(i13(j,ii)) = 1.0d0
            end do
            do j = 1, n14(ii)
               mscale(i14(j,ii)) = 1.0d0
            end do
            do j = 1, n15(ii)
               mscale(i15(j,ii)) = 1.0d0
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine emoeba3b  --  neighbor list multipole analysis  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "emoeba3b" calculates the atomic multipole interaction energy
c     using a neighbor list, and partitions the energy among the atoms
c
c
      subroutine emoeba3b
      use sizes
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use chgpot
      use couple
      use energi
      use group
      use inform
      use inter
      use iounit
      use math
      use molcul
      use mplpot
      use mpole
      use neigh
      use potent
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      integer ix,iy,iz
      integer kx,ky,kz
      real*8 e,f,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 qrix,qriy,qriz
      real*8 qrkx,qrky,qrkz
      real*8 dri,drk,dik
      real*8 qrri,qrrk
      real*8 qrrik,qik
      real*8 diqrk,dkqri
      real*8 term1,term2,term3
      real*8 term4,term5
      real*8, allocatable :: mscale(:)
      logical proceed
      logical header,huge
      logical usei,usek
      character*6 mode
c
c
c     zero out total atomic multipole energy and partitioning
c
      nem = 0
      em = 0.0d0
      do i = 1, n
         aem(i) = 0.0d0
      end do
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. npole.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Atomic Multipole Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         mscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(npole,ipole,x,y,z,xaxis,yaxis,zaxis,rpole,use,
!$OMP& n12,i12,n13,i13,n14,i14,n15,i15,m2scale,m3scale,m4scale,
!$OMP& m5scale,nelst,elst,use_group,use_intra,use_bounds,off2,
!$OMP& f,molcule,name,verbose,debug,header,iout)
!$OMP& firstprivate(mscale) shared (em,einter,nem,aem)
!$OMP DO reduction(+:em,einter,nem,aem) schedule(guided)
c
c     calculate the multipole interaction energy term
c
      do i = 1, npole
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         usei = (use(ii) .or. use(iz) .or. use(ix) .or. use(iy))
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do kkk = 1, nelst(i)
            k = elst(kkk,i)
            kk = ipole(k)
            kz = zaxis(k)
            kx = xaxis(k)
            ky = yaxis(k)
            usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
            if (proceed) then
               xr = x(kk) - xi
               yr = y(kk) - yi
               zr = z(kk) - zi
               if (use_bounds)  call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  ck = rpole(1,k)
                  dkx = rpole(2,k)
                  dky = rpole(3,k)
                  dkz = rpole(4,k)
                  qkxx = rpole(5,k)
                  qkxy = rpole(6,k)
                  qkxz = rpole(7,k)
                  qkyy = rpole(9,k)
                  qkyz = rpole(10,k)
                  qkzz = rpole(13,k)
c
c     get reciprocal distance terms for this interaction
c
                  rr1 = f * mscale(kk) / r
                  rr3 = rr1 / r2
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr9 = 7.0d0 * rr7 / r2
c
c     intermediates involving moments and distance separation
c
                  dri = dix*xr + diy*yr + diz*zr
                  drk = dkx*xr + dky*yr + dkz*zr
                  dik = dix*dkx + diy*dky + diz*dkz
                  qrix = qixx*xr + qixy*yr + qixz*zr
                  qriy = qixy*xr + qiyy*yr + qiyz*zr
                  qriz = qixz*xr + qiyz*yr + qizz*zr
                  qrkx = qkxx*xr + qkxy*yr + qkxz*zr
                  qrky = qkxy*xr + qkyy*yr + qkyz*zr
                  qrkz = qkxz*xr + qkyz*yr + qkzz*zr
                  qrri = qrix*xr + qriy*yr + qriz*zr
                  qrrk = qrkx*xr + qrky*yr + qrkz*zr
                  qrrik = qrix*qrkx + qriy*qrky + qriz*qrkz
                  qik = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                     + qixx*qkxx + qiyy*qkyy + qizz*qkzz
                  diqrk = dix*qrkx + diy*qrky + diz*qrkz
                  dkqri = dkx*qrix + dky*qriy + dkz*qriz
c
c     calculate intermediate terms for multipole interaction
c
                  term1 = ci*ck
                  term2 = ck*dri - ci*drk + dik
                  term3 = ci*qrrk + ck*qrri - dri*drk
     &                       + 2.0d0*(dkqri-diqrk+qik)
                  term4 = dri*qrrk - drk*qrri - 4.0d0*qrrik
                  term5 = qrri*qrrk
c
c     compute the energy contribution for this interaction
c
                  e = term1*rr1 + term2*rr3 + term3*rr5
     &                   + term4*rr7 + term5*rr9
                  if (use_group)  e = e * fgrp
c
c     increment the overall multipole energy components
c
                  if (e .ne. 0.0d0) then
                     nem = nem + 1
                     em = em + e
                     aem(ii) = aem(ii) + 0.5d0*e
                     aem(kk) = aem(kk) + 0.5d0*e
                     if (molcule(ii) .ne. molcule(kk))
     &                  einter = einter + e
                  end if
c
c     print message if the energy of this interaction is large
c
                  huge = (abs(e) .gt. 100.0d0)
                  if ((debug.and.e.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,20)
   20                   format (/,' Individual Atomic Multipole',
     &                             ' Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             15x,'Distance',8x,'Energy',/)
                     end if
                     write (iout,30)  ii,name(ii),kk,name(kk),r,e
   30                format (' M-Pole',4x,2(i7,'-',a3),9x,
     &                          f10.4,2x,f12.4)
                  end if
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
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
      deallocate (mscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine emoeba3c  --  Ewald multipole analysis via loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "emoeba3c" calculates the atomic multipole interaction energy
c     using a particle mesh Ewald summation and double loop, and
c     partitions the energy among the atoms
c
c
      subroutine emoeba3c
      use sizes
      use action
      use analyz
      use atoms
      use boxes
      use chgpot
      use energi
      use ewald
      use math
      use mpole
      implicit none
      integer i,ii
      real*8 e,f
      real*8 term,fterm
      real*8 cii,dii,qii
      real*8 xd,yd,zd
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
c
c
c     zero out the multipole and polarization energies
c
      nem = 0
      em = 0.0d0
      do i = 1, n
         aem(i) = 0.0d0
      end do
      if (npole .eq. 0)  return
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the real space part of the Ewald summation
c
      call emreal3c
c
c     compute the reciprocal space part of the Ewald summation
c
      call emrecip
c
c     compute the self-energy part of the Ewald summation
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      do i = 1, npole
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         cii = ci*ci
         dii = dix*dix + diy*diy + diz*diz
         qii = 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
     &            + qixx*qixx + qiyy*qiyy + qizz*qizz
         e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
         em = em + e
         nem = nem + 1
         aem(i) = aem(i) + e
      end do
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            xd = xd + dix + rpole(1,i)*x(ii)
            yd = yd + diy + rpole(1,i)*y(ii)
            zd = zd + diz + rpole(1,i)*z(ii)
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         e = term * (xd*xd+yd*yd+zd*zd)
         em = em + e
         nem = nem + 1
         do i = 1, npole
            ii = ipole(i)
            aem(ii) = aem(ii) + e/dble(npole)
         end do
      end if
      return
      end
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine emoebareal3c  --  real space mpole analysis via loop  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "emoebareal3c" evaluates the real space portion of the Ewald sum
c     energy due to atomic multipole interactions and partitions
c     the energy among the atoms
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
      subroutine emoebareal3c
      use sizes
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use cell
      use chgpot
      use couple
      use energi
      use ewald
      use inform
      use inter
      use iounit
      use math
      use molcul
      use mplpot
      use mpole
      use potent
      use shunt
      implicit none
      integer i,j,k
      integer ii,kk
      integer jcell
      real*8 e,efull,f
      real*8 bfac,erfc
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 scalekk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 qrix,qriy,qriz
      real*8 qrkx,qrky,qrkz
      real*8 dri,drk,dik
      real*8 qrri,qrrk
      real*8 qrrik,qik
      real*8 diqrk,dkqri
      real*8 term1,term2,term3
      real*8 term4,term5
      real*8 bn(0:4)
      real*8, allocatable :: mscale(:)
      logical header,huge
      character*6 mode
      external erfc
c
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. npole.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Atomic Multipole Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         mscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     compute the real space portion of the Ewald summation
c
      do i = 1, npole-1
         ii = ipole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do k = i+1, npole
            kk = ipole(k)
            xr = x(kk) - xi
            yr = y(kk) - yi
            zr = z(kk) - zi
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 4
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
               do j = 0, 4
                  bn(j) = f * bn(j)
               end do
c
c     intermediates involving moments and distance separation
c
               dri = dix*xr + diy*yr + diz*zr
               drk = dkx*xr + dky*yr + dkz*zr
               dik = dix*dkx + diy*dky + diz*dkz
               qrix = qixx*xr + qixy*yr + qixz*zr
               qriy = qixy*xr + qiyy*yr + qiyz*zr
               qriz = qixz*xr + qiyz*yr + qizz*zr
               qrkx = qkxx*xr + qkxy*yr + qkxz*zr
               qrky = qkxy*xr + qkyy*yr + qkyz*zr
               qrkz = qkxz*xr + qkyz*yr + qkzz*zr
               qrri = qrix*xr + qriy*yr + qriz*zr
               qrrk = qrkx*xr + qrky*yr + qrkz*zr
               qrrik = qrix*qrkx + qriy*qrky + qriz*qrkz
               qik = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                  + qixx*qkxx + qiyy*qkyy + qizz*qkzz
               diqrk = dix*qrkx + diy*qrky + diz*qrkz
               dkqri = dkx*qrix + dky*qriy + dkz*qriz
c
c     calculate intermediate terms for multipole interaction
c
               term1 = ci*ck
               term2 = ck*dri - ci*drk + dik
               term3 = ci*qrrk + ck*qrri - dri*drk
     &                    + 2.0d0*(dkqri-diqrk+qik)
               term4 = dri*qrrk - drk*qrri - 4.0d0*qrrik
               term5 = qrri*qrrk
c
c     compute the full undamped energy for this interaction
c
               efull = term1*rr1 + term2*rr3 + term3*rr5
     &                    + term4*rr7 + term5*rr9
               efull = mscale(kk) * efull
               if (efull .ne. 0.0d0) then
                  nem = nem + 1
                  aem(ii) = aem(ii) + 0.5d0*efull
                  aem(kk) = aem(kk) + 0.5d0*efull
                  if (molcule(ii) .ne. molcule(kk))
     &               einter = einter + efull
               end if
c
c     modify error function terms to account for scaling
c
               scalekk = 1.0d0 - mscale(kk)
               rr1 = bn(0) - scalekk*rr1
               rr3 = bn(1) - scalekk*rr3
               rr5 = bn(2) - scalekk*rr5
               rr7 = bn(3) - scalekk*rr7
               rr9 = bn(4) - scalekk*rr9
c
c     compute the energy contribution for this interaction
c
               e = term1*rr1 + term2*rr3 + term3*rr5
     &                + term4*rr7 + term5*rr9
               em = em + e
c
c     print a message if the energy of this interaction is large
c
               huge = (abs(efull) .gt. 100.0d0)
               if ((debug.and.efull.ne.0.0d0)
     &               .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,20)
   20                format (/,' Individual Atomic Multipole',
     &                          ' Interactions :',
     &                       //,' Type',14x,'Atom Names',
     &                          15x,'Distance',8x,'Energy',/)
                  end if
                  write (iout,30)  ii,name(ii),kk,name(kk),r,efull
   30             format (' M-Pole',4x,2(i7,'-',a3),9x,
     &                       f10.4,2x,f12.4)
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (use_replica) then
c
c     calculate interaction energy with other unit cells
c
         do i = 1, npole
            ii = ipole(i)
            xi = x(ii)
            yi = y(ii)
            zi = z(ii)
            ci = rpole(1,i)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            qixx = rpole(5,i)
            qixy = rpole(6,i)
            qixz = rpole(7,i)
            qiyy = rpole(9,i)
            qiyz = rpole(10,i)
            qizz = rpole(13,i)
            do j = 1, n12(ii)
               mscale(i12(j,ii)) = m2scale
            end do
            do j = 1, n13(ii)
               mscale(i13(j,ii)) = m3scale
            end do
            do j = 1, n14(ii)
               mscale(i14(j,ii)) = m4scale
            end do
            do j = 1, n15(ii)
               mscale(i15(j,ii)) = m5scale
            end do
c
c     evaluate all sites within the cutoff distance
c
            do k = i, npole
               kk = ipole(k)
               do jcell = 1, ncell
                  xr = x(kk) - xi
                  yr = y(kk) - yi
                  zr = z(kk) - zi
                  call imager (xr,yr,zr,jcell)
                  r2 = xr*xr + yr* yr + zr*zr
                  if (.not. (use_polymer .and. r2.le.polycut2))
     &               mscale(kk) = 1.0d0
                  if (r2 .le. off2) then
                     r = sqrt(r2)
                     ck = rpole(1,k)
                     dkx = rpole(2,k)
                     dky = rpole(3,k)
                     dkz = rpole(4,k)
                     qkxx = rpole(5,k)
                     qkxy = rpole(6,k)
                     qkxz = rpole(7,k)
                     qkyy = rpole(9,k)
                     qkyz = rpole(10,k)
                     qkzz = rpole(13,k)
c
c     get reciprocal distance terms for this interaction
c
                     rr1 = f / r
                     rr3 = rr1 / r2
                     rr5 = 3.0d0 * rr3 / r2
                     rr7 = 5.0d0 * rr5 / r2
                     rr9 = 7.0d0 * rr7 / r2
c
c     calculate the real space Ewald error function terms
c
                     ralpha = aewald * r
                     bn(0) = erfc(ralpha) / r
                     alsq2 = 2.0d0 * aewald**2
                     alsq2n = 0.0d0
                     if (aewald .gt. 0.0d0)
     &                  alsq2n = 1.0d0 / (sqrtpi*aewald)
                     exp2a = exp(-ralpha**2)
                     do j = 1, 4
                        bfac = dble(j+j-1)
                        alsq2n = alsq2 * alsq2n
                        bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
                     end do
                     do j = 0, 4
                        bn(j) = f * bn(j)
                     end do
c
c     intermediates involving moments and distance separation
c
                     dri = dix*xr + diy*yr + diz*zr
                     drk = dkx*xr + dky*yr + dkz*zr
                     dik = dix*dkx + diy*dky + diz*dkz
                     qrix = qixx*xr + qixy*yr + qixz*zr
                     qriy = qixy*xr + qiyy*yr + qiyz*zr
                     qriz = qixz*xr + qiyz*yr + qizz*zr
                     qrkx = qkxx*xr + qkxy*yr + qkxz*zr
                     qrky = qkxy*xr + qkyy*yr + qkyz*zr
                     qrkz = qkxz*xr + qkyz*yr + qkzz*zr
                     qrri = qrix*xr + qriy*yr + qriz*zr
                     qrrk = qrkx*xr + qrky*yr + qrkz*zr
                     qrrik = qrix*qrkx + qriy*qrky + qriz*qrkz
                     qik = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                        + qixx*qkxx + qiyy*qkyy + qizz*qkzz
                     diqrk = dix*qrkx + diy*qrky + diz*qrkz
                     dkqri = dkx*qrix + dky*qriy + dkz*qriz
c
c     calculate intermediate terms for multipole interaction
c
                     term1 = ci*ck
                     term2 = ck*dri - ci*drk + dik
                     term3 = ci*qrrk + ck*qrri - dri*drk
     &                          + 2.0d0*(dkqri-diqrk+qik)
                     term4 = dri*qrrk - drk*qrri - 4.0d0*qrrik
                     term5 = qrri*qrrk
c
c     compute the full undamped energy for this interaction
c
                     efull = term1*rr1 + term2*rr3 + term3*rr5
     &                          + term4*rr7 + term5*rr9
                     efull = mscale(kk) * efull
                     if (ii .eq. kk)  efull = 0.5d0 * efull
                     if (efull .ne. 0.0d0) then
                        nem = nem + 1
                        aem(ii) = aem(ii) + 0.5d0*efull
                        aem(kk) = aem(kk) + 0.5d0*efull
                        einter = einter + efull
                     end if
c
c     modify distances to account for Ewald and exclusions
c
                     scalekk = 1.0d0 - mscale(kk)
                     rr1 = bn(0) - scalekk*rr1
                     rr3 = bn(1) - scalekk*rr3
                     rr5 = bn(2) - scalekk*rr5
                     rr7 = bn(3) - scalekk*rr7
                     rr9 = bn(4) - scalekk*rr9
c
c     compute the energy contribution for this interaction
c
                     e = term1*rr1 + term2*rr3 + term3*rr5
     &                      + term4*rr7 + term5*rr9
                     if (ii .eq. kk)  e = 0.5d0 * e
                     em = em + e
c
c     print message if the energy of this interaction is large
c
                     huge = (abs(efull) .gt. 100.0d0)
                     if ((debug.and.efull.ne.0.0d0)
     &                     .or. (verbose.and.huge)) then
                        if (header) then
                           header = .false.
                           write (iout,40)
   40                      format (/,' Individual Atomic Multipole',
     &                                ' Interactions :',
     &                             //,' Type',14x,'Atom Names',
     &                                15x,'Distance',8x,'Energy',/)
                        end if
                        write (iout,50)  ii,name(ii),kk,name(kk),r,efull
   50                   format (' M-Pole',4x,2(i7,'-',a3),1x,
     &                             '(XTAL)',2x,f10.4,2x,f12.4)
                     end if
                  end if
               end do
            end do
c
c     reset exclusion coefficients for connected atoms
c
            do j = 1, n12(ii)
               mscale(i12(j,ii)) = 1.0d0
            end do
            do j = 1, n13(ii)
               mscale(i13(j,ii)) = 1.0d0
            end do
            do j = 1, n14(ii)
               mscale(i14(j,ii)) = 1.0d0
            end do
            do j = 1, n15(ii)
               mscale(i15(j,ii)) = 1.0d0
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine emoeba3d  --  Ewald multipole analysis via list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "emoeba3d" calculates the atomic multipole interaction energy
c     using particle mesh Ewald summation and a neighbor list, and
c     partitions the energy among the atoms
c
c
      subroutine emoeba3d
      use sizes
      use action
      use analyz
      use atoms
      use boxes
      use chgpot
      use energi
      use ewald
      use math
      use mpole
      implicit none
      integer i,ii
      real*8 e,f
      real*8 term,fterm
      real*8 cii,dii,qii
      real*8 xd,yd,zd
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
c
c
c     zero out the multipole and polarization energies
c
      nem = 0
      em = 0.0d0
      do i = 1, n
         aem(i) = 0.0d0
      end do
      if (npole .eq. 0)  return
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the real space part of the Ewald summation
c
      call emreal3d
c
c     compute the reciprocal space part of the Ewald summation
c
      call emrecip
c
c     compute the self-energy part of the Ewald summation
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      do i = 1, npole
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         cii = ci*ci
         dii = dix*dix + diy*diy + diz*diz
         qii = 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
     &            + qixx*qixx + qiyy*qiyy + qizz*qizz
         e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
         em = em + e
         nem = nem + 1
         aem(i) = aem(i) + e
      end do
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            xd = xd + dix + rpole(1,i)*x(ii)
            yd = yd + diy + rpole(1,i)*y(ii)
            zd = zd + diz + rpole(1,i)*z(ii)
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         e = term * (xd*xd+yd*yd+zd*zd)
         em = em + e
         nem = nem + 1
         do i = 1, npole
            ii = ipole(i)
            aem(ii) = aem(ii) + e/dble(npole)
         end do
      end if
      return
      end
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine emoebareal3d  --  real space mpole analysis via list  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "emoebareal3d" evaluates the real space portion of the Ewald sum
c     energy due to atomic multipole interactions, and partitions
c     the energy among the atoms using a pairwise neighbor list
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
      subroutine emoebareal3d
      use sizes
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use chgpot
      use couple
      use energi
      use ewald
      use inform
      use inter
      use iounit
      use math
      use molcul
      use mplpot
      use mpole
      use neigh
      use potent
      use shunt
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      real*8 e,efull,f
      real*8 bfac,erfc
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 scalekk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 qrix,qriy,qriz
      real*8 qrkx,qrky,qrkz
      real*8 dri,drk,dik
      real*8 qrri,qrrk
      real*8 qrrik,qik
      real*8 diqrk,dkqri
      real*8 term1,term2,term3
      real*8 term4,term5
      real*8 bn(0:4)
      real*8, allocatable :: mscale(:)
      logical header,huge
      character*6 mode
      external erfc
c
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. npole.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Atomic Multipole Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         mscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(npole,ipole,x,y,z,rpole,n12,i12,n13,i13,n14,i14,n15,
!$OMP& i15,m2scale,m3scale,m4scale,m5scale,nelst,elst,use_bounds,
!$OMP& f,off2,aewald,molcule,name,verbose,debug,header,iout)
!$OMP& firstprivate(mscale) shared (em,einter,nem,aem)
!$OMP DO reduction(+:em,einter,nem,aem) schedule(guided)
c
c     compute the real space portion of the Ewald summation
c
      do i = 1, npole
         ii = ipole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do kkk = 1, nelst(i)
            k = elst(kkk,i)
            kk = ipole(k)
            xr = x(kk) - xi
            yr = y(kk) - yi
            zr = z(kk) - zi
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 4
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
               do j = 0, 4
                  bn(j) = f * bn(j)
               end do
c
c     intermediates involving moments and distance separation
c
               dri = dix*xr + diy*yr + diz*zr
               drk = dkx*xr + dky*yr + dkz*zr
               dik = dix*dkx + diy*dky + diz*dkz
               qrix = qixx*xr + qixy*yr + qixz*zr
               qriy = qixy*xr + qiyy*yr + qiyz*zr
               qriz = qixz*xr + qiyz*yr + qizz*zr
               qrkx = qkxx*xr + qkxy*yr + qkxz*zr
               qrky = qkxy*xr + qkyy*yr + qkyz*zr
               qrkz = qkxz*xr + qkyz*yr + qkzz*zr
               qrri = qrix*xr + qriy*yr + qriz*zr
               qrrk = qrkx*xr + qrky*yr + qrkz*zr
               qrrik = qrix*qrkx + qriy*qrky + qriz*qrkz
               qik = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                  + qixx*qkxx + qiyy*qkyy + qizz*qkzz
               diqrk = dix*qrkx + diy*qrky + diz*qrkz
               dkqri = dkx*qrix + dky*qriy + dkz*qriz
c
c     calculate intermediate terms for multipole interaction
c
               term1 = ci*ck
               term2 = ck*dri - ci*drk + dik
               term3 = ci*qrrk + ck*qrri - dri*drk
     &                    + 2.0d0*(dkqri-diqrk+qik)
               term4 = dri*qrrk - drk*qrri - 4.0d0*qrrik
               term5 = qrri*qrrk
c
c     compute the full undamped energy for this interaction
c
               efull = term1*rr1 + term2*rr3 + term3*rr5
     &                    + term4*rr7 + term5*rr9
               efull = mscale(kk) * efull
               if (efull .ne. 0.0d0) then
                  nem = nem + 1
                  aem(ii) = aem(ii) + 0.5d0*efull
                  aem(kk) = aem(kk) + 0.5d0*efull
                  if (molcule(ii) .ne. molcule(kk))
     &               einter = einter + efull
               end if
c
c     modify error function terms to account for scaling
c
               scalekk = 1.0d0 - mscale(kk)
               rr1 = bn(0) - scalekk*rr1
               rr3 = bn(1) - scalekk*rr3
               rr5 = bn(2) - scalekk*rr5
               rr7 = bn(3) - scalekk*rr7
               rr9 = bn(4) - scalekk*rr9
c
c     compute the energy contribution for this interaction
c
               e = term1*rr1 + term2*rr3 + term3*rr5
     &                + term4*rr7 + term5*rr9
               em = em + e
c
c     print a message if the energy of this interaction is large
c
               huge = (abs(efull) .gt. 100.0d0)
               if ((debug.and.efull.ne.0.0d0)
     &               .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,20)
   20                format (/,' Individual Atomic Multipole',
     &                          ' Interactions :',
     &                       //,' Type',14x,'Atom Names',
     &                          15x,'Distance',8x,'Energy',/)
                  end if
                  write (iout,30)  ii,name(ii),kk,name(kk),r,efull
   30             format (' M-Pole',4x,2(i7,'-',a3),9x,
     &                       f10.4,2x,f12.4)
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
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
      deallocate (mscale)
      return
      end
