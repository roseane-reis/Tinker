c
c
c     #############################################################
c     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ejosh1  --  mpole/polar energy & derivatives  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ejosh1" calculates the multipole and dipole polarization
c     energy and derivatives with respect to Cartesian coordinates
c
c
      subroutine ejosh1
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
            call ejosh1d
         else
            call ejosh1c
         end if
      else
         if (use_mlist) then
            call ejosh1b
         else
            call ejosh1a
         end if
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine ejosh1a  --  double loop multipole derivatives  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "ejosh1a" calculates the multipole energy and derivatives with 
c     respect to Cartesian coordinates using a pairwise double loop
c
c
      subroutine ejosh1a
      use sizes
      use atoms
      use bound
      use cell
      use chgpot
      use couple
      use deriv
      use energi
      use group
      use inter
      use molcul
      use mplpot
      use mpole
      use polgrp
      use polpot
      use shunt
      use usage
      use virial
      use chgpen
      use disp
      use pauli
      implicit none
      integer i,j,k
      integer ii,kk,jcell
      integer ix,iy,iz
      integer kx,ky,kz
      integer iax,iay,iaz
      integer rorder
      real*8 e,de,f,fgrp
      real*8 ecc,ecv,evc,evv
      real*8 dei,dek,deik
      real*8 e_ele,e_disp,e_pauli
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 r7,r8
      real*8 alphai,alphak
      real*8 corei,vali
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 corek,valk
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dikx,diky,dikz
      real*8 dirx,diry,dirz
      real*8 dkrx,dkry,dkrz
      real*8 qrix,qriy,qriz
      real*8 qrkx,qrky,qrkz
      real*8 qrixr,qriyr,qrizr
      real*8 qrkxr,qrkyr,qrkzr
      real*8 qrrx,qrry,qrrz
      real*8 qikrx,qikry,qikrz
      real*8 qkirx,qkiry,qkirz
      real*8 qikrxr,qikryr,qikrzr
      real*8 qkirxr,qkiryr,qkirzr
      real*8 diqkx,diqky,diqkz
      real*8 dkqix,dkqiy,dkqiz
      real*8 diqkxr,diqkyr,diqkzr
      real*8 dkqixr,dkqiyr,dkqizr
      real*8 dqiqkx,dqiqky,dqiqkz
      real*8 dri,drk,qrri,qrrk
      real*8 diqrk,dkqri
      real*8 dik,qik,qrrik
      real*8 term1,term2,term3
      real*8 term4,term5,term6
      real*8 term1ik,term2ik,term3ik
      real*8 term4ik,term5ik,term6ik
      real*8 dterm1ik,dterm2ik,dterm3ik
      real*8 dterm4ik,dterm5ik,dterm6ik
      real*8 term1i,term2i,term3i,term4i
      real*8 term1k,term2k,term3k,term5k
      real*8 dterm1i,dterm2i,dterm3i,dterm4i
      real*8 dterm1k,dterm2k,dterm3k,dterm5k
      real*8 rr6
      real*8 c6i,c6k,c6ik
      real*8 displam
      real*8 pvali,pvalk
      real*8 overlapi,overlapk,oik
      real*8 apauli,apaulk
      real*8 frcx,frcy,frcz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 ttmi(3),ttmk(3)
      real*8 fid(3),fkd(3)
      real*8 fix(3),fiy(3),fiz(3)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: tem(:,:)
      real*8, allocatable :: tepr(:,:)
      real*8, allocatable :: lambdai(:)
      real*8, allocatable :: lambdak(:)
      real*8, allocatable :: lambdaik(:)
      logical proceed,usei,usek
      character*6 mode
c
c
c     zero out the atomic multipole energy and derivatives
c
      em = 0.0d0
      edis = 0.0d0
      epr = 0.0d0
      do i = 1, n
         do j = 1, 3
            dem(j,i) = 0.0d0
            dedis(j,i) = 0.0d0
            depr(j,i) = 0.0d0
            permfield(j,i) = 0.0d0
         end do
      end do
      if (npole .eq. 0)  return
c
c     set maximum power of 1/r for damping 
c     (11 for quadrupole-quadrupole force) 
c
      rorder = 11
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (dscale(n))
      allocate (tem(3,n))
      allocate (tepr(3,n))
      allocate (lambdai(rorder))
      allocate (lambdak(rorder))
      allocate (lambdaik(rorder))
c
c     initialize connected atom scaling and torque arrays
c
      do i = 1, n
         mscale(i) = 1.0d0
         dscale(i) = 1.0d0
         do j = 1, 3
            tem(j,i) = 0.0d0
            tepr(j,i) = 0.0d0
         end do
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
c     compute the multipole interaction energy and gradient
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
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
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
            if (.not. proceed)  goto 10
            xr = x(kk) - xi
            yr = y(kk) - yi
            zr = z(kk) - zi
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
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
c               rr1 = f * mscale(kk) / r
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               rr11 = 9.0d0 * rr9 / r2
c
c     intermediates involving moments and distance separation
c
               dikx = diy*dkz - diz*dky
               diky = diz*dkx - dix*dkz
               dikz = dix*dky - diy*dkx
               dirx = diy*zr - diz*yr
               diry = diz*xr - dix*zr
               dirz = dix*yr - diy*xr
               dkrx = dky*zr - dkz*yr
               dkry = dkz*xr - dkx*zr
               dkrz = dkx*yr - dky*xr
c     same as dir in dfield0a
               dri = dix*xr + diy*yr + diz*zr
c     same as dkr in dfield0a
               drk = dkx*xr + dky*yr + dkz*zr
c
               dik = dix*dkx + diy*dky + diz*dkz
c     same as qi* in dfield0a
               qrix = qixx*xr + qixy*yr + qixz*zr
               qriy = qixy*xr + qiyy*yr + qiyz*zr
               qriz = qixz*xr + qiyz*yr + qizz*zr
c     same as qk* in dfield0a
               qrkx = qkxx*xr + qkxy*yr + qkxz*zr
               qrky = qkxy*xr + qkyy*yr + qkyz*zr
               qrkz = qkxz*xr + qkyz*yr + qkzz*zr
c     same as qir in dfield0a
               qrri = qrix*xr + qriy*yr + qriz*zr
c     same as qkr in dfield0a
               qrrk = qrkx*xr + qrky*yr + qrkz*zr
c
               qrrik = qrix*qrkx + qriy*qrky + qriz*qrkz
               qik = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                  + qixx*qkxx + qiyy*qkyy + qizz*qkzz
               qrixr = qriz*yr - qriy*zr
               qriyr = qrix*zr - qriz*xr
               qrizr = qriy*xr - qrix*yr
               qrkxr = qrkz*yr - qrky*zr
               qrkyr = qrkx*zr - qrkz*xr
               qrkzr = qrky*xr - qrkx*yr
               qrrx = qrky*qriz - qrkz*qriy
               qrry = qrkz*qrix - qrkx*qriz
               qrrz = qrkx*qriy - qrky*qrix
               qikrx = qixx*qrkx + qixy*qrky + qixz*qrkz
               qikry = qixy*qrkx + qiyy*qrky + qiyz*qrkz
               qikrz = qixz*qrkx + qiyz*qrky + qizz*qrkz
               qkirx = qkxx*qrix + qkxy*qriy + qkxz*qriz
               qkiry = qkxy*qrix + qkyy*qriy + qkyz*qriz
               qkirz = qkxz*qrix + qkyz*qriy + qkzz*qriz
               qikrxr = qikrz*yr - qikry*zr
               qikryr = qikrx*zr - qikrz*xr
               qikrzr = qikry*xr - qikrx*yr
               qkirxr = qkirz*yr - qkiry*zr
               qkiryr = qkirx*zr - qkirz*xr
               qkirzr = qkiry*xr - qkirx*yr
               diqkx = dix*qkxx + diy*qkxy + diz*qkxz
               diqky = dix*qkxy + diy*qkyy + diz*qkyz
               diqkz = dix*qkxz + diy*qkyz + diz*qkzz
               dkqix = dkx*qixx + dky*qixy + dkz*qixz
               dkqiy = dkx*qixy + dky*qiyy + dkz*qiyz
               dkqiz = dkx*qixz + dky*qiyz + dkz*qizz
               diqrk = dix*qrkx + diy*qrky + diz*qrkz
               dkqri = dkx*qrix + dky*qriy + dkz*qriz
               diqkxr = diqkz*yr - diqky*zr
               diqkyr = diqkx*zr - diqkz*xr
               diqkzr = diqky*xr - diqkx*yr
               dkqixr = dkqiz*yr - dkqiy*zr
               dkqiyr = dkqix*zr - dkqiz*xr
               dkqizr = dkqiy*xr - dkqix*yr
               dqiqkx = diy*qrkz - diz*qrky + dky*qriz - dkz*qriy
     &                     - 2.0d0*(qixy*qkxz+qiyy*qkyz+qiyz*qkzz
     &                             -qixz*qkxy-qiyz*qkyy-qizz*qkyz)
               dqiqky = diz*qrkx - dix*qrkz + dkz*qrix - dkx*qriz
     &                     - 2.0d0*(qixz*qkxx+qiyz*qkxy+qizz*qkxz
     &                             -qixx*qkxz-qixy*qkyz-qixz*qkzz)
               dqiqkz = dix*qrky - diy*qrkx + dkx*qriy - dky*qrix
     &                     - 2.0d0*(qixx*qkxy+qixy*qkyy+qixz*qkyz
     &                             -qixy*qkxx-qiyy*qkxy-qiyz*qkxz)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     electrostatic energy and gradient!
c
c     calculate valence - valence interaction intermediate terms
c
               term1ik = vali*valk
               term2ik = valk*dri - vali*drk + dik
               term3ik = vali*qrrk + valk*qrri - dri*drk
     &              + 2.0d0*(dkqri-diqrk+qik)
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
c     compute damping function
c
               call damphlike(r,rorder,alphai,alphak,
     &              lambdai,lambdak,lambdaik)
c
c     compute the valence - valence energy contribution for this interaction
c
               evv = term1ik*rr1*lambdaik(1) + 
     &              term2ik*rr3*lambdaik(3) + 
     &              term3ik*rr5*lambdaik(5) +
     &              term4ik*rr7*lambdaik(7) + 
     &              term5ik*rr9*lambdaik(9)
c     
c     compute the core - valence energy contribution for this interaction
c     
               ecv = term1i*rr1*lambdai(1) +
     &              term2i*rr3*lambdai(3) +
     &              term3i*rr5*lambdai(5)
c     
               evc = term1k*rr1*lambdak(1) +
     &              term2k*rr3*lambdak(3) +
     &              term3k*rr5*lambdak(5)
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
               em = em + e_ele
               if (molcule(ii) .ne. molcule(kk))
     &            einter = einter + e_ele
c
c     calculate intermediate terms for force and torque
c
c     valence i - valence k
c
               deik = term1ik*rr3*lambdaik(3) + 
     &              term2ik*rr5*lambdaik(5) + 
     &              term3ik*rr7*lambdaik(7) +
     &              term4ik*rr9*lambdaik(9) +
     &              term5ik*rr11*lambdaik(11)
               dterm1ik = -valk*rr3*lambdaik(3) + 
     &              drk*rr5*lambdaik(5) - 
     &              qrrk*rr7*lambdaik(7)
               dterm2ik = vali*rr3*lambdaik(3) + 
     &              dri*rr5*lambdaik(5) + 
     &              qrri*rr7*lambdaik(7)
               dterm3ik = 2.0d0 * rr5 * lambdaik(5)
               dterm4ik = 2.0d0 * (-valk*rr5*lambdaik(5) +
     &              drk*rr7*lambdaik(7) - 
     &              qrrk*rr9*lambdaik(9))
               dterm5ik = 2.0d0 * (-vali*rr5*lambdaik(5) -
     &              dri*rr7*lambdaik(7) - 
     &              qrri*rr9*lambdaik(9))
               dterm6ik = 4.0d0 * rr7*lambdaik(7)
c
c     compute the force components for this interaction
c
               frcx = deik*xr + dterm1ik*dix + dterm2ik*dkx
     &                   + dterm3ik*(diqkx-dkqix) + dterm4ik*qrix
     &                   + dterm5ik*qrkx + dterm6ik*(qikrx+qkirx)
               frcy = deik*yr + dterm1ik*diy + dterm2ik*dky
     &                   + dterm3ik*(diqky-dkqiy) + dterm4ik*qriy
     &                   + dterm5ik*qrky + dterm6ik*(qikry+qkiry)
               frcz = deik*zr + dterm1ik*diz + dterm2ik*dkz
     &                   + dterm3ik*(diqkz-dkqiz) + dterm4ik*qriz
     &                   + dterm5ik*qrkz + dterm6ik*(qikrz+qkirz)
c
c     save permanent electric field for induced dipole calculation
c     note: this already has the core contribution
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
c     compute the torque components for this interaction
c
               ttmi(1) = -rr3*lambdaik(3)*dikx + dterm1ik*dirx + 
     &              dterm3ik*(dqiqkx+dkqixr)
     &                      - dterm4ik*qrixr - dterm6ik*(qikrxr+qrrx)
               ttmi(2) = -rr3*lambdaik(3)*diky + dterm1ik*diry + 
     &              dterm3ik*(dqiqky+dkqiyr)
     &                      - dterm4ik*qriyr - dterm6ik*(qikryr+qrry)
               ttmi(3) = -rr3*lambdaik(3)*dikz + dterm1ik*dirz + 
     &              dterm3ik*(dqiqkz+dkqizr)
     &                      - dterm4ik*qrizr - dterm6ik*(qikrzr+qrrz)
               ttmk(1) = rr3*lambdaik(3)*dikx + dterm2ik*dkrx - 
     &              dterm3ik*(dqiqkx+diqkxr)
     &                      - dterm5ik*qrkxr - dterm6ik*(qkirxr-qrrx)
               ttmk(2) = rr3*lambdaik(3)*diky + dterm2ik*dkry - 
     &              dterm3ik*(dqiqky+diqkyr)
     &                      - dterm5ik*qrkyr - dterm6ik*(qkiryr-qrry)
               ttmk(3) = rr3*lambdaik(3)*dikz + dterm2ik*dkrz - 
     &              dterm3ik*(dqiqkz+diqkzr)
     &                      - dterm5ik*qrkzr - dterm6ik*(qkirzr-qrrz)
c
c     valence - core force and torque
c
               dei = term1i*rr3*lambdai(3) +
     &              term2i*rr5*lambdai(5) +
     &              term3i*rr7*lambdai(7) 
               dek = term1k*rr3*lambdak(3) +
     &              term2k*rr5*lambdak(5) +
     &              term3k*rr7*lambdak(7) 
               dterm1i = -corek*rr3*lambdai(3)
               dterm2k = corei*rr3*lambdak(3) 
               dterm4i = 2.0d0 * (-corek*rr5*lambdai(5))
               dterm5k = 2.0d0 * (-corei*rr5*lambdak(5)) 
c
               frcx = frcx + (dei + dek)*xr + dterm1i*dix + dterm2k*dkx
     &                   + dterm4i*qrix
     &                   + dterm5k*qrkx  
               frcy = frcy + (dei + dek)*yr + dterm1i*diy + dterm2k*dky
     &                   + dterm4i*qriy
     &                   + dterm5k*qrky  
               frcz = frcz + (dei + dek)*zr + dterm1i*diz + dterm2k*dkz
     &                   + dterm4i*qriz
     &                   + dterm5k*qrkz  
c
               ttmi(1) = ttmi(1) + dterm1i*dirx 
     &              - dterm4i*qrixr 
               ttmi(2) = ttmi(2) + dterm1i*diry 
     &              - dterm4i*qriyr 
               ttmi(3) = ttmi(3) + dterm1i*dirz 
     &              - dterm4i*qrizr 
               ttmk(1) = ttmk(1) + dterm2k*dkrx 
     &              - dterm5k*qrkxr 
               ttmk(2) = ttmk(2) + dterm2k*dkry 
     &              - dterm5k*qrkyr 
               ttmk(3) = ttmk(3) + dterm2k*dkrz 
     &              - dterm5k*qrkzr 
c
c     core - core force (no torque)
c
               de = term1*rr3

               frcx = frcx + de*xr 
               frcy = frcy + de*yr 
               frcz = frcz + de*zr
c
c     add in permitivity
c
               frcx = frcx * f * mscale(kk)
               frcy = frcy * f * mscale(kk)
               frcz = frcz * f * mscale(kk)
               ttmi(1) = ttmi(1) * f * mscale(kk)
               ttmi(2) = ttmi(2) * f * mscale(kk)
               ttmi(3) = ttmi(3) * f * mscale(kk)
               ttmk(1) = ttmk(1) * f * mscale(kk)
               ttmk(2) = ttmk(2) * f * mscale(kk)
               ttmk(3) = ttmk(3) * f * mscale(kk)
c
c     force and torque components scaled by group membership
c
               if (use_group) then
                  frcx = fgrp * frcx
                  frcy = fgrp * frcy
                  frcz = fgrp * frcz
                  do j = 1, 3
                     ttmi(j) = fgrp * ttmi(j)
                     ttmk(j) = fgrp * ttmk(j)
                  end do
               end if
c
c     increment force-based gradient and torque on first site
c
               dem(1,ii) = dem(1,ii) + frcx
               dem(2,ii) = dem(2,ii) + frcy
               dem(3,ii) = dem(3,ii) + frcz
               tem(1,i) = tem(1,i) + ttmi(1)
               tem(2,i) = tem(2,i) + ttmi(2)
               tem(3,i) = tem(3,i) + ttmi(3)
c
c     increment force-based gradient and torque on second site
c
               dem(1,kk) = dem(1,kk) - frcx
               dem(2,kk) = dem(2,kk) - frcy
               dem(3,kk) = dem(3,kk) - frcz
               tem(1,k) = tem(1,k) + ttmk(1)
               tem(2,k) = tem(2,k) + ttmk(2)
               tem(3,k) = tem(3,k) + ttmk(3)
c
c     increment electric field on both sites
c
               do j = 1, 3
                  permfield(j,i) = permfield(j,i) + fid(j)*dscale(kk)
                  permfield(j,k) = permfield(j,k) + fkd(j)*dscale(kk)
c
c                  permfield(j,i) = permfield(j,i) + fid(j)*dscale(kk)
c                  permfield(j,k) = permfield(j,k) + fkd(j)*dscale(kk)
c                  permfieldp(j,i) = permfieldp(j,i) + fid(j)*pscale(kk)
c                  permfieldp(j,k) = permfieldp(j,k) + fkd(j)*pscale(kk)
               end do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     dispersion energy and gradient!
c
               rr6 = rr3**2
               r7 = r**7
               r8 = r**8
               c6ik = c6i*c6k
               displam = 0.5d0*(3.0d0*lambdaik(5) - lambdaik(3))
               e_disp = -c6ik * rr6 * mscale(kk) * displam**2
c
               if (use_group)  e_disp = e_disp * fgrp
               edis = edis + e_disp
               if (molcule(ii) .ne. molcule(kk))
     &            einter = einter + e_disp
c
c     get specific terms for damped dispersion gradient
c
c     semi-hack: i'm using lambda(10) to stuff in the neccessary terms
c     see DispersionDamping.mw notebook
c
               de = c6ik*6.0d0*(displam**2)/r8  -
     &              c6ik*2.0d0*displam*lambdaik(10)/r7
c
               frcx = de * xr * mscale(kk)
               frcy = de * yr * mscale(kk)
               frcz = de * zr * mscale(kk)
c
c     force components scaled by group membership
c
               if (use_group) then
                  frcx = fgrp * frcx
                  frcy = fgrp * frcy
                  frcz = fgrp * frcz
               end if
c
c     increment force-based gradient on first site 
c
               dedis(1,ii) = dedis(1,ii) - frcx
               dedis(2,ii) = dedis(2,ii) - frcy
               dedis(3,ii) = dedis(3,ii) - frcz
c
c     increment force-based gradient on second site 
c
               dedis(1,kk) = dedis(1,kk) + frcx
               dedis(2,kk) = dedis(2,kk) + frcy
               dedis(3,kk) = dedis(3,kk) + frcz
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     pauli repulsion energy and gradient !
c
               call damppauli(r,r2,rr1,rr3,rr5,rr7,rr9,rr11,rorder,
     &              apauli,apaulk,lambdaik)
c
c     recompute terms with number of pauli valence electrons 
c
               term1ik = pvali*pvalk
               term2ik = pvalk*dri - pvali*drk + dik
               term3ik = pvali*qrrk + pvalk*qrri - dri*drk
     &              + 2.0d0*(dkqri-diqrk+qik)
c     
c     compute valence - valence energy contribution for this interaction
c     (pauli repulsion has no terms involving the core)
c
               evv = term1ik*lambdaik(1) +
     &              term2ik*lambdaik(3) +
     &              term3ik*lambdaik(5) +
     &              term4ik*lambdaik(7) +
     &              term5ik*lambdaik(9)
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
               epr = epr + e_pauli
               if (molcule(ii) .ne. molcule(kk))
     &              einter = einter + e_pauli
c     
c     now compute gradient in same way as electrostatics
c
c
c     calculate intermediate terms for force and torque
c
c     valence i - valence k
c
               deik = term1ik*lambdaik(3) +
     &              term2ik*lambdaik(5) +
     &              term3ik*lambdaik(7) +
     &              term4ik*lambdaik(9) +
     &              term5ik*lambdaik(11)
               dterm1ik = -pvalk*lambdaik(3) +
     &              drk*lambdaik(5) -
     &              qrrk*lambdaik(7)
               dterm2ik = pvali*lambdaik(3) +
     &              dri*lambdaik(5) +
     &              qrri*lambdaik(7)
               dterm3ik = 2.0d0 * lambdaik(5)
               dterm4ik = 2.0d0 * (-pvalk*lambdaik(5) +
     &              drk*lambdaik(7) -
     &              qrrk*lambdaik(9))
               dterm5ik = 2.0d0 * (-pvali*lambdaik(5) -
     &              dri*lambdaik(7) -
     &              qrri*lambdaik(9))
               dterm6ik = 4.0d0 * lambdaik(7)
c     
c     compute the force components for this interaction
c     
               frcx = deik*xr + dterm1ik*dix + dterm2ik*dkx
     &              + dterm3ik*(diqkx-dkqix) + dterm4ik*qrix
     &              + dterm5ik*qrkx + dterm6ik*(qikrx+qkirx)
               frcy = deik*yr + dterm1ik*diy + dterm2ik*dky
     &              + dterm3ik*(diqky-dkqiy) + dterm4ik*qriy
     &              + dterm5ik*qrky + dterm6ik*(qikry+qkiry)
               frcz = deik*zr + dterm1ik*diz + dterm2ik*dkz
     &              + dterm3ik*(diqkz-dkqiz) + dterm4ik*qriz
     &              + dterm5ik*qrkz + dterm6ik*(qikrz+qkirz)
c
               frcx = frcx*rr1
               frcy = frcy*rr1
               frcz = frcz*rr1
c
c     compute the torque components for this interaction
c
               ttmi(1) = -lambdaik(3)*dikx + dterm1ik*dirx + 
     &              dterm3ik*(dqiqkx+dkqixr)
     &              - dterm4ik*qrixr - dterm6ik*(qikrxr+qrrx)
               ttmi(2) = -lambdaik(3)*diky + dterm1ik*diry + 
     &              dterm3ik*(dqiqky+dkqiyr)
     &              - dterm4ik*qriyr - dterm6ik*(qikryr+qrry)
               ttmi(3) = -lambdaik(3)*dikz + dterm1ik*dirz + 
     &              dterm3ik*(dqiqkz+dkqizr)
     &              - dterm4ik*qrizr - dterm6ik*(qikrzr+qrrz)
               ttmk(1) = lambdaik(3)*dikx + dterm2ik*dkrx - 
     &              dterm3ik*(dqiqkx+diqkxr)
     &              - dterm5ik*qrkxr - dterm6ik*(qkirxr-qrrx)
               ttmk(2) = lambdaik(3)*diky + dterm2ik*dkry - 
     &              dterm3ik*(dqiqky+diqkyr)
     &              - dterm5ik*qrkyr - dterm6ik*(qkiryr-qrry)
               ttmk(3) = lambdaik(3)*dikz + dterm2ik*dkrz - 
     &              dterm3ik*(dqiqkz+diqkzr)
     &              - dterm5ik*qrkzr - dterm6ik*(qkirzr-qrrz)
c
               ttmi(1) = ttmi(1)*rr1*oik*mscale(kk)
               ttmi(2) = ttmi(2)*rr1*oik*mscale(kk)
               ttmi(3) = ttmi(3)*rr1*oik*mscale(kk)
               ttmk(1) = ttmk(1)*rr1*oik*mscale(kk)
               ttmk(2) = ttmk(2)*rr1*oik*mscale(kk)
               ttmk(3) = ttmk(3)*rr1*oik*mscale(kk)
c
c     now take chain rule terms of 1/r
c
               frcx = frcx + evv*rr3*xr
               frcy = frcy + evv*rr3*yr
               frcz = frcz + evv*rr3*zr
c
c     multiply by overlap prefactor
c
               frcx = oik*mscale(kk)*frcx
               frcy = oik*mscale(kk)*frcy
               frcz = oik*mscale(kk)*frcz
c
c     force and torque components scaled by group membership
c
               if (use_group) then
                  frcx = fgrp * frcx
                  frcy = fgrp * frcy
                  frcz = fgrp * frcz
                  do j = 1, 3
                     ttmi(j) = fgrp * ttmi(j)
                     ttmk(j) = fgrp * ttmk(j)
                  end do
               end if
c
c     increment force-based gradient and torque on first site
c
               depr(1,ii) = depr(1,ii) + frcx
               depr(2,ii) = depr(2,ii) + frcy
               depr(3,ii) = depr(3,ii) + frcz
               tepr(1,i) = tepr(1,i) + ttmi(1)
               tepr(2,i) = tepr(2,i) + ttmi(2)
               tepr(3,i) = tepr(3,i) + ttmi(3)
c
c     increment force-based gradient and torque on second site
c
               depr(1,kk) = depr(1,kk) - frcx
               depr(2,kk) = depr(2,kk) - frcy
               depr(3,kk) = depr(3,kk) - frcz
               tepr(1,k) = tepr(1,k) + ttmk(1)
               tepr(2,k) = tepr(2,k) + ttmk(2)
               tepr(3,k) = tepr(3,k) + ttmk(3)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     increment the virial due to pairwise Cartesian forces
c
               vxx = -xr * frcx
               vxy = -yr * frcx
               vxz = -zr * frcx
               vyy = -yr * frcy
               vyz = -zr * frcy
               vzz = -zr * frcz
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vxy
               vir(3,1) = vir(3,1) + vxz
               vir(1,2) = vir(1,2) + vxy
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vyz
               vir(1,3) = vir(1,3) + vxz
               vir(2,3) = vir(2,3) + vyz
               vir(3,3) = vir(3,3) + vzz
            end if
   10       continue
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
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
         end do
      end do
c
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (use_replica) then
c
c     calculate interaction with other unit cells
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
         do k = i, npole
            kk = ipole(k)
            kz = zaxis(k)
            kx = xaxis(k)
            ky = yaxis(k)
            usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
            if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
            if (.not. proceed)  goto 20
            do jcell = 1, ncell
            xr = x(kk) - xi
            yr = y(kk) - yi
            zr = z(kk) - zi
            call imager (xr,yr,zr,jcell)
            r2 = xr*xr + yr*yr + zr*zr
            if (.not. (use_polymer .and. r2.le.polycut2)) then
               mscale(kk) = 1.0d0
            end if
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
               rr11 = 9.0d0 * rr9 / r2
c
c     intermediates involving moments and distance separation
c
               dikx = diy*dkz - diz*dky
               diky = diz*dkx - dix*dkz
               dikz = dix*dky - diy*dkx
               dirx = diy*zr - diz*yr
               diry = diz*xr - dix*zr
               dirz = dix*yr - diy*xr
               dkrx = dky*zr - dkz*yr
               dkry = dkz*xr - dkx*zr
               dkrz = dkx*yr - dky*xr
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
               qrixr = qriz*yr - qriy*zr
               qriyr = qrix*zr - qriz*xr
               qrizr = qriy*xr - qrix*yr
               qrkxr = qrkz*yr - qrky*zr
               qrkyr = qrkx*zr - qrkz*xr
               qrkzr = qrky*xr - qrkx*yr
               qrrx = qrky*qriz - qrkz*qriy
               qrry = qrkz*qrix - qrkx*qriz
               qrrz = qrkx*qriy - qrky*qrix
               qikrx = qixx*qrkx + qixy*qrky + qixz*qrkz
               qikry = qixy*qrkx + qiyy*qrky + qiyz*qrkz
               qikrz = qixz*qrkx + qiyz*qrky + qizz*qrkz
               qkirx = qkxx*qrix + qkxy*qriy + qkxz*qriz
               qkiry = qkxy*qrix + qkyy*qriy + qkyz*qriz
               qkirz = qkxz*qrix + qkyz*qriy + qkzz*qriz
               qikrxr = qikrz*yr - qikry*zr
               qikryr = qikrx*zr - qikrz*xr
               qikrzr = qikry*xr - qikrx*yr
               qkirxr = qkirz*yr - qkiry*zr
               qkiryr = qkirx*zr - qkirz*xr
               qkirzr = qkiry*xr - qkirx*yr
               diqkx = dix*qkxx + diy*qkxy + diz*qkxz
               diqky = dix*qkxy + diy*qkyy + diz*qkyz
               diqkz = dix*qkxz + diy*qkyz + diz*qkzz
               dkqix = dkx*qixx + dky*qixy + dkz*qixz
               dkqiy = dkx*qixy + dky*qiyy + dkz*qiyz
               dkqiz = dkx*qixz + dky*qiyz + dkz*qizz
               diqrk = dix*qrkx + diy*qrky + diz*qrkz
               dkqri = dkx*qrix + dky*qriy + dkz*qriz
               diqkxr = diqkz*yr - diqky*zr
               diqkyr = diqkx*zr - diqkz*xr
               diqkzr = diqky*xr - diqkx*yr
               dkqixr = dkqiz*yr - dkqiy*zr
               dkqiyr = dkqix*zr - dkqiz*xr
               dkqizr = dkqiy*xr - dkqix*yr
               dqiqkx = diy*qrkz - diz*qrky + dky*qriz - dkz*qriy
     &                     - 2.0d0*(qixy*qkxz+qiyy*qkyz+qiyz*qkzz
     &                             -qixz*qkxy-qiyz*qkyy-qizz*qkyz)
               dqiqky = diz*qrkx - dix*qrkz + dkz*qrix - dkx*qriz
     &                     - 2.0d0*(qixz*qkxx+qiyz*qkxy+qizz*qkxz
     &                             -qixx*qkxz-qixy*qkyz-qixz*qkzz)
               dqiqkz = dix*qrky - diy*qrkx + dkx*qriy - dky*qrix
     &                     - 2.0d0*(qixx*qkxy+qixy*qkyy+qixz*qkyz
     &                             -qixy*qkxx-qiyy*qkxy-qiyz*qkxz)
c
c     calculate intermediate terms for multipole energy
c
               term1 = ci*ck
               term2 = ck*dri - ci*drk + dik
               term3 = ci*qrrk + ck*qrri - dri*drk
     &                    + 2.0d0*(dkqri-diqrk+qik)
               term4 = dri*qrrk - drk*qrri - 4.0d0*qrrik
               term5 = qrri*qrrk
c
c     compute the energy contribution for this interaction
c
               e = term1*rr1 + term2*rr3 + term3*rr5
     &                + term4*rr7 + term5*rr9
               if (ii .eq. kk)  e = 0.5d0 * e
               em = em + e
               einter = einter + e
c
c     calculate intermediate terms for force and torque
c
               de = term1*rr3 + term2*rr5 + term3*rr7
     &                 + term4*rr9 + term5*rr11
               term1 = -ck*rr3 + drk*rr5 - qrrk*rr7
               term2 = ci*rr3 + dri*rr5 + qrri*rr7
               term3 = 2.0d0 * rr5
               term4 = 2.0d0 * (-ck*rr5+drk*rr7-qrrk*rr9)
               term5 = 2.0d0 * (-ci*rr5-dri*rr7-qrri*rr9)
               term6 = 4.0d0 * rr7
c
c     compute the force components for this interaction
c
               frcx = de*xr + term1*dix + term2*dkx
     &                   + term3*(diqkx-dkqix) + term4*qrix
     &                   + term5*qrkx + term6*(qikrx+qkirx)
               frcy = de*yr + term1*diy + term2*dky
     &                   + term3*(diqky-dkqiy) + term4*qriy
     &                   + term5*qrky + term6*(qikry+qkiry)
               frcz = de*zr + term1*diz + term2*dkz
     &                   + term3*(diqkz-dkqiz) + term4*qriz
     &                   + term5*qrkz + term6*(qikrz+qkirz)
c
c     compute the torque components for this interaction
c
               ttmi(1) = -rr3*dikx + term1*dirx + term3*(dqiqkx+dkqixr)
     &                      - term4*qrixr - term6*(qikrxr+qrrx)
               ttmi(2) = -rr3*diky + term1*diry + term3*(dqiqky+dkqiyr)
     &                      - term4*qriyr - term6*(qikryr+qrry)
               ttmi(3) = -rr3*dikz + term1*dirz + term3*(dqiqkz+dkqizr)
     &                      - term4*qrizr - term6*(qikrzr+qrrz)
               ttmk(1) = rr3*dikx + term2*dkrx - term3*(dqiqkx+diqkxr)
     &                      - term5*qrkxr - term6*(qkirxr-qrrx)
               ttmk(2) = rr3*diky + term2*dkry - term3*(dqiqky+diqkyr)
     &                      - term5*qrkyr - term6*(qkiryr-qrry)
               ttmk(3) = rr3*dikz + term2*dkrz - term3*(dqiqkz+diqkzr)
     &                      - term5*qrkzr - term6*(qkirzr-qrrz)
c
c     force and torque scaled for self-interactions and groups
c
               if (ii .eq. kk) then
                  frcx = 0.5d0 * frcx
                  frcy = 0.5d0 * frcy
                  frcz = 0.5d0 * frcz
                  do j = 1, 3
                     ttmi(j) = 0.5d0 * ttmi(j)
                     ttmk(j) = 0.5d0 * ttmk(j)
                  end do
               end if
               if (use_group) then
                  frcx = fgrp * frcx
                  frcy = fgrp * frcy
                  frcz = fgrp * frcz
                  do j = 1, 3
                     ttmi(j) = fgrp * ttmi(j)
                     ttmk(j) = fgrp * ttmk(j)
                  end do
               end if
c
c     increment force-based gradient and torque on first site
c
               dem(1,ii) = dem(1,ii) + frcx
               dem(2,ii) = dem(2,ii) + frcy
               dem(3,ii) = dem(3,ii) + frcz
               tem(1,i) = tem(1,i) + ttmi(1)
               tem(2,i) = tem(2,i) + ttmi(2)
               tem(3,i) = tem(3,i) + ttmi(3)
c
c     increment force-based gradient and torque on second site
c
               dem(1,kk) = dem(1,kk) - frcx
               dem(2,kk) = dem(2,kk) - frcy
               dem(3,kk) = dem(3,kk) - frcz
               tem(1,k) = tem(1,k) + ttmk(1)
               tem(2,k) = tem(2,k) + ttmk(2)
               tem(3,k) = tem(3,k) + ttmk(3)
c
c     increment the virial due to pairwise Cartesian forces
c
               vxx = -xr * frcx
               vxy = -yr * frcx
               vxz = -zr * frcx
               vyy = -yr * frcy
               vyz = -zr * frcy
               vzz = -zr * frcz
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vxy
               vir(3,1) = vir(3,1) + vxz
               vir(1,2) = vir(1,2) + vxy
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vyz
               vir(1,3) = vir(1,3) + vxz
               vir(2,3) = vir(2,3) + vyz
               vir(3,3) = vir(3,3) + vzz
            end if
            end do
   20       continue
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
c     resolve site torques then increment forces and virial
c
      do i = 1, npole
c
c     resolve pauli repulsion torques
c
         call torque (i,tepr(1,i),fix,fiy,fiz,depr)
c
         call torque (i,tem(1,i),fix,fiy,fiz,dem)
         ii = ipole(i)
         iaz = zaxis(i)
         iax = xaxis(i)
         iay = yaxis(i)
         if (iaz .eq. 0)  iaz = ii
         if (iax .eq. 0)  iax = ii
         if (iay .eq. 0)  iay = ii
         xiz = x(iaz) - x(ii)
         yiz = y(iaz) - y(ii)
         ziz = z(iaz) - z(ii)
         xix = x(iax) - x(ii)
         yix = y(iax) - y(ii)
         zix = z(iax) - z(ii)
         xiy = x(iay) - x(ii)
         yiy = y(iay) - y(ii)
         ziy = z(iay) - z(ii)
         vxx = xix*fix(1) + xiy*fiy(1) + xiz*fiz(1)
         vxy = 0.5d0 * (yix*fix(1) + yiy*fiy(1) + yiz*fiz(1)
     &                    + xix*fix(2) + xiy*fiy(2) + xiz*fiz(2))
         vxz = 0.5d0 * (zix*fix(1) + ziy*fiy(1) + ziz*fiz(1)
     &                    + xix*fix(3) + xiy*fiy(3) + xiz*fiz(3)) 
         vyy = yix*fix(2) + yiy*fiy(2) + yiz*fiz(2)
         vyz = 0.5d0 * (zix*fix(2) + ziy*fiy(2) + ziz*fiz(2)
     &                    + yix*fix(3) + yiy*fiy(3) + yiz*fiz(3))
         vzz = zix*fix(3) + ziy*fiy(3) + ziz*fiz(3)
         vir(1,1) = vir(1,1) + vxx
         vir(2,1) = vir(2,1) + vxy
         vir(3,1) = vir(3,1) + vxz
         vir(1,2) = vir(1,2) + vxy
         vir(2,2) = vir(2,2) + vyy
         vir(3,2) = vir(3,2) + vyz
         vir(1,3) = vir(1,3) + vxz
         vir(2,3) = vir(2,3) + vyz
         vir(3,3) = vir(3,3) + vzz
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      deallocate (tem)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ejosh1b  --  neighbor list multipole derivs  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ejosh1b" calculates the multipole energy and derivatives
c     with respect to Cartesian coordinates using a neighbor list
c
c
      subroutine ejosh1b
      use sizes
      use atoms
      use bound
      use chgpot
      use couple
      use deriv
      use energi
      use group
      use inter
      use molcul
      use mplpot
      use mpole
      use neigh
      use shunt
      use usage
      use virial
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      integer ix,iy,iz
      integer kx,ky,kz
      integer iax,iay,iaz
      real*8 e,de,f,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dikx,diky,dikz
      real*8 dirx,diry,dirz
      real*8 dkrx,dkry,dkrz
      real*8 qrix,qriy,qriz
      real*8 qrkx,qrky,qrkz
      real*8 qrixr,qriyr,qrizr
      real*8 qrkxr,qrkyr,qrkzr
      real*8 qrrx,qrry,qrrz
      real*8 qikrx,qikry,qikrz
      real*8 qkirx,qkiry,qkirz
      real*8 qikrxr,qikryr,qikrzr
      real*8 qkirxr,qkiryr,qkirzr
      real*8 diqkx,diqky,diqkz
      real*8 dkqix,dkqiy,dkqiz
      real*8 diqkxr,diqkyr,diqkzr
      real*8 dkqixr,dkqiyr,dkqizr
      real*8 dqiqkx,dqiqky,dqiqkz
      real*8 dri,drk,qrri,qrrk
      real*8 diqrk,dkqri
      real*8 dik,qik,qrrik
      real*8 term1,term2,term3
      real*8 term4,term5,term6
      real*8 frcx,frcy,frcz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 ttmi(3),ttmk(3)
      real*8 fix(3),fiy(3),fiz(3)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: tem(:,:)
      logical proceed,usei,usek
      character*6 mode
c
c
c     zero out the atomic multipole energy and derivatives
c
      em = 0.0d0
      do i = 1, n
         do j = 1, 3
            dem(j,i) = 0.0d0
         end do
      end do
      if (npole .eq. 0)  return
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (tem(3,n))
c
c     initialize connected atom scaling and torque arrays
c
      do i = 1, n
         mscale(i) = 1.0d0
         do j = 1, 3
            tem(j,i) = 0.0d0
         end do
      end do
c
c     set conversion factor, cutoff and scaling coefficients
c
      f = electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(npole,ipole,x,y,z,xaxis,yaxis,zaxis,rpole,use,n12,
!$OMP& i12,n13,i13,n14,i14,n15,i15,m2scale,m3scale,m4scale,m5scale,
!$OMP& nelst,elst,use_group,use_intra,use_bounds,off2,f,molcule)
!$OMP& firstprivate(mscale) shared (em,einter,dem,tem,vir)
!$OMP DO reduction(+:em,einter,dem,tem,vir) schedule(guided)
c
c     compute the multipole interaction energy and gradient
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
            if (.not. proceed)  goto 10
            xr = x(kk) - xi
            yr = y(kk) - yi
            zr = z(kk) - zi
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
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
               rr11 = 9.0d0 * rr9 / r2
c
c     intermediates involving moments and distance separation
c
               dikx = diy*dkz - diz*dky
               diky = diz*dkx - dix*dkz
               dikz = dix*dky - diy*dkx
               dirx = diy*zr - diz*yr
               diry = diz*xr - dix*zr
               dirz = dix*yr - diy*xr
               dkrx = dky*zr - dkz*yr
               dkry = dkz*xr - dkx*zr
               dkrz = dkx*yr - dky*xr
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
               qrixr = qriz*yr - qriy*zr
               qriyr = qrix*zr - qriz*xr
               qrizr = qriy*xr - qrix*yr
               qrkxr = qrkz*yr - qrky*zr
               qrkyr = qrkx*zr - qrkz*xr
               qrkzr = qrky*xr - qrkx*yr
               qrrx = qrky*qriz - qrkz*qriy
               qrry = qrkz*qrix - qrkx*qriz
               qrrz = qrkx*qriy - qrky*qrix
               qikrx = qixx*qrkx + qixy*qrky + qixz*qrkz
               qikry = qixy*qrkx + qiyy*qrky + qiyz*qrkz
               qikrz = qixz*qrkx + qiyz*qrky + qizz*qrkz
               qkirx = qkxx*qrix + qkxy*qriy + qkxz*qriz
               qkiry = qkxy*qrix + qkyy*qriy + qkyz*qriz
               qkirz = qkxz*qrix + qkyz*qriy + qkzz*qriz
               qikrxr = qikrz*yr - qikry*zr
               qikryr = qikrx*zr - qikrz*xr
               qikrzr = qikry*xr - qikrx*yr
               qkirxr = qkirz*yr - qkiry*zr
               qkiryr = qkirx*zr - qkirz*xr
               qkirzr = qkiry*xr - qkirx*yr
               diqkx = dix*qkxx + diy*qkxy + diz*qkxz
               diqky = dix*qkxy + diy*qkyy + diz*qkyz
               diqkz = dix*qkxz + diy*qkyz + diz*qkzz
               dkqix = dkx*qixx + dky*qixy + dkz*qixz
               dkqiy = dkx*qixy + dky*qiyy + dkz*qiyz
               dkqiz = dkx*qixz + dky*qiyz + dkz*qizz
               diqrk = dix*qrkx + diy*qrky + diz*qrkz
               dkqri = dkx*qrix + dky*qriy + dkz*qriz
               diqkxr = diqkz*yr - diqky*zr
               diqkyr = diqkx*zr - diqkz*xr
               diqkzr = diqky*xr - diqkx*yr
               dkqixr = dkqiz*yr - dkqiy*zr
               dkqiyr = dkqix*zr - dkqiz*xr
               dkqizr = dkqiy*xr - dkqix*yr
               dqiqkx = diy*qrkz - diz*qrky + dky*qriz - dkz*qriy
     &                     - 2.0d0*(qixy*qkxz+qiyy*qkyz+qiyz*qkzz
     &                             -qixz*qkxy-qiyz*qkyy-qizz*qkyz)
               dqiqky = diz*qrkx - dix*qrkz + dkz*qrix - dkx*qriz
     &                     - 2.0d0*(qixz*qkxx+qiyz*qkxy+qizz*qkxz
     &                             -qixx*qkxz-qixy*qkyz-qixz*qkzz)
               dqiqkz = dix*qrky - diy*qrkx + dkx*qriy - dky*qrix
     &                     - 2.0d0*(qixx*qkxy+qixy*qkyy+qixz*qkyz
     &                             -qixy*qkxx-qiyy*qkxy-qiyz*qkxz)
c
c     calculate intermediate terms for multipole energy
c
               term1 = ci*ck
               term2 = ck*dri - ci*drk + dik
               term3 = ci*qrrk + ck*qrri - dri*drk
     &                    + 2.0d0*(dkqri-diqrk+qik)
               term4 = dri*qrrk - drk*qrri - 4.0d0*qrrik
               term5 = qrri*qrrk
c
c     compute the energy contribution for this interaction
c
               e = term1*rr1 + term2*rr3 + term3*rr5
     &                + term4*rr7 + term5*rr9
               if (use_group)  e = e * fgrp
               em = em + e
               if (molcule(ii) .ne. molcule(kk))
     &            einter = einter + e
c
c     calculate intermediate terms for force and torque
c
               de = term1*rr3 + term2*rr5 + term3*rr7
     &                 + term4*rr9 + term5*rr11
               term1 = -ck*rr3 + drk*rr5 - qrrk*rr7
               term2 = ci*rr3 + dri*rr5 + qrri*rr7
               term3 = 2.0d0 * rr5
               term4 = 2.0d0 * (-ck*rr5+drk*rr7-qrrk*rr9)
               term5 = 2.0d0 * (-ci*rr5-dri*rr7-qrri*rr9)
               term6 = 4.0d0 * rr7
c
c     compute the force components for this interaction
c
               frcx = de*xr + term1*dix + term2*dkx
     &                   + term3*(diqkx-dkqix) + term4*qrix
     &                   + term5*qrkx + term6*(qikrx+qkirx)
               frcy = de*yr + term1*diy + term2*dky
     &                   + term3*(diqky-dkqiy) + term4*qriy
     &                   + term5*qrky + term6*(qikry+qkiry)
               frcz = de*zr + term1*diz + term2*dkz
     &                   + term3*(diqkz-dkqiz) + term4*qriz
     &                   + term5*qrkz + term6*(qikrz+qkirz)
c
c     compute the torque components for this interaction
c
               ttmi(1) = -rr3*dikx + term1*dirx + term3*(dqiqkx+dkqixr)
     &                      - term4*qrixr - term6*(qikrxr+qrrx)
               ttmi(2) = -rr3*diky + term1*diry + term3*(dqiqky+dkqiyr)
     &                      - term4*qriyr - term6*(qikryr+qrry)
               ttmi(3) = -rr3*dikz + term1*dirz + term3*(dqiqkz+dkqizr)
     &                      - term4*qrizr - term6*(qikrzr+qrrz)
               ttmk(1) = rr3*dikx + term2*dkrx - term3*(dqiqkx+diqkxr)
     &                      - term5*qrkxr - term6*(qkirxr-qrrx)
               ttmk(2) = rr3*diky + term2*dkry - term3*(dqiqky+diqkyr)
     &                      - term5*qrkyr - term6*(qkiryr-qrry)
               ttmk(3) = rr3*dikz + term2*dkrz - term3*(dqiqkz+diqkzr)
     &                      - term5*qrkzr - term6*(qkirzr-qrrz)
c
c     increment force-based gradient and torque on first site
c
               dem(1,ii) = dem(1,ii) + frcx
               dem(2,ii) = dem(2,ii) + frcy
               dem(3,ii) = dem(3,ii) + frcz
               tem(1,i) = tem(1,i) + ttmi(1)
               tem(2,i) = tem(2,i) + ttmi(2)
               tem(3,i) = tem(3,i) + ttmi(3)
c
c     increment force-based gradient and torque on second site
c
               dem(1,kk) = dem(1,kk) - frcx
               dem(2,kk) = dem(2,kk) - frcy
               dem(3,kk) = dem(3,kk) - frcz
               tem(1,k) = tem(1,k) + ttmk(1)
               tem(2,k) = tem(2,k) + ttmk(2)
               tem(3,k) = tem(3,k) + ttmk(3)
c
c     increment the virial due to pairwise Cartesian forces
c
               vxx = -xr * frcx
               vxy = -yr * frcx
               vxz = -zr * frcx
               vyy = -yr * frcy
               vyz = -zr * frcy
               vzz = -zr * frcz
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vxy
               vir(3,1) = vir(3,1) + vxz
               vir(1,2) = vir(1,2) + vxy
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vyz
               vir(1,3) = vir(1,3) + vxz
               vir(2,3) = vir(2,3) + vyz
               vir(3,3) = vir(3,3) + vzz
            end if
   10       continue
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
!$OMP DO reduction(+:dem,vir) schedule(guided)
c
c     resolve site torques then increment forces and virial
c
      do i = 1, npole
         call torque (i,tem(1,i),fix,fiy,fiz,dem)
         ii = ipole(i)
         iaz = zaxis(i)
         iax = xaxis(i)
         iay = yaxis(i)
         if (iaz .eq. 0)  iaz = ii
         if (iax .eq. 0)  iax = ii
         if (iay .eq. 0)  iay = ii
         xiz = x(iaz) - x(ii)
         yiz = y(iaz) - y(ii)
         ziz = z(iaz) - z(ii)
         xix = x(iax) - x(ii)
         yix = y(iax) - y(ii)
         zix = z(iax) - z(ii)
         xiy = x(iay) - x(ii)
         yiy = y(iay) - y(ii)
         ziy = z(iay) - z(ii)
         vxx = xix*fix(1) + xiy*fiy(1) + xiz*fiz(1)
         vxy = 0.5d0 * (yix*fix(1) + yiy*fiy(1) + yiz*fiz(1)
     &                    + xix*fix(2) + xiy*fiy(2) + xiz*fiz(2))
         vxz = 0.5d0 * (zix*fix(1) + ziy*fiy(1) + ziz*fiz(1)
     &                    + xix*fix(3) + xiy*fiy(3) + xiz*fiz(3)) 
         vyy = yix*fix(2) + yiy*fiy(2) + yiz*fiz(2)
         vyz = 0.5d0 * (zix*fix(2) + ziy*fiy(2) + ziz*fiz(2)
     &                    + yix*fix(3) + yiy*fiy(3) + yiz*fiz(3))
         vzz = zix*fix(3) + ziy*fiy(3) + ziz*fiz(3)
         vir(1,1) = vir(1,1) + vxx
         vir(2,1) = vir(2,1) + vxy
         vir(3,1) = vir(3,1) + vxz
         vir(1,2) = vir(1,2) + vxy
         vir(2,2) = vir(2,2) + vyy
         vir(3,2) = vir(3,2) + vyz
         vir(1,3) = vir(1,3) + vxz
         vir(2,3) = vir(2,3) + vyz
         vir(3,3) = vir(3,3) + vzz
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
      deallocate (tem)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ejosh1c  --  Ewald multipole derivs via loop  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ejosh1c" calculates the multipole energy and derivatives
c     with respect to Cartesian coordinates using particle mesh
c     Ewald summation and a double loop
c
c
      subroutine ejosh1c
      use sizes
      use atoms
      use boxes
      use chgpot
      use deriv
      use energi
      use ewald
      use math
      use mpole
      use virial
      implicit none
      integer i,j,ii
      real*8 e,f
      real*8 term,fterm
      real*8 cii,dii,qii
      real*8 xd,yd,zd
      real*8 xq,yq,zq
      real*8 xv,yv,zv,vterm
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 xdfield,ydfield
      real*8 zdfield
      real*8 trq(3),frcx(3)
      real*8 frcy(3),frcz(3)
c
c
c     zero out the atomic multipole energy and derivatives
c
      em = 0.0d0
      do i = 1, n
         do j = 1, 3
            dem(j,i) = 0.0d0
         end do
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
      call ejoshreal1c
c
c     compute the reciprocal space part of the Ewald summation
c
      call ejoshrecip1
c
c     compute the Ewald self-energy term over all the atoms
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
            xd = xd + rpole(2,i) + rpole(1,i)*x(ii)
            yd = yd + rpole(3,i) + rpole(1,i)*y(ii)
            zd = zd + rpole(4,i) + rpole(1,i)*z(ii)
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         em = em + term*(xd*xd+yd*yd+zd*zd)
         do i = 1, npole
            ii = ipole(i)
            dem(1,ii) = dem(1,ii) + 2.0d0*term*rpole(1,i)*xd
            dem(2,ii) = dem(2,ii) + 2.0d0*term*rpole(1,i)*yd
            dem(3,ii) = dem(3,ii) + 2.0d0*term*rpole(1,i)*zd
         end do
         xdfield = -2.0d0 * term * xd
         ydfield = -2.0d0 * term * yd
         zdfield = -2.0d0 * term * zd
         do i = 1, npole
            trq(1) = rpole(3,i)*zdfield - rpole(4,i)*ydfield
            trq(2) = rpole(4,i)*xdfield - rpole(2,i)*zdfield
            trq(3) = rpole(2,i)*ydfield - rpole(3,i)*xdfield
            call torque (i,trq,frcx,frcy,frcz,dem)
         end do
c
c     boundary correction to virial due to overall cell dipole
c
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         xq = 0.0d0
         yq = 0.0d0
         zq = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            xd = xd + rpole(2,i)
            yd = yd + rpole(3,i)
            zd = zd + rpole(4,i)
            xq = xq + rpole(1,i)*x(ii)
            yq = yq + rpole(1,i)*y(ii)
            zq = zq + rpole(1,i)*z(ii)
         end do
         xv = xd * xq
         yv = yd * yq
         zv = zd * zq
         vterm = term * (xd*xd + yd*yd + zd*zd + 2.0d0*(xv+yv+zv)
     &                      + xq*xq + yq*yq + zq*zq)
         vir(1,1) = vir(1,1) + 2.0d0*term*(xq*xq+xv) + vterm
         vir(2,1) = vir(2,1) + 2.0d0*term*(xq*yq+xv)
         vir(3,1) = vir(3,1) + 2.0d0*term*(xq*zq+xv)
         vir(1,2) = vir(1,2) + 2.0d0*term*(yq*xq+yv)
         vir(2,2) = vir(2,2) + 2.0d0*term*(yq*yq+yv) + vterm
         vir(3,2) = vir(3,2) + 2.0d0*term*(yq*zq+yv)
         vir(1,3) = vir(1,3) + 2.0d0*term*(zq*xq+zv)
         vir(2,3) = vir(2,3) + 2.0d0*term*(zq*yq+zv)
         vir(3,3) = vir(3,3) + 2.0d0*term*(zq*zq+zv) + vterm
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ejoshreal1c  --  Ewald real space derivs via loop  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ejoshreal1c" evaluates the real space portion of the Ewald
c     summation energy and gradient due to multipole interactions
c     via a double loop
c
c
      subroutine ejoshreal1c
      use sizes
      use atoms
      use bound
      use cell
      use chgpot
      use couple
      use deriv
      use energi
      use ewald
      use inter
      use math
      use molcul
      use mplpot
      use mpole
      use shunt
      use virial
      implicit none
      integer i,j,k
      integer ii,kk,jcell
      integer iax,iay,iaz
      real*8 e,efull,de,f
      real*8 bfac,erfc
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 scalekk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dikx,diky,dikz
      real*8 dirx,diry,dirz
      real*8 dkrx,dkry,dkrz
      real*8 qrix,qriy,qriz
      real*8 qrkx,qrky,qrkz
      real*8 qrixr,qriyr,qrizr
      real*8 qrkxr,qrkyr,qrkzr
      real*8 qrrx,qrry,qrrz
      real*8 qikrx,qikry,qikrz
      real*8 qkirx,qkiry,qkirz
      real*8 qikrxr,qikryr,qikrzr
      real*8 qkirxr,qkiryr,qkirzr
      real*8 diqkx,diqky,diqkz
      real*8 dkqix,dkqiy,dkqiz
      real*8 diqkxr,diqkyr,diqkzr
      real*8 dkqixr,dkqiyr,dkqizr
      real*8 dqiqkx,dqiqky,dqiqkz
      real*8 dri,drk,qrri,qrrk
      real*8 diqrk,dkqri
      real*8 dik,qik,qrrik
      real*8 term1,term2,term3
      real*8 term4,term5,term6
      real*8 frcx,frcy,frcz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 ttmi(3),ttmk(3)
      real*8 fix(3),fiy(3),fiz(3)
      real*8 bn(0:5)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: tem(:,:)
      character*6 mode
      external erfc
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (tem(3,n))
c
c     initialize connected atom scaling and torque arrays
c
      do i = 1, n
         mscale(i) = 1.0d0
         do j = 1, 3
            tem(j,i) = 0.0d0
         end do
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
            r2 = xr*xr + yr*yr + zr*zr
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
               rr11 = 9.0d0 * rr9 / r2
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 5
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
               do j = 0, 5
                  bn(j) = f * bn(j)
               end do
c
c     intermediates involving moments and distance separation
c
               dikx = diy*dkz - diz*dky
               diky = diz*dkx - dix*dkz
               dikz = dix*dky - diy*dkx
               dirx = diy*zr - diz*yr
               diry = diz*xr - dix*zr
               dirz = dix*yr - diy*xr
               dkrx = dky*zr - dkz*yr
               dkry = dkz*xr - dkx*zr
               dkrz = dkx*yr - dky*xr
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
               qrixr = qriz*yr - qriy*zr
               qriyr = qrix*zr - qriz*xr
               qrizr = qriy*xr - qrix*yr
               qrkxr = qrkz*yr - qrky*zr
               qrkyr = qrkx*zr - qrkz*xr
               qrkzr = qrky*xr - qrkx*yr
               qrrx = qrky*qriz - qrkz*qriy
               qrry = qrkz*qrix - qrkx*qriz
               qrrz = qrkx*qriy - qrky*qrix
               qikrx = qixx*qrkx + qixy*qrky + qixz*qrkz
               qikry = qixy*qrkx + qiyy*qrky + qiyz*qrkz
               qikrz = qixz*qrkx + qiyz*qrky + qizz*qrkz
               qkirx = qkxx*qrix + qkxy*qriy + qkxz*qriz
               qkiry = qkxy*qrix + qkyy*qriy + qkyz*qriz
               qkirz = qkxz*qrix + qkyz*qriy + qkzz*qriz
               qikrxr = qikrz*yr - qikry*zr
               qikryr = qikrx*zr - qikrz*xr
               qikrzr = qikry*xr - qikrx*yr
               qkirxr = qkirz*yr - qkiry*zr
               qkiryr = qkirx*zr - qkirz*xr
               qkirzr = qkiry*xr - qkirx*yr
               diqkx = dix*qkxx + diy*qkxy + diz*qkxz
               diqky = dix*qkxy + diy*qkyy + diz*qkyz
               diqkz = dix*qkxz + diy*qkyz + diz*qkzz
               dkqix = dkx*qixx + dky*qixy + dkz*qixz
               dkqiy = dkx*qixy + dky*qiyy + dkz*qiyz
               dkqiz = dkx*qixz + dky*qiyz + dkz*qizz
               diqrk = dix*qrkx + diy*qrky + diz*qrkz
               dkqri = dkx*qrix + dky*qriy + dkz*qriz
               diqkxr = diqkz*yr - diqky*zr
               diqkyr = diqkx*zr - diqkz*xr
               diqkzr = diqky*xr - diqkx*yr
               dkqixr = dkqiz*yr - dkqiy*zr
               dkqiyr = dkqix*zr - dkqiz*xr
               dkqizr = dkqiy*xr - dkqix*yr
               dqiqkx = diy*qrkz - diz*qrky + dky*qriz - dkz*qriy
     &                     - 2.0d0*(qixy*qkxz+qiyy*qkyz+qiyz*qkzz
     &                             -qixz*qkxy-qiyz*qkyy-qizz*qkyz)
               dqiqky = diz*qrkx - dix*qrkz + dkz*qrix - dkx*qriz
     &                     - 2.0d0*(qixz*qkxx+qiyz*qkxy+qizz*qkxz
     &                             -qixx*qkxz-qixy*qkyz-qixz*qkzz)
               dqiqkz = dix*qrky - diy*qrkx + dkx*qriy - dky*qrix
     &                     - 2.0d0*(qixx*qkxy+qixy*qkyy+qixz*qkyz
     &                             -qixy*qkxx-qiyy*qkxy-qiyz*qkxz)
c
c     calculate intermediate terms for multipole energy
c
               term1 = ci*ck
               term2 = ck*dri - ci*drk + dik
               term3 = ci*qrrk + ck*qrri - dri*drk
     &                    + 2.0d0*(dkqri-diqrk+qik)
               term4 = dri*qrrk - drk*qrri - 4.0d0*qrrik
               term5 = qrri*qrrk
c
c     compute the full energy without any Ewald scaling
c
               efull = term1*rr1 + term2*rr3 + term3*rr5
     &                    + term4*rr7 + term5*rr9
               efull = efull * mscale(kk)
               if (molcule(ii) .ne. molcule(kk))
     &            einter = einter + efull
c
c     modify distances to account for Ewald and exclusions
c
               scalekk = 1.0d0 - mscale(kk)
               rr1 = bn(0) - scalekk*rr1
               rr3 = bn(1) - scalekk*rr3
               rr5 = bn(2) - scalekk*rr5
               rr7 = bn(3) - scalekk*rr7
               rr9 = bn(4) - scalekk*rr9
               rr11 = bn(5) - scalekk*rr11
c
c     compute the energy contribution for this interaction
c
               e = term1*rr1 + term2*rr3 + term3*rr5
     &                + term4*rr7 + term5*rr9
               em = em + e
c
c     calculate intermediate terms for force and torque
c
               de = term1*rr3 + term2*rr5 + term3*rr7
     &                 + term4*rr9 + term5*rr11
               term1 = -ck*rr3 + drk*rr5 - qrrk*rr7
               term2 = ci*rr3 + dri*rr5 + qrri*rr7
               term3 = 2.0d0 * rr5
               term4 = 2.0d0 * (-ck*rr5+drk*rr7-qrrk*rr9)
               term5 = 2.0d0 * (-ci*rr5-dri*rr7-qrri*rr9)
               term6 = 4.0d0 * rr7
c
c     compute the force components for this interaction
c
               frcx = de*xr + term1*dix + term2*dkx
     &                   + term3*(diqkx-dkqix) + term4*qrix
     &                   + term5*qrkx + term6*(qikrx+qkirx)
               frcy = de*yr + term1*diy + term2*dky
     &                   + term3*(diqky-dkqiy) + term4*qriy
     &                   + term5*qrky + term6*(qikry+qkiry)
               frcz = de*zr + term1*diz + term2*dkz
     &                   + term3*(diqkz-dkqiz) + term4*qriz
     &                   + term5*qrkz + term6*(qikrz+qkirz)
c
c     compute the torque components for this interaction
c
               ttmi(1) = -rr3*dikx + term1*dirx + term3*(dqiqkx+dkqixr)
     &                      - term4*qrixr - term6*(qikrxr+qrrx)
               ttmi(2) = -rr3*diky + term1*diry + term3*(dqiqky+dkqiyr)
     &                      - term4*qriyr - term6*(qikryr+qrry)
               ttmi(3) = -rr3*dikz + term1*dirz + term3*(dqiqkz+dkqizr)
     &                      - term4*qrizr - term6*(qikrzr+qrrz)
               ttmk(1) = rr3*dikx + term2*dkrx - term3*(dqiqkx+diqkxr)
     &                      - term5*qrkxr - term6*(qkirxr-qrrx)
               ttmk(2) = rr3*diky + term2*dkry - term3*(dqiqky+diqkyr)
     &                      - term5*qrkyr - term6*(qkiryr-qrry)
               ttmk(3) = rr3*dikz + term2*dkrz - term3*(dqiqkz+diqkzr)
     &                      - term5*qrkzr - term6*(qkirzr-qrrz)
c
c     increment force-based gradient and torque on first site
c
               dem(1,ii) = dem(1,ii) + frcx
               dem(2,ii) = dem(2,ii) + frcy
               dem(3,ii) = dem(3,ii) + frcz
               tem(1,i) = tem(1,i) + ttmi(1)
               tem(2,i) = tem(2,i) + ttmi(2)
               tem(3,i) = tem(3,i) + ttmi(3)
c
c     increment force-based gradient and torque on second site
c
               dem(1,kk) = dem(1,kk) - frcx
               dem(2,kk) = dem(2,kk) - frcy
               dem(3,kk) = dem(3,kk) - frcz
               tem(1,k) = tem(1,k) + ttmk(1)
               tem(2,k) = tem(2,k) + ttmk(2)
               tem(3,k) = tem(3,k) + ttmk(3)
c
c     increment the virial due to pairwise Cartesian forces
c
               vxx = -xr * frcx
               vxy = -yr * frcx
               vxz = -zr * frcx
               vyy = -yr * frcy
               vyz = -zr * frcy
               vzz = -zr * frcz
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vxy
               vir(3,1) = vir(3,1) + vxz
               vir(1,2) = vir(1,2) + vxy
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vyz
               vir(1,3) = vir(1,3) + vxz
               vir(2,3) = vir(2,3) + vyz
               vir(3,3) = vir(3,3) + vzz
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
c     calculate interaction with other unit cells
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
            r2 = xr*xr + yr*yr + zr*zr
            if (.not. (use_polymer .and. r2.le.polycut2)) then
               mscale(kk) = 1.0d0
            end if
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
               rr11 = 9.0d0 * rr9 / r2
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 5
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
               do j = 0, 5
                  bn(j) = f * bn(j)
               end do
c
c     intermediates involving moments and distance separation
c
               dikx = diy*dkz - diz*dky
               diky = diz*dkx - dix*dkz
               dikz = dix*dky - diy*dkx
               dirx = diy*zr - diz*yr
               diry = diz*xr - dix*zr
               dirz = dix*yr - diy*xr
               dkrx = dky*zr - dkz*yr
               dkry = dkz*xr - dkx*zr
               dkrz = dkx*yr - dky*xr
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
               qrixr = qriz*yr - qriy*zr
               qriyr = qrix*zr - qriz*xr
               qrizr = qriy*xr - qrix*yr
               qrkxr = qrkz*yr - qrky*zr
               qrkyr = qrkx*zr - qrkz*xr
               qrkzr = qrky*xr - qrkx*yr
               qrrx = qrky*qriz - qrkz*qriy
               qrry = qrkz*qrix - qrkx*qriz
               qrrz = qrkx*qriy - qrky*qrix
               qikrx = qixx*qrkx + qixy*qrky + qixz*qrkz
               qikry = qixy*qrkx + qiyy*qrky + qiyz*qrkz
               qikrz = qixz*qrkx + qiyz*qrky + qizz*qrkz
               qkirx = qkxx*qrix + qkxy*qriy + qkxz*qriz
               qkiry = qkxy*qrix + qkyy*qriy + qkyz*qriz
               qkirz = qkxz*qrix + qkyz*qriy + qkzz*qriz
               qikrxr = qikrz*yr - qikry*zr
               qikryr = qikrx*zr - qikrz*xr
               qikrzr = qikry*xr - qikrx*yr
               qkirxr = qkirz*yr - qkiry*zr
               qkiryr = qkirx*zr - qkirz*xr
               qkirzr = qkiry*xr - qkirx*yr
               diqkx = dix*qkxx + diy*qkxy + diz*qkxz
               diqky = dix*qkxy + diy*qkyy + diz*qkyz
               diqkz = dix*qkxz + diy*qkyz + diz*qkzz
               dkqix = dkx*qixx + dky*qixy + dkz*qixz
               dkqiy = dkx*qixy + dky*qiyy + dkz*qiyz
               dkqiz = dkx*qixz + dky*qiyz + dkz*qizz
               diqrk = dix*qrkx + diy*qrky + diz*qrkz
               dkqri = dkx*qrix + dky*qriy + dkz*qriz
               diqkxr = diqkz*yr - diqky*zr
               diqkyr = diqkx*zr - diqkz*xr
               diqkzr = diqky*xr - diqkx*yr
               dkqixr = dkqiz*yr - dkqiy*zr
               dkqiyr = dkqix*zr - dkqiz*xr
               dkqizr = dkqiy*xr - dkqix*yr
               dqiqkx = diy*qrkz - diz*qrky + dky*qriz - dkz*qriy
     &                     - 2.0d0*(qixy*qkxz+qiyy*qkyz+qiyz*qkzz
     &                             -qixz*qkxy-qiyz*qkyy-qizz*qkyz)
               dqiqky = diz*qrkx - dix*qrkz + dkz*qrix - dkx*qriz
     &                     - 2.0d0*(qixz*qkxx+qiyz*qkxy+qizz*qkxz
     &                             -qixx*qkxz-qixy*qkyz-qixz*qkzz)
               dqiqkz = dix*qrky - diy*qrkx + dkx*qriy - dky*qrix
     &                     - 2.0d0*(qixx*qkxy+qixy*qkyy+qixz*qkyz
     &                             -qixy*qkxx-qiyy*qkxy-qiyz*qkxz)
c
c     calculate intermediate terms for multipole energy
c
               term1 = ci*ck
               term2 = ck*dri - ci*drk + dik
               term3 = ci*qrrk + ck*qrri - dri*drk
     &                    + 2.0d0*(dkqri-diqrk+qik)
               term4 = dri*qrrk - drk*qrri - 4.0d0*qrrik
               term5 = qrri*qrrk
c
c     compute the full energy without any Ewald scaling
c
               efull = term1*rr1 + term2*rr3 + term3*rr5
     &                    + term4*rr7 + term5*rr9
               efull = efull * mscale(kk)
               if (ii .eq. kk)  efull = 0.5d0 * e
               einter = einter + efull
c
c     modify distances to account for Ewald and exclusions
c
               scalekk = 1.0d0 - mscale(kk)
               rr1 = bn(0) - scalekk*rr1
               rr3 = bn(1) - scalekk*rr3
               rr5 = bn(2) - scalekk*rr5
               rr7 = bn(3) - scalekk*rr7
               rr9 = bn(4) - scalekk*rr9
               rr11 = bn(5) - scalekk*rr11
c
c     compute the energy contribution for this interaction
c
               e = term1*rr1 + term2*rr3 + term3*rr5
     &                + term4*rr7 + term5*rr9
               if (ii .eq. kk)  e = 0.5d0 * e
               em = em + e
c
c     calculate intermediate terms for force and torque
c
               de = term1*rr3 + term2*rr5 + term3*rr7
     &                 + term4*rr9 + term5*rr11
               term1 = -ck*rr3 + drk*rr5 - qrrk*rr7
               term2 = ci*rr3 + dri*rr5 + qrri*rr7
               term3 = 2.0d0 * rr5
               term4 = 2.0d0 * (-ck*rr5+drk*rr7-qrrk*rr9)
               term5 = 2.0d0 * (-ci*rr5-dri*rr7-qrri*rr9)
               term6 = 4.0d0 * rr7
c
c     compute the force components for this interaction
c
               frcx = de*xr + term1*dix + term2*dkx
     &                   + term3*(diqkx-dkqix) + term4*qrix
     &                   + term5*qrkx + term6*(qikrx+qkirx)
               frcy = de*yr + term1*diy + term2*dky
     &                   + term3*(diqky-dkqiy) + term4*qriy
     &                   + term5*qrky + term6*(qikry+qkiry)
               frcz = de*zr + term1*diz + term2*dkz
     &                   + term3*(diqkz-dkqiz) + term4*qriz
     &                   + term5*qrkz + term6*(qikrz+qkirz)
c
c     compute the torque components for this interaction
c
               ttmi(1) = -rr3*dikx + term1*dirx + term3*(dqiqkx+dkqixr)
     &                      - term4*qrixr - term6*(qikrxr+qrrx)
               ttmi(2) = -rr3*diky + term1*diry + term3*(dqiqky+dkqiyr)
     &                      - term4*qriyr - term6*(qikryr+qrry)
               ttmi(3) = -rr3*dikz + term1*dirz + term3*(dqiqkz+dkqizr)
     &                      - term4*qrizr - term6*(qikrzr+qrrz)
               ttmk(1) = rr3*dikx + term2*dkrx - term3*(dqiqkx+diqkxr)
     &                      - term5*qrkxr - term6*(qkirxr-qrrx)
               ttmk(2) = rr3*diky + term2*dkry - term3*(dqiqky+diqkyr)
     &                      - term5*qrkyr - term6*(qkiryr-qrry)
               ttmk(3) = rr3*dikz + term2*dkrz - term3*(dqiqkz+diqkzr)
     &                      - term5*qrkzr - term6*(qkirzr-qrrz)
c
c     force and torque components scaled for self-interactions
c
               if (ii .eq. kk) then
                  frcx = 0.5d0 * frcx
                  frcy = 0.5d0 * frcy
                  frcz = 0.5d0 * frcz
                  do j = 1, 3
                     ttmi(j) = 0.5d0 * ttmi(j)
                     ttmk(j) = 0.5d0 * ttmk(j)
                  end do
               end if
c
c     increment force-based gradient and torque on first site
c
               dem(1,ii) = dem(1,ii) + frcx
               dem(2,ii) = dem(2,ii) + frcy
               dem(3,ii) = dem(3,ii) + frcz
               tem(1,i) = tem(1,i) + ttmi(1)
               tem(2,i) = tem(2,i) + ttmi(2)
               tem(3,i) = tem(3,i) + ttmi(3)
c
c     increment force-based gradient and torque on second site
c
               dem(1,kk) = dem(1,kk) - frcx
               dem(2,kk) = dem(2,kk) - frcy
               dem(3,kk) = dem(3,kk) - frcz
               tem(1,k) = tem(1,k) + ttmk(1)
               tem(2,k) = tem(2,k) + ttmk(2)
               tem(3,k) = tem(3,k) + ttmk(3)
c
c     increment the virial due to pairwise Cartesian forces
c
               vxx = -xr * frcx
               vxy = -yr * frcx
               vxz = -zr * frcx
               vyy = -yr * frcy
               vyz = -zr * frcy
               vzz = -zr * frcz
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vxy
               vir(3,1) = vir(3,1) + vxz
               vir(1,2) = vir(1,2) + vxy
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vyz
               vir(1,3) = vir(1,3) + vxz
               vir(2,3) = vir(2,3) + vyz
               vir(3,3) = vir(3,3) + vzz
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
c     resolve site torques then increment forces and virial
c
      do i = 1, npole
         call torque (i,tem(1,i),fix,fiy,fiz,dem)
         ii = ipole(i)
         iaz = zaxis(i)
         iax = xaxis(i)
         iay = yaxis(i)
         if (iaz .eq. 0)  iaz = ii
         if (iax .eq. 0)  iax = ii
         if (iay .eq. 0)  iay = ii
         xiz = x(iaz) - x(ii)
         yiz = y(iaz) - y(ii)
         ziz = z(iaz) - z(ii)
         xix = x(iax) - x(ii)
         yix = y(iax) - y(ii)
         zix = z(iax) - z(ii)
         xiy = x(iay) - x(ii)
         yiy = y(iay) - y(ii)
         ziy = z(iay) - z(ii)
         vxx = xix*fix(1) + xiy*fiy(1) + xiz*fiz(1)
         vxy = 0.5d0 * (yix*fix(1) + yiy*fiy(1) + yiz*fiz(1)
     &                    + xix*fix(2) + xiy*fiy(2) + xiz*fiz(2))
         vxz = 0.5d0 * (zix*fix(1) + ziy*fiy(1) + ziz*fiz(1)
     &                    + xix*fix(3) + xiy*fiy(3) + xiz*fiz(3)) 
         vyy = yix*fix(2) + yiy*fiy(2) + yiz*fiz(2)
         vyz = 0.5d0 * (zix*fix(2) + ziy*fiy(2) + ziz*fiz(2)
     &                    + yix*fix(3) + yiy*fiy(3) + yiz*fiz(3))
         vzz = zix*fix(3) + ziy*fiy(3) + ziz*fiz(3)
         vir(1,1) = vir(1,1) + vxx
         vir(2,1) = vir(2,1) + vxy
         vir(3,1) = vir(3,1) + vxz
         vir(1,2) = vir(1,2) + vxy
         vir(2,2) = vir(2,2) + vyy
         vir(3,2) = vir(3,2) + vyz
         vir(1,3) = vir(1,3) + vxz
         vir(2,3) = vir(2,3) + vyz
         vir(3,3) = vir(3,3) + vzz
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      deallocate (tem)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ejosh1d  --  Ewald multipole derivs via list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ejosh1d" calculates the multipole energy and derivatives
c     with respect to Cartesian coordinates using particle mesh Ewald
c     summation and a neighbor list
c
c
      subroutine ejosh1d
      use sizes
      use atoms
      use boxes
      use chgpot
      use deriv
      use energi
      use ewald
      use math
      use mpole
      use virial
      use disp
      implicit none
      integer i,j,ii
      real*8 e,f
      real*8 pre
      real*8 term,fterm
      real*8 cii,dii,qii
      real*8 xd,yd,zd
      real*8 xq,yq,zq
      real*8 xv,yv,zv,vterm
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 xdfield,ydfield
      real*8 zdfield
      real*8 trq(3),frcx(3)
      real*8 frcy(3),frcz(3)
c
c
c     zero out the atomic multipole energy and derivatives
c
      em = 0.0d0
      edis = 0.0d0
      epr = 0.0d0
      do i = 1, n
         do j = 1, 3
            dem(j,i) = 0.0d0
            dedis(j,i) = 0.0d0
            depr(j,i) = 0.0d0
            permfield(j,i) = 0.0d0
         end do
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
      call ejoshreal1d
c
c     compute the reciprocal space part of the electrostatic Ewald summation
c
      call emrecip1
c
c     compute the reciprocal space part of the dispersion Ewald summation
c
      call edisprecip1
c
c     compute the Ewald self-energy term over all the atoms
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
      end do
c
c     compute the self term for the dispersion Ewald summation
c
      do i = 1, npole
         pre = adewald**6/12.0d0
         edis = edis + pre*csix(i)*csix(i)
      end do
c
c     compute self energy portion of electrostatic field
c
      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      do i = 1, npole
         do j = 1, 3
            permfield(j,i) = permfield(j,i) + term*rpole(j+1,i)
         end do
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
            xd = xd + rpole(2,i) + rpole(1,i)*x(ii)
            yd = yd + rpole(3,i) + rpole(1,i)*y(ii)
            zd = zd + rpole(4,i) + rpole(1,i)*z(ii)
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         em = em + term*(xd*xd+yd*yd+zd*zd)
         do i = 1, npole
            ii = ipole(i)
            dem(1,ii) = dem(1,ii) + 2.0d0*term*rpole(1,i)*xd
            dem(2,ii) = dem(2,ii) + 2.0d0*term*rpole(1,i)*yd
            dem(3,ii) = dem(3,ii) + 2.0d0*term*rpole(1,i)*zd
         end do
         xdfield = -2.0d0 * term * xd
         ydfield = -2.0d0 * term * yd
         zdfield = -2.0d0 * term * zd
         do i = 1, npole
            trq(1) = rpole(3,i)*zdfield - rpole(4,i)*ydfield
            trq(2) = rpole(4,i)*xdfield - rpole(2,i)*zdfield
            trq(3) = rpole(2,i)*ydfield - rpole(3,i)*xdfield
            call torque (i,trq,frcx,frcy,frcz,dem)
         end do
c
c     boundary correction to virial due to overall cell dipole
c
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         xq = 0.0d0
         yq = 0.0d0
         zq = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            xd = xd + rpole(2,i)
            yd = yd + rpole(3,i)
            zd = zd + rpole(4,i)
            xq = xq + rpole(1,i)*x(ii)
            yq = yq + rpole(1,i)*y(ii)
            zq = zq + rpole(1,i)*z(ii)
         end do
         xv = xd * xq
         yv = yd * yq
         zv = zd * zq
         vterm = term * (xd*xd + yd*yd + zd*zd + 2.0d0*(xv+yv+zv)
     &                      + xq*xq + yq*yq + zq*zq)
         vir(1,1) = vir(1,1) + 2.0d0*term*(xq*xq+xv) + vterm
         vir(2,1) = vir(2,1) + 2.0d0*term*(xq*yq+xv)
         vir(3,1) = vir(3,1) + 2.0d0*term*(xq*zq+xv)
         vir(1,2) = vir(1,2) + 2.0d0*term*(yq*xq+yv)
         vir(2,2) = vir(2,2) + 2.0d0*term*(yq*yq+yv) + vterm
         vir(3,2) = vir(3,2) + 2.0d0*term*(yq*zq+yv)
         vir(1,3) = vir(1,3) + 2.0d0*term*(zq*xq+zv)
         vir(2,3) = vir(2,3) + 2.0d0*term*(zq*yq+zv)
         vir(3,3) = vir(3,3) + 2.0d0*term*(zq*zq+zv) + vterm
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ejoshreal1d  --  Ewald real space derivs via list  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ejoshreal1d" evaluates the real space portion of the Ewald
c     summation energy and gradient due to multipole interactions
c     via a neighbor list
c
c
      subroutine ejoshreal1d
      use sizes
      use atoms
      use atomid
      use bound
      use chgpot
      use couple
      use deriv
      use energi
      use ewald
      use inter
      use math
      use molcul
      use mplpot
      use mpole
      use neigh
      use polgrp
      use polpot
      use shunt
      use virial
      use chgpen
      use disp
      use pauli
      use tarray
      use openmp
      implicit none
      integer i,j,k,m
      integer ii,kk,kkk
      integer iax,iay,iaz
      integer nlocal,nchunk
      integer tid,maxlocal
!$    integer omp_get_thread_num
      integer, allocatable :: toffset(:)
      integer, allocatable :: ilocal(:,:)
      real*8 e,efull,de,f
      real*8 dei,dek,deik
      real*8 ecc,ecv,evc,evv
      real*8 e_ele,e_disp,e_pauli
      real*8 alphai,alphak
      real*8 bfac,erfc
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 scalekk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 r,r2,r6,r7,r8
      real*8 rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 rr1core,rr3core
      real*8 rr1i,rr3i,rr5i
      real*8 rr7i,rr9i,rr11i
      real*8 rr1k,rr3k,rr5k
      real*8 rr7k,rr9k,rr11k
      real*8 rr1ik,rr3ik,rr5ik
      real*8 rr7ik,rr9ik,rr11ik
      real*8 urr3ik,urr5ik
      real*8 corei,vali
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 corek,valk
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dikx,diky,dikz
      real*8 dirx,diry,dirz
      real*8 dkrx,dkry,dkrz
      real*8 qrix,qriy,qriz
      real*8 qrkx,qrky,qrkz
      real*8 qrixr,qriyr,qrizr
      real*8 qrkxr,qrkyr,qrkzr
      real*8 qrrx,qrry,qrrz
      real*8 qikrx,qikry,qikrz
      real*8 qkirx,qkiry,qkirz
      real*8 qikrxr,qikryr,qikrzr
      real*8 qkirxr,qkiryr,qkirzr
      real*8 diqkx,diqky,diqkz
      real*8 dkqix,dkqiy,dkqiz
      real*8 diqkxr,diqkyr,diqkzr
      real*8 dkqixr,dkqiyr,dkqizr
      real*8 dqiqkx,dqiqky,dqiqkz
      real*8 dri,drk,qrri,qrrk
      real*8 diqrk,dkqri
      real*8 dik,qik,qrrik
      real*8 term1,term2,term3
      real*8 term4,term5,term6
      real*8 term1ik,term2ik,term3ik
      real*8 term4ik,term5ik,term6ik
      real*8 dterm1ik,dterm2ik,dterm3ik
      real*8 dterm4ik,dterm5ik,dterm6ik
      real*8 term1i,term2i,term3i,term4i
      real*8 term1k,term2k,term3k,term5k
      real*8 dterm1i,dterm2i,dterm3i,dterm4i
      real*8 dterm1k,dterm2k,dterm3k,dterm5k
      real*8 frcx,frcy,frcz
      real*8 frc_elex,frc_eley,frc_elez
      real*8 frc_dispx,frc_dispy,frc_dispz
      real*8 frc_prx,frc_pry,frc_prz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 vxx_ele,vyy_ele,vzz_ele
      real*8 vxy_ele,vxz_ele,vyz_ele
      real*8 vxx_disp,vyy_disp,vzz_disp
      real*8 vxy_disp,vxz_disp,vyz_disp
      real*8 vxx_pr,vyy_pr,vzz_pr
      real*8 vxy_pr,vxz_pr,vyz_pr
      real*8 rr6
      real*8 c6i,c6k,c6ik
      real*8 displam
      real*8 damp,ddamp,term,expterm
      real*8 ralpha2
      real*8 pvali,pvalk
      real*8 overlapi,overlapk,oik
      real*8 apauli,apaulk
      real*8 fid(3),fkd(3)
      real*8 ttmi(3),ttmk(3)
      real*8 fix(3),fiy(3),fiz(3)
      real*8 fixem(3),fiyem(3),fizem(3)
      real*8 fixpr(3),fiypr(3),fizpr(3)
      real*8 bn(0:5)
      real*8 lambdai(11),lambdak(11),lambdaik(11)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: muscale(:)
      real*8, allocatable :: tem(:,:)
      real*8, allocatable :: tepr(:,:)
      real*8, allocatable :: dlocal(:,:)
      character*6 mode
      external erfc
c
c
c     values for storage of mutual polarization intermediates
c
      nchunk = int(0.5d0*dble(npole)/dble(nthread)) + 1
      maxlocal = int(dble(npole)*dble(maxelst)/dble(nthread))
      nlocal = 0
      ntpair = 0
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (dscale(n))
      allocate (muscale(n))
      allocate (tem(3,n))
      allocate (tepr(3,n))
      allocate (toffset(0:nthread-1))
c
c     initialize connected atom scaling and torque arrays
c
      do i = 1, n
         mscale(i) = 1.0d0
         dscale(i) = 1.0d0
         muscale(i) = 1.0d0
         do j = 1, 3
            tem(j,i) = 0.0d0
            tepr(j,i) = 0.0d0
         end do
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
!$OMP& shared(npole,ipole,x,y,z,rpole,n12,i12,n13,i13,n14,i14,
!$OMP& np11,np12,np13,np14,ip11,ip12,ip13,ip14,
!$OMP& n15,i15,m2scale,m3scale,m4scale,m5scale,atomic,
!$OMP& d1scale,d2scale,d3scale,d4scale,
!$OMP& mu2scale,mu3scale,mu4scale,mu5scale,nelst,elst,
!$OMP& use_bounds,f,off2,aewald,molcule,xaxis,yaxis,zaxis,
!$OMP& monopole,alphaele,csix,overpauli,alphapauli,monopauli,
!$OMP& adewald,
!$OMP& ntpair,tindex,                         
!$OMP& tdipdip,toffset,maxlocal,maxelst,                     
!$OMP& nthread,nchunk)
!$OMP& firstprivate(mscale,dscale,muscale,nlocal) shared (em,
!$OMP& dem,tem,vir,permfield,edis,dedis,epr,depr,tepr)
c
c     perform dynamic allocation of some local arrays
c
      allocate (ilocal(2,maxlocal))
      allocate (dlocal(6,maxlocal))
!$OMP DO reduction(+:em,dem,tem,vir,
!$OMP& permfield,edis,dedis,epr,depr,tepr) schedule(guided)
c
c     compute the real space portion of the Ewald summation
c
      do i = 1, npole
         ii = ipole(i)
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
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
            muscale(i12(j,ii)) = mu2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
            muscale(i13(j,ii)) = mu3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
            muscale(i14(j,ii)) = mu4scale
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
            muscale(i15(j,ii)) = mu5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
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
            r2 = xr*xr + yr*yr + zr*zr
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
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               rr11 = 9.0d0 * rr9 / r2
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 5
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
c               do j = 0, 5
c                  bn(j) = f * bn(j)
c               end do
c
c     intermediates involving moments and distance separation
c
               dikx = diy*dkz - diz*dky
               diky = diz*dkx - dix*dkz
               dikz = dix*dky - diy*dkx
               dirx = diy*zr - diz*yr
               diry = diz*xr - dix*zr
               dirz = dix*yr - diy*xr
               dkrx = dky*zr - dkz*yr
               dkry = dkz*xr - dkx*zr
               dkrz = dkx*yr - dky*xr
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
               qrixr = qriz*yr - qriy*zr
               qriyr = qrix*zr - qriz*xr
               qrizr = qriy*xr - qrix*yr
               qrkxr = qrkz*yr - qrky*zr
               qrkyr = qrkx*zr - qrkz*xr
               qrkzr = qrky*xr - qrkx*yr
               qrrx = qrky*qriz - qrkz*qriy
               qrry = qrkz*qrix - qrkx*qriz
               qrrz = qrkx*qriy - qrky*qrix
               qikrx = qixx*qrkx + qixy*qrky + qixz*qrkz
               qikry = qixy*qrkx + qiyy*qrky + qiyz*qrkz
               qikrz = qixz*qrkx + qiyz*qrky + qizz*qrkz
               qkirx = qkxx*qrix + qkxy*qriy + qkxz*qriz
               qkiry = qkxy*qrix + qkyy*qriy + qkyz*qriz
               qkirz = qkxz*qrix + qkyz*qriy + qkzz*qriz
               qikrxr = qikrz*yr - qikry*zr
               qikryr = qikrx*zr - qikrz*xr
               qikrzr = qikry*xr - qikrx*yr
               qkirxr = qkirz*yr - qkiry*zr
               qkiryr = qkirx*zr - qkirz*xr
               qkirzr = qkiry*xr - qkirx*yr
               diqkx = dix*qkxx + diy*qkxy + diz*qkxz
               diqky = dix*qkxy + diy*qkyy + diz*qkyz
               diqkz = dix*qkxz + diy*qkyz + diz*qkzz
               dkqix = dkx*qixx + dky*qixy + dkz*qixz
               dkqiy = dkx*qixy + dky*qiyy + dkz*qiyz
               dkqiz = dkx*qixz + dky*qiyz + dkz*qizz
               diqrk = dix*qrkx + diy*qrky + diz*qrkz
               dkqri = dkx*qrix + dky*qriy + dkz*qriz
               diqkxr = diqkz*yr - diqky*zr
               diqkyr = diqkx*zr - diqkz*xr
               diqkzr = diqky*xr - diqkx*yr
               dkqixr = dkqiz*yr - dkqiy*zr
               dkqiyr = dkqix*zr - dkqiz*xr
               dkqizr = dkqiy*xr - dkqix*yr
               dqiqkx = diy*qrkz - diz*qrky + dky*qriz - dkz*qriy
     &                     - 2.0d0*(qixy*qkxz+qiyy*qkyz+qiyz*qkzz
     &                             -qixz*qkxy-qiyz*qkyy-qizz*qkyz)
               dqiqky = diz*qrkx - dix*qrkz + dkz*qrix - dkx*qriz
     &                     - 2.0d0*(qixz*qkxx+qiyz*qkxy+qizz*qkxz
     &                             -qixx*qkxz-qixy*qkyz-qixz*qkzz)
               dqiqkz = dix*qrky - diy*qrkx + dkx*qriy - dky*qrix
     &                     - 2.0d0*(qixx*qkxy+qixy*qkyy+qixz*qkyz
     &                             -qixy*qkxx-qiyy*qkxy-qiyz*qkxz)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     electrostatic energy and gradient!
c
c
c     calculate valence - valence interaction intermediate terms
c
               term1ik = vali*valk
               term2ik = valk*dri - vali*drk + dik
               term3ik = vali*qrrk + valk*qrri - dri*drk
     &              + 2.0d0*(dkqri-diqrk+qik)
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
c     compute damping function
c
               call damphlike(r,11,alphai,alphak,
     &              lambdai,lambdak,lambdaik)
c
c     note: removed einter accumulation
c
c
c     modify distances to account for Ewald and exclusions
c
c               scalekk = 1.0d0 - mscale(kk)
c               rr1 = bn(0) - scalekk*rr1
c               rr3 = bn(1) - scalekk*rr3
c               rr5 = bn(2) - scalekk*rr5
c               rr7 = bn(3) - scalekk*rr7
c               rr9 = bn(4) - scalekk*rr9
c               rr11 = bn(5) - scalekk*rr11
c
               rr1core = bn(0) - (1.0d0 - mscale(kk))*rr1
               rr3core = bn(1) - (1.0d0 - mscale(kk))*rr3
c               rr5core = bn(2) - (1.0d0 - mscale(kk))*rr5
c
               rr1i = bn(0) - (1.0d0 - mscale(kk)*lambdai(1))*rr1
               rr3i = bn(1) - (1.0d0 - mscale(kk)*lambdai(3))*rr3
               rr5i = bn(2) - (1.0d0 - mscale(kk)*lambdai(5))*rr5
               rr7i = bn(3) - (1.0d0 - mscale(kk)*lambdai(7))*rr7
c     
               rr1k = bn(0) - (1.0d0 - mscale(kk)*lambdak(1))*rr1
               rr3k = bn(1) - (1.0d0 - mscale(kk)*lambdak(3))*rr3
               rr5k = bn(2) - (1.0d0 - mscale(kk)*lambdak(5))*rr5
               rr7k = bn(3) - (1.0d0 - mscale(kk)*lambdak(7))*rr7
c
               rr1ik = bn(0) - (1.0d0 - mscale(kk)*lambdaik(1))*rr1
               rr3ik = bn(1) - (1.0d0 - mscale(kk)*lambdaik(3))*rr3
               rr5ik = bn(2) - (1.0d0 - mscale(kk)*lambdaik(5))*rr5
               rr7ik = bn(3) - (1.0d0 - mscale(kk)*lambdaik(7))*rr7
               rr9ik = bn(4) - (1.0d0 - mscale(kk)*lambdaik(9))*rr9
               rr11ik = bn(5) - (1.0d0 - mscale(kk)*lambdaik(11))*rr11
c
c     compute the valence - valence energy contribution for this interaction
c
               evv = term1ik*rr1ik + 
     &              term2ik*rr3ik + 
     &              term3ik*rr5ik +
     &              term4ik*rr7ik + 
     &              term5ik*rr9ik
c     
c     compute the core - valence energy contribution for this interaction
c     
               ecv = term1i*rr1i +
     &              term2i*rr3i +
     &              term3i*rr5i
c     
               evc = term1k*rr1k +
     &              term2k*rr3k +
     &              term3k*rr5k
c     
c     compute the core - core energy contribution for this interaction
c     
               ecc = term1*rr1core
c
c     compute the energy contribution for this interaction
c
               e = evv + ecv + evc + ecc
               em = em + f*e
c
c     calculate intermediate terms for force and torque
c
c
c     valence i - valence k
c
               deik = term1ik*rr3ik + 
     &              term2ik*rr5ik + 
     &              term3ik*rr7ik +
     &              term4ik*rr9ik +
     &              term5ik*rr11ik
               dterm1ik = -valk*rr3ik + 
     &              drk*rr5ik - 
     &              qrrk*rr7ik
               dterm2ik = vali*rr3ik + 
     &              dri*rr5ik + 
     &              qrri*rr7ik
               dterm3ik = 2.0d0 * rr5ik
               dterm4ik = 2.0d0 * (-valk*rr5ik +
     &              drk*rr7ik - 
     &              qrrk*rr9ik)
               dterm5ik = 2.0d0 * (-vali*rr5ik -
     &              dri*rr7ik - 
     &              qrri*rr9ik)
               dterm6ik = 4.0d0 * rr7ik
c
c     compute the force components for this interaction
c
               frc_elex = deik*xr + dterm1ik*dix + dterm2ik*dkx
     &                   + dterm3ik*(diqkx-dkqix) + dterm4ik*qrix
     &                   + dterm5ik*qrkx + dterm6ik*(qikrx+qkirx)
               frc_eley = deik*yr + dterm1ik*diy + dterm2ik*dky
     &                   + dterm3ik*(diqky-dkqiy) + dterm4ik*qriy
     &                   + dterm5ik*qrky + dterm6ik*(qikry+qkiry)
               frc_elez = deik*zr + dterm1ik*diz + dterm2ik*dkz
     &                   + dterm3ik*(diqkz-dkqiz) + dterm4ik*qriz
     &                   + dterm5ik*qrkz + dterm6ik*(qikrz+qkirz)
c
c     compute the torque components for this interaction
c
               ttmi(1) = -rr3ik*dikx + dterm1ik*dirx + 
     &              dterm3ik*(dqiqkx+dkqixr)
     &                      - dterm4ik*qrixr - dterm6ik*(qikrxr+qrrx)
               ttmi(2) = -rr3ik*diky + dterm1ik*diry + 
     &              dterm3ik*(dqiqky+dkqiyr)
     &                      - dterm4ik*qriyr - dterm6ik*(qikryr+qrry)
               ttmi(3) = -rr3ik*dikz + dterm1ik*dirz + 
     &              dterm3ik*(dqiqkz+dkqizr)
     &                      - dterm4ik*qrizr - dterm6ik*(qikrzr+qrrz)
               ttmk(1) = rr3ik*dikx + dterm2ik*dkrx - 
     &              dterm3ik*(dqiqkx+diqkxr)
     &                      - dterm5ik*qrkxr - dterm6ik*(qkirxr-qrrx)
               ttmk(2) = rr3ik*diky + dterm2ik*dkry - 
     &              dterm3ik*(dqiqky+diqkyr)
     &                      - dterm5ik*qrkyr - dterm6ik*(qkiryr-qrry)
               ttmk(3) = rr3ik*dikz + dterm2ik*dkrz - 
     &              dterm3ik*(dqiqkz+diqkzr)
     &                      - dterm5ik*qrkzr - dterm6ik*(qkirzr-qrrz)
c
c     valence - core force and torque
c
               dei = term1i*rr3i +
     &              term2i*rr5i +
     &              term3i*rr7i 
               dek = term1k*rr3k +
     &              term2k*rr5k +
     &              term3k*rr7k 
               dterm1i = -corek*rr3i
               dterm2k = corei*rr3k 
               dterm4i = 2.0d0 * (-corek*rr5i)
               dterm5k = 2.0d0 * (-corei*rr5k) 
c
               frc_elex = frc_elex + (dei + dek)*xr + dterm1i*dix 
     &              + dterm2k*dkx
     &              + dterm4i*qrix
     &              + dterm5k*qrkx  
               frc_eley = frc_eley + (dei + dek)*yr + dterm1i*diy 
     &              + dterm2k*dky
     &              + dterm4i*qriy
     &              + dterm5k*qrky  
               frc_elez = frc_elez + (dei + dek)*zr + dterm1i*diz 
     &              + dterm2k*dkz
     &              + dterm4i*qriz
     &              + dterm5k*qrkz  
c
               ttmi(1) = ttmi(1) + dterm1i*dirx 
     &              - dterm4i*qrixr 
               ttmi(2) = ttmi(2) + dterm1i*diry 
     &              - dterm4i*qriyr 
               ttmi(3) = ttmi(3) + dterm1i*dirz 
     &              - dterm4i*qrizr 
               ttmk(1) = ttmk(1) + dterm2k*dkrx 
     &              - dterm5k*qrkxr 
               ttmk(2) = ttmk(2) + dterm2k*dkry 
     &              - dterm5k*qrkyr 
               ttmk(3) = ttmk(3) + dterm2k*dkrz 
     &              - dterm5k*qrkzr 
c
c     core - core force (no torque)
c
               de = term1*rr3core
c
               frc_elex = frc_elex + de*xr 
               frc_eley = frc_eley + de*yr 
               frc_elez = frc_elez + de*zr
c
c     increment force-based gradient and torque on first site
c
               dem(1,ii) = dem(1,ii) + frc_elex * f 
               dem(2,ii) = dem(2,ii) + frc_eley * f 
               dem(3,ii) = dem(3,ii) + frc_elez * f 
               tem(1,i) = tem(1,i) + ttmi(1) * f 
               tem(2,i) = tem(2,i) + ttmi(2) * f 
               tem(3,i) = tem(3,i) + ttmi(3) * f 
c
c     increment force-based gradient and torque on second site
c
               dem(1,kk) = dem(1,kk) - frc_elex * f 
               dem(2,kk) = dem(2,kk) - frc_eley * f 
               dem(3,kk) = dem(3,kk) - frc_elez * f 
               tem(1,k) = tem(1,k) + ttmk(1) * f 
               tem(2,k) = tem(2,k) + ttmk(2) * f 
               tem(3,k) = tem(3,k) + ttmk(3) * f
c
c
c     save permanent electric field for induced dipole calculation
c     note: this already has the core contribution
c
               rr1core = bn(0) - (1.0d0 - dscale(kk))*rr1
               rr3core = bn(1) - (1.0d0 - dscale(kk))*rr3
c
               rr1i = bn(0) - (1.0d0 - dscale(kk)*lambdai(1))*rr1
               rr3i = bn(1) - (1.0d0 - dscale(kk)*lambdai(3))*rr3
               rr5i = bn(2) - (1.0d0 - dscale(kk)*lambdai(5))*rr5
               rr7i = bn(3) - (1.0d0 - dscale(kk)*lambdai(7))*rr7
c
               rr1k = bn(0) - (1.0d0 - dscale(kk)*lambdak(1))*rr1
               rr3k = bn(1) - (1.0d0 - dscale(kk)*lambdak(3))*rr3
               rr5k = bn(2) - (1.0d0 - dscale(kk)*lambdak(5))*rr5
               rr7k = bn(3) - (1.0d0 - dscale(kk)*lambdak(7))*rr7
c
               fid(1) = -xr*(rr3core*corek + rr3k*valk -
     &              rr5k*drk + rr7k*qrrk)
     &              - rr3k*dkx + 2.0d0*rr5k*qrkx
               fid(2) = -yr*(rr3core*corek + rr3k*valk -
     &              rr5k*drk+rr7k*qrrk)
     &              - rr3k*dky + 2.0d0*rr5k*qrky
               fid(3) = -zr*(rr3core*corek + rr3k*valk -
     &              rr5k*drk+rr7k*qrrk)
     &              - rr3k*dkz + 2.0d0*rr5k*qrkz
               fkd(1) = xr*(rr3core*corei + rr3i*vali +
     &              rr5i*dri + rr7i*qrri)
     &              - rr3i*dix - 2.0d0*rr5i*qrix
               fkd(2) = yr*(rr3core*corei + rr3i*vali +
     &              rr5i*dri + rr7i*qrri)
     &              - rr3i*diy - 2.0d0*rr5i*qriy
               fkd(3) = zr*(rr3core*corei + rr3i*vali +
     &              rr5i*dri + rr7i*qrri)
     &              - rr3i*diz - 2.0d0*rr5i*qriz
c
c     increment electric field on both sites
c
               do j = 1, 3
                  permfield(j,i) = permfield(j,i) + fid(j)
                  permfield(j,k) = permfield(j,k) + fkd(j)
               end do
c
c     save diple - dipole t matrix for mutual induction
c
c     INSERT MUTUAL EXCLUSION RULES HERE!!!
c
               urr3ik = bn(1) - (1.0d0 - muscale(kk)*lambdaik(3))*rr3
               urr5ik = bn(2) - (1.0d0 - muscale(kk)*lambdaik(5))*rr5
               nlocal = nlocal + 1
               ilocal(1,nlocal) = i
               ilocal(2,nlocal) = k
               dlocal(1,nlocal) = -urr3ik + urr5ik*xr*xr
               dlocal(2,nlocal) = urr5ik*xr*yr
               dlocal(3,nlocal) = urr5ik*xr*zr
               dlocal(4,nlocal) = -urr3ik + urr5ik*yr*yr
               dlocal(5,nlocal) = urr5ik*yr*zr
               dlocal(6,nlocal) = -urr3ik + urr5ik*zr*zr
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     dispersion energy and gradient!
c
               r6 = r2**3
               r7 = r6*r
               r8 = r6*r2
               c6ik = c6i*c6k
               displam = 0.5d0*(3.0d0*lambdaik(5) - lambdaik(3))
c
c     get dispersion ewald coefficients
c
               ralpha2 = r2 * adewald**2
               damp = 1.0d0
               ddamp = 0.0d0
               if (ralpha2 .lt. 50.0d0) then
                  expterm = exp(-ralpha2)
                  term = 1.0d0 + ralpha2 + 0.5d0*ralpha2**2
                  damp = term*expterm
                  ddamp = -(ralpha2**3)*expterm/r
               end if
c
c     compute real space dispersion energy
c
               e_disp = -c6ik*(damp +
     &              (mscale(kk)*displam**2 - 1.0d0))/r6
               edis = edis + e_disp
c
c     semi-hack: i'm using lambda(10) to stuff in the neccessary terms              
c     see DispersionDamping.mw notebook
c
               de = c6ik*6.0d0*(damp + mscale(kk)*displam**2 - 
     &              1.0d0)/r8 - c6ik*ddamp/r7 - 
     &              c6ik*mscale(kk)*2.0d0*displam*lambdaik(10)/r7
c
               frc_dispx = de * xr
               frc_dispy = de * yr
               frc_dispz = de * zr
c
c     increment force-based gradient on first site 
c
               dedis(1,ii) = dedis(1,ii) - frc_dispx
               dedis(2,ii) = dedis(2,ii) - frc_dispy
               dedis(3,ii) = dedis(3,ii) - frc_dispz
c
c     increment force-based gradient on second site 
c
               dedis(1,kk) = dedis(1,kk) + frc_dispx
               dedis(2,kk) = dedis(2,kk) + frc_dispy
               dedis(3,kk) = dedis(3,kk) + frc_dispz
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     pauli repulsion energy and gradient!
c
               call damppauli(r,r2,rr1,rr3,rr5,rr7,rr9,rr11,11,
     &              apauli,apaulk,lambdaik)
c
c     recompute terms with number of pauli valence electrons 
c
               term1ik = pvali*pvalk
               term2ik = pvalk*dri - pvali*drk + dik
               term3ik = pvali*qrrk + pvalk*qrri - dri*drk
     &              + 2.0d0*(dkqri-diqrk+qik)
c     
c     compute valence - valence energy contribution for this interaction
c     (pauli repulsion has no terms involving the core)
c
               evv = term1ik*lambdaik(1) +
     &              term2ik*lambdaik(3) +
     &              term3ik*lambdaik(5) +
     &              term4ik*lambdaik(7) +
     &              term5ik*lambdaik(9)
c
c     combining rule for pauli repulsion prefactor 
c
               oik = overlapi*overlapk
c
c     total pauli repulsion energy 
c
               e_pauli = oik * mscale(kk) * evv * rr1 
               epr = epr + e_pauli
c     
c     now compute gradient in same way as electrostatics
c
c
c     calculate intermediate terms for force and torque
c
c     valence i - valence k
c
               deik = term1ik*lambdaik(3) +
     &              term2ik*lambdaik(5) +
     &              term3ik*lambdaik(7) +
     &              term4ik*lambdaik(9) +
     &              term5ik*lambdaik(11)
               dterm1ik = -pvalk*lambdaik(3) +
     &              drk*lambdaik(5) -
     &              qrrk*lambdaik(7)
               dterm2ik = pvali*lambdaik(3) +
     &              dri*lambdaik(5) +
     &              qrri*lambdaik(7)
               dterm3ik = 2.0d0 * lambdaik(5)
               dterm4ik = 2.0d0 * (-pvalk*lambdaik(5) +
     &              drk*lambdaik(7) -
     &              qrrk*lambdaik(9))
               dterm5ik = 2.0d0 * (-pvali*lambdaik(5) -
     &              dri*lambdaik(7) -
     &              qrri*lambdaik(9))
               dterm6ik = 4.0d0 * lambdaik(7)
c     
c     compute the force components for this interaction
c     
               frc_prx = deik*xr + dterm1ik*dix + dterm2ik*dkx
     &              + dterm3ik*(diqkx-dkqix) + dterm4ik*qrix
     &              + dterm5ik*qrkx + dterm6ik*(qikrx+qkirx)
               frc_pry = deik*yr + dterm1ik*diy + dterm2ik*dky
     &              + dterm3ik*(diqky-dkqiy) + dterm4ik*qriy
     &              + dterm5ik*qrky + dterm6ik*(qikry+qkiry)
               frc_prz = deik*zr + dterm1ik*diz + dterm2ik*dkz
     &              + dterm3ik*(diqkz-dkqiz) + dterm4ik*qriz
     &              + dterm5ik*qrkz + dterm6ik*(qikrz+qkirz)
c
               frc_prx = frc_prx*rr1
               frc_pry = frc_pry*rr1
               frc_prz = frc_prz*rr1
c
c     compute the torque components for this interaction
c
               ttmi(1) = -lambdaik(3)*dikx + dterm1ik*dirx + 
     &              dterm3ik*(dqiqkx+dkqixr)
     &              - dterm4ik*qrixr - dterm6ik*(qikrxr+qrrx)
               ttmi(2) = -lambdaik(3)*diky + dterm1ik*diry + 
     &              dterm3ik*(dqiqky+dkqiyr)
     &              - dterm4ik*qriyr - dterm6ik*(qikryr+qrry)
               ttmi(3) = -lambdaik(3)*dikz + dterm1ik*dirz + 
     &              dterm3ik*(dqiqkz+dkqizr)
     &              - dterm4ik*qrizr - dterm6ik*(qikrzr+qrrz)
               ttmk(1) = lambdaik(3)*dikx + dterm2ik*dkrx - 
     &              dterm3ik*(dqiqkx+diqkxr)
     &              - dterm5ik*qrkxr - dterm6ik*(qkirxr-qrrx)
               ttmk(2) = lambdaik(3)*diky + dterm2ik*dkry - 
     &              dterm3ik*(dqiqky+diqkyr)
     &              - dterm5ik*qrkyr - dterm6ik*(qkiryr-qrry)
               ttmk(3) = lambdaik(3)*dikz + dterm2ik*dkrz - 
     &              dterm3ik*(dqiqkz+diqkzr)
     &              - dterm5ik*qrkzr - dterm6ik*(qkirzr-qrrz)
c
               ttmi(1) = ttmi(1)*rr1*oik*mscale(kk)
               ttmi(2) = ttmi(2)*rr1*oik*mscale(kk)
               ttmi(3) = ttmi(3)*rr1*oik*mscale(kk)
               ttmk(1) = ttmk(1)*rr1*oik*mscale(kk)
               ttmk(2) = ttmk(2)*rr1*oik*mscale(kk)
               ttmk(3) = ttmk(3)*rr1*oik*mscale(kk)
c
c     now take chain rule terms of 1/r
c     note: there are no torques for chain rule terms?
c
               frc_prx = frc_prx + evv*rr3*xr
               frc_pry = frc_pry + evv*rr3*yr
               frc_prz = frc_prz + evv*rr3*zr
c
c     multiply by overlap prefactor
c
               frc_prx = oik*mscale(kk)*frc_prx
               frc_pry = oik*mscale(kk)*frc_pry
               frc_prz = oik*mscale(kk)*frc_prz
c
c     increment force-based gradient and torque on first site
c
               depr(1,ii) = depr(1,ii) + frc_prx
               depr(2,ii) = depr(2,ii) + frc_pry
               depr(3,ii) = depr(3,ii) + frc_prz
               tepr(1,i) = tepr(1,i) + ttmi(1)
               tepr(2,i) = tepr(2,i) + ttmi(2)
               tepr(3,i) = tepr(3,i) + ttmi(3)
c
c     increment force-based gradient and torque on second site
c
               depr(1,kk) = depr(1,kk) - frc_prx
               depr(2,kk) = depr(2,kk) - frc_pry
               depr(3,kk) = depr(3,kk) - frc_prz
               tepr(1,k) = tepr(1,k) + ttmk(1)
               tepr(2,k) = tepr(2,k) + ttmk(2)
               tepr(3,k) = tepr(3,k) + ttmk(3)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     increment the virial due to pairwise Cartesian forces
c
               vxx_ele = -xr * f * frc_elex
               vxy_ele = -yr * f * frc_elex
               vxz_ele = -zr * f * frc_elex
               vyy_ele = -yr * f * frc_eley
               vyz_ele = -zr * f * frc_eley
               vzz_ele = -zr * f * frc_elez
c
               vxx_disp = xr * frc_dispx
               vxy_disp = yr * frc_dispx
               vxz_disp = zr * frc_dispx
               vyy_disp = yr * frc_dispy
               vyz_disp = zr * frc_dispy
               vzz_disp = zr * frc_dispz
c
               vxx_pr = -xr * frc_prx
               vxy_pr = -yr * frc_prx
               vxz_pr = -zr * frc_prx
               vyy_pr = -yr * frc_pry
               vyz_pr = -zr * frc_pry
               vzz_pr = -zr * frc_prz
c
               vir(1,1) = vir(1,1) + vxx_ele + vxx_disp + vxx_pr 
               vir(2,1) = vir(2,1) + vxy_ele + vxy_disp + vxy_pr
               vir(3,1) = vir(3,1) + vxz_ele + vxz_disp + vxz_pr
               vir(1,2) = vir(1,2) + vxy_ele + vxy_disp + vxy_pr
               vir(2,2) = vir(2,2) + vyy_ele + vyy_disp + vyy_pr
               vir(3,2) = vir(3,2) + vyz_ele + vyz_disp + vyz_pr
               vir(1,3) = vir(1,3) + vxz_ele + vxz_disp + vxz_pr
               vir(2,3) = vir(2,3) + vyz_ele + vyz_disp + vyz_pr
               vir(3,3) = vir(3,3) + vzz_ele + vzz_disp + vzz_pr
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
            muscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
            muscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
            muscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
            muscale(i15(j,ii)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
c
c     find offset into global arrays for the current thread
c
!$OMP CRITICAL
      tid = 0
!$    tid = omp_get_thread_num ()
      toffset(tid) = ntpair
      ntpair = ntpair + nlocal
!$OMP END CRITICAL
c
c     store terms used later for mutual polarization
c
      k = toffset(tid)
      do i = 1, nlocal
         m = k + i
         tindex(1,m) = ilocal(1,i)
         tindex(2,m) = ilocal(2,i)
         do j = 1, 6
            tdipdip(j,m) = dlocal(j,i)
         end do
      end do
      deallocate (ilocal)
      deallocate (dlocal)
c
c
!$OMP DO reduction(+:dem,depr,vir) schedule(guided)
c
c     resolve site torques then increment forces and virial
c
      do i = 1, npole
c
c     torques from paulir repulsion
c
         call torque (i,tepr(1,i),fixpr,fiypr,fizpr,depr)
c
c     torques from electrostatics
c
         call torque (i,tem(1,i),fixem,fiyem,fizem,dem)
c
c     add together electrostatic and pauli repulsion torques
c
         fix = fixpr + fixem
         fiy = fiypr + fiyem
         fiz = fizpr + fizem
c
         ii = ipole(i)
         iaz = zaxis(i)
         iax = xaxis(i)
         iay = yaxis(i)
         if (iaz .eq. 0)  iaz = ii
         if (iax .eq. 0)  iax = ii
         if (iay .eq. 0)  iay = ii
         xiz = x(iaz) - x(ii)
         yiz = y(iaz) - y(ii)
         ziz = z(iaz) - z(ii)
         xix = x(iax) - x(ii)
         yix = y(iax) - y(ii)
         zix = z(iax) - z(ii)
         xiy = x(iay) - x(ii)
         yiy = y(iay) - y(ii)
         ziy = z(iay) - z(ii)
         vxx = xix*fix(1) + xiy*fiy(1) + xiz*fiz(1)
         vxy = 0.5d0 * (yix*fix(1) + yiy*fiy(1) + yiz*fiz(1)
     &                    + xix*fix(2) + xiy*fiy(2) + xiz*fiz(2))
         vxz = 0.5d0 * (zix*fix(1) + ziy*fiy(1) + ziz*fiz(1)
     &                    + xix*fix(3) + xiy*fiy(3) + xiz*fiz(3)) 
         vyy = yix*fix(2) + yiy*fiy(2) + yiz*fiz(2)
         vyz = 0.5d0 * (zix*fix(2) + ziy*fiy(2) + ziz*fiz(2)
     &                    + yix*fix(3) + yiy*fiy(3) + yiz*fiz(3))
         vzz = zix*fix(3) + ziy*fiy(3) + ziz*fiz(3)
         vir(1,1) = vir(1,1) + vxx
         vir(2,1) = vir(2,1) + vxy
         vir(3,1) = vir(3,1) + vxz
         vir(1,2) = vir(1,2) + vxy
         vir(2,2) = vir(2,2) + vyy
         vir(3,2) = vir(3,2) + vyz
         vir(1,3) = vir(1,3) + vxz
         vir(2,3) = vir(2,3) + vyz
         vir(3,3) = vir(3,3) + vzz
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
      deallocate (muscale)
      deallocate (tem)
      deallocate (tepr)
c
c      do i = 1, 6
c         do j =1, n*maxelst
c            if (tdipdip(i,j) .ne. 0.0d0) then
c               print *,"ejosh1 tdipdip",tdipdip(i,j)
c            end if
c         end do
c      end do
c
      return
      end
c
c
c     ####################################################################
c     ##                                                                ##
c     ##  subroutine ejoshrecip1  --  PME recip multipole energy & derivs  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     "ejoshrecip1" evaluates the reciprocal space portion of the particle
c     mesh Ewald summation energy and gradient due to multipoles
c
c     literature reference:
c
c     C. Sagui, L. G. Pedersen and T. A. Darden, "Towards an Accurate
c     Representation of Electrostatics in Classical Force Fields:
c     Efficient Implementation of Multipolar Interactions in
c     Biomolecular Simulations", Journal of Chemical Physics, 120,
c     73-87 (2004)
c
c     modifications for nonperiodic systems suggested by Tom Darden
c     during May 2007
c
c
      subroutine ejoshrecip1
      use sizes
      use atoms
      use bound
      use boxes
      use chgpot
      use deriv
      use energi
      use ewald
      use math
      use mpole
      use mrecip
      use pme
      use virial
      implicit none
      integer i,j,k,ii
      integer k1,k2,k3
      integer m1,m2,m3
      integer iax,iay,iaz
      integer ntot,nff
      integer nf1,nf2,nf3
      integer deriv1(10)
      integer deriv2(10)
      integer deriv3(10)
      real*8 e,eterm,f
      real*8 r1,r2,r3
      real*8 h1,h2,h3
      real*8 f1,f2,f3
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 volterm,denom
      real*8 hsq,expterm
      real*8 term,pterm
      real*8 vterm,struc2
      real*8 trq(3),fix(3)
      real*8 fiy(3),fiz(3)
c
c     indices into the electrostatic field array
c
      data deriv1  / 2, 5,  8,  9, 11, 16, 18, 14, 15, 20 /
      data deriv2  / 3, 8,  6, 10, 14, 12, 19, 16, 20, 17 /
      data deriv3  / 4, 9, 10,  7, 15, 17, 13, 20, 18, 19 /
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
      f = electric / dielec
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(cmp)) then
         if (size(cmp) .lt. 10*npole) then
            deallocate (cmp)
            deallocate (fmp)
            deallocate (cphi)
            deallocate (fphi)
         end if
      end if
      if (.not. allocated(cmp)) then
         allocate (cmp(10,npole))
         allocate (fmp(10,npole))
         allocate (cphi(10,npole))
         allocate (fphi(20,npole))
      end if
c
c     zero out the temporary virial accumulation variables
c
      vxx = 0.0d0
      vxy = 0.0d0
      vxz = 0.0d0
      vyy = 0.0d0
      vyz = 0.0d0
      vzz = 0.0d0
c
c     copy multipole moments and coordinates to local storage
c
      do i = 1, npole
         cmp(1,i) = rpole(1,i)
         cmp(2,i) = rpole(2,i)
         cmp(3,i) = rpole(3,i)
         cmp(4,i) = rpole(4,i)
         cmp(5,i) = rpole(5,i)
         cmp(6,i) = rpole(9,i)
         cmp(7,i) = rpole(13,i)
         cmp(8,i) = 2.0d0 * rpole(6,i)
         cmp(9,i) = 2.0d0 * rpole(7,i)
         cmp(10,i) = 2.0d0 * rpole(10,i)
      end do
c
c     compute the arrays of B-spline coefficients
c
      call bspline_fill
      call table_fill
c
c     assign permanent multipoles to PME grid and perform
c     the 3-D FFT forward transformation
c
      call cmp_to_fmp (cmp,fmp)
      call grid_mpole (fmp)
      call fftfront
c
c     initialize variables required for the scalar summation
c
      ntot = nfft1 * nfft2 * nfft3
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nff = nfft1 * nfft2
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
c
c     make the scalar summation over reciprocal lattice
c
      do i = 1, ntot-1
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
            eterm = 0.5d0 * electric * expterm * struc2
            vterm = (2.0d0/hsq) * (1.0d0-term) * eterm
            vxx = vxx + h1*h1*vterm - eterm
            vxy = vxy + h1*h2*vterm
            vxz = vxz + h1*h3*vterm
            vyy = vyy + h2*h2*vterm - eterm
            vyz = vyz + h2*h3*vterm
            vzz = vzz + h3*h3*vterm - eterm
         end if
         qfac(k1,k2,k3) = expterm
      end do
c
c     save the virial for use in polarization computation
c
      vmxx = vxx
      vmxy = vxy
      vmxz = vxz
      vmyy = vyy
      vmyz = vyz
      vmzz = vzz
c
c     account for zeroth grid point for nonperiodic system
c
      qfac(1,1,1) = 0.0d0
      if (.not. use_bounds) then
         expterm = 0.5d0 * pi / xbox
         struc2 = qgrid(1,1,1,1)**2 + qgrid(2,1,1,1)**2
         e = f * expterm * struc2
         em = em + e
         qfac(1,1,1) = expterm
      end if
c
c     complete the transformation of the PME grid
c
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
               term = qfac(i,j,k)
               qgrid(1,i,j,k) = term * qgrid(1,i,j,k)
               qgrid(2,i,j,k) = term * qgrid(2,i,j,k)
            end do
         end do
      end do
c
c     perform 3-D FFT backward transform and get potential
c
      call fftback
      call fphi_mpole (fphi)
      do i = 1, npole
         do j = 1, 20
            fphi(j,i) = f * fphi(j,i)
         end do
      end do
      call fphi_to_cphi (fphi,cphi)
c
c     increment the permanent multipole energy and gradient
c
      e = 0.0d0
      do i = 1, npole
         f1 = 0.0d0
         f2 = 0.0d0
         f3 = 0.0d0
         do k = 1, 10
            e = e + fmp(k,i)*fphi(k,i)
            f1 = f1 + fmp(k,i)*fphi(deriv1(k),i)
            f2 = f2 + fmp(k,i)*fphi(deriv2(k),i)
            f3 = f3 + fmp(k,i)*fphi(deriv3(k),i)
         end do
         f1 = dble(nfft1) * f1
         f2 = dble(nfft2) * f2
         f3 = dble(nfft3) * f3
         h1 = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
         h2 = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
         h3 = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
         ii = ipole(i)
         dem(1,ii) = dem(1,ii) + h1
         dem(2,ii) = dem(2,ii) + h2
         dem(3,ii) = dem(3,ii) + h3
      end do
      e = 0.5d0 * e
      em = em + e
c
c     increment the permanent multipole virial contributions
c
      do i = 1, npole
         vxx = vxx - cmp(2,i)*cphi(2,i) - 2.0d0*cmp(5,i)*cphi(5,i)
     &            - cmp(8,i)*cphi(8,i) - cmp(9,i)*cphi(9,i)
         vxy = vxy - 0.5d0*(cmp(3,i)*cphi(2,i)+cmp(2,i)*cphi(3,i))
     &            - (cmp(5,i)+cmp(6,i))*cphi(8,i)
     &            - 0.5d0*cmp(8,i)*(cphi(5,i)+cphi(6,i))
     &            - 0.5d0*(cmp(9,i)*cphi(10,i)+cmp(10,i)*cphi(9,i))
         vxz = vxz - 0.5d0*(cmp(4,i)*cphi(2,i)+cmp(2,i)*cphi(4,i))
     &            - (cmp(5,i)+cmp(7,i))*cphi(9,i)
     &            - 0.5d0*cmp(9,i)*(cphi(5,i)+cphi(7,i))
     &            - 0.5d0*(cmp(8,i)*cphi(10,i)+cmp(10,i)*cphi(8,i))
         vyy = vyy - cmp(3,i)*cphi(3,i) - 2.0d0*cmp(6,i)*cphi(6,i)
     &            - cmp(8,i)*cphi(8,i) - cmp(10,i)*cphi(10,i)
         vyz = vyz - 0.5d0*(cmp(4,i)*cphi(3,i)+cmp(3,i)*cphi(4,i))
     &            - (cmp(6,i)+cmp(7,i))*cphi(10,i)
     &            - 0.5d0*cmp(10,i)*(cphi(6,i)+cphi(7,i))
     &            - 0.5d0*(cmp(8,i)*cphi(9,i)+cmp(9,i)*cphi(8,i))
         vzz = vzz - cmp(4,i)*cphi(4,i) - 2.0d0*cmp(7,i)*cphi(7,i)
     &            - cmp(9,i)*cphi(9,i) - cmp(10,i)*cphi(10,i)
      end do
c
c     resolve site torques then increment forces and virial
c
      do i = 1, npole
         trq(1) = cmp(4,i)*cphi(3,i) - cmp(3,i)*cphi(4,i)
     &               + 2.0d0*(cmp(7,i)-cmp(6,i))*cphi(10,i)
     &               + cmp(9,i)*cphi(8,i) + cmp(10,i)*cphi(6,i)
     &               - cmp(8,i)*cphi(9,i) - cmp(10,i)*cphi(7,i)
         trq(2) = cmp(2,i)*cphi(4,i) - cmp(4,i)*cphi(2,i)
     &               + 2.0d0*(cmp(5,i)-cmp(7,i))*cphi(9,i)
     &               + cmp(8,i)*cphi(10,i) + cmp(9,i)*cphi(7,i)
     &               - cmp(9,i)*cphi(5,i) - cmp(10,i)*cphi(8,i)
         trq(3) = cmp(3,i)*cphi(2,i) - cmp(2,i)*cphi(3,i)
     &               + 2.0d0*(cmp(6,i)-cmp(5,i))*cphi(8,i)
     &               + cmp(8,i)*cphi(5,i) + cmp(10,i)*cphi(9,i)
     &               - cmp(8,i)*cphi(6,i) - cmp(9,i)*cphi(10,i)
         call torque (i,trq,fix,fiy,fiz,dem)
         ii = ipole(i)
         iaz = zaxis(i)
         iax = xaxis(i)
         iay = yaxis(i)
         if (iaz .eq. 0)  iaz = ii
         if (iax .eq. 0)  iax = ii
         if (iay .eq. 0)  iay = ii
         xiz = x(iaz) - x(ii)
         yiz = y(iaz) - y(ii)
         ziz = z(iaz) - z(ii)
         xix = x(iax) - x(ii)
         yix = y(iax) - y(ii)
         zix = z(iax) - z(ii)
         xiy = x(iay) - x(ii)
         yiy = y(iay) - y(ii)
         ziy = z(iay) - z(ii)
         vxx = vxx + xix*fix(1) + xiy*fiy(1) + xiz*fiz(1)
         vxy = vxy + 0.5d0*(yix*fix(1) + yiy*fiy(1) + yiz*fiz(1)
     &                        + xix*fix(2) + xiy*fiy(2) + xiz*fiz(2))
         vxz = vxz + 0.5d0*(zix*fix(1) + ziy*fiy(1) + ziz*fiz(1)
     &                        + xix*fix(3) + xiy*fiy(3) + xiz*fiz(3)) 
         vyy = vyy + yix*fix(2) + yiy*fiy(2) + yiz*fiz(2)
         vyz = vyz + 0.5d0*(zix*fix(2) + ziy*fiy(2) + ziz*fiz(2)
     &                        + yix*fix(3) + yiy*fiy(3) + yiz*fiz(3))
         vzz = vzz + zix*fix(3) + ziy*fiy(3) + ziz*fiz(3)
      end do
c
c     increment the total internal virial tensor components
c
      vir(1,1) = vir(1,1) + vxx
      vir(2,1) = vir(2,1) + vxy
      vir(3,1) = vir(3,1) + vxz
      vir(1,2) = vir(1,2) + vxy
      vir(2,2) = vir(2,2) + vyy
      vir(3,2) = vir(3,2) + vyz
      vir(1,3) = vir(1,3) + vxz
      vir(2,3) = vir(2,3) + vyz
      vir(3,3) = vir(3,3) + vzz
      return
      end
