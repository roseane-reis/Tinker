c
c     #############################################################
c     ## COPYRIGHT (C) 2016 by Josh Rackers & Jay William Ponder ##
c     ##                   All Rights Reserved                   ## 
c     ############################################################# 
c
c     ####################################################################
c     ##                                                                ##
c     ##  damping  --  set of routines to generate damping coefficients ##
c     ##                                                                ##
c     ####################################################################
c
c
c     "damping" generates the damping coefficients that go with 
c     corresponding powers of r for various damping functions
c
c
c
c     ####################################################################
c     ##                                                                ##
c     ##  subroutine dampewald  --  generate ewald damping coefficients ##
c     ##                                                                ##
c     ####################################################################
c
c
c     "dampewald" generates the damping coefficients for the ewald error
c     function damping that go with corresponding powers of r
c
      subroutine dampewald(i,k,rorder,r,r2,scale)
      use sizes
      use ewald
      use math
      implicit none
      integer i,k,j,maxj
      integer rorder
      real*8 r,r2
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 bfac
      real*8 scale(*)
      real*8, allocatable :: bn(:)
c
c     set damping factors to one
c
      do j = 1, rorder
         scale(j) = 1.0d0
      end do
c
c     set max order needed for ewald damping factors
c
      maxj = (rorder - 1)/2
      allocate (bn(0:maxj))
c     
c     compute ewald damping factors
c
      ralpha = aewald * r
      bn(0) = erfc(ralpha) / r
      scale(1) = bn(0)
      alsq2 = 2.0d0 * aewald**2
      alsq2n = 0.0d0
      if (aewald .gt. 0.0d0) alsq2n = 1.0d0 / (sqrtpi*aewald)
      exp2a = exp(-ralpha**2)
      do j = 1, maxj
         bfac = dble(j+j-1)
         alsq2n = alsq2 * alsq2n
         bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
         scale(2*j+1) = bn(j)
      end do
c
c     deallocate local array
c
      deallocate (bn)
      return
      end
c
c
c     ####################################################################
c     ##                                                                ##
c     ##  subroutine dampthole  --  generate thole damping coefficients ##
c     ##                                                                ##
c     ####################################################################
c
c
c     "dampthole" generates the damping coefficients for the thole 
c     damping functional form that go with corresponding powers of r
c
      subroutine dampthole(i,k,rorder,r,scale)
      use sizes
      use polar
      implicit none
      integer i,k,j
      integer rorder
      real*8 r
      real*8 damp,expdamp
      real*8 pdi,pti,pdk,ptk
      real*8 pgamma
      real*8 scale(*)
c     
c     set damping factors to one
c
      do j = 1, rorder
         scale(j) = 1.0d0
      end do
c
c     read in damping parameters
c
      pdi = pdamp(i)
      pti = thole(i)
      pdk = pdamp(k)
      ptk = thole(k)
c
c     assign thole damping scale factors
c
      damp = pdi * pdk
      if (damp .ne. 0.0d0) then
         pgamma = min(pti,ptk)
         damp = pgamma * (r/damp)**3
         if (damp .lt. 50.0d0) then
            expdamp = exp(-damp)
            scale(3) = 1.0d0 - expdamp
            scale(5) = 1.0d0 - (1.0d0 + damp)*expdamp
            if (rorder.ge.7) then
               scale(7) = 1.0d0 - (1.0d0 + damp + 0.6d0*damp**2)*expdamp
            end if
            if (rorder.ge.9) then
               scale(9) = 1.0d0 - (1.0d0 + damp + 
     &              (18.0d0/35.0d0)*damp**2 + (9.0d0/35.0d0)*damp**3)*
     &              expdamp
            end if
         end if
      end if
      return
      end
c
c
c     ######################################################################
c     ##                                                                  ##
c     ##  subroutine damphlike  --  generate gordon damping coefficents   ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     "damphlike" generates the damping coefficients for hydrogen-like
c     damping functional form that go with corresponding powers of r
c
c     one-site scale factors are used for nuclear-electron interactions
c
      subroutine damphlike(r,rorder,alphai,alphak,
     &                 lambdai,lambdak,lambdaik)
      implicit none
      integer i
      integer rorder
      real*8 r
      real*8 alphai,alphak
      real*8 dampi,dampk
      real*8 expdampi,expdampk
      real*8 termi,termk
      real*8 lambdai(*),lambdak(*)
      real*8 lambdaik(*)
c
c
c     compute common factors for damping
c
      dampi = alphai*r
      dampk = alphak*r
      expdampi = exp(-dampi)
      expdampk = exp(-dampk)
c
c     hydrogen-like damping model
c     for charge-charge form see Gordon damping model 1
c
      lambdai(1) = 1.0d0 - (1.0d0 + dampi/2.0d0)*expdampi
      lambdak(1) = 1.0d0 - (1.0d0 + dampk/2.0d0)*expdampk
      lambdai(3) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi**2)*expdampi
      lambdak(3) = 1.0d0 - (1.0d0 + dampk + 0.5d0*dampk**2)*expdampk
      lambdai(5) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi**2 + 
     &     (1.0d0/6.0d0)*dampi**3)*expdampi
      lambdak(5) = 1.0d0 - (1.0d0 + dampk + 0.5d0*dampk**2 + 
     &     (1.0d0/6.0d0)*dampk**3)*expdampk
      if (rorder .ge. 9) then
         lambdai(7) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi**2 +
     &        (1.0d0/6.0d0)*dampi**3 + (1.0d0/30.0d0)*dampi**4)*
     &        expdampi
         lambdak(7) = 1.0d0 -(1.0d0 + dampk + 0.5d0*dampk**2 +
     &        (1.0d0/6.0d0)*dampk**3 + (1.0d0/30.0d0)*dampk**4)*
     &        expdampk
         lambdai(9) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi**2 +
     &        (1.0d0/6.0d0)*dampi**3 + (4.0d0/105.0d0)*dampi**4 +
     &        (1.0d0/210.0d0)*dampi**5)*expdampi
         lambdak(9) = 1.0d0 - (1.0d0 + dampk + 0.5d0*dampk**2 +
     &        (1.0d0/6.0d0)*dampk**3 + (4.0d0/105.0d0)*dampk**4 +
     &        (1.0d0/210.0d0)*dampk**5)*expdampk
      end if
c      if (abs(alphai - alphak) .ge. 0.0001d0) then
      if (alphai .ne. alphak) then
         termi = alphak**2/(alphak**2 - alphai**2)
         termk = alphai**2/(alphai**2 - alphak**2)
         lambdaik(1) = 1.0d0 - (termi**2)*
     &        (1.0d0 + 2.0d0*termk + 0.5d0*dampi)*expdampi -
     &        (termk**2)*(1.0d0 + 2.0d0*termi + 0.5d0*dampk)*
     &        expdampk
         lambdaik(3) = 1.0d0 - (termi**2)*(1.0d0 + dampi + 
     &        0.5d0*dampi**2)*expdampi - 
     &        (termk**2)*(1.0d0 + dampk + 0.5d0*dampk**2)*expdampk
     &        - 2.0d0*(termi**2)*termk*(1.0d0 + dampi)*expdampi
     &        - 2.0d0*(termk**2)*termi*(1.0d0 + dampk)*expdampk
         lambdaik(5) = 1.0d0 - (termi**2)*
     &        (1.0d0 + dampi + 0.5d0*dampi**2 + 
     &        (1.0d0/6.0d0)*dampi**3)*expdampi - (termk**2)*
     &        (1.0d0 + dampk + 0.5d0*dampk**2 + 
     &        (1.0d0/6.0d0)*dampk**3)*expdampk - 
     &        2.0d0*(termi**2)*termk*(1.0 + dampi + 
     &        (1.0d0/3.0d0)*dampi**2)*expdampi - 
     &        2.0d0*(termk**2)*termi*(1.0 + dampk +
     &        (1.0d0/3.0d0)*dampk**2)*expdampk
         if (rorder .ge. 7) then
            lambdaik(7) = 1.0d0 - (termi**2)*
     &           (1.0d0 + dampi + 0.5d0*dampi**2 + 
     &           (1.0d0/6.0d0)*dampi**3 + (1.0d0/30.0d0)*dampi**4)*
     &           expdampi - 
     &           (termk**2)*(1.0d0 + dampk + 0.5d0*dampk**2 + 
     &           (1.0d0/6.0d0)*dampk**3 + (1.0d0/30.0d0)*dampk**4)*
     &           expdampk - 
     &           2.0d0*(termi**2)*termk*(1.0d0 + dampi +
     &           (2.0d0/5.0d0)*dampi**2 + (1.0d0/15.0d0)*dampi**3)*
     &           expdampi - 
     &           2.0d0*(termk**2)*termi*(1.0d0 + dampk + 
     &           (2.0d0/5.0d0)*dampk**2 + (1.0d0/15.0d0)*dampk**3)*
     &           expdampk
         end if
         if (rorder .ge. 9) then
            lambdaik(9) = 1.0d0 - (termi**2)*
     &           (1.0d0 + dampi + 0.5d0*dampi**2 + 
     &           (1.0d0/6.0d0)*dampi**3 + (4.0d0/105.0d0)*dampi**4 +
     &           (1.0d0/210.0d0)*dampi**5)*expdampi - (termk**2)*
     &           (1.0d0 + dampk + 0.5d0*dampk**2 +
     &           (1.0d0/6.0d0)*dampk**3 + (4.0d0/105.0d0)*dampk**4 +
     &           (1.0d0/210.0d0)*dampk**5)*expdampk -
     &           2.0d0*(termi**2)*termk*
     &           (1.0d0 + dampi + (3.0d0/7.0d0)*dampi**2 + 
     &           (2.0d0/21.0d0)*dampi**3 + (1.0d0/105.0d0)*dampi**4)*
     &           expdampi - 
     &           2.0d0*(termk**2)*termi*
     &           (1.0d0 + dampk + (3.0d0/7.0d0)*dampk**2 +
     &           (2.0d0/21.0d0)*dampk**3 + (1.0d0/105.0d0)*dampk**4)*
     &           expdampk
         end if
         if (rorder .ge. 11) then
            lambdaik(11) = 1.0d0 - (termi**2)*
     &           (1.0d0 + dampi + 0.5d0*dampi**2 +
     &           (1.0d0/6.0d0)*dampi**3 + (5.0d0/126.0d0)*dampi**4 +
     &           (2.0d0/315.0d0)*dampi**5 + (1.0d0/1890.0d0)*dampi**6)*
     &           expdampi -
     &           (termk**2)*
     &           (1.0d0+ dampk+ 0.5d0*dampk**2 +
     &           (1.0d0/6.0d0)*dampk**3+ (5.0d0/126.0d0)*dampk**4 +
     &           (2.0d0/315.0d0)*dampk**5 + (1.0d0/1890.0d0)*dampk**6)*
     &           expdampk -
     &           2.0d0*(termi**2)*termk*
     &           (1.0d0 + dampi + (4.0d0/9.0d0)*dampi**2 +
     &           (1.0d0/9.0d0)*dampi**3 + (1.0d0/63.0d0)*dampi**4 +
     &           (1.0d0/945.0d0)*dampi**5)*expdampi - 
     &           2.0d0*(termk**2)*termi*
     &           (1.0d0+ dampk+ (4.0d0/9.0d0)*dampk**2 +
     &           (1.0d0/9.0d0)*dampk**3 + (1.0d0/63.0d0)*dampk**4 +
     &           (1.0d0/945.0d0)*dampk**5)*expdampk 
         end if
      else
         lambdaik(1) = 1.0d0 - (1.0d0 + (11.0d0/16.0d0)*dampi + 
     &        (3.0d0/16.0d0)*dampi**2 + (1.0d0/48.0d0)*dampi**3)
     &        *expdampi
         lambdaik(3) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi**2 + 
     &        (7.0d0/48.0d0)*dampi**3 + (1.0d0/48.0d0)*dampi**4)
     &        *expdampi
         lambdaik(5) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi**2 + 
     &        (1.0d0/6.0d0)*dampi**3 + (1.0d0/24.0d0)*dampi**4 +
     &        (1.0d0/144.0d0)*dampi**5)*expdampi
         if (rorder .ge. 7) then
            lambdaik(7) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi**2 + 
     &           (1.0d0/6.0d0)*dampi**3 + (1.0d0/24.0d0)*dampi**4 + 
     &           (1.0d0/120.0d0)*dampi**5 + (1.0d0/720.0d0)*dampi**6)
     &           *expdampi
         end if
         if (rorder .ge. 9) then
            lambdaik(9) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi**2 + 
     &           (1.0d0/6.0d0)*dampi**3 + (1.0d0/24.0d0)*dampi**4 + 
     &           (1.0d0/120.0d0)*dampi**5 + (1.0d0/720.0d0)*dampi**6+ 
     &           (1.0d0/5040.0d0)*dampi**7)*expdampi
         end if
         if (rorder .ge. 11) then
            lambdaik(11) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi**2 +
     &           (1.0d0/6.0d0)*dampi**3 + (1.0d0/24.0d0)*dampi**4 + 
     &           (1.0d0/120.0d0)*dampi**5 + (1.0d0/720.0d0)*dampi**6 +
     &           (1.0d0/5040.0d0)*dampi**7 + (1.0d0/45360.0d0)*dampi**8)
     &           *expdampi
         end if
      end if
c
c     get dispersion damping terms for gradient
c
      if (rorder .ge. 11) then
         if (alphai .ne. alphak) then
            lambdaik(10) = 0.25d0*(dampi**2)*(r*(alphai**2)*(termi**2) +
     &           4.0d0*alphai*(termi**2)*termk - alphai*termi**2)*
     &           expdampi +
     &           0.25d0*(dampk**2)*(r*(alphak**2)*(termk**2) +
     &           4.0d0*alphak*(termk**2)*termi - alphak*termk**2)*
     &           expdampk 

c            lambdaik(10) = (r**2)*(termi*(alphai**3)*expdampi +
c     &           termk*(alphak**3)*expdampk)/2.0d0
         else
            lambdaik(10) = (dampi**2)*alphai*
     &           (dampi**3 - 3.0d0*dampi - 3.0d0)*
     &           expdampi/96.0d0
         end if
      end if
c     
      return
      end
c
c
c     ######################################################################
c     ##                                                                  ##
c     ##  subroutine damppauli  --  generate gordon damping coefficents   ##
c     ##                                                                  ##
c     ######################################################################
c
c
c     "damppauli" generates the damping coefficients for hydrogen-like
c     damping functional form that go with corresponding powers of r
c
c     one-site scale factors are used for nuclear-electron interactions
c
      subroutine damppauli(r,r2,rr1,rr3,rr5,rr7,rr9,rr11,rorder,
     &     apauli,apaulk,lambdaik)
      implicit none
      integer rorder
      real*8 r,r2
      real*8 rr1,rr3,rr5
      real*8 rr7,rr9,rr11
      real*8 apauli,apaulk
      real*8 apauli2,apaulk2
      real*8 dampi,dampk
      real*8 expdampi,expdampk
      real*8 term,termi,termk
      real*8 pre
      real*8 diff
      real*8 s,ds,dds,ddds,dddds,ddddds
      real*8 lambdaik(*)
c
      apauli2 = 0.5d0*apauli
      apaulk2 = 0.5d0*apaulk
      dampi = apauli2*r
      dampk = apaulk2*r
      expdampi = exp(-dampi)
      expdampk = exp(-dampk)
c
      diff = abs(apauli - apaulk)
      if (diff .ge. 0.001d0) then
         term = apauli2**2 - apaulk2**2
         pre = 64.0d0*(apauli**3)*(apaulk**3)/(term**4)
c
      s = (dampi - 0.4D1 / term * apaulk2 * apauli2) * expdampk + (dampk
     & + 0.4D1 / term * apaulk2 * apauli2) * expdampi
c
      ds = (apauli2 * apaulk2 * r2 - 0.4D1 * apauli2 * apaulk2 ** 2 / te
     &rm * r - 0.4D1 / term * apaulk2 * apauli2) * expdampk + (apauli2 *
     & apaulk2 * r2 + 0.4D1 / term * apaulk2 * apauli2 + 0.4D1 * apauli2
     & ** 2 * apaulk2 / term * r) * expdampi
c
      dds = (apauli2 * apaulk2 * r2 / 0.3D1 - 0.4D1 * apauli2 * apaulk2 
     &** 2 / term * r + apauli2 * apaulk2 ** 2 * r ** 3 / 0.3D1 - 0.4D1 
     &/ 0.3D1 * apauli2 * apaulk2 ** 3 / term * r2 - 0.4D1 / term * apau
     &lk2 * apauli2) * expdampk + (apauli2 * apaulk2 * r2 / 0.3D1 + 0.4D
     &1 / term * apaulk2 * apauli2 + apauli2 ** 2 * apaulk2 * r ** 3 / 0
     &.3D1 + 0.4D1 * apauli2 ** 2 * apaulk2 / term * r + 0.4D1 / 0.3D1 *
     & apauli2 ** 3 * apaulk2 / term * r2) * expdampi

c
      ddds = (-0.4D1 / term * apaulk2 * apauli2 - 0.4D1 * apauli2 * apau
     &lk2 ** 2 / term * r - 0.8D1 / 0.5D1 * apauli2 * apaulk2 ** 3 / ter
     &m * r2 + apauli2 * apaulk2 * r2 / 0.5D1 + apauli2 * apaulk2 ** 2 *
     & r ** 3 / 0.5D1 + apauli2 * apaulk2 ** 3 * r ** 4 / 0.15D2 - 0.4D1
     & / 0.15D2 * apauli2 * apaulk2 ** 4 / term * r ** 3) * expdampk + (
     &apauli2 ** 3 * apaulk2 * r ** 4 / 0.15D2 + 0.4D1 / 0.15D2 * apauli
     &2 ** 4 * apaulk2 / term * r ** 3 + 0.4D1 / term * apaulk2 * apauli
     &2 + 0.4D1 * apauli2 ** 2 * apaulk2 / term * r + 0.8D1 / 0.5D1 * ap
     &auli2 ** 3 * apaulk2 / term * r2 + apauli2 * apaulk2 * r2 / 0.5D1 
     &+ apauli2 ** 2 * apaulk2 * r ** 3 / 0.5D1) * expdampi
c
      dddds = (-0.12D2 / 0.7D1 * apauli2 * apaulk2 ** 3 / term * r2 - 0.
     &8D1 / 0.21D2 * apauli2 * apaulk2 ** 4 / term * r ** 3 + apauli2 * 
     &apaulk2 * r2 / 0.7D1 + apauli2 * apaulk2 ** 2 * r ** 3 / 0.7D1 - 0
     &.4D1 / term * apaulk2 * apauli2 + 0.2D1 / 0.35D2 * apauli2 * apaul
     &k2 ** 3 * r ** 4 - 0.4D1 / 0.105D3 * apauli2 * apaulk2 ** 5 / term
     & * r ** 4 + apauli2 * apaulk2 ** 4 * r ** 5 / 0.105D3 - 0.4D1 * ap
     &auli2 * apaulk2 ** 2 / term * r) * expdampk + (apauli2 ** 2 * apau
     &lk2 * r ** 3 / 0.7D1 + 0.4D1 / term * apaulk2 * apauli2 + 0.2D1 / 
     &0.35D2 * apauli2 ** 3 * apaulk2 * r ** 4 + 0.4D1 / 0.105D3 * apaul
     &i2 ** 5 * apaulk2 / term * r ** 4 + apauli2 ** 4 * apaulk2 * r ** 
     &5 / 0.105D3 + 0.4D1 * apauli2 ** 2 * apaulk2 / term * r + 0.12D2 /
     & 0.7D1 * apauli2 ** 3 * apaulk2 / term * r2 + 0.8D1 / 0.21D2 * apa
     &uli2 ** 4 * apaulk2 / term * r ** 3 + apauli2 * apaulk2 * r2 / 0.7
     &D1) * expdampi
c
         if (rorder .ge. 11) then
      ddddds = (0.2D1 / 0.189D3 * apauli2 * apaulk2 ** 4 * r ** 5 + apau
     &li2 * apaulk2 ** 5 * r ** 6 / 0.945D3 + apaulk2 ** 3 * apauli2 * r
     & ** 4 / 0.21D2 - 0.4D1 / 0.9D1 / term * apaulk2 ** 4 * apauli2 * r
     & ** 3 - 0.4D1 / 0.945D3 * apauli2 * apaulk2 ** 6 / term * r ** 5 -
     & 0.16D2 / 0.9D1 / term * apaulk2 ** 3 * apauli2 * r2 - 0.4D1 / 0.6
     &3D2 * apauli2 * apaulk2 ** 5 / term * r ** 4 + apaulk2 * apauli2 *
     & r2 / 0.9D1 + apaulk2 ** 2 * apauli2 * r ** 3 / 0.9D1 - 0.4D1 / te
     &rm * apaulk2 * apauli2 - 0.4D1 / term * apaulk2 ** 2 * apauli2 * r
     &) * expdampk + (apaulk2 * apauli2 ** 3 * r ** 4 / 0.21D2 + 0.4D1 /
     & term * apaulk2 * apauli2 ** 2 * r + 0.2D1 / 0.189D3 * apauli2 ** 
     &4 * apaulk2 * r ** 5 + apauli2 ** 5 * apaulk2 * r ** 6 / 0.945D3 +
     & apaulk2 * apauli2 ** 2 * r ** 3 / 0.9D1 + 0.4D1 / term * apaulk2 
     &* apauli2 + 0.4D1 / 0.9D1 / term * apaulk2 * apauli2 ** 4 * r ** 3
     & + apaulk2 * apauli2 * r2 / 0.9D1 + 0.4D1 / 0.63D2 * apauli2 ** 5 
     &* apaulk2 / term * r ** 4 + 0.4D1 / 0.945D3 * apauli2 ** 6 * apaul
     &k2 / term * r ** 5 + 0.16D2 / 0.9D1 / term * apaulk2 * apauli2 ** 
     &3 * r2) * expdampi
         end if
c
      else
         apauli = min(apauli,apaulk)
         apauli2 = min(apauli2,apaulk2)
         dampi = apauli2*r
         expdampi = exp(-dampi)
         pre = (apauli**6)/(apauli2**6)        
c
      s = (r + r2 * apauli2 + r ** 3 * apauli2 ** 2 / 0.3D1) * expdampi
      ds = (r ** 3 * apauli2 ** 2 / 0.3D1 + r ** 4 * apauli2 ** 3 / 0.3D
     &1) * expdampi
      dds = apauli2 ** 4 * expdampi * r ** 5 / 0.9D1
      ddds = apauli2 ** 5 * expdampi * r ** 6 / 0.45D2
      dddds = (apauli2 ** 5 * r ** 6 / 0.315D3 + apauli2 ** 6 * r ** 7 /
     & 0.315D3) * expdampi
         if (rorder .ge. 11) then
      ddddds = (apauli2 ** 5 * r ** 6 / 0.945D3 + apauli2 ** 6 * r ** 7
     & / 0.945D3 + apauli2 ** 7 * r ** 8 / 0.2835D4) * expdampi
         end if
c
      end if
c
c     turn partial derivatives into full derivatives
c
      s = s*rr1
      ds = ds*rr3
      dds = dds*rr5
      ddds = ddds*rr7
      dddds = dddds*rr9
      ddddds = ddddds*rr11
c     
      lambdaik(1) = pre*s**2
      lambdaik(3) = 2.0d0*pre*s*ds
      lambdaik(5) = 2.0d0*pre*(ds*ds + s*dds)
      lambdaik(7) = 2.0d0*pre*(dds*ds + ds*dds + ds*dds + s*ddds)
      lambdaik(9) = 2.0d0*pre*(ddds*ds + dds*dds + dds*dds + ds*ddds + 
     &        dds*dds + ds*ddds + ds*ddds + s*dddds)
      if (rorder .ge. 11) then
         lambdaik(11) = 2.0d0*pre*(dddds*ds + ddds*dds + 
     &        3.0d0*(ddds*dds + 
     &        dds*ddds) + dds*ddds + ds*dddds + 
     &        2.0d0*(dds*ddds + ds*dddds)+
     &        ds*dddds + s*ddddds)
      end if
c
      return
      end
