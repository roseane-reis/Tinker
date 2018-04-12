c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2018  by  Joshua A Rackers    ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  chgpen.f  --  specifics for charge penetration funtion    ##
c     ##                                                            ##
c     ################################################################
c
c
c     maxmono          maximum number of monopole division (2)
c
c     alpha_perm       charge penetration damping parameter
c     val_ele          number of valence electrons
c     monopole         value of core charge and valence charge
c     penetration      charge penetration function
c     num_ele          global rule for setting core and valence charge
c     alphaele         charge penetration damping parameter by atom
c
c
c
      module chgpen
      use sizes
      implicit none
      integer maxmono
      parameter (maxmono=2)
      real*8 alpha_perm(maxclass)
      real*8 val_ele(maxclass)
      real*8, allocatable :: monopole(:,:)
      real*8, allocatable :: alphaele(:)
      character*20 penetration
      character*20 num_ele
      save
      end
