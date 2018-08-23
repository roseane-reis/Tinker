c
c     ctindex   indexing mode (atom type or class) for dispersion parameters
c
c
      module xtrapot
      implicit none
      character*5 ctindex
      character*5 transtyp
      integer nct
      integer maxct
      integer, allocatable :: ict(:)
      parameter (maxct = 2000)
      real*8 ct_chg(maxct)
      real*8 alpha_ct(maxct)
      real*8, allocatable :: chgct(:)
      real*8, allocatable :: alphact(:)
      save
      end
