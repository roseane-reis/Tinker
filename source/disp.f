c
c
c
      module disp
      implicit none
      logical dispdamp
      integer ndisp
      integer maxdisp
      integer, allocatable :: idisp(:)
      parameter (maxdisp = 2000)
      real*8 dispsix(maxdisp)
      real*8 dispalpha(maxdisp)
      real*8, allocatable :: csix(:)
      real*8, allocatable :: adisp(:)
      save
      end
