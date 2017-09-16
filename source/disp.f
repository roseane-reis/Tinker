c
c
c
      module disp
      implicit none
      integer ndisp
      integer maxdisp
      integer, allocatable :: idisp(:)
      parameter (maxdisp = 2000)
      real*8 dispsix(maxdisp)
      real*8, allocatable :: csix(:)
      save
      end
