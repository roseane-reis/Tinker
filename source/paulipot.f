c
c     pr2scale     factor by which 1-2 dispersion interactions are scaled
c     pr3scale     factor by which 1-3 dispersion interactions are scaled
c     pr4scale     factor by which 1-4 dispersion interactions are scaled
c     pauliindex   indexing mode (atom type or class) for dispersion parameters
c
c
      module paulipot
      implicit none
      real*8 pr2scale,pr3scale
      real*8 pr4scale,pr5scale
      character*5 pauliindex
      save
      end
