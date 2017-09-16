c
c     dis2scale     factor by which 1-2 dispersion interactions are scaled
c     dis3scale     factor by which 1-3 dispersion interactions are scaled
c     dis4scale     factor by which 1-4 dispersion interactions are scaled
c     dispindex   indexing mode (atom type or class) for dispersion parameters
c
c
      module disppot
      implicit none
      real*8 dis2scale,dis3scale
      real*8 dis4scale,dis5scale
      character*5 dispindex
      save
      end
