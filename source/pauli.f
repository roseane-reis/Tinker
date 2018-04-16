c
c
c
      module pauli
      implicit none
      integer npauli
      integer maxpauli
      integer, allocatable :: ipauli(:)
      parameter (maxpauli = 2000)
      real*8 prepauli(maxpauli)
      real*8 alpha_pauli(maxpauli)
      real*8 pauli_valence(maxpauli)
      real*8, allocatable :: overpauli(:)
      real*8, allocatable :: alphapauli(:)
      real*8, allocatable :: monopauli(:)
      save
      end
