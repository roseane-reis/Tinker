ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
      subroutine kpauli
      use sizes
      use atomid
      use atoms
      use pauli
      use mpole
      use paulipot
      use inform
      use iounit
      use keys
      implicit none
      integer i,k
      integer id
      integer ic,it
      integer next
      real*8 c
      real*8 a
      real*8 v
      real*8 tmp
      real*8 atmp
      real*8 vtmp
      real*8 zi,ci
      logical header
      logical aheader
      logical vheader
      character*16 blank
      character*20 keyword
      character*240 record
      character*240 string
c
c     process keywords containing pauli repulsion parameters
c
      blank = '        '
      header = .true.
      aheader = .true.
      vheader = .true.
      id = 0
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:11) .eq. 'PAULI-SIZE ') then
            call getnumb (record,k,next)
            if (k.ge.1 .and. k.le.maxclass) then
               c = prepauli(k)
               string = record(next:240)
               read (string,*,err=10,end=10)  c
 10            continue
               if (header .and. .not.silent) then
                  header = .false.
                  if (pauliindex .eq. 'CLASS') then
                     write (iout,20)
 20                  format (/,' Additional Pauli Size Parameters :',
     &                    //,5x,'Atom Class',11x,'C6')
                  else
                     write (iout,30)
 30                  format (/,' Additional Pauli Size Parameters :',
     &                    //,5x,'Atom Type',12x,'C6')
                  end if
               end if
               prepauli(k) = c
               if (.not. silent) then
                  write (iout,40)  k,c
 40               format (4x,i6,8x,f12.4)
               end if
               id = id + 1
            else if (k .gt. maxclass) then
               write (iout,50)  maxclass
 50            format (/,' KVDW  --  Only Atom Classes through',i4,
     &              ' are Allowed')
               abort = .true.
            end if
         end if
c
c     read in pauli damping parameter
c
         if (keyword(1:12) .eq. 'ALPHA-PAULI ') then
            call getnumb (record,k,next)
            if (k.ge.1 .and. k.le.maxclass) then
               a = alpha_pauli(k)
               string = record(next:240)
               read (string,*,err=60,end=60)  a
 60            continue
               if (aheader .and. .not.silent) then
                  aheader = .false.
                  if (pauliindex .eq. 'CLASS') then
                     write (iout,70)
 70                  format (/,' Additional Pauli Damping '
     &                   'Parameters: ', //,5x,'Atom Class',10x,'Alpha')
                  else
                     write (iout,80)
 80                  format (/,' Additional Pauli Damping ',
     &                   'Parameters: ', //,5x,'Atom Type',11x,'Alpha')
                  end if
               end if
               alpha_pauli(k) = a
               if (.not. silent) then
                  write (iout,90)  k,a
 90               format (4x,i6,8x,f12.4)
               end if
               id = id + 1
            else if (k .gt. maxclass) then
               write (iout,100)  maxclass
 100            format (/,' KVDW  --  Only Atom Classes through',i4,
     &              ' are Allowed')
               abort = .true.
            end if
         end if
c
c     read in number of valence electrons for pauli repulsion
c
         if (keyword(1:18) .eq. 'PAULI-VALENCE-ELE ') then
            call getnumb (record,k,next)
            if (k.ge.1 .and. k.le.maxclass) then
               v = pauli_valence(k)
               string = record(next:240)
               read (string,*,err=110,end=110)  v
 110           continue
               if (vheader .and. .not.silent) then
                  vheader = .false.
                  if (pauliindex .eq. 'CLASS') then
                     write (iout,120)
 120                 format (/,' Additional Pauli Valence Electron '
     &                    'Parameters: ',//,5x,'Atom Class',10x,'Val')
                  else
                     write (iout,130)
 130                 format (/,' Additional Pauli Valence Electron ',
     &                    'Parameters: ', //,5x,'Atom Type',11x,'Val')
                  end if
               end if
               pauli_valence(k) = v
               if (.not. silent) then
                  write (iout,140)  k,v
 140              format (4x,i6,8x,f12.4)
               end if
               id = id + 1
            else if (k .gt. maxclass) then
               write (iout,150)  maxclass
 150           format (/,' KVDW  --  Only Atom Classes through',i4,
     &              ' are Allowed')
               abort = .true.
            end if
         end if
      end do
c
c     check for keywords indicating a global rule for assigning number of electrons
c
c     for now this feature is left out
c
c      do i = 1, nkey
c         next = 1
c         record = keyline(i)
c         call gettext (record,keyword,next)
c         call upcase (keyword)
c         string = record(next:120)
c         if (keyword(1:9) .eq. 'PAULI-CORE-ELE ') then
c            read (string,*,err=160,end=160)  pauli_num_ele
c         end if
c 160     continue
c      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(ipauli))  deallocate (ipauli)
      if (allocated(overpauli)) deallocate (overpauli)
      if (allocated(alphapauli)) deallocate (alphapauli)
      if (allocated(monopauli)) deallocate (monopauli)
      allocate (ipauli(n))
      allocate (overpauli(n))
      allocate (alphapauli(n))
      allocate (monopauli(n))
c
c     zero out pauli parameters
c
      do i = 1, n
         overpauli(i) = 0.0d0
         alphapauli(i) = 0.0d0
         monopauli(i) = 0.0d0
      end do
c
c     set default number of valence electrons to valence + partial charge
c
      do i = 1, n
         ci = pole(1,i)
         zi = dble(atomic(i))
         if (atomic(i) .gt. 2)  zi = zi - 2.0d0
         if (atomic(i) .gt. 10)  zi = zi - 8.0d0
         if (atomic(i) .gt. 18)  zi = zi - 8.0d0
         if (atomic(i) .gt. 20)  zi = zi - 10.0d0
         if (zi .eq. ci) then
            if (atomic(i) .lt. 10) then
               zi = 2.0d0 + ci
            else
               zi = 8.0d0 + ci
            end if
         end if
c
c     positive convention
c
         if (pauli_valence(class(i)) .eq. 0.0d0) then 
            pauli_valence(class(i)) = -(ci - zi)
         end if
      end do
c
c     assign prefactor and damping (alpha) parameters 
c     
      npauli = 0
      do i = 1, n
         ic = class(i)
         it = type(i)
         if (pauliindex .eq. 'CLASS') then
            tmp = prepauli(ic)
            atmp = alpha_pauli(ic)
            vtmp = pauli_valence(ic)
         else
            tmp = prepauli(it)
            atmp = alpha_pauli(it)
            vtmp = pauli_valence(it)
         end if
         if (tmp .ne. 0.0d0) then 
            npauli = npauli + 1
            if (pauliindex .eq. 'CLASS') then
               overpauli(npauli) = tmp
               alphapauli(npauli) = atmp
c     make sure monopauli is negative
               monopauli(npauli) = -vtmp
            else
               overpauli(npauli) = tmp
               alphapauli(npauli) = atmp
               monopauli(npauli) = -vtmp
            end if
            ipauli(npauli) = i
         end if
c         print *,"monopauli",i,atomic(i),monopauli(i),alpha_pauli(35)
      end do
      return
      end
