ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
      subroutine kextra
      use sizes
      use atomid
      use atoms
      use mpole
      use xtrapot
      use inform
      use iounit
      use keys
      use xtrapot
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
c     process keywords containing charge transfer parameters
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
         if (keyword(1:10) .eq. 'CT-CHARGE ') then
            call getnumb (record,k,next)
            if (k.ge.1 .and. k.le.maxclass) then
               c = ct_chg(k)
               string = record(next:240)
               read (string,*,err=10,end=10)  c
 10            continue
               if (header .and. .not.silent) then
                  header = .false.
                  if (ctindex .eq. 'CLASS') then
                     write (iout,20)
 20                  format (/,' Additional CT-Charge Parameters :',
     &                    //,5x,'Atom Class',11x,'C6')
                  else
                     write (iout,30)
 30                  format (/,' Additional CT-Charge Parameters :',
     &                    //,5x,'Atom Type',12x,'C6')
                  end if
               end if
               ct_chg(k) = c
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
c     read in ct damping parameter
c
         if (keyword(1:9) .eq. 'ALPHA-CT ') then
            call getnumb (record,k,next)
            if (k.ge.1 .and. k.le.maxclass) then
               a = alpha_ct(k)
               string = record(next:240)
               read (string,*,err=60,end=60)  a
 60            continue
               if (aheader .and. .not.silent) then
                  aheader = .false.
                  if (ctindex .eq. 'CLASS') then
                     write (iout,70)
 70                  format (/,' Additional CT Damping '
     &                   'Parameters: ', //,5x,'Atom Class',10x,'Alpha')
                  else
                     write (iout,80)
 80                  format (/,' Additional CT Damping ',
     &                   'Parameters: ', //,5x,'Atom Type',11x,'Alpha')
                  end if
               end if
               alpha_ct(k) = a
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
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(ict))  deallocate (ict)
      if (allocated(chgct)) deallocate (chgct)
      if (allocated(alphact)) deallocate (alphact)
      allocate (ict(n))
      allocate (chgct(n))
      allocate (alphact(n))
c
c     assign prefactor and damping (alpha) parameters 
c     
      nct = 0
      do i = 1, n
         ic = class(i)
         it = type(i)
         if (ctindex .eq. 'CLASS') then
            tmp = ct_chg(ic)
            atmp = alpha_ct(ic)
         else
            tmp = ct_chg(it)
            atmp = alpha_ct(it)
         end if
         if (tmp .ne. 0.0d0) then 
            nct = nct + 1
            if (ctindex .eq. 'CLASS') then
               chgct(nct) = tmp
               alphact(nct) = atmp
            else
               chgct(nct) = tmp
               alphact(nct) = atmp
            end if
            ict(nct) = i
         end if
      end do
c
c
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
         if (keyword(1:14) .eq. 'TRANSFER-TYPE ') then
c            call getword (record,penetration,next)
            read (string,*,err=110,end=110)  transtyp
         end if
 110     continue
      end do
c
      return
      end
