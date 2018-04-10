ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
      subroutine kdisp
      use sizes
      use atomid
      use atoms
      use disp
      use disppot
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
      real*8 tmp
      real*8 atmp
      logical header
      logical aheader
      character*16 blank
      character*20 keyword
      character*240 record
      character*240 string
c
c     process keywords containing dispersion parameters
c
      blank = '        '
      header = .true.
      aheader = .true.
      id = 0
      dispdamp = .false.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:15) .eq. 'DISPERSION-DAMP') then
            dispdamp = .true.
         end if
         if (keyword(1:5) .eq. 'CSIX ') then
            call getnumb (record,k,next)
            if (k.ge.1 .and. k.le.maxclass) then
               c = dispsix(k)
               string = record(next:240)
               read (string,*,err=10,end=10)  c
 10            continue
               if (header .and. .not.silent) then
                  header = .false.
                  if (dispindex .eq. 'CLASS') then
                     write (iout,20)
 20                  format (/,' Additional C6 Dispersion Parameters :',
     &                    //,5x,'Atom Class',11x,'C6')
                  else
                     write (iout,30)
 30                  format (/,' Additional C6 Dispersion Parameters :',
     &                    //,5x,'Atom Type',12x,'C6')
                  end if
               end if
               dispsix(k) = c
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
c     read in dispersion damping parameter
c
         if (keyword(1:17) .eq. 'ALPHA-DISPERSION ') then
            call getnumb (record,k,next)
            if (k.ge.1 .and. k.le.maxclass) then
               a = dispalpha(k)
               string = record(next:240)
               read (string,*,err=60,end=60)  a
 60            continue
               if (aheader .and. .not.silent) then
                  aheader = .false.
                  if (dispindex .eq. 'CLASS') then
                     write (iout,70)
 70                  format (/,' Additional Dispersion Alpha '
     &                   'Parameters: ', //,5x,'Atom Class',10x,'Alpha')
                  else
                     write (iout,80)
 80                  format (/,' Additional Dispersion Alpha ',
     &                   'Parameters: ', //,5x,'Atom Type',11x,'Alpha')
                  end if
               end if
               dispalpha(k) = a
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
      if (allocated(idisp))  deallocate (idisp)
      if (allocated(csix)) deallocate (csix)
      if (allocated(adisp)) deallocate (adisp)
      allocate (idisp(n))
      allocate (csix(n))
      allocate (adisp(n))
c
c     zero out dispersion
c
      do i = 1, n
         csix(i) = 0.0d0
         adisp(i) = 0.0d0
      end do
c
c     assign c6 and dispersion damping (alpha) parameters 
c     
      ndisp = 0
      do i = 1, n
         ic = class(i)
         it = type(i)
         if (dispindex .eq. 'CLASS') then
            tmp = dispsix(ic)
            atmp = dispalpha(ic)
         else
            tmp = dispsix(it)
            atmp = dispalpha(it)
         end if
         if (tmp .ne. 0.0d0) then 
            ndisp = ndisp + 1
            if (dispindex .eq. 'CLASS') then
               csix(ndisp) = tmp
               adisp(ndisp) = atmp
            else
               csix(ndisp) = tmp
               adisp(ndisp) = atmp
            end if
            idisp(ndisp) = i
         end if
      end do
      return
      end
         
