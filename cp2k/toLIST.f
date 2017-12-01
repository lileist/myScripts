      program main
      implicit none
      character*2 elem1
      integer x(1000)
      character*70 a,inputFile,non0,non1,non2,non3
      integer i,natom,elemO,elemCe

      elemO = 0
      elemCe = 0


!      speci1=120
!      speci2=273

      print*,"Please enter the name of the input file"
      read(*,*),inputFile
      open(3,file=inputFile,readonly)
      open(40,file='list.input')

!calculation part
      print*,"Please enter number of atoms"
      read(*,*),natom

!      read(3,'()')
!      read(3,'()')
!      read(3,'()')
!      read(3,'()')
!      read(3,'()')
      read(3,'()')
      read(3,'()')
      read(3,'()')
      read(3,'()')
      read(3,'()')
      read(3,'()')

      do i=1,natom
         read(3,*),x(i)
  100    format(A6,4x,I3)
         write(40,100),'  LIST',x(i)
      enddo
      END
