      program main
      implicit none
      character*2 elem1
      real*8 x(1000),y(1000),z(1000)
      character*70 a,inputFile,non0,non1,non2,non3
      integer i,natom,elemO,elemCe

      elemO = 0
      elemCe = 0


!      speci1=120
!      speci2=273

      print*,"Please enter the name of the input file"
      read(*,*),inputFile
      open(3,file=inputFile,readonly)
      open(40,file='geometry.xyz')

!calculation part
      print*,"Please enter number of atoms"
      read(*,*),natom

      read(3,'()')
      read(3,'()')
      read(3,'()')
      read(3,'()')
      read(3,'()')

      write(40,'(I4)'),natom
      write(40,'()')
      do i=1,natom
         read(3,*),non0,x(i),y(i),z(i),non1,non2,non3,elem1
  100    format(A2,4x,f15.9,2x,f15.9,2x,f15.9,3x)
         write(40,100),elem1,x(i),y(i),z(i)
      enddo
      END
