      program main
      implicit none
      character*2 atom_I,atom_II,atom_III,elem1,atom_IIII
      character*2 atom_VI
      character*70 state(1000)
      real*8 x(1000),y(1000),z(1000),z_max,z_min
      real*8 a_x,a_y,a_z
      real*8 b_x,b_y,b_z
      real*8 c_x,c_y,c_z
      character*70 a,inputFile
      integer natom,i,natom_VI
      integer natom_I,natom_II,natom_III,natom_IIII

      print*,"Please enter the name of the input file"
      read(*,*),inputFile
      print*,"Please enter min and max of z in cartision coord"
      read(*,*),z_min,z_max
      open(3,file=inputFile)
      open(40,file='POSCAR')

!calculation part

!read and write lattice parameters
      read(3,'()')
      read(3,'()')

      read(3,*),a_x,a_y,a_z
      read(3,*),b_x,b_y,b_z
      read(3,*),c_x,c_y,c_z

      read(3,*),atom_I,atom_II,atom_III,atom_IIII,atom_VI
      read(3,*),natom_I,natom_II,natom_III,natom_IIII,natom_VI

      print*,"test"
      write(40,'(5(2x,A4))'),atom_I,atom_II,atom_III,atom_IIII,atom_VI
      write(40,*),"    1.00000000"
      write(40,'(3(2x,f14.10))'),a_x,a_y,a_z
      write(40,'(3(2x,f14.10))'),b_x,b_y,b_z
      write(40,'(3(2x,f14.10))'),c_x,c_y,c_z
      write(40,'(5(2x,A4))'),atom_I,atom_II,atom_III,atom_IIII,atom_VI
      write(40,200),natom_I,natom_II,natom_III,natom_IIII,natom_VI
  200  format(5(2x,I4))
!converte to fractional coordinates
      read(3,'()')
!fix the atoms with z between z_min and z_min
      natom = natom_I + natom_II + natom_III + natom_IIII+natom_VI
      do i=1,natom
         read(3,*),x(i),y(i),z(i)
         if(z(i).lt.z_max.and.z(i).gt.z_min)then
            state(i)='F  F  F'
         else
            state(i)='T  T  T'
         endif
      enddo
      
      write(40,'(A18)'),"Selective dynamics"
      write(40,'(A6)'),"Direct"
      do i=1,natom
  100    format(4x,f19.16,2x,f19.16,2x,f19.16,3x,A7)
          write(40,100),x(i),y(i),z(i),state(i)
      enddo
      END
