      program main
      real*8 dis,sump
      character*2 empt1,empt2
      integer*8 numb,k
      real*8 v(500)
      real*8 pot(50000,10000)
      
      nx=48
      ny=72
      nxy=10
      nz=400
      nline=INT(nxy/5)
      remain=MOD(nxy,5)
      do i=1,nz
        v(i)=0d0
      enddo
!      speci1=120
!      speci2=273

      open(3,file='cal.cube',readonly)
      open(40,file='pot.dat')
      do k=1,nz      
       i=1
       do j=1,nline
         if(remain>=1.AND.j==nline)then
            read(3,*),(pot(i,k),i=i,i+remain-1)
         else
         read(3,*),pot(i,k),pot(i+1,k),pot(i+2,k),pot(i+3,k),pot(i+4,k)
         i=i+5
         endif
       enddo
      enddo
      do k=1,nz
        do i=1,nxy
         v(k)=v(k)+pot(i,k)
        enddo
         v(k)=v(k)/nxy
        write(40,'(f8.5,4x,f10.5)'),k*0.0715,v(k)
      enddo
      close (3)
 
      end

