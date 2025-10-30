      program test
      implicit double precision (a-h,o-z)
      real*8 x(3,5000),F(3,5000),E(5000)
      real*8 AL(3,3)
      integer iat(5000)
      
      open(10,file="MOVEMENT.orig")
      rewind(10)
      open(12,file="MOVEMENT")
      rewind(12)

      write(6,*) "input num_step"
      read(5,*) num_step

      do istep=1,num_step
      read(10,*) num
      read(10,*)
      read(10,*) AL(1,1),AL(2,1),AL(3,1)
      read(10,*) AL(1,2),AL(2,2),AL(3,2)
      read(10,*) AL(1,3),AL(2,3),AL(3,3)
      read(10,*) 
      Etot=0.d0
      do i=1,num
      read(10,*) iat(i),x(1,i),x(2,i),x(3,i),m1,m2,m3,E(i),
     &   F(1,i),F(2,i),F(3,i)
      Etot=Etot+E(i)
      enddo
      read(10,*) 
    
      write(12,200) num,Etot,Etot,0.d0 
200   format(i5," atoms,Iteration (f) = 0.1E+01, Etot,Ep,Ek (eV)=",
     &  3(E14.6,1x)) 
      write(12,*) "Lattice vector (Angstrom)"
      write(12,"(3(E18.10,1x))") AL(1,1),AL(2,1),AL(3,1)
      write(12,"(3(E18.10,1x))") AL(1,2),AL(2,2),AL(3,2)
      write(12,"(3(E18.10,1x))") AL(1,3),AL(2,3),AL(3,3)
      write(12,*) "Position (normalized) "
      do i=1,num
      write(12,"(i4,2x,3(f14.10,1x),' 1 1 1')") 
     &     iat(i),x(1,i),x(2,i),x(3,i)
      enddo
      write(12,*) "Force (eV/Angstrom)"
      do i=1,num
      write(12,"(i4,2x,3(E14.7,1x))") iat(i), -F(1,i),-F(2,i),-F(3,i)
      enddo
      write(12,"('Atomic-Energy, Etot(eV),E_nonloc(eV),Q_atom:dE(eV)=',
     &    E14.7)") 0.d0
      do i=1,num
      write(12,"(i6,2x,3(E15.8,1x))") iat(i),E(i),0.d0,1.d0 
      enddo
      write(12,*) "-------------------------"

      enddo
      close(10)
      close(12)

      stop
      end

