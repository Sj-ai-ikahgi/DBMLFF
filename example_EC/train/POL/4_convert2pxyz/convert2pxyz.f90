program fit_direct

      implicit double precision (a-h,o-z)


      real*8 AL(3,3),vol,vol_n
      real*8 ALm(3,3)
      real*8 ALtt(3,3)
      real*8 xc_cent(3)
      integer icorner(3)
      real*8 Etot
      real*8 pi

      real*8 Epp(3,100)

      real*8 Ecut2,fact_kin2,Rbox
      integer n1,n2,n3       ! box for the whole AL
      integer nm1_all(1000),nm2_all(1000),nm3_all(1000)    !

      integer iatom(200),iatom2(200)
      real*8 xatom(3,200)
      integer ind(2,200)
      integer indp(6)
      integer indpp(3,3,200)


      real*8 proj(3,3,1000),projI(3,3,1000)
      integer ipol(4,1000)
      integer indb(3,200)

      real*8 pxyz(3,2,1000)
      real*8 dxp1(3),dxp2(3),dxp3(3)
      real*8 pp11(3),pp22(3)

      integer iflag_mol(1000)
      character*40 f_xatom,char4,char2,txt1,txt2

      real*8,allocatable,dimension (:,:) :: AA
      real*8,allocatable,dimension (:) :: BB,ppp,SS,work

      integer iatom_type(20)
      real*8 param(4,20)

      integer tnum

      !cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      open(10,file="xatom0.config")
      rewind(10)
      read(10,*) natom
      read(10,*)
      read(10,*) AL(1,1),AL(2,1),AL(3,1)
      read(10,*) AL(1,2),AL(2,2),AL(3,2)
      read(10,*) AL(1,3),AL(2,3),AL(3,3)
      read(10,*)
      do i=1,natom
      read(10,*) iatom(i),xatom(1,i),xatom(2,i),xatom(3,i)
      enddo
      close(10)


      open(23,file="Etot0")
      rewind(23)
      read(23,*) txt1,txt2,E0,dE0
      close(23)
      !cccccccccccccccccccccccccccccccccccccccccccccccccccc
      open(10,file="polarization.out",action="read")
      open(12,file="pxyz.outC",action="write")
      rewind(10)
      rewind(12)

      do 1000 iii=1,100000
            read(10,*,iostat=ierr)
            if(ierr.ne.0) exit
            read(10,*) npoint,strength,wcore,xx1,xx2,xx3
            read(10,*) char4,QV
            read(10,*) natom2
            if(natom2.ne.natom*2) then
                  write(6,*) "natom2.ne.natom*2,stop",natom2,natom
                  stop
            endif
            do i=1,natom
                  do ii=1,2
                        read(10,*) iatom2(i),pxyz(1,ii,i),pxyz(2,ii,i),pxyz(3,ii,i)
                  enddo
            enddo
            read(10,*) char4,char2,E1,dE


            write(12,"(2(f10.5,1x),3(f15.10,1x))") strength,wcore,xx1,xx2,xx3
            if(ierr.ne.0) exit
            write(12,"('(rho+I)*vr(eV)=', f25.14)") QV

            write(12,*) natom
            write(12,*) "Lattice vector"
            write(12,"(3(f25.14,1x))") AL(1,1),AL(2,1),AL(3,1)
            write(12,"(3(f25.14,1x))") AL(1,2),AL(2,2),AL(3,2)
            write(12,"(3(f25.14,1x))") AL(1,3),AL(2,3),AL(3,3)
            write(12,*) "Position"
            do i=1,natom
                  write(12,"(i5,3(f12.9,1x),' 1 1 1')") iatom(i),xatom(1,i),xatom(2,i),xatom(3,i)
            enddo
            write(12,*) "----------------------"
            !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            write(12,*) natom2
            do i=1,natom
                  do ii=1,2
                  write(12,"(i4,3(E20.10,1x))") iatom2(i),pxyz(1,ii,i),pxyz(2,ii,i),pxyz(3,ii,i)
                  enddo
            enddo
            write(12,"('E_tot(eV)  =', 2(E25.15,1x))") E0,dE0
            write(12,"('E_tot(eV)  =', 2(E25.15,1x))") E1,dE


      1000    continue

      close(10)
      close(12)

      stop
end
