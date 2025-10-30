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

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccc
   open(9,file="polar_param.input")
   rewind(9)
   read(9,*) natom_m
   natom=natom_m
   nparam=0
   do ii=1,natom_m
      read(9,*)i1,ipol(1,ii),ipol(2,ii),ipol(3,ii),ipol(4,ii)
      read(9,*) (indp(i),i=1,6)
      do i=1,6
         if(indp(i).gt.nparam) nparam=indp(i)
      enddo
      indpp(1,1,ii)=indp(1)
      indpp(2,2,ii)=indp(2)
      indpp(3,3,ii)=indp(3)
      indpp(1,2,ii)=indp(4)
      indpp(2,1,ii)=indp(4)
      indpp(1,3,ii)=indp(5)
      indpp(3,1,ii)=indp(5)
      indpp(2,3,ii)=indp(6)
      indpp(3,2,ii)=indp(6)
   enddo
   read(9,*) numb
   do ii=1,numb
      read(9,*) indb(1,ii),indb(2,ii),indb(3,ii)  ! indb(3,ii) is theindex of parameter
      if(indb(3,ii).gt.nparam) nparam=indb(3,ii)
   enddo
   read(9,*) ntype_atom
   do ii=1,ntype_atom
      read(9,*) iatom_type(ii),param(1,ii),param(2,ii),param(3,ii),param(4,ii)
   enddo
   close(9)

   write(6,*) "nparam=",nparam
   nump=4*nparam
   allocate(AA(nump,nump))
   allocate(BB(nump))
   allocate(ppp(nump))
   allocate(SS(nump))
   lwork=30*nump+60
   allocate(work(lwork))


   !cccccccccccccccccccccccccccccccccccccccccccccccccccc
   open(10,file="pxyz.out",action="read")
   rewind(10)


   AA=0.d0
   BB=0.d0

   do 1000 iii=1,100000
      read(10,*,iostat=ierr)
      if(ierr.ne.0) exit
      read(10,*) char4,QV
      read(10,*) natom
      read(10,*)
      read(10,*) AL(1,1),AL(2,1),AL(3,1)
      read(10,*) AL(1,2),AL(2,2),AL(3,2)
      read(10,*) AL(1,3),AL(2,3),AL(3,3)
      read(10,*)
      do i=1,natom
         read(10,*) iatom(i),xatom(1,i),xatom(2,i),xatom(3,i)
      enddo

      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do iat=1,natom

         iat1=ipol(1,iat)
         iat2=ipol(2,iat)
         iat3=ipol(3,iat)
         iat4=ipol(4,iat)

         do i=1,3
            dxp1(i)=xatom(i,iat2)-xatom(i,iat1)
            dxp2(i)=xatom(i,iat4)-xatom(i,iat3)
            if(abs(dxp1(i)+1).lt.abs(dxp1(i))) dxp1(i)=dxp1(i)+1
            if(abs(dxp1(i)-1).lt.abs(dxp1(i))) dxp1(i)=dxp1(i)-1
            if(abs(dxp2(i)+1).lt.abs(dxp2(i))) dxp2(i)=dxp2(i)+1
            if(abs(dxp2(i)-1).lt.abs(dxp2(i))) dxp2(i)=dxp2(i)-1
         enddo


         dxpp1=AL(1,1)*dxp1(1)+AL(1,2)*dxp1(2)+AL(1,3)*dxp1(3)
         dypp1=AL(2,1)*dxp1(1)+AL(2,2)*dxp1(2)+AL(2,3)*dxp1(3)
         dzpp1=AL(3,1)*dxp1(1)+AL(3,2)*dxp1(2)+AL(3,3)*dxp1(3)

         dxpp2=AL(1,1)*dxp2(1)+AL(1,2)*dxp2(2)+AL(1,3)*dxp2(3)
         dypp2=AL(2,1)*dxp2(1)+AL(2,2)*dxp2(2)+AL(2,3)*dxp2(3)
         dzpp2=AL(3,1)*dxp2(1)+AL(3,2)*dxp2(2)+AL(3,3)*dxp2(3)

         dxpp3=dypp1*dzpp2-dzpp1*dypp2   ! dpp3=dpp1 x dpp2
         dypp3=dzpp1*dxpp2-dxpp1*dzpp2
         dzpp3=dxpp1*dypp2-dypp1*dxpp2

         d=dsqrt(dxpp1**2+dypp1**2+dzpp1**2)
         dxpp1=dxpp1/d
         dypp1=dypp1/d
         dzpp1=dzpp1/d

         d=dsqrt(dxpp2**2+dypp2**2+dzpp2**2)
         dxpp2=dxpp2/d
         dypp2=dypp2/d
         dzpp2=dzpp2/d

         d=dsqrt(dxpp3**2+dypp3**2+dzpp3**2)
         dxpp3=dxpp3/d
         dypp3=dypp3/d
         dzpp3=dzpp3/d

         !------------------------------------------------

         proj(1,1,iat)=dxpp1
         proj(1,2,iat)=dypp1
         proj(1,3,iat)=dzpp1

         proj(2,1,iat)=dxpp2
         proj(2,2,iat)=dypp2
         proj(2,3,iat)=dzpp2

         proj(3,1,iat)=dxpp3
         proj(3,2,iat)=dypp3
         proj(3,3,iat)=dzpp3

      enddo
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do jj=1,numb
         jj1=indb(1,jj)
         jj2=indb(2,jj)
         dx1=xatom(1,jj1)-xatom(1,jj2)
         dx2=xatom(2,jj1)-xatom(2,jj2)
         dx3=xatom(3,jj1)-xatom(3,jj2)
         dx=AL(1,1)*dx1+AL(1,2)*dx2+AL(1,3)*dx3
         dy=AL(2,1)*dx1+AL(2,2)*dx2+AL(2,3)*dx3
         dz=AL(3,1)*dx1+AL(3,2)*dx2+AL(3,3)*dx3
         d=dsqrt(dx**2+dy**2+dz**2)
         Epp(1,jj)=dx/d
         Epp(2,jj)=dy/d
         Epp(3,jj)=dz/d
      enddo

      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      read(10,*)
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
      read(10,*) char4,char2,E0,dE
      read(10,*) char4,char2,E1,dE
      !      write(6,*) "test",E1,E1-(E0+QV)
      dE=E1-(E0+QV)


      ppp=0.d0

      do iat=1,natom


         do jj1=1,3    ! the 3-direction of a given atom
            do jj2=1,3

               if(indpp(jj1,jj2,iat).ne.0) then

                  ip=(indpp(jj1,jj2,iat)-1)*4

                  do ii1=1,2
                     do ii2=1,2

                        ip=ip+1
                        ! note, for jj1.eq.jj2, there are only three independent parameter.
                        !  We are relying on the regulation to automatically make them the same

                        pp1=proj(jj1,1,iat)*pxyz(1,ii1,iat)+proj(jj1,2,iat)*pxyz(2,ii1,iat)+proj(jj1,3,iat)*pxyz(3,ii1,iat)
                        pp2=proj(jj2,1,iat)*pxyz(1,ii2,iat)+proj(jj2,2,iat)*pxyz(2,ii2,iat)+proj(jj2,3,iat)*pxyz(3,ii2,iat)

                        ppp(ip)=ppp(ip)+pp1*pp2

                     enddo
                  enddo

               endif


            enddo
         enddo

      enddo

      !cccccccccccccccccccccccccccccccccccccccccc

      do jj=1,numb
         jj1=indb(1,jj)
         jj2=indb(2,jj)

         ip=(indb(3,jj)-1)*4
         do ii1=1,2
            pp1=pxyz(1,ii1,jj1)*Epp(1,jj)+pxyz(2,ii1,jj1)*Epp(2,jj)+pxyz(3,ii1,jj1)*Epp(3,jj)
            do ii2=1,2
               pp2=pxyz(1,ii2,jj2)*Epp(1,jj)+pxyz(2,ii2,jj2)*Epp(2,jj)+pxyz(3,ii2,jj2)*Epp(3,jj)
               ip=ip+1
               ppp(ip)=ppp(ip)+pp1*pp2
            enddo
         enddo
      enddo


      do ip1=1,nump
         do ip2=1,nump
            AA(ip1,ip2)=AA(ip1,ip2)+ppp(ip1)*ppp(ip2)
         enddo
      enddo

      do ip1=1,nump
         BB(ip1)=BB(ip1)+dE*ppp(ip1)
      enddo


1000 continue
   write(6,*) "num of cases=",iii
   close(10)

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   do ip1=1,nump
      AA(ip1,ip1)=AA(ip1,ip1)+1.D-8
   enddo

   call dgelss(nump,nump,1,AA,nump,BB, &
      nump,SS,-0.1,irank,work,lwork,info)

   open(25,file="fit.BB")
   rewind(25)
   do ip=1,nump
      write(25,*) ip,BB(ip)
   enddo
   close(25)


   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   open(10,file="pxyz.out",action="read")
   rewind(10)
   open(12,file="test.fit")
   rewind(12)

   do 1001 iii=1,100000
      read(10,*,iostat=ierr)
      if(ierr.ne.0) exit
      read(10,*) char4,QV
      read(10,*) natom
      read(10,*)
      read(10,*) AL(1,1),AL(2,1),AL(3,1)
      read(10,*) AL(1,2),AL(2,2),AL(3,2)
      read(10,*) AL(1,3),AL(2,3),AL(3,3)
      read(10,*)
      do i=1,natom
         read(10,*) iatom(i),xatom(1,i),xatom(2,i),xatom(3,i)
      enddo

      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do iat=1,natom

         iat1=ipol(1,iat)
         iat2=ipol(2,iat)
         iat3=ipol(3,iat)
         iat4=ipol(4,iat)

         do i=1,3
            dxp1(i)=xatom(i,iat2)-xatom(i,iat1)
            dxp2(i)=xatom(i,iat4)-xatom(i,iat3)
            if(abs(dxp1(i)+1).lt.abs(dxp1(i))) dxp1(i)=dxp1(i)+1
            if(abs(dxp1(i)-1).lt.abs(dxp1(i))) dxp1(i)=dxp1(i)-1
            if(abs(dxp2(i)+1).lt.abs(dxp2(i))) dxp2(i)=dxp2(i)+1
            if(abs(dxp2(i)-1).lt.abs(dxp2(i))) dxp2(i)=dxp2(i)-1
         enddo


         dxpp1=AL(1,1)*dxp1(1)+AL(1,2)*dxp1(2)+AL(1,3)*dxp1(3)
         dypp1=AL(2,1)*dxp1(1)+AL(2,2)*dxp1(2)+AL(2,3)*dxp1(3)
         dzpp1=AL(3,1)*dxp1(1)+AL(3,2)*dxp1(2)+AL(3,3)*dxp1(3)

         dxpp2=AL(1,1)*dxp2(1)+AL(1,2)*dxp2(2)+AL(1,3)*dxp2(3)
         dypp2=AL(2,1)*dxp2(1)+AL(2,2)*dxp2(2)+AL(2,3)*dxp2(3)
         dzpp2=AL(3,1)*dxp2(1)+AL(3,2)*dxp2(2)+AL(3,3)*dxp2(3)

         dxpp3=dypp1*dzpp2-dzpp1*dypp2   ! dpp3=dpp1 x dpp2
         dypp3=dzpp1*dxpp2-dxpp1*dzpp2
         dzpp3=dxpp1*dypp2-dypp1*dxpp2

         d=dsqrt(dxpp1**2+dypp1**2+dzpp1**2)
         dxpp1=dxpp1/d
         dypp1=dypp1/d
         dzpp1=dzpp1/d

         d=dsqrt(dxpp2**2+dypp2**2+dzpp2**2)
         dxpp2=dxpp2/d
         dypp2=dypp2/d
         dzpp2=dzpp2/d

         d=dsqrt(dxpp3**2+dypp3**2+dzpp3**2)
         dxpp3=dxpp3/d
         dypp3=dypp3/d
         dzpp3=dzpp3/d

         !------------------------------------------------

         proj(1,1,iat)=dxpp1
         proj(1,2,iat)=dypp1
         proj(1,3,iat)=dzpp1

         proj(2,1,iat)=dxpp2
         proj(2,2,iat)=dypp2
         proj(2,3,iat)=dzpp2

         proj(3,1,iat)=dxpp3
         proj(3,2,iat)=dypp3
         proj(3,3,iat)=dzpp3

      enddo
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do jj=1,numb
         jj1=indb(1,jj)
         jj2=indb(2,jj)
         dx1=xatom(1,jj1)-xatom(1,jj2)
         dx2=xatom(2,jj1)-xatom(2,jj2)
         dx3=xatom(3,jj1)-xatom(3,jj2)
         dx=AL(1,1)*dx1+AL(1,2)*dx2+AL(1,3)*dx3
         dy=AL(2,1)*dx1+AL(2,2)*dx2+AL(2,3)*dx3
         dz=AL(3,1)*dx1+AL(3,2)*dx2+AL(3,3)*dx3
         d=dsqrt(dx**2+dy**2+dz**2)
         Epp(1,jj)=dx/d
         Epp(2,jj)=dy/d
         Epp(3,jj)=dz/d
      enddo

      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      read(10,*)
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
      read(10,*) char4,char2,E0,dE
      read(10,*) char4,char2,E1,dE
      !      write(6,*) "test",E1,E1-(E0+QV)
      dE=E1-(E0+QV)


      ppp=0.d0

      Etest=0.d0

      do iat=1,natom


         do jj1=1,3    ! the 3-direction of a given atom
            do jj2=1,3

               if(indpp(jj1,jj2,iat).ne.0) then

                  ip=(indpp(jj1,jj2,iat)-1)*4

                  do ii1=1,2
                     do ii2=1,2

                        ip=ip+1
                        ! note, for jj1.eq.jj2, there are only three independent parameter.
                        !  We are relying on the regulation to automatically make them the same

                        pp1=proj(jj1,1,iat)*pxyz(1,ii1,iat)+proj(jj1,2,iat)*pxyz(2,ii1,iat)+proj(jj1,3,iat)*pxyz(3,ii1,iat)
                        pp2=proj(jj2,1,iat)*pxyz(1,ii2,iat)+proj(jj2,2,iat)*pxyz(2,ii2,iat)+proj(jj2,3,iat)*pxyz(3,ii2,iat)

                        ppp(ip)=ppp(ip)+pp1*pp2
                        Etest=Etest+pp1*pp2*BB(ip)

                     enddo
                  enddo

               endif


            enddo
         enddo

      enddo

      !cccccccccccccccccccccccccccccccccccccccccc

      do jj=1,numb
         jj1=indb(1,jj)
         jj2=indb(2,jj)

         ip=(indb(3,jj)-1)*4
         do ii1=1,2
            pp1=pxyz(1,ii1,jj1)*Epp(1,jj)+pxyz(2,ii1,jj1)*Epp(2,jj)+pxyz(3,ii1,jj1)*Epp(3,jj)
            do ii2=1,2
               pp2=pxyz(1,ii2,jj2)*Epp(1,jj)+pxyz(2,ii2,jj2)*Epp(2,jj)+pxyz(3,ii2,jj2)*Epp(3,jj)
               ip=ip+1
               Etest=Etest+pp1*pp2*BB(ip)
            enddo
         enddo
      enddo

      write(12,*) dE,Etest


1001 continue
   write(6,*) "num of cases=",iii
   close(10)
   close(12)

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


   stop
end

