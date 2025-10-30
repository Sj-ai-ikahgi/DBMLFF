program gen_config2
   implicit double precision (a-h,o-z)
   real*8 AL(3,3),AL2(3,3)
   real*8 xatom(3,200)
   real*8 xp(3,10),strength(10),wcore(10)
   integer iatom(200),iatom2(200)
   real*8 px_at(200),py_at(200),pz_at(200)
   real*8,allocatable,dimension (:,:,:) :: vr,rho
   real*8,allocatable,dimension (:) :: vr_tmp
   integer iatom_type(40)
   integer itype_atom(200)
   real*8 param(4,40)

   character*40 filename

   real*8 dw, wfactor


   call getarg(1,filename)
   read(filename,*) num_pt

   write(6,*) "num_pt=",num_pt

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

   !cccccccccccccccccccccccccccccccccccccccccccccccccc
   open(9,file="polar_param.input")
   rewind(9)
   read(9,*) natom_m
   if(natom_m.ne.natom) then
      write(6,*) "natom_m.ne.natom,stop"
      stop
   endif

   do ii=1,natom_m
      read(9,*)
      read(9,*)
   enddo
   read(9,*) numb
   do ii=1,numb
      read(9,*)
   enddo
   read(9,*) ntype_atom
   do ii=1,ntype_atom
      read(9,*) iatom_type(ii),param(1,ii),param(2,ii),param(3,ii),param(4,ii)
   enddo
   close(9)

   do i=1,natom
      itype=0
      do ii=1,ntype_atom
         if(iatom(i).eq.iatom_type(ii)) itype=ii
      enddo
      if(itype.eq.0) then
         write(6,*) "iatom type not found in polar_param", iatom(i)
         stop
      endif
      itype_atom(i)=itype
   enddo

   !cccccccccccccccccccccccccccccccccccccccccccccccccc

   open(13,file="OUT.RHO0",form="unformatted")
   rewind(13)
   read(13) n1,n2,n3,nnodes
   read(13) AL2

   nr=n1*n2*n3
   nr_n=nr/nnodes
   allocate(vr_tmp(nr_n))
   allocate(rho(n1,n2,n3))

   do iread=1,nnodes
      read(13) vr_tmp

      do ii=1,nr_n
         jj=ii+(iread-1)*nr_n
         i=(jj-1)/(n2*n3)+1
         j=(jj-1-(i-1)*n2*n3)/n3+1
         k=jj-(i-1)*n2*n3-(j-1)*n3
         rho(i,j,k)=vr_tmp(ii)
      enddo
   enddo
   close(13)
!cccccccccccccccccccccccccccccccccccccccccccccccccc

   open(12,file="point.all")
   rewind(12)
   do ii=1,num_pt
      read(12,*) npt
      do jj=1,npt
         read(12,*) strength(jj),wcore(jj),xp(1,jj),xp(2,jj),xp(3,jj)
      enddo
   enddo
   close(12)

   ! open(24,file="polarization.out",position="append")
   ! write(24,*) "-----------------------------"
   ! !       if(npt.eq.1) then
   ! write(24,"(i4,5(2(f10.4,1x),1x,3(f14.10,1x),1x))")  &
   !    npt,(strength(jj),wcore(jj),xp(1,jj),xp(2,jj),xp(3,jj),jj=1,npt)
   ! !       endif
   ! !       if(npt.eq.2) then
   ! !       write(24,"(i4,3(f14.10,1x))") npt,xp(1,1),xp(2,1),xp(3,1)
   ! !       endif

   allocate(vr(n1,n2,n3))
   vr=0.d0

   AL=AL/0.529177

   do jj=1,npt
      do k=1,n3
         do j=1,n2
            do i=1,n1
               dx1=(i-1.d0)/n1-xp(1,jj)
               dx2=(j-1.d0)/n2-xp(2,jj)
               dx3=(k-1.d0)/n3-xp(3,jj)
               dx=AL(1,1)*dx1+AL(1,2)*dx2+AL(1,3)*dx3
               dy=AL(2,1)*dx1+AL(2,2)*dx2+AL(2,3)*dx3
               dz=AL(3,1)*dx1+AL(3,2)*dx2+AL(3,3)*dx3
               d=dsqrt(dx**2+dy**2+dz**2)
               fact=erf(d)
              !       vr(i,j,k)=vr(i,j,k)-(fact/(d+0.01)-4.d0*exp(-d**2/0.75d0**2))*strength(jj)
              !   (3.d0,0.95) is the Li+ potential fitting.
               vr(i,j,k)=vr(i,j,k)+(-fact/(d+0.01)+wcore(jj)*exp(-d**2/0.95d0**2))*strength(jj)
            enddo
         enddo
      enddo
   enddo
       !       vr=0.5*vr
       !       vr=1.d0*vr

   !cccccccccccccccccccccccccccccccccccccccccccccccccc
   open(13,file="IN.VEXT",form="unformatted")
   rewind(13)
   write(13) n1,n2,n3,nnodes
   write(13) AL2

   nr=n1*n2*n3
   nr_n=nr/nnodes

   do iread=1,nnodes
      do ii=1,nr_n
         jj=ii+(iread-1)*nr_n
         i=(jj-1)/(n2*n3)+1
         j=(jj-1-(i-1)*n2*n3)/n3+1
         k=jj-(i-1)*n2*n3-(j-1)*n3
         vr_tmp(ii)=vr(i,j,k)
      enddo
      write(13) vr_tmp
   enddo
   close(13)
   !ccccccccccccccccccccccccccccccccccccccccccccccccccccc
   vol=AL(1,1)*(AL(2,2)*AL(3,3)-AL(3,2)*AL(2,3))+ &
      AL(2,1)*(AL(3,2)*AL(1,3)-AL(1,2)*AL(3,3))+ &
      AL(3,1)*(AL(1,2)*AL(2,3)-AL(2,2)*AL(1,3))
   vol=dabs(vol)
   vol_n=vol/(n1*n2*n3)
   write(6,*) "vol_n=",vol_n

   Q=0.d0
   QV=0.d0
   do k=1,n3
      do j=1,n2
         do i=1,n1
            Q=Q+rho(i,j,k)
            QV=QV+rho(i,j,k)*vr(i,j,k)
         enddo
      enddo
   enddo
   Q=Q*vol/(n1*n2*n3)
   QV=QV*vol/(n1*n2*n3)
   write(6,*) "Q,QV(eV)=",Q,QV*27.21138602
   !ccccccccccccccccccccccccccccccccccccccccccccccccccccc
   E_IV=0.d0
   do iat=1,natom
      if(iatom(iat).eq.8) ch=6
      if(iatom(iat).eq.6) ch=4
      if(iatom(iat).eq.1) ch=1
      if(iatom(iat).eq.15) ch=5
      if(iatom(iat).eq.9) ch=7
      x1=xatom(1,iat)
      x2=xatom(2,iat)
      x3=xatom(3,iat)
      i1=x1*n1+1
      i2=x2*n2+1
      i3=x3*n3+1
      f1=i1-x1*n1
      f2=i2-x2*n2
      f3=i3-x3*n3
      vt000=vr(i1,i2,i3)
      vt100=vr(i1+1,i2,i3)
      vt010=vr(i1,i2+1,i3)
      vt001=vr(i1,i2,i3+1)
      vt110=vr(i1+1,i2+1,i3)
      vt101=vr(i1+1,i2,i3+1)
      vt011=vr(i1,i2+1,i3+1)
      vt111=vr(i1+1,i2+1,i3+1)
      v_av=vt000*f1*f2*f3+vt100*(1-f1)*f2*f3+  &
         vt010*f1*(1-f2)*f3+vt001*f1*f2*(1-f3)+ &
         vt110*(1-f1)*(1-f2)*f3+vt101*(1-f1)*f2*(1-f3)+  &
         vt011*f1*(1-f2)*(1-f3)+vt111*(1-f1)*(1-f2)*(1-f3)
      E_IV=E_IV-ch*v_av
   enddo

   ! Note, in PWmat, the ion-Vext energy is represented by the point charge
   ! in Vr, not the v_ion(pseudo) and Vr. So, the use of E_IV is correct.
   ! Rho is the DFT charge density
   write(6,*) "E_IV=",E_IV*27.21138602
   write(6,*) "(rho+I)*vr(eV)=",(QV+E_IV)*27.21138602
   ! write(24,*) "(rho+I)*vr(eV)=",(QV+E_IV)*27.21138602
   !ccccccccccccccccccccccccccccccccccccccccccccccccccccc

   wfactor = 1.5
   dw = 0.5*wfactor

   iat2=0
   do iat=1,natom
      x10=xatom(1,iat)
      x20=xatom(2,iat)
      x30=xatom(3,iat)

      !cccc the charge density fitting result from neutral atoms
      itype=itype_atom(iat)
      a1=param(1,itype)
      dw1=param(2,itype)
      a2=param(3,itype)
      dw2=param(4,itype)

      dw1 = dw1*wfactor
      dw2 = dw2*wfactor

      px=0.d0
      py=0.d0
      pz=0.d0
      px2=0.d0
      py2=0.d0
      pz2=0.d0

      do k=1,n3
         do j=1,n2
            do i=1,n1

               dx1=(i-1.d0)/n1-x10
               dx2=(j-1.d0)/n2-x20
               dx3=(k-1.d0)/n3-x30
               dx=AL(1,1)*dx1+AL(1,2)*dx2+AL(1,3)*dx3
               dy=AL(2,1)*dx1+AL(2,2)*dx2+AL(2,3)*dx3
               dz=AL(3,1)*dx1+AL(3,2)*dx2+AL(3,3)*dx3
               d=dsqrt(dx**2+dy**2+dz**2)
               w=exp(-d**2/dw**2)

               d_tmp=dsqrt(d**2+0.3**2)

               w2=a1*exp(-d_tmp/dw1)*(1.d0-a2*exp(-(d/dw2)**2))

               px=px+w*dx/d_tmp*vr(i,j,k)
               py=py+w*dy/d_tmp*vr(i,j,k)
               pz=pz+w*dz/d_tmp*vr(i,j,k)
               px2=px2+w2*dx/d_tmp*vr(i,j,k)
               py2=py2+w2*dy/d_tmp*vr(i,j,k)
               pz2=pz2+w2*dz/d_tmp*vr(i,j,k)
            enddo
         enddo
      enddo

      iat2=iat2+1
      px_at(iat2)=px*vol_n
      py_at(iat2)=py*vol_n
      pz_at(iat2)=pz*vol_n
      iatom2(iat2)=iatom(iat)
      iat2=iat2+1
      px_at(iat2)=px2*vol_n
      py_at(iat2)=py2*vol_n
      pz_at(iat2)=pz2*vol_n
      iatom2(iat2)=iatom(iat)
   enddo
   !ccccccccccccccccccccccccccccccccccccccccccccccccccccc

   open(10,file="pxyz.out.test")
   rewind(10)
   write(10,*) 2*natom
   ! write(24,*) 2*natom
   do iat=1,2*natom
      write(10,"(i4,2x,3(E16.8,1x))") iatom2(iat),px_at(iat),py_at(iat),pz_at(iat)
      ! write(24,"(i4,2x,3(E16.8,1x))") iatom2(iat),px_at(iat),py_at(iat),pz_at(iat)
   enddo
   close(10)
   ! close(24)

   !ccccccccccccccccccccccccccccccccccccccccccccccccccccc

   stop
end














