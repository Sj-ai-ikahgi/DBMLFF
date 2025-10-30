program  vion_Coulomb_fitV33
   implicit double precision (a-h,o-z)

   parameter (matom=2000)

   real*8 ALI(3,3),AL(3,3),vol,vins,AL_f(3,3)
   complex*16 cci,cc,c_zero,c_one
   real*8 xatom(3,matom)
   real*8 x11(10*matom),y11(10*matom),z11(10*matom)
   integer iatom(matom)
   integer icent(10,10*matom),nat_cent(10*matom)
   integer itype_cent(matom)
   real*8 w_cent(10,10*matom)
   real*8  Zcent(10*matom)
   real*8 xatom_cent(3,10*matom)
   real*8 xcent(3),sum_at(matom)


   real*8,allocatable,dimension(:,:) :: AA
   real*8,allocatable,dimension(:) :: BB
   real*8,allocatable,dimension(:) :: SS,work
   real*8,allocatable,dimension(:,:,:) :: box,box2
   real*8,allocatable,dimension(:) :: Q_type
   real*8,allocatable,dimension(:,:,:,:) :: rho_basis
   real*8,allocatable,dimension(:,:,:) :: rho,rho_gen,rho_orig
   real*8,allocatable,dimension(:,:,:) :: rho_in
   real*8,allocatable,dimension(:,:,:) :: amask
   real*8,allocatable,dimension(:) :: R_grid
   real*8,allocatable,dimension(:) :: vr_tmp
   real*8,allocatable,dimension(:,:) :: funcr2,funcr2_in,funcr_basis
   real*8,allocatable,dimension(:,:) :: funcr2_tmp,funcr2_tmp2
   real*8,allocatable,dimension(:,:) :: funcr2_tmp3
   real*8,allocatable,dimension(:,:,:) :: rho_gauss
   real*8,allocatable,dimension(:,:) :: amoment,amoment_subset
   real*8,allocatable,dimension(:) :: ww,amoment0
   real*8,allocatable,dimension(:) :: charge_type,charge_type_in
   integer,allocatable,dimension(:) :: itype_subset


   character*40 filename_xatom(100),filename_rho(100)
   character(len=100) fileId

   open(11,file="rho_fitV33.input")
   rewind(11)
   read(11,*) nsystem
   ! do ii=1,nsystem
   !    read(11,*) filename_xatom(ii),filename_rho(ii)
   ! enddo
   read(11,*) iflag_in,iflag_mask,dGauss
   read(11,*) num_G,Rmax,Rmin,Ztot,fact_min
   Rmax=Rmax/0.529177
   Rmin=Rmin/0.529177
   dGauss=dGauss/0.529177
   read(11,*) ncent,ntype
   do ii=1,ncent
      read(11,*) itype_cent(ii),zcent(ii),num,(icent(jj,ii),jj=1,num),(w_cent(jj,ii),jj=1,num)
      nat_cent(ii)=num
   enddo
   close(11)

   allocate(Q_type(ntype))

   !ccccccccccccccccccccccccccccccccccccccccccccccccc

   allocate(R_grid(num_G+1))

   iflag=1

   if(iflag.eq.1) then
      alpha=(Rmax/Rmin)**(1.d0/num_G)
      do ii=1,num_G+1
         R_grid(ii)=Rmin*alpha**(ii-1)-Rmin
         R_grid(ii)=R_grid(ii)+Rmin*(ii-1)
         if(ii.gt.1) R_grid(ii)=R_grid(ii)+Rmin*fact_min
      enddo
      R_grid(1)=-1.D-10
      fact=Rmax/R_grid(num_G+1)
      !      fact=Rmax/R_grid(num_G)
      R_grid=R_grid*fact
   endif

   !cccccccccccccccccccccccccccccccccccccccccccccc
   if(iflag.eq.2) then
      do ii=2,num_G+1
         R_grid(ii)=(ii-1)*Rmax/num_G
      enddo
      R_grid(1)=-1.D-10
   endif
   !cccccccccccccccccccccccccccccccccccccccccc

   write(6,*) "R_grid",R_grid(1:3)

   num_Gtot=num_G*ntype
   !ccccccccccccccccccccccccccccccccccccccccccccccccc
   Rm2=2*Rmax
   nr=50000
   yy2_max=0.d0
   do ir=1,nr
      r=Rm2*(ir-1.d0)/nr
      if(ir.eq.1) r=1.E-5
      yy=r**2*exp(-(r-0.5*R_grid(3)/1.8)**2/(R_grid(3)/1.8)**2*2.772)
      if(yy.lt.1.D-20) yy=0.d0
      if(yy.gt.yy2_max) yy2_max=yy
   enddo
   !ccccccccccccccccccccccccccccccccccccccccccccccccc
   if(iflag_in.eq.1) then
      open(13,file="funcr_atom.fit.bin.in",form="unformatted")
      rewind(13)
      read(13) nr_in,Rm2_in,ntype_in
      if(ntype.ne.ntype_in) then
         write(6,*) "ntype,ntype_in not the same,stop",ntype,ntype_in
         stop
      endif
      allocate(funcr2_in(nr,ntype))
      allocate(charge_type_in(ntype))
      read(13) charge_type_in
      read(13) funcr2_in    ! don't do smoothing
      close(13)
   endif

   !ccccccccccccccccccccccccccccccccccccccccccccccccc
   allocate(BB(num_Gtot))
   allocate(AA(num_Gtot,num_Gtot))
   allocate(SS(num_Gtot))
   lwork=30*num_Gtot+60
   allocate(work(lwork))
   !ccccccccccccccccccccccccccccccccccccccccccccccccc

   BB=0.d0
   AA=0.d0


   do 5000 iii=1,nsystem


      ! different systems can have different lattice, but the number of atoms
      ! and their orders must the be same
      write(fileId,"(I5)") iii
      open(10,file="./data/xatom"//trim(adjustl(fileId))//".config")
      write(6,*) "@@@ read in config ", iii, "./data/xatom"//trim(adjustl(fileId))//".config"
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

      AL=AL/0.529177d0
      vol=al(3,1)*(al(1,2)*al(2,3)-al(1,3)*al(2,2))   &
         +al(3,2)*(al(1,3)*al(2,1)-al(1,1)*al(2,3))  &
         +al(3,3)*(al(1,1)*al(2,2)-al(1,2)*al(2,1))

      vol=abs(vol)

      do ii=1,ncent
         do ixyz=1,3
            x1=xatom(ixyz,icent(1,ii))
            xc=x1*w_cent(1,ii)
            w_sum=w_cent(1,ii)
            do jj=2,nat_cent(ii)
               x2=xatom(ixyz,icent(jj,ii))
               if(abs(x2+1-x1).lt.abs(x2-x1)) x2=x2+1
               if(abs(x2-1-x1).lt.abs(x2-x1)) x2=x2-1
               xc=xc+x2*w_cent(jj,ii)
               w_sum=w_sum+w_cent(jj,ii)
            enddo
            xc=xc/w_sum
            xatom_cent(ixyz,ii)=xc
         enddo
      enddo


      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      open(11,file="./data/OUT.RHO"//trim(adjustl(fileId)),form="unformatted")
      write(6,*) "@@@ read in rho ", iii, "./data/OUT.RHO"//trim(adjustl(fileId))
      rewind(11)
      read(11) n1,n2,n3,nnodes
      read(11) AL_f

      write(6,*) "n1, n2, n3 ", n1, n2, n3

      nr=n1*n2*n3
      nr_n=nr/nnodes
      allocate(vr_tmp(nr_n))
      allocate(rho(n1,n2,n3))
      do iread=1,nnodes
         read(11) vr_tmp
         do ii=1,nr_n
            jj=ii+(iread-1)*nr_n
            i=(jj-1)/(n2*n3)+1
            j=(jj-1-(i-1)*n2*n3)/n3+1
            k=jj-(i-1)*n2*n3-(j-1)*n3
            rho(i,j,k)=vr_tmp(ii)
         enddo
      enddo
      deallocate(vr_tmp)
      close(11)

      allocate(amask(n1,n2,n3))
      amask=0.d0
      charge0=0.d0
      do k=1,n3
         do j=1,n2
            do i=1,n1
               charge0=charge0+rho(i,j,k)
               if(iflag_mask.eq.1) then
                  amask(i,j,k)=0.001**2/(rho(i,j,k)**2+0.001**2)
               else
                  amask(i,j,k)=1.d0
               endif

            enddo
         enddo
      enddo
      charge0=charge0*vol/(n1*n2*n3)
      write(6,*) "total charge=",charge0,Ztot
      if(abs(Ztot).gt.0.001) then   ! not neutral one
         fact=Ztot/charge0
         rho=rho*fact
         charge0=Ztot
      endif


      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(iflag_in.eq.1) then

         Rcut_in=Rm2_in/1.7
         allocate(rho_in(n1,n2,n3))
         rho_in=0.d0
         do iat=1,ncent
            itype=itype_cent(iat)
            do k=1,n3
               dx3=(k-1.d0)/n3-xatom_cent(3,iat)
               if(abs(dx3+1).lt.abs(dx3)) dx3=dx3+1
               if(abs(dx3-1).lt.abs(dx3)) dx3=dx3-1
               do j=1,n2
                  dx2=(j-1.d0)/n2-xatom_cent(2,iat)
                  if(abs(dx2+1).lt.abs(dx2)) dx2=dx2+1
                  if(abs(dx2-1).lt.abs(dx2)) dx2=dx2-1
                  do i=1,n1
                     dx1=(i-1.d0)/n1-xatom_cent(1,iat)
                     if(abs(dx1+1).lt.abs(dx1)) dx1=dx1+1
                     if(abs(dx1-1).lt.abs(dx1)) dx1=dx1-1

                     dx=AL(1,1)*dx1+AL(1,2)*dx2+AL(1,3)*dx3
                     dy=AL(2,1)*dx1+AL(2,2)*dx2+AL(2,3)*dx3
                     dz=AL(3,1)*dx1+AL(3,2)*dx2+AL(3,3)*dx3
                     d=sqrt(dx**2+dy**2+dz**2)

                     if(d.lt.Rcut_in) then
                        fact11=nr_in/Rm2_in
                        yy=d*fact11
                        ir=yy
                        x=yy-ir
                        f1=1-x-0.5d0*x*(1-x)
                        f2=x+x*(1-x)
                        f3=-0.5d0*x*(1-x)
                        ir=ir+1
                        rho_in(i,j,k)=rho_in(i,j,k)+ &
                           funcr2_in(ir,itype)*f1+funcr2_in(ir+1,itype)*f2+ &
                           funcr2_in(ir+2,itype)*f3
                     endif
                  enddo
               enddo
            enddo
         enddo

         rho=rho-rho_in
         deallocate(rho_in)
      endif
      !ccccccccccccccccccccccccccccccccccccccccccccc
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      allocate(rho_basis(n1,n2,n3,num_Gtot))
      rho_basis=0.d0

      pi=4*datan(1.d0)


      do iat=1,ncent
         itype=itype_cent(iat)
         x1=mod(xatom_cent(1,iat)+1.d0,1.d0)
         x2=mod(xatom_cent(2,iat)+1.d0,1.d0)
         x3=mod(xatom_cent(3,iat)+1.d0,1.d0)

         do k=1,n3
            dx3=(k-1.d0)/n3-x3
            if(abs(dx3+1).lt.abs(dx3)) dx3=dx3+1
            if(abs(dx3-1).lt.abs(dx3)) dx3=dx3-1
            do j=1,n2
               dx2=(j-1.d0)/n2-x2
               if(abs(dx2+1).lt.abs(dx2)) dx2=dx2+1
               if(abs(dx2-1).lt.abs(dx2)) dx2=dx2-1
               do i=1,n1
                  dx1=(i-1.d0)/n1-x1
                  if(abs(dx1+1).lt.abs(dx1)) dx1=dx1+1
                  if(abs(dx1-1).lt.abs(dx1)) dx1=dx1-1

                  dx=AL(1,1)*dx1+AL(1,2)*dx2+AL(1,3)*dx3
                  dy=AL(2,1)*dx1+AL(2,2)*dx2+AL(2,3)*dx3
                  dz=AL(3,1)*dx1+AL(3,2)*dx2+AL(3,3)*dx3

                  d=dsqrt(dx**2+dy**2+dz**2)



                  !        if(d.lt.R_grid(3)) then
                  !         rho_basis(i,j,k,1+(itype-1)*num_G)=rho_basis(i,j,k,1+(itype-1)*num_G)+ &
                  !                   cos(d/R_grid(3)*pi/2)**2
                  ! each basis might not be just one ringe!
                  !        endif
                  !       if(d.lt.R_grid(3)) then
                  !         rho_basis(i,j,k,2+(itype-1)*num_G)=rho_basis(i,j,k,2+(itype-1)*num_G)+ &
                  !              sin(d/R_grid(3)*pi)**2
                  !        endif
                  !        if(d.gt.R_grid(2)/2.and.d.lt.R_grid(4)) then
                  !         rho_basis(i,j,k,3+(itype-1)*num_G)=rho_basis(i,j,k,3+(itype-1)*num_G)+ &
                  !              sin((d-R_grid(2)/2)/(R_grid(4)-R_grid(2)/2)*pi)**2
                  !        endif
                  !        do ii=4,num_G
                  !        if(d.gt.R_grid(ii-2).and.d.lt.R_grid(ii+1)) then
                  !         rho_basis(i,j,k,ii+(itype-1)*num_G)=rho_basis(i,j,k,ii+(itype-1)*num_G)+&
                  !              sin((d-R_grid(ii-2))/(R_grid(ii+1)-R_grid(ii-2))*pi)**2
                  ! each basis might not be just one ringe!
                  !        endif


                  yy=exp(-d**2/R_grid(3)**2*2.772)
                  if(yy.lt.1.D-20) yy=0.d0
                  rho_basis(i,j,k,1+(itype-1)*num_G)=rho_basis(i,j,k,1+(itype-1)*num_G)+ yy


                  yy=d**2*exp(-(d-0.5*R_grid(3)/1.8)**2/(R_grid(3)/1.8)**2*2.772)/yy2_max
                  if(yy.lt.1.D-20) yy=0.d0
                  rho_basis(i,j,k,2+(itype-1)*num_G)=rho_basis(i,j,k,2+(itype-1)*num_G)+ yy


                  yy=exp(-(d-(0.25*R_grid(2)+0.5*R_grid(4)))**2/((R_grid(4)-R_grid(2)/2))**2*10.0)
                  if(yy.lt.1.D-20) yy=0.d0
                  rho_basis(i,j,k,3+(itype-1)*num_G)=rho_basis(i,j,k,3+(itype-1)*num_G)+ yy


                  do ii=4,num_G
                     yy=exp(-(d-0.5*(R_grid(ii+1)+R_grid(ii-2)))**2/(R_grid(ii+1)-R_grid(ii-2))**2*10.0)
                     if(yy.lt.1.D-20) yy=0.d0
                     rho_basis(i,j,k,ii+(itype-1)*num_G)=rho_basis(i,j,k,ii+(itype-1)*num_G)+yy
                  enddo

               enddo
            enddo
         enddo
      enddo
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do ii1=1,num_Gtot

         sum=0.d0
         do k=1,n3
            do j=1,n2
               do i=1,n1
                  sum=sum+rho(i,j,k)*rho_basis(i,j,k,ii1)*amask(i,j,k)
               enddo
            enddo
         enddo
         BB(ii1)=BB(ii1)+sum
      enddo
      !ccccccccccccccccccccccccccccc

      do ii1=1,num_Gtot
         do ii2=1,ii1

            !ccccccccc I can use blas to get AA
            sum=0.d0
            do k=1,n3
               do j=1,n2
                  do i=1,n1
                     sum=sum+rho_basis(i,j,k,ii1)*rho_basis(i,j,k,ii2)*amask(i,j,k)
                  enddo
               enddo
            enddo
            AA(ii1,ii2)=AA(ii1,ii2)+sum
            if(ii1.ne.ii2) then
               AA(ii2,ii1)=AA(ii2,ii1)+sum
            endif
         enddo
      enddo

      if(iii.ne.nsystem) then
            !   we keep the last one for a test
         deallocate(rho_basis)
         deallocate(rho)
         deallocate(amask)
      endif


5000 continue   ! system

   write(6,*) "finished prepare A and B"

   delta=0.d0
   do ii=1,num_Gtot
      AA(ii,ii)=AA(ii,ii)+delta
   enddo

   call dgelss(num_Gtot,num_Gtot,1,AA,num_Gtot, BB, &
      num_Gtot,SS,-0.1,irank,work,lwork,info)

   write(6,*) "end solve linear equation"
      !       write(6,*) "test,BB",BB

   !ccccccccccccccccccccccccccccccc
   open(10,file="BB.store",form="unformatted")
   rewind(10)
   write(10) num_G,ntype,num_G*ntype
   write(10) BB
   close(10)

   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   allocate(rho_orig(n1,n2,n3))
   rho_orig=rho


   sum=0.d0
   diff=0.d0
   sum2=0.d0
   sum0=0.d0

   do k=1,n3
      do j=1,n2
         do i=1,n1

            sum_rho=0.d0
            do ii=1,num_Gtot
               sum_rho=sum_rho+BB(ii)*rho_basis(i,j,k,ii)
            enddo

            sum=sum+rho(i,j,k)**2*amask(i,j,k)
            diff=diff+abs(rho(i,j,k)-sum_rho)**2*amask(i,j,k)
            rho(i,j,k)=sum_rho
            sum2=sum2+rho(i,j,k)
            sum0=sum0+rho_orig(i,j,k)
         enddo
      enddo
   enddo
   sum=sum*vol/(n1*n2*n3)
   sum2=sum2*vol/(n1*n2*n3)
   sum0=sum0*vol/(n1*n2*n3)
   diff=diff*vol/(n1*n2*n3)
   write(6,*) "(last image) sum**2,diff**2",sum,diff
   write(6,*) "(last image) tot charge=",sum0,sum2

   dQtot=sum0-sum2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc


1000 continue


   Rm2=2*Rmax
   nr=50000

   if(iflag_in.eq.1) then
      Rm2=Rm2_in
      nr=nr_in
   endif


   allocate(funcr_basis(nr,num_G))
   allocate(funcr2(nr,ntype))
   funcr_basis=0.d0
   do ir=1,nr
      r=Rm2*(ir-1)/nr
      if(ir.eq.1) r=1.E-5

      !ccccccccccccccccccccccccccccccccccccccccccccccccccc
      !        if(r.lt.R_grid(3)) then
      !          funcr_basis(ir,1)=cos(r/R_grid(3)*pi/2)**2
      !        endif


      !        if(r.lt.R_grid(3)) then
      !          funcr_basis(ir,2)=sin(r/R_grid(3)*pi)**2
      !        endif

      !        if(r.gt.R_grid(2)/2.and.r.lt.R_grid(4)) then
      !          funcr_basis(ir,3)=  &
      !              sin((r-R_grid(2)/2)/(R_grid(4)-R_grid(2)/2)*pi)**2
      !        endif

      !        do ii=4,num_G
      !        if(r.gt.R_grid(ii-2).and.r.lt.R_grid(ii+1)) then
      !          funcr_basis(ir,i)=  &
      !              sin((r-R_grid(ii-2))/(R_grid(ii+1)-R_grid(ii-2))*pi)**2
      !        endif
      !         enddo
      !ccccccccccccccccccccccccccccccccccccccccccccccccccc
      yy=exp(-r**2/R_grid(3)**2*2.772)
      if(yy.lt.1.D-20) yy=0.d0
      funcr_basis(ir,1)=yy


      yy=r**2*exp(-(r-0.5*R_grid(3)/1.8)**2/(R_grid(3)/1.8)**2*2.772)/yy2_max
      if(yy.lt.1.D-20) yy=0.d0
      funcr_basis(ir,2)=yy

      yy=exp(-(r-(0.25*R_grid(2)+0.5*R_grid(4)))**2/((R_grid(4)-R_grid(2)/2))**2*10.0)
      if(yy.lt.1.D-20) yy=0.d0
      funcr_basis(ir,3)=yy

      do ii=4,num_G
         yy=exp(-(r-0.5*(R_grid(ii+1)+R_grid(ii-2)))**2/(R_grid(ii+1)-R_grid(ii-2))**2*10.0)
         if(yy.lt.1.D-20) yy=0.d0
         funcr_basis(ir,ii)=yy
      enddo
      !ccccccccccccccccccccccccccccccccccccccccccccccccccc
   enddo


   Qtot=0.d0
   do iat=1,ncent
      Qtot=Qtot+Zcent(iat)
   enddo


   do itype=1,ntype
      dQ=-1000.d0
      do iat=1,ncent
         if(itype_cent(iat).eq.itype) then
            dQ=Zcent(iat)/Qtot*dQtot
         endif
      enddo
      if(dQ.lt.-999.d0) then
         write(6,*) "something wrong, dQ", itype,dQ
      endif
      write(6,*) "itype,dQ",itype,dQ


      do ir=1,nr
         r=Rm2*(ir-1)/nr
         sum=0.d0
         do ii=1,num_G
            sum=sum+funcr_basis(ir,ii)*BB(ii+(itype-1)*num_G)
         enddo
         funcr2(ir,itype)=sum
         ! dQ is from the dQtot, not for each individual Zcent(iat)
         funcr2(ir,itype)=funcr2(ir,itype)+exp(-(r/dGauss)**2)/(dGauss**3*pi**1.5)*dQ
         if(abs(funcr2(ir,itype)).lt.1.D-30) funcr2(ir,itype)=0.d0
      enddo
   enddo

   allocate(charge_type(ntype))


   if(iflag_in.eq.1) then
      funcr2=funcr2+funcr2_in
   endif


   do itype=1,ntype
      sum=0.d0
      do ir=1,nr
         r=Rm2*(ir-1)/nr
         sum=sum+funcr2(ir,itype)*r**2
      enddo
      sum=sum*4*pi*Rm2/nr
      charge_type(itype)=sum
   enddo

   sum=0.d0
   do iat=1,ncent
      sum=sum+charge_type(itype_cent(iat))
   enddo
   write(6,*) "rescale: tot charge Z, sum(funcr)",ztot,sum
   funcr2=funcr2*ztot/sum
   charge_type=charge_type*ztot/sum

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   open(13,file="funcr_atom.fit")
   rewind(13)
   do ir=1,nr
      r=Rm2*(ir-1)/nr
      write(13,"(50(E11.4,1x))") r*0.5219177,(funcr2(ir,ii),ii=1,ntype)
   enddo
   close(13)

   open(13,file="funcr_atom.fit.bin",form="unformatted")
   rewind(13)
   write(13) nr,Rm2,ntype
   write(13) charge_type
   write(13) funcr2
   close(13)

   open(13,file="funcr_basis")
   rewind(13)
   do ir=1,nr,10
      r=Rm2*(ir-1)/nr
      write(13,"(50(E11.4,1x))") r*0.5219177,(funcr_basis(ir,ii),ii=1,10)
   enddo
   close(13)



   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   Rcut=Rm2/1.7
   allocate(rho_in(n1,n2,n3))
   rho_in=0.d0
   do iat=1,ncent
      itype=itype_cent(iat)
      do k=1,n3
         dx3=(k-1.d0)/n3-xatom_cent(3,iat)
         if(abs(dx3+1).lt.abs(dx3)) dx3=dx3+1
         if(abs(dx3-1).lt.abs(dx3)) dx3=dx3-1
         do j=1,n2
            dx2=(j-1.d0)/n2-xatom_cent(2,iat)
            if(abs(dx2+1).lt.abs(dx2)) dx2=dx2+1
            if(abs(dx2-1).lt.abs(dx2)) dx2=dx2-1
            do i=1,n1
               dx1=(i-1.d0)/n1-xatom_cent(1,iat)
               if(abs(dx1+1).lt.abs(dx1)) dx1=dx1+1
               if(abs(dx1-1).lt.abs(dx1)) dx1=dx1-1

               dx=AL(1,1)*dx1+AL(1,2)*dx2+AL(1,3)*dx3
               dy=AL(2,1)*dx1+AL(2,2)*dx2+AL(2,3)*dx3
               dz=AL(3,1)*dx1+AL(3,2)*dx2+AL(3,3)*dx3
               d=sqrt(dx**2+dy**2+dz**2)

               if(d.lt.Rcut) then
                  fact11=nr/Rm2
                  yy=d*fact11
                  ir=yy
                  x=yy-ir
                  f1=1-x-0.5d0*x*(1-x)
                  f2=x+x*(1-x)
                  f3=-0.5d0*x*(1-x)
                  ir=ir+1
                  rho_in(i,j,k)=rho_in(i,j,k)+ &
                     funcr2(ir,itype)*f1+funcr2(ir+1,itype)*f2+ &
                     funcr2(ir+2,itype)*f3
               endif
            enddo
         enddo
      enddo
   enddo

   open(11,file="OUT.FIT_rho",form="unformatted")
   rewind(11)
   write(11) n1,n2,n3,nnodes
   write(11) AL_f

   nr=n1*n2*n3
   nr_n=nr/nnodes
   allocate(vr_tmp(nr_n))
   do iread=1,nnodes
      do ii=1,nr_n
         jj=ii+(iread-1)*nr_n
         i=(jj-1)/(n2*n3)+1
         j=(jj-1-(i-1)*n2*n3)/n3+1
         k=jj-(i-1)*n2*n3-(j-1)*n3
         vr_tmp(ii)=rho_in(i,j,k)
      enddo
      write(11) vr_tmp
   enddo
   deallocate(vr_tmp)
   close(11)

   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   stop
end
