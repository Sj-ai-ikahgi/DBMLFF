subroutine find_feature_3b_type_ppp(natom,itype_atom,Rc,Rc2,n3b1,n3b2, &
   num_neigh,list_neigh, &
   dR_neigh,iat_neigh,ntype,grid31_2,grid32_2, &
   feat_all,dfeat_all,nfeat0m,m_neigh,n3b1m,n3b2m,nfeat_atom, &
   xatom,AL,iat11,iat22,iat33,iat44,idir)
   use mod_mpi
   implicit none
   integer ntype
   integer natom,n2b,n3b1(ntype),n3b2(ntype)
   integer nfeat0m,m_neigh,itype0
   integer n3b1m,n3b2m
   integer nfeat_atom(natom)
   integer nfeat_atom_tmp(natom)
   integer itype_atom(natom)
   real*8 Rc(ntype),Rc2(ntype)
   real*8 dR_neigh(3,m_neigh,ntype,natom)
   real*8 dR_neigh_alltype(3,m_neigh,natom)
   integer iat_neigh(m_neigh,ntype,natom),list_neigh(m_neigh,ntype,natom)
   integer num_neigh(ntype,natom)
   integer num_neigh_alltype(natom)
   integer nperiod(3)
   integer iflag,i,j,num,iat,itype
   integer i1,i2,i3,itype1,itype2,j1,j2,iat1,iat2
   real*8 d,dx1,dx2,dx3,dx,dy,dz,dd
   real*8 grid31_2(2,n3b1m,ntype),grid32_2(2,n3b2m,ntype)
   real*8 pi,pi2,x,f1
   integer iflag_grid


   real*8 feat3_tmp(n3b1m,m_neigh,ntype)
   real*8 dfeat3_tmp(n3b1m,m_neigh,ntype,3)
   integer ind_f(n3b1m,m_neigh,ntype,natom)
   integer num_k(m_neigh,ntype),num12
   real*8 f32(n3b2m),df32(n3b2m,2,3)
   integer inf_f32(n3b2m),k,k1,k2,k12,j12,ii_f,jj,jj1,jj2,nneigh,ii
   real*8 y,y2
   integer itype12,ind_f32(n3b2m)
   integer ind_all_neigh(m_neigh,ntype,natom),list_neigh_alltype(m_neigh,natom)

   real*8 feat3(n3b1m*n3b1m*n3b2m,ntype*(ntype+1)/2,natom_n)
   real*8 dfeat3(n3b1m*n3b1m*n3b2m,ntype*(ntype+1)/2,natom_n,m_neigh,3)
   real*8 feat_all(nfeat0m,natom_n)
   real*8 dfeat_all(nfeat0m,natom_n,m_neigh,3)
   integer ierr

   ! ---------- new variables for type_ppp ----------
   real*8 ppp_dir(3,3)
   real*8 xatom(3,natom)
   real*8 AL(3,3)
   integer iat11,iat22,iat33,iat44
   integer idir
   real*8 dx11_frac(3),dx22_frac(3)
   real*8 dx11(3),dx22(3),dx33(3)
   real*8 d11,d22,d33
   real*8 ddx11_diat22(3,3),ddx22_diat44(3,3),ddx33_diat22(3,3),ddx33_diat44(3,3)
   integer jj11,jj22,jj33,jj44
   real*8 drange
   real*8 temp1(3),temp2(n3b1m,m_neigh,ntype),temp3,temp4,temp5,temp6,temp7,temp8
   real*8 y2_saving(n3b1m,m_neigh,ntype)

   ! ---------- calc ppp_dir and dppp_dir ----------
   ppp_dir=0.d0
   ddx11_diat22=0.d0
   ddx22_diat44=0.d0
   ddx33_diat22=0.d0
   ddx33_diat44=0.d0

   do i=1,3
      dx11_frac(i)=xatom(i,iat22)-xatom(i,iat11)
      dx22_frac(i)=xatom(i,iat44)-xatom(i,iat33)
      do while(abs(dx11_frac(i)+1).lt.abs(dx11_frac(i)))
         dx11_frac(i)=dx11_frac(i)+1
      enddo
      do while(abs(dx11_frac(i)-1).lt.abs(dx11_frac(i)))
         dx11_frac(i)=dx11_frac(i)-1
      enddo
      do while(abs(dx22_frac(i)+1).lt.abs(dx22_frac(i)))
         dx22_frac(i)=dx22_frac(i)+1
      enddo
      do while(abs(dx22_frac(i)-1).lt.abs(dx22_frac(i)))
         dx22_frac(i)=dx22_frac(i)-1
      enddo
   enddo

   dx11(:)=AL(:,1)*dx11_frac(1)+AL(:,2)*dx11_frac(2)+AL(:,3)*dx11_frac(3)
   dx22(:)=AL(:,1)*dx22_frac(1)+AL(:,2)*dx22_frac(2)+AL(:,3)*dx22_frac(3)

   dx33(1)=dx11(2)*dx22(3)-dx11(3)*dx22(2)     ! dx33=dx11 x dx22
   dx33(2)=dx11(3)*dx22(1)-dx11(1)*dx22(3)
   dx33(3)=dx11(1)*dx22(2)-dx11(2)*dx22(1)

   ! ddx11_diat22(3,3) and ddx22_diat44(3,3) are symmetric matrix
   ! you can consider both ddx11_diat22(1,2) and ddx11_diat22(2,1) as ddx11(1)/dx22
   ! however, to be consistent with ddx33_d, we should use ddx11_diat22(1,2) as ddx11(1)/dx22
   d11=dsqrt(dx11(1)**2+dx11(2)**2+dx11(3)**2)
   do i1=1,3
      do i2=1,3
         ddx11_diat22(i1,i2)=-dx11(i1)*dx11(i2)/d11**3
      enddo
      ddx11_diat22(i1,i1)=ddx11_diat22(i1,i1)+1/d11
   enddo

   d22=dsqrt(dx22(1)**2+dx22(2)**2+dx22(3)**2)
   do i1=1,3
      do i2=1,3
         ddx22_diat44(i1,i2)=-dx22(i1)*dx22(i2)/d22**3
      enddo
      ddx22_diat44(i1,i1)=ddx22_diat44(i1,i1)+1/d22
   enddo

   ! ddx33_diat22(1,2) is ddx33(1)/dx22
   ! ddx33_diat44(2,3) is ddx33(2)/dx43
   d33=dsqrt(dx33(1)**2+dx33(2)**2+dx33(3)**2)
   do i1=1,3
      ! derivative for the 1/d part
      ddx33_diat22(i1,1)=-dx33(i1)/d33**3*(-dx33(2)*dx22(3)+dx33(3)*dx22(2))
      ddx33_diat22(i1,2)=-dx33(i1)/d33**3*(dx33(1)*dx22(3)-dx33(3)*dx22(1))
      ddx33_diat22(i1,3)=-dx33(i1)/d33**3*(-dx33(1)*dx22(2)+dx33(2)*dx22(1))

      ddx33_diat44(i1,1)=-dx33(i1)/d33**3*(dx33(2)*dx11(3)-dx33(3)*dx11(2))
      ddx33_diat44(i1,2)=-dx33(i1)/d33**3*(-dx33(1)*dx11(3)+dx33(3)*dx11(1))
      ddx33_diat44(i1,3)=-dx33(i1)/d33**3*(dx33(1)*dx11(2)-dx33(2)*dx11(1))
   enddo

   ddx33_diat22(1,2)=ddx33_diat22(1,2)+dx22(3)/d33
   ddx33_diat22(1,3)=ddx33_diat22(1,3)-dx22(2)/d33
   ddx33_diat22(2,1)=ddx33_diat22(2,1)-dx22(3)/d33
   ddx33_diat22(2,3)=ddx33_diat22(2,3)+dx22(1)/d33
   ddx33_diat22(3,1)=ddx33_diat22(3,1)+dx22(2)/d33
   ddx33_diat22(3,2)=ddx33_diat22(3,2)-dx22(1)/d33

   ddx33_diat44(1,2)=ddx33_diat44(1,2)-dx11(3)/d33
   ddx33_diat44(1,3)=ddx33_diat44(1,3)+dx11(2)/d33
   ddx33_diat44(2,1)=ddx33_diat44(2,1)+dx11(3)/d33
   ddx33_diat44(2,3)=ddx33_diat44(2,3)-dx11(1)/d33
   ddx33_diat44(3,1)=ddx33_diat44(3,1)-dx11(2)/d33
   ddx33_diat44(3,2)=ddx33_diat44(3,2)+dx11(1)/d33

   dx11=dx11/d11
   dx22=dx22/d22
   dx33=dx33/d33

   ppp_dir(:,1)=dx11
   ppp_dir(:,2)=dx22
   ppp_dir(:,3)=dx33

   num_neigh_alltype=0
   do iat=1,natom
      num=1
      list_neigh_alltype(1,iat)=iat   ! the first neighbore is itself
      dR_neigh_alltype(:,1,iat)=0.d0

      do  itype=1,ntype
         do   j=1,num_neigh(itype,iat)
            num=num+1
            if(num.gt.m_neigh) then
               write(6,*) "total num_neigh.gt.m_neigh,stop",m_neigh
               stop
            endif
            ind_all_neigh(j,itype,iat)=num
            list_neigh_alltype(num,iat)=list_neigh(j,itype,iat)
            dR_neigh_alltype(:,num,iat)=dR_neigh(:,j,itype,iat)
         enddo
      enddo
      num_neigh_alltype(iat)=num
   enddo

   !ccccccccccccccccccccccccccccccccccccccccc

   pi=4*datan(1.d0)
   pi2=2*pi

   feat3=0.d0
   dfeat3=0.d0

   iat1=0
   do 3000 iat=1,natom
    !   if(mod(iat-1,nnodes).eq.inode-1) then
         iat1=iat1+1

         itype0=itype_atom(iat)

         if(num_neigh_alltype(iat).gt.1) then
            jj11=0
            jj22=0
            jj33=0
            jj44=0
            do i=1,num_neigh_alltype(iat)
               if(list_neigh_alltype(i,iat).eq.iat11) jj11=i
               if(list_neigh_alltype(i,iat).eq.iat22) jj22=i
               if(list_neigh_alltype(i,iat).eq.iat33) jj33=i
               if(list_neigh_alltype(i,iat).eq.iat44) jj44=i
            enddo
            if((idir.eq.1).and.((jj11.eq.0).or.(jj22.eq.0))) then
               write(6,*) "(idir.eq.1).and.((jj11.eq.0).or.(jj22.eq.0)),stop"
               stop
            endif
            if((idir.eq.2).and.((jj33.eq.0).or.(jj44.eq.0))) then
               write(6,*) "(idir.eq.2).and.((jj33.eq.0).or.(jj44.eq.0)),stop"
               stop
            endif
            if((idir.eq.3).and.((jj11.eq.0).or.(jj22.eq.0).or.(jj33.eq.0).or.(jj44.eq.0))) then
               write(6,*) "(idir.eq.3).and.((jj11.eq.0).or.(jj22.eq.0).or.(jj33.eq.0).or.(jj44.eq.0)),stop"
               stop
            endif
         endif

         do 1000 itype=1,ntype
            do 1000 j=1,num_neigh(itype,iat) ! for itype, the self is not in the neighbore list

               jj=ind_all_neigh(j,itype,iat)   ! but itself is in jj=1

               dd=dR_neigh(1,j,itype,iat)**2+dR_neigh(2,j,itype,iat)**2+dR_neigh(3,j,itype,iat)**2
               d=dsqrt(dd)

               !ccccccccccccccccccccccccccccccccccccccccccccccccccccc

               num=0
               do k=1,n3b1(itype0)

                  if(d.ge.grid31_2(1,k,itype0).and.d.lt.grid31_2(2,k,itype0)) then
                     num=num+1

                     drange=grid31_2(2,k,itype0)-grid31_2(1,k,itype0)
                     x=(d-grid31_2(1,k,itype0))/drange
                     y=(x-0.5d0)*pi2
                     f1=0.5d0*(cos(y)+1)
                     feat3_tmp(num,j,itype)=f1
                     ind_f(num,j,itype,iat)=k
                     y2=-pi*sin(y)/(d*drange)
                     y2_saving(num,j,itype)=y2
                     temp2(num,j,itype)=(drange*sin(y)-pi2*d*cos(y))/(dd*d*drange**2)
                     dfeat3_tmp(num,j,itype,:)=y2*dR_neigh(:,j,itype,iat)
                  endif
               enddo
               num_k(j,itype)=num

            !cccccccccccc So, one Rij will always have two features k, k+1  (1,2)
            !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
1000     continue


        !   Now, the three body feature
        !ccccccccccccccccccccccccccccccccccccc



         do 2000 itype2=1,ntype
            do 2000 itype1=1,itype2

               itype12=itype1+((itype2-1)*itype2)/2

               do 2000 j2=1,num_neigh(itype2,iat)
                  do 2000 j1=1,num_neigh(itype1,iat)

                     !      if(itype1.eq.itype2.and.j1.ge.j2) goto 2000
                     if(itype1.eq.itype2.and.j1.eq.j2) goto 2000

                     jj1=ind_all_neigh(j1,itype1,iat)
                     jj2=ind_all_neigh(j2,itype2,iat)

                     dd=(dR_neigh(1,j1,itype1,iat)-dR_neigh(1,j2,itype2,iat))**2+ &
                        (dR_neigh(2,j1,itype1,iat)-dR_neigh(2,j2,itype2,iat))**2+ &
                        (dR_neigh(3,j1,itype1,iat)-dR_neigh(3,j2,itype2,iat))**2

                     d=dsqrt(dd)


                     if(d.gt.Rc2(itype0).or.d.lt.1.D-4) goto 2000


                     num=0
                     do k=1,n3b2(itype0)
                        if(d.ge.grid32_2(1,k,itype0).and.d.lt.grid32_2(2,k,itype0)) then
                           num=num+1

                           x=(d-grid32_2(1,k,itype0))/(grid32_2(2,k,itype0)-grid32_2(1,k,itype0))
                           y=(x-0.5d0)*pi2
                           f1=0.5d0*(cos(y)+1)
                           f32(num)=f1
                           ind_f32(num)=k
                           y2=-pi*sin(y)/(d*(grid32_2(2,k,itype0)-grid32_2(1,k,itype0)))
                           df32(num,1,:)=y2*(dR_neigh(:,j1,itype1,iat)-dR_neigh(:,j2,itype2,iat))
                           df32(num,2,:)=-df32(num,1,:)
                        endif

                     enddo
                     num12=num



                    !cccccccccccccccccccccccc
                    !   Each R has two k features, so for the three R, we have the following
                     do i1=1,num_k(j1,itype1)
                        do i2=1,num_k(j2,itype2)
                           do j12=1,num12
                              k1=ind_f(i1,j1,itype1,iat)
                              k2=ind_f(i2,j2,itype2,iat)
                              k12=ind_f32(j12)

                              ii_f=0
                              if(itype1.ne.itype2) then
                                 ii_f=k1+(k2-1)*n3b1(itype0)+(k12-1)*n3b1(itype0)**2
                              endif
                              if(itype1.eq.itype2.and.k1.le.k2) then
                                 ii_f=k1+((k2-1)*k2)/2+(k12-1)*(n3b1(itype0)*(n3b1(itype0)+1))/2
                              endif

                              if(ii_f.ne.0) then
                                 ! ---------- new feat for type_ppp ----------
                                 ! feat3(ii_f,itype12,iat1)=feat3(ii_f,itype12,iat1)+ &
                                 ! feat3_tmp(i1,j1,itype1)*feat3_tmp(i2,j2,itype2)*f32(j12)
                                 temp1(:)=dfeat3_tmp(i1,j1,itype1,:)*feat3_tmp(i2,j2,itype2)*f32(j12)+ &
                                    feat3_tmp(i1,j1,itype1)*dfeat3_tmp(i2,j2,itype2,:)*f32(j12)
                                 feat3(ii_f,itype12,iat1)=feat3(ii_f,itype12,iat1)- &
                                    temp1(1)*ppp_dir(1,idir)- &
                                    temp1(2)*ppp_dir(2,idir)- &
                                    temp1(3)*ppp_dir(3,idir)

                                 ! ---------- new dfeat for type_ppp ----------
                                 ! dfeat3(ii_f,itype12,iat1,jj1,:)=dfeat3(ii_f,itype12,iat1,jj1,:)+ &
                                 ! dfeat3_tmp(i1,j1,itype1,:)*feat3_tmp(i2,j2,itype2)*f32(j12)+ &
                                 ! feat3_tmp(i1,j1,itype1)*feat3_tmp(i2,j2,itype2)*df32(j12,1,:)
                                 ! dfeat3(ii_f,itype12,iat1,jj2,:)=dfeat3(ii_f,itype12,iat1,jj2,:)+ &
                                 ! feat3_tmp(i1,j1,itype1)*dfeat3_tmp(i2,j2,itype2,:)*f32(j12)+ &
                                 ! feat3_tmp(i1,j1,itype1)*feat3_tmp(i2,j2,itype2)*df32(j12,2,:)
                                 ! dfeat3(ii_f,itype12,iat1,1,:)=dfeat3(ii_f,itype12,iat1,1,:)- &
                                 ! dfeat3_tmp(i1,j1,itype1,:)*feat3_tmp(i2,j2,itype2)*f32(j12)- &
                                 ! feat3_tmp(i1,j1,itype1)*dfeat3_tmp(i2,j2,itype2,:)*f32(j12)
                                 temp3=dR_neigh(1,j1,itype1,iat)*feat3_tmp(i2,j2,itype2)*f32(j12)*ppp_dir(1,idir)+ &
                                    dR_neigh(2,j1,itype1,iat)*feat3_tmp(i2,j2,itype2)*f32(j12)*ppp_dir(2,idir)+ &
                                    dR_neigh(3,j1,itype1,iat)*feat3_tmp(i2,j2,itype2)*f32(j12)*ppp_dir(3,idir)
                                 temp4=dfeat3_tmp(i1,j1,itype1,1)*f32(j12)*ppp_dir(1,idir)+ &
                                    dfeat3_tmp(i1,j1,itype1,2)*f32(j12)*ppp_dir(2,idir)+ &
                                    dfeat3_tmp(i1,j1,itype1,3)*f32(j12)*ppp_dir(3,idir)
                                 temp5=dfeat3_tmp(i1,j1,itype1,1)*feat3_tmp(i2,j2,itype2)*ppp_dir(1,idir)+ &
                                    dfeat3_tmp(i1,j1,itype1,2)*feat3_tmp(i2,j2,itype2)*ppp_dir(2,idir)+ &
                                    dfeat3_tmp(i1,j1,itype1,3)*feat3_tmp(i2,j2,itype2)*ppp_dir(3,idir)
                                 temp6=dR_neigh(1,j2,itype2,iat)*feat3_tmp(i1,j1,itype1)*f32(j12)*ppp_dir(1,idir)+ &
                                    dR_neigh(2,j2,itype2,iat)*feat3_tmp(i1,j1,itype1)*f32(j12)*ppp_dir(2,idir)+ &
                                    dR_neigh(3,j2,itype2,iat)*feat3_tmp(i1,j1,itype1)*f32(j12)*ppp_dir(3,idir)
                                 temp7=dfeat3_tmp(i2,j2,itype2,1)*f32(j12)*ppp_dir(1,idir)+ &
                                    dfeat3_tmp(i2,j2,itype2,2)*f32(j12)*ppp_dir(2,idir)+ &
                                    dfeat3_tmp(i2,j2,itype2,3)*f32(j12)*ppp_dir(3,idir)
                                 temp8=dfeat3_tmp(i2,j2,itype2,1)*feat3_tmp(i1,j1,itype1)*ppp_dir(1,idir)+ &
                                    dfeat3_tmp(i2,j2,itype2,2)*feat3_tmp(i1,j1,itype1)*ppp_dir(2,idir)+ &
                                    dfeat3_tmp(i2,j2,itype2,3)*feat3_tmp(i1,j1,itype1)*ppp_dir(3,idir)
                                 dfeat3(ii_f,itype12,iat1,jj1,:)=dfeat3(ii_f,itype12,iat1,jj1,:)- &
                                    pi*dR_neigh(:,j1,itype1,iat)*temp2(i1,j1,itype1)*temp3-y2_saving(i1,j1,itype1)*feat3_tmp(i2,j2,itype2)*f32(j12)*ppp_dir(:,idir)- &
                                    df32(j12,1,:)*temp5- &
                                    dfeat3_tmp(i1,j1,itype1,:)*temp7- &
                                    df32(j12,1,:)*temp8
                                 dfeat3(ii_f,itype12,iat1,jj2,:)=dfeat3(ii_f,itype12,iat1,jj2,:)- &
                                    dfeat3_tmp(i2,j2,itype2,:)*temp4- &
                                    df32(j12,2,:)*temp5- &
                                    pi*dR_neigh(:,j2,itype2,iat)*temp2(i2,j2,itype2)*temp6-y2_saving(i2,j2,itype2)*feat3_tmp(i1,j1,itype1)*f32(j12)*ppp_dir(:,idir)- &
                                    df32(j12,2,:)*temp8
                                 dfeat3(ii_f,itype12,iat1,1,:)=dfeat3(ii_f,itype12,iat1,1,:)+ &
                                    pi*dR_neigh(:,j1,itype1,iat)*temp2(i1,j1,itype1)*temp3+y2_saving(i1,j1,itype1)*feat3_tmp(i2,j2,itype2)*f32(j12)*ppp_dir(:,idir)+ &
                                    dfeat3_tmp(i2,j2,itype2,:)*temp4+ &
                                    pi*dR_neigh(:,j2,itype2,iat)*temp2(i2,j2,itype2)*temp6+y2_saving(i2,j2,itype2)*feat3_tmp(i1,j1,itype1)*f32(j12)*ppp_dir(:,idir)+ &
                                    dfeat3_tmp(i1,j1,itype1,:)*temp7
                                 if(idir.eq.1) then
                                    dfeat3(ii_f,itype12,iat1,jj11,:)=dfeat3(ii_f,itype12,iat1,jj11,:)+ &
                                       temp1(1)*ddx11_diat22(1,:)+ &
                                       temp1(2)*ddx11_diat22(2,:)+ &
                                       temp1(3)*ddx11_diat22(3,:)
                                    dfeat3(ii_f,itype12,iat1,jj22,:)=dfeat3(ii_f,itype12,iat1,jj22,:)- &
                                       temp1(1)*ddx11_diat22(1,:)- &
                                       temp1(2)*ddx11_diat22(2,:)- &
                                       temp1(3)*ddx11_diat22(3,:)
                                 endif
                                 if(idir.eq.2) then
                                    dfeat3(ii_f,itype12,iat1,jj33,:)=dfeat3(ii_f,itype12,iat1,jj33,:)+ &
                                       temp1(1)*ddx22_diat44(1,:)+ &
                                       temp1(2)*ddx22_diat44(2,:)+ &
                                       temp1(3)*ddx22_diat44(3,:)
                                    dfeat3(ii_f,itype12,iat1,jj44,:)=dfeat3(ii_f,itype12,iat1,jj44,:)- &
                                       temp1(1)*ddx22_diat44(1,:)- &
                                       temp1(2)*ddx22_diat44(2,:)- &
                                       temp1(3)*ddx22_diat44(3,:)
                                 endif
                                 if(idir.eq.3) then
                                    dfeat3(ii_f,itype12,iat1,jj11,:)=dfeat3(ii_f,itype12,iat1,jj11,:)+ &
                                       temp1(1)*ddx33_diat22(1,:)+ &
                                       temp1(2)*ddx33_diat22(2,:)+ &
                                       temp1(3)*ddx33_diat22(3,:)
                                    dfeat3(ii_f,itype12,iat1,jj22,:)=dfeat3(ii_f,itype12,iat1,jj22,:)- &
                                       temp1(1)*ddx33_diat22(1,:)- &
                                       temp1(2)*ddx33_diat22(2,:)- &
                                       temp1(3)*ddx33_diat22(3,:)
                                    dfeat3(ii_f,itype12,iat1,jj33,:)=dfeat3(ii_f,itype12,iat1,jj33,:)+ &
                                       temp1(1)*ddx33_diat44(1,:)+ &
                                       temp1(2)*ddx33_diat44(2,:)+ &
                                       temp1(3)*ddx33_diat44(3,:)
                                    dfeat3(ii_f,itype12,iat1,jj44,:)=dfeat3(ii_f,itype12,iat1,jj44,:)- &
                                       temp1(1)*ddx33_diat44(1,:)- &
                                       temp1(2)*ddx33_diat44(2,:)- &
                                       temp1(3)*ddx33_diat44(3,:)
                                 endif
                                 !cccc (ii_f,itype12) is the feature index
                              endif

                           enddo
                        enddo
                     enddo
2000     continue

    !   endif
3000 continue

    !cccccccccccccccccccccccccccccccccccccccccccccccc
    !cccccccccccccccccccccccccccccccccccccccccccccccc
    !   Now, we collect everything together, collapse the index (k,itype)
    !   and feat2,feat3, into a single feature.

   nfeat_atom_tmp=0

   iat1=0
   do 5000 iat=1,natom
    !   if(mod(iat-1,nnodes).eq.inode-1) then
         iat1=iat1+1

         itype0=itype_atom(iat)
         nneigh=num_neigh_alltype(iat)
         num=0

         do itype2=1,ntype
            do itype1=1,itype2

               itype12=itype1+((itype2-1)*itype2)/2

               do k1=1,n3b1(itype0)
                  do k2=1,n3b1(itype0)
                     do k12=1,n3b2(itype0)

                        ii_f=0
                        if(itype1.ne.itype2) then
                           ii_f=k1+(k2-1)*n3b1(itype0)+(k12-1)*n3b1(itype0)**2
                        endif
                        if(itype1.eq.itype2.and.k1.le.k2) then
                           ii_f=k1+((k2-1)*k2)/2+(k12-1)*(n3b1(itype0)*(n3b1(itype0)+1))/2
                        endif

                        if(ii_f.gt.0) then
                           num=num+1
                           feat_all(num,iat1)=feat3(ii_f,itype12,iat1)
                           dfeat_all(num,iat1,1:nneigh,:)=dfeat3(ii_f,itype12,iat1,1:nneigh,:)
                        endif

                     enddo
                  enddo
               enddo
            enddo
         enddo
         nfeat_atom_tmp(iat)=num
         if(num.gt.nfeat0m) then
            write(6,*) "num.gt.nfeat0m,stop",num,nfeat0m
            stop
         endif

    !   endif
5000 continue

   !   call mpi_allreduce(nfeat_atom_tmp,nfeat_atom,natom,MPI_INTEGER,MPI_SUM,MPI_COMM_MOL,ierr)

   nfeat_atom = nfeat_atom_tmp

    !ccccccccccccccccccccccccccccccccccc
    !  Now, we have to redefine the dfeat_all in another way.
    !  dfeat_all(nfeat,iat,jneigh,3) means:
    !  d_jth_feat_of_iat/d_R(jth_neigh_of_iat)
    !  dfeat_allR(nfeat,iat,jneigh,3) means:
    !  d_jth_feat_of_jth_neigh/d_R(iat)
    !  Now, just output dfeat_all
    !cccccccccccccccccccccccccccccccccccccc

   return
end subroutine find_feature_3b_type_ppp
