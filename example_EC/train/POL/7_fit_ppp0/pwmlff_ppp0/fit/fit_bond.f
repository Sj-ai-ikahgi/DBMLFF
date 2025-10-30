       program linear_forceMM
       implicit double precision (a-h,o-z)


       integer lwork
       integer,allocatable,dimension(:) :: iatom,iatom_type,itype_atom
       real*8,allocatable,dimension(:) :: Energy,Energyt
       real*8,allocatable,dimension(:,:) :: feat,feat2,feat22_type
       real*8,allocatable,dimension(:,:,:) :: feat_type,feat2_type
       real*8,allocatable,dimension(:,:) :: feat2_group
       real*8,allocatable,dimension(:) :: energy_group
       integer,allocatable,dimension(:) :: num_neigh,num,num_atomtype
       integer,allocatable,dimension(:) :: num_neight
       integer,allocatable,dimension(:,:) :: list_neigh,ind_type

       real*8,allocatable,dimension(:,:,:,:) :: dfeat,dfeat2
       real*8,allocatable,dimension(:,:,:) :: dfeat_type,dfeat2_type

       real*8,allocatable,dimension(:,:) :: AA,AA_tmp
       real*8,allocatable,dimension(:) :: BB,BB_tmp

       real*8,allocatable,dimension(:,:,:) :: Gfeat_type
       real*8,allocatable,dimension(:,:) :: Gfeat_tmp

       real*8,allocatable,dimension(:,:,:) :: AA_type
       real*8,allocatable,dimension(:,:) :: BB_type

       real*8,allocatable,dimension(:,:) :: SS_tmp,SS_tmp2

       integer,allocatable,dimension(:) :: ipiv

       real*8,allocatable,dimension(:,:) :: w_feat
       real*8,allocatable,dimension(:,:,:) :: feat2_ref

       real*8,allocatable,dimension(:,:,:) :: PV
       real*8,allocatable,dimension(:,:) :: feat2_shift,feat2_scale


       real*8,allocatable,dimension(:,:) :: WW,VV,QQ
       real*8,allocatable,dimension(:,:,:,:) :: SS

       real*8,allocatable,dimension(:,:) :: Gfeat2,dGfeat2

       real*8,allocatable,dimension(:,:) :: force

     
       real*8,allocatable,dimension(:,:) :: xatom
       real*8,allocatable,dimension(:) :: rad_atom
       real*8,allocatable,dimension(:,:,:) :: wp_atom
       real*8 AL(3,3),pi,dE,dFx,dFy,dFz,AL_tmp(3,3)

       real*8,allocatable,dimension(:,:) :: xatom_tmp

 
       integer,allocatable,dimension(:) :: num_inv
       integer,allocatable,dimension(:,:) :: index_inv,index_inv2

       integer,allocatable,dimension(:) :: nfeat1,nfeat2,
     &    nfeat2i
       integer,allocatable,dimension(:,:) :: nfeat,ipos_feat

       real*8, allocatable, dimension (:,:) :: dfeat_tmp
       real*8, allocatable, dimension (:,:) :: feat_ftype
       integer,allocatable, dimension (:) :: iat_tmp,jneigh_tmp,
     &   ifeat_tmp
       integer num_tmp,jj
       ! character(len=80) dfeat_n(400)
       character(len=1000) trainSetFileDir(400)
       character(len=1000) trainSetDir
       character(len=1000) MOVEMENTDir,dfeatDir,infoDir,trainDataDir,
     &   MOVEMENTallDir
       integer sys_num,sys
       integer nfeat1tm(100),ifeat_type(100),nfeat1t(100)

        real*8, allocatable,dimension (:) :: bond_length,bond_alpha
        real*8, allocatable,dimension (:) :: Ebond
        real*8, allocatable,dimension (:) :: Ebond2,bond_alpha2
        integer, allocatable,dimension (:) :: itype_bond
        integer, allocatable,dimension (:,:) :: iat_bond
        integer nbond,ibond,nbond_type
        real*8 weight_system(100),weight_tmp
        integer num_system
        
        character(len=20) str

  

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



       open(10,file="fit_linearMM.input")
       rewind(10)
       read(10,*) ntype,m_neigh
       allocate(itype_atom(ntype))
       allocate(nfeat1(ntype))
       allocate(nfeat2(ntype))
       allocate(nfeat2i(ntype))
       allocate(rad_atom(ntype))
       allocate(wp_atom(ntype,ntype,2))
       wp_atom=0.d0
       do i=1,ntype
       read(10,*) itype_atom(i)
!     &  rad_atom(i),wp_atom(i)
       enddo
       read(10,*) weight_E,weight_E0,weight_F,delta
       read(10,*) dwidth
       read(10,*) num_system
       do ii=1,num_system
       read(10,*) weight_system(ii)
       enddo
       close(10)

       open(10,file="vdw_fitB.ntype")
       rewind(10)
       read(10,*) ntype_t,nterm
       if(ntype_t.ne.ntype) then
       write(6,*) "ntype not same in vwd_fitB.ntype,something wrong"
       stop
       endif
       do itype1=1,ntype
       read(10,*) itype_t,rad_atom(itype1),E_ave_vdw,
     &            ((wp_atom(i,itype1,j1),i=1,ntype),j1=1,nterm)
       enddo
       close(10)

        open(10,file="bond.molecule")
        rewind(10)
        read(10,*) nbond,nbond_type,Etot0M
        allocate(bond_length(nbond))
        allocate(bond_alpha(nbond))
        allocate(iat_bond(2,nbond))
        allocate(Ebond(nbond))
        allocate(itype_bond(nbond))
        do i=1,nbond
        read(10,*) itype_bond(i),iat_bond(1,i),iat_bond(2,i),
     &   bond_length(i), bond_alpha(i)
        enddo
        close(10)


       open(10,file="feat.info")
       rewind(10)
       read(10,*) iflag_PCA   ! this can be used to turn off degmm part
       read(10,*) nfeat_type
       do kkk=1,nfeat_type
         read(10,*) ifeat_type(kkk)   ! the index (1,2,3) of the feature type
       enddo
       read(10,*) ntype_tmp
        if(ntype_tmp.ne.ntype) then
          write(6,*) 
     & "ntype of atom not same, fit_linearMM.input, feat.info, stop"
          write(6,*) ntype,ntype_tmp
          stop
         endif
       allocate(nfeat(ntype,nfeat_type))
       allocate(ipos_feat(ntype,nfeat_type))
        do i=1,ntype
         read(10,*) iatom_tmp,nfeat1(i),nfeat2(i)   ! these nfeat1,nfeat2 include all ftype
          if(iatom_tmp.ne.itype_atom(i)) then
          write(6,*) "iatom not same, fit_linearMM.input, feat.info"
          write(6,*) iatom_tmp,itype_atom(i)
          stop
          endif
        enddo
        
        do ii=1,ntype
        read(10,*) (nfeat(ii,kkk),kkk=1,nfeat_type)
        enddo
        close(10)

!   nfeat1(ii) the total (all iftype) num of feature for iatom type ii (sum_kk nfeat(ii,kk))
!   nfeat2(ii) the total num of PCA feature for iatom type ii
        
        do ii=1,ntype
        ipos_feat(ii,1)=0
        do kkk=2,nfeat_type
        ipos_feat(ii,kkk)=ipos_feat(ii,kkk-1)+nfeat(ii,kkk-1)
        enddo
        enddo
 


cccccccc Right now, nfeat1,nfeat2,for different types
cccccccc must be the same. We will change that later, allow them 
cccccccc to be different
       nfeat1m=0   ! the original feature
       nfeat2m=0   ! the new PCA, PV feature
       nfeat2tot=0 ! tht total feature of diff atom type
       nfeat2i=0   ! the starting point
       nfeat2i(1)=0
       do i=1,ntype
       if(nfeat1(i).gt.nfeat1m) nfeat1m=nfeat1(i)
       if(nfeat2(i).gt.nfeat2m) nfeat2m=nfeat2(i)
       nfeat2tot=nfeat2tot+nfeat2(i)
       if(i.gt.1) then
       nfeat2i(i)=nfeat2i(i-1)+nfeat2(i-1)
       endif

       enddo


       allocate(w_feat(nfeat2m,ntype))
       do itype=1,ntype
       
       write(str,*) itype
       str=adjustl(str)
       
       open(10,file="weight_feat."//trim(str))
       rewind(10)
       do j=1,nfeat2(itype)
       read(10,*) j1,w_feat(j,itype)
       w_feat(j,itype)=w_feat(j,itype)**2
       enddo
       close(10)
       enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       allocate(PV(nfeat1m,nfeat2m,ntype))
       allocate(feat2_shift(nfeat2m,ntype))
       allocate(feat2_scale(nfeat2m,ntype))
       do itype=1,ntype
       
       write(str,*) itype
       str=adjustl(str)
       
       open(11,file="feat_PV."//trim(str),form="unformatted")
       rewind(11)
       read(11) nfeat1_tmp,nfeat2_tmp
       if(nfeat2_tmp.ne.nfeat2(itype)) then
       write(6,*) "nfeat2.not.same,feat2_ref",itype,nfeat2_tmp,
     &   nfeat2(itype)
       stop
       endif
       if(nfeat1_tmp.ne.nfeat1(itype)) then
       write(6,*) "nfeat1.not.same,feat2_ref",itype,nfeat1_tmp,
     &   nfeat1(itype)
       stop
       endif
       read(11) PV(1:nfeat1(itype),1:nfeat2(itype),itype)
       read(11) feat2_shift(1:nfeat2(itype),itype)
       read(11) feat2_scale(1:nfeat2(itype),itype)
       close(11)
       enddo
!cccccccccccccccccccccccccccccccccccccccccccccccccc

       allocate(AA(nbond_type,nbond_type))
       allocate(BB(nbond_type))
       allocate(ipiv(nbond_type))

       allocate(num(ntype))
       allocate(num_atomtype(ntype))
       allocate(AA_tmp(nfeat2tot,nfeat2tot))
       allocate(BB_tmp(nfeat2tot))
       allocate(AA_type(nfeat2m,nfeat2m,ntype))
       allocate(BB_type(nfeat2m,ntype))

       ! sys_num=400
       open(13,file="location")
       rewind(13)
       read(13,*) sys_num  !,trainSetDir
       read(13,'(a1000)') trainSetDir
       ! allocate(trainSetFileDir(sys_num))
       do i=1,sys_num
       read(13,'(a1000)') trainSetFileDir(i)    
       enddo
       close(13)
       ! MOVEMENTallDir=trim(trainSetDir)//"/MOVEMENT"
       ! trainDataDir=trim(trainSetDir)//"/trainData.txt"

       if(num_system.ne.sys_num) then
       write(6,*) "n_sys from fit_linearMM.input,location different",
     &     num_system,sys_num
       stop
       endif

       AA=0.d0
       BB=0.d0

       open(99,file="Ebond.pred")
       rewind(99)

       do 900 sys=1,sys_num

       do 777 kkk=1,nfeat_type
       ! MOVEMENTDir=trim(trainSetFileDir(sys))//"/MOVEMENT"
       
       write(str,*) ifeat_type(kkk)
       str=adjustl(str)
       
       dfeatDir=trim(trainSetFileDir(sys))//"/dfeat.fbin.Ftype"
     &          //trim(str)
       open(1000+kkk,file=dfeatDir,action="read",access="stream",
     &     form="unformatted")
       rewind(1000+kkk)
       read(1000+kkk) nimaget,natomt,nfeat1tm(kkk),
     &    m_neight

!      nfeat1tm(kkk) is the max(nfeat(ii,kkk)) for all ii(iatype)
  
       if(kkk.eq.1) then
       nimage=nimaget
       natom=natomt
       m_neigh=m_neight
       else
       if(nimaget.ne.nimage.or.natomt.ne.natom.or.
     &  m_neight.ne.m_neigh) then
       write(6,*) "param changed in diff ifeat_type"
       write(6,*) nimage,natom,m_neigh
       write(6,*) nimaget,natomt,m_neight
       stop
       endif
       endif


       read(1000+kkk) ntype_tmp,(nfeat1t(ii),ii=1,ntype_tmp)
!    This is one etra line, perhaps we don't need it

!! for this kkk_ftype, for each atom type, ii, the num of feature is nfeat1t(ii)
!cccccccccccccccccccccccccccccccccccccccccccccccc
       if(ntype_tmp.ne.ntype) then
       write(6,*) "ntype_tmp.ne.ntype,dfeat.fbin,stop"
       write(6,*) ntype_tmp,ntype
       stop
       endif

        do ii=1,ntype
        if(nfeat1t(ii).ne.nfeat(ii,kkk)) then   ! the num of feat for ii_th iatype, and kkk_th feat type
        write(6,*) "nfeat1t not the same, dfeat.fbin,stop"
        write(6,*) nfeat1t(ii),nfeat(ii,kkk),ii,kkk
        stop
        endif
        enddo

       if(kkk.eq.1) then
         if(sys.ne.1) then
         deallocate(iatom)
         endif
       allocate(iatom(natom))
       endif
       read(1000+kkk) iatom      ! The same for different kkk

777    continue
!ccccccccccccccccccccccccccccccccccccccccccccccccc
       
           
       if (sys.ne.1) then      

       deallocate(iatom_type)
       deallocate(Energy)
       deallocate(Energyt)
       deallocate(num_neight)
       deallocate(feat)
       deallocate(feat2)
       deallocate(feat_type)
       deallocate(feat2_type)
       deallocate(feat22_type)
       deallocate(num_neigh)
       deallocate(list_neigh)
       deallocate(ind_type)
       deallocate(dfeat)
       deallocate(dfeat_type)
       deallocate(dfeat2_type)
       deallocate(dfeat2)
       deallocate(xatom)
       deallocate(feat2_group)
       deallocate(energy_group)

       deallocate(num_inv)
       deallocate(index_inv)
       deallocate(index_inv2)
       deallocate(force)
       deallocate(VV)
       deallocate(SS)

       endif

!cccccccccccccccccccccccccccccccccccccccccccccccc

       allocate(iatom_type(natom))
       allocate(Energy(natom))
       allocate(Energyt(natom))
       allocate(num_neight(natom))
       allocate(feat(nfeat1m,natom))  
! nfeat1m is the max(nfeat1(ii)) for ii(iatype), nfeat1(ii)=sum_kkk nfeat(ii,kkk)
! nfeat1m is the max num of total feature (sum over all feature type)
       allocate(feat2(nfeat2m,natom))
       allocate(feat_type(nfeat1m,natom,ntype))
       allocate(feat2_type(nfeat2m,natom,ntype))
       allocate(feat22_type(nfeat2m,ntype))
       allocate(num_neigh(natom))
       allocate(list_neigh(m_neigh,natom))
       allocate(ind_type(natom,ntype))
       allocate(dfeat(nfeat1m,natom,m_neigh,3))
       allocate(dfeat_type(nfeat1m,natom*m_neigh*3,ntype))
       allocate(dfeat2_type(nfeat2m,natom*m_neigh*3,ntype))
       allocate(dfeat2(nfeat2m,natom,m_neigh,3))
       allocate(xatom(3,natom))
       allocate(feat2_group(nfeat2tot,natom))
       allocate(energy_group(natom))

       dfeat=0.d0
       dfeat_type=0.d0

       allocate(num_inv(natom))
       allocate(index_inv(3*m_neigh,natom))
       allocate(index_inv2(3*m_neigh,natom))
       allocate(force(3,natom))
       allocate(VV(nfeat2tot,3*natom))
       allocate(SS(nfeat2m,natom,3,ntype))



       pi=4*datan(1.d0)


       do i=1,natom
        iitype=0
        do itype=1,ntype
        if(itype_atom(itype).eq.iatom(i)) then
        iitype=itype
        endif
        enddo
        if(iitype.eq.0) then
        write(6,*) "this type not found", iatom(i)
        endif
        iatom_type(i)=iitype
      enddo


       num_atomtype=0
       do i=1,natom
       itype=iatom_type(i)
       num_atomtype(itype)=num_atomtype(itype)+1
       enddo



!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


       num_tmp=0


       do 3000 image=1,nimage

       write(6,*) "image=",image,sys

       AA_type=0.d0
       BB_type=0.d0


!-----------------------------------------------------------------
!-----------------------------------------------------------------
!---- read in the feature from different kkk, and put them together
!-----------------------------------------------------------------
       dfeat(:,:,:,:)=0.0
       feat(:,:)=0.0
       do 778 kkk=1,nfeat_type

       allocate(feat_ftype(nfeat1tm(kkk),natom))
       read(1000+kkk) energy   ! repeated
       read(1000+kkk) force    ! repeated
       read(1000+kkk) feat_ftype

       if(kkk.eq.1) then
       energyt=energy
       else
        diff=0.d0
        do ii=1,natom
        diff=diff+abs(energyt(ii)-energy(ii))
        enddo
        if(diff.gt.1.E-9) then
        write(6,*) "energy Ei not the same for diff ifeature type, stop"
        stop
        endif
       endif


       do iat=1,natom
       itype=iatom_type(iat)
       do ii=1,nfeat(itype,kkk)
       feat(ii+ipos_feat(itype,kkk),iat)=feat_ftype(ii,iat)   ! put different kkk together
       enddo
       enddo
       deallocate(feat_ftype)
       
       read(1000+kkk) num_neigh     ! this is actually the num_neighM (of Rc_M)
       read(1000+kkk) list_neigh    ! this is actually the list_neighM (of Rc_M)
!    the above should be the same for different kkk. 
!    Perhaps we should check it later. Here we proceed without checking 
        if(kkk.eq.1) then
        num_neight=num_neigh
        else
         diff=0.d0
         do ii=1,natom
         diff=diff+abs(num_neight(ii)-num_neigh(ii))
         enddo
         if(diff.gt.1.E-9) then
         write(6,*) "num_neigh not the same for diff ifeature type,stop"
         stop
         endif
        endif


!TODO:
       ! read(10) dfeat
       read(1000+kkk) num_tmp
       allocate(dfeat_tmp(3,num_tmp))
       allocate(iat_tmp(num_tmp))
       allocate(jneigh_tmp(num_tmp))
       allocate(ifeat_tmp(num_tmp))
       read(1000+kkk) iat_tmp
       read(1000+kkk) jneigh_tmp
       read(1000+kkk) ifeat_tmp
       read(1000+kkk) dfeat_tmp
       
       read(1000+kkk) xatom    ! xatom(3,natom), repeated for diff kkk
       read(1000+kkk) AL       ! AL(3,3), repeated for diff kkk




       do jj=1,num_tmp


       itype2=iatom_type(list_neigh(jneigh_tmp(jj),iat_tmp(jj))) ! itype2: the type of the neighbor
       dfeat(ifeat_tmp(jj)+ipos_feat(itype2,kkk),
     &                             iat_tmp(jj),jneigh_tmp(jj),:)
     &          =dfeat_tmp(:,jj)
!  Place dfeat from different iftype into the same dfeat
       enddo
       deallocate(dfeat_tmp)
       deallocate(iat_tmp)
       deallocate(jneigh_tmp)
       deallocate(ifeat_tmp)

778    continue
cccccccccccccccccccccccccccccccccccccccccccccccccc
       
       num=0
       do i=1,natom
       itype=iatom_type(i)
       num(itype)=num(itype)+1
       ind_type(num(itype),itype)=i
       feat_type(:,num(itype),itype)=feat(:,i)
!  we have to seperate the feature into different iatype, since they have different PV
!  The num of total feature for different iatype is nfeat1(iatype)
       enddo

ccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccc

        Ebond_tot=0.d0
        do ibond=1,nbond

        i=iat_bond(1,ibond)
        j=iat_bond(2,ibond)

        dx1=mod(xatom(1,j)-xatom(1,i)+100.d0,1.d0)
        if(abs(dx1-1).lt.abs(dx1)) dx1=dx1-1
        dx2=mod(xatom(2,j)-xatom(2,i)+100.d0,1.d0)
        if(abs(dx2-1).lt.abs(dx2)) dx2=dx2-1
        dx3=mod(xatom(3,j)-xatom(3,i)+100.d0,1.d0)
        if(abs(dx3-1).lt.abs(dx3)) dx3=dx3-1
        dx=AL(1,1)*dx1+AL(1,2)*dx2+AL(1,3)*dx3
        dy=AL(2,1)*dx1+AL(2,2)*dx2+AL(2,3)*dx3
        dz=AL(3,1)*dx1+AL(3,2)*dx2+AL(3,3)*dx3
        dd=dsqrt(dx**2+dy**2+dz**2)

        Ebond(ibond)=(dd-bond_length(ibond))**2
        Ebond_tot=Ebond_tot+Ebond(ibond)*bond_alpha(ibond)

        enddo  ! ibond


!ccccccccccccccccccccccccccccccccccccccccccc
        Etot=0.d0
        do i=1,natom
        Etot=Etot+energy(i)
        enddo

        Etot=Etot-Etot0M

        write(99,*) Etot,Ebond_tot

        do ib1=1,nbond
        itype=itype_bond(ib1)
        BB(itype)=BB(itype)+Etot*Ebond(ib1)*weight_system(sys)
        enddo

        do ib1=1,nbond
        do ib2=1,nbond
        itype1=itype_bond(ib1)
        itype2=itype_bond(ib2)

        AA(itype1,itype2)=AA(itype1,itype2)+Ebond(ib1)*Ebond(ib2)*
     &     weight_system(sys)
        enddo
        enddo
      
ccccccccccccccccccccccccccccccccccccccccccc

3000   continue




       close(10)



900    continue

       close(99)


       do k1=1,nbond_type
       AA(k1,k1)=AA(k1,k1)+delta
       enddo

       do kkk=1,nfeat_type
       close(1000+kkk)
       enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       
       call dgesv(nbond_type,1,AA,nbond_type,ipiv,BB,
     &     nbond_type,info)  
!cccccccccccccccccccccccccccccccccccccccccccccccccc
       open(10,file="bond.fitting")
       rewind(10)
       write(10,*) nbond,nbond_type,Etot0M
       do i=1,nbond
       itype=itype_bond(i)
       write(10,"(3(i5,1x),2(E17.10,1x))") itype_bond(i),iat_bond(1,i),
     &   iat_bond(2,i), bond_length(i),BB(itype)
       enddo
       close(10)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       

       stop
       end

       
