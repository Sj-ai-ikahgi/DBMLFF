module calc_ftype1_ppp
    
    ! ---------- USAGE!!! ----------
    ! call load_model_type1_ppp_from_mem (mod_ppp.f90)
    !     -> call set_image_info_type1_ppp
    ! call calc_pol_dir (mod_calc_pol_dir.f90)
    ! call load_neighbore (find_neighbore_ppp.f90) or load_neighbore_from_mem (mod_ppp.f90)
    !     -> call calc_dR_neigh(find_neighbore_ppp.f90)
    
    use mod_mpi
    use mod_md, only: A_AU_1
    use mod_param_densityFF, only: indb,ipol
    use mod_calc_pol_dir, only: dxpp,dxpp1_diat2,dxpp2_diat4,dxpp3_diat2,dxpp3_diat4
    use find_neighbore_ppp, only: num_neigh,list_neigh,dR_neigh,num_neigh_alltype, &
        list_neigh_alltype,ind_all_neigh
    
    implicit none
    
    ! ---------- load module ----------
    ! loaded from mem in load_model_type1_ppp_from_mem (mod_ppp.f90)
    ! one should call load_model_type1_ppp_from_mem (mod_ppp.f90) and then use this module
    real*8 Rc_M
    integer m_neigh
    integer ntype
    integer iat_type(100)
    real*8 Rc_type(100)
    integer n2b_type(100)
    integer iflag_ftype
    integer n2bm
    integer nfeat0m
    real*8,allocatable,dimension(:,:,:) :: grid2_2
    
    ! passed from set_image_info_type1_ppp
    integer natom
    integer,allocatable,dimension(:) :: iatom
    
    ! ---------- public ----------
    real*8,allocatable,dimension(:,:) :: feat_M1
    real*8,allocatable,dimension(:,:,:,:) :: dfeat_M1
    integer nfeat0M1
    
    real*8,allocatable,dimension(:,:,:,:) :: dfeat_iat0             ! dfeat/dR_self
    real*8,allocatable,dimension(:,:,:,:,:,:) :: ddfeat_iat1        ! d(dfeat/dR_self)/dR_neigh, d(dfeat/dR_self)/dR_self=-d(dfeat/dR_self)/dR_neigh
    
    ! ---------- private ----------
    integer,allocatable,dimension(:) :: itype_atom

    private :: itype_atom
    
contains
    subroutine set_image_info_type1_ppp(iatom_in,natom_in)
        integer,intent(in) :: iatom_in(natom_in)
        integer,intent(in) :: natom_in
        
        integer i,j
        
        ! ----------
        
        if(allocated(iatom)) then
            deallocate(iatom)
            deallocate(itype_atom)
        endif
        allocate(iatom(natom_in))
        allocate(itype_atom(natom_in))
        
        ! ----------
        
        natom=natom_in
        iatom(1:natom)=iatom_in(1:natom)
        
        ! ----------
        
        do i=1,natom
            do j=1,ntype
                if(iatom(i).eq.iat_type(j)) then
                    itype_atom(i)=j
                endif
            enddo
            if (itype_atom(i).eq.0) then
                write(6,*) "this atom type didn't found",itype_atom(i)
                stop
            endif
        enddo
    end subroutine set_image_info_type1_ppp
    
    
    subroutine calc_orig_dfeat_and_ddfeat_2b_ppp()
        ! --------------------------------------------------
        ! added_by_fx
        ! --------------------------------------------------
        ! O,M   dfeat_iat0      ! dfeat/dR_self
        ! O,M   ddfeat_iat1     ! d(dfeat/dR_self)/dR_neigh, d(dfeat/dR_self)/dR_self=-d(dfeat/dR_self)/dR_neigh
        ! --------------------------------------------------
        implicit none
        
        integer i,j,k
        integer iat,itype
        integer itype0
        integer num
        integer jj
        
        real*8 pi,pi2
        
        real*8 d,dd
        real*8 drange,x,y,f1,y2
        real*8 temp
        
        ! ----------
        
        if(allocated(dfeat_iat0)) then
            deallocate(dfeat_iat0)
            deallocate(ddfeat_iat1)
        endif
        allocate(dfeat_iat0(3,n2bm,ntype,natom))
        allocate(ddfeat_iat1(3,m_neigh,3,n2bm,ntype,natom))
        
        ! ----------
        
        pi=4*datan(1.d0)
        pi2=2*pi
        
        dfeat_iat0=0.d0
        do 3000 iat=1,natom
            itype0=itype_atom(iat)
            
            do 1000 itype=1,ntype
                do 1000 j=1,num_neigh(itype,iat)
                    jj=ind_all_neigh(j,itype,iat)
                    
                    dd=dR_neigh(1,j,itype,iat)**2+dR_neigh(2,j,itype,iat)**2+dR_neigh(3,j,itype,iat)**2
                    d=dsqrt(dd)
                    
                    do k=1,n2b_type(itype0)
                        if(d.ge.grid2_2(1,k,itype0).and.d.lt.grid2_2(2,k,itype0)) then
                            drange=grid2_2(2,k,itype0)-grid2_2(1,k,itype0)
                            x=(d-grid2_2(1,k,itype0))/drange
                            y=(x-0.5d0)*pi2
                            f1=0.5d0*(cos(y)+1)
                            y2=-pi*sin(y)/(d*drange)
                            
                            dfeat_iat0(:,k,itype,iat)=dfeat_iat0(:,k,itype,iat)-y2*dR_neigh(:,j,itype,iat)
                            
                            temp=pi*(drange*sin(y)-pi2*d*cos(y))/(dd*d*drange**2)
                            ddfeat_iat1(1,jj,1,k,itype,iat)=-temp*dR_neigh(1,j,itype,iat)*dR_neigh(1,j,itype,iat)-y2
                            ddfeat_iat1(2,jj,1,k,itype,iat)=-temp*dR_neigh(2,j,itype,iat)*dR_neigh(1,j,itype,iat)
                            ddfeat_iat1(3,jj,1,k,itype,iat)=-temp*dR_neigh(3,j,itype,iat)*dR_neigh(1,j,itype,iat)
                            ddfeat_iat1(1,jj,2,k,itype,iat)=-temp*dR_neigh(1,j,itype,iat)*dR_neigh(2,j,itype,iat)
                            ddfeat_iat1(2,jj,2,k,itype,iat)=-temp*dR_neigh(2,j,itype,iat)*dR_neigh(2,j,itype,iat)-y2
                            ddfeat_iat1(3,jj,2,k,itype,iat)=-temp*dR_neigh(3,j,itype,iat)*dR_neigh(2,j,itype,iat)
                            ddfeat_iat1(1,jj,3,k,itype,iat)=-temp*dR_neigh(1,j,itype,iat)*dR_neigh(3,j,itype,iat)
                            ddfeat_iat1(2,jj,3,k,itype,iat)=-temp*dR_neigh(2,j,itype,iat)*dR_neigh(3,j,itype,iat)
                            ddfeat_iat1(3,jj,3,k,itype,iat)=-temp*dR_neigh(3,j,itype,iat)*dR_neigh(3,j,itype,iat)-y2
                        else
                            ddfeat_iat1(1:3,jj,1:3,k,itype,iat)=0.d0
                        endif
                    enddo               
            1000 continue
        3000 continue
    end subroutine calc_orig_dfeat_and_ddfeat_2b_ppp
    
    subroutine calc_new_feat_and_dfeat_2b_ppp(itype_mol,proj_type,iatom_or_ibond,ippp)
        implicit none
        
        integer,intent(in) :: itype_mol
        integer,intent(in) :: proj_type  ! 1, atom; 2, bond
        integer,intent(in) :: iatom_or_ibond
        integer,intent(in) :: ippp
        
        integer i,j,k
        integer iat,itype,i2b,ineigh
        integer iatom_nonzero_feat,idir
        integer jj
        integer iat1,iat2,iat3,iat4
        integer jj1,jj2,jj3,jj4
        integer num
        integer itype0
        
        real*8 real_bin(3)
        
        ! ----------
        
        if(allocated(feat_M1)) then
            deallocate(feat_M1)
            deallocate(dfeat_M1)
        endif
        allocate(feat_M1(nfeat0m,natom))
        allocate(dfeat_M1(nfeat0m,natom,m_neigh,3))
        
        ! ----------
        
        iat1=0
        iat2=0
        iat3=0
        iat4=0
        if(proj_type.eq.1) then
            iat1=ipol(1,iatom_or_ibond,itype_mol)
            iat2=ipol(2,iatom_or_ibond,itype_mol)
            iat3=ipol(3,iatom_or_ibond,itype_mol)
            iat4=ipol(4,iatom_or_ibond,itype_mol)
            
            iatom_nonzero_feat=iatom_or_ibond
            idir=(ippp-1)/2+1
        elseif(proj_type.eq.2) then
            iat1=indb(1,iatom_or_ibond,itype_mol)
            iat2=indb(2,iatom_or_ibond,itype_mol)
            
            if((ippp.eq.1).or.(ippp.eq.2)) iatom_nonzero_feat=iat1
            if((ippp.eq.3).or.(ippp.eq.4)) iatom_nonzero_feat=iat2
            idir=1
        endif
        
        jj1=0
        jj2=0
        jj3=0
        jj4=0
        do i=1,num_neigh_alltype(iatom_nonzero_feat)
            if(list_neigh_alltype(i,iatom_nonzero_feat).eq.iat1) jj1=i
            if(list_neigh_alltype(i,iatom_nonzero_feat).eq.iat2) jj2=i
            if(list_neigh_alltype(i,iatom_nonzero_feat).eq.iat3) jj3=i
            if(list_neigh_alltype(i,iatom_nonzero_feat).eq.iat4) jj4=i
        enddo
        
        if((idir.eq.1).and.((jj1.eq.0).or.(jj2.eq.0))) then
            write(6,*) "(idir.eq.1).and.((jj1.eq.0).or.(jj2.eq.0)),stop"
            stop
        endif
        if((idir.eq.2).and.((jj3.eq.0).or.(jj4.eq.0))) then
            write(6,*) "(idir.eq.2).and.((jj3.eq.0).or.(jj4.eq.0)),stop"
            stop
        endif
        if((idir.eq.3).and.((jj1.eq.0).or.(jj2.eq.0).or.(jj3.eq.0).or.(jj4.eq.0))) then
            write(6,*) "(idir.eq.3).and.((jj1.eq.0).or.(jj2.eq.0).or.(jj3.eq.0).or.(jj4.eq.0)),stop"
            stop
        endif
        
        ! ----------
        
        itype0=itype_atom(iatom_nonzero_feat)
        
        num=0
        do itype=1,ntype
            do i2b=1,n2b_type(itype0)
                num=num+1
                
                ! feat part
                feat_M1(num,iatom_nonzero_feat)=dfeat_iat0(1,i2b,itype,iatom_nonzero_feat)*dxpp(1,idir,iatom_or_ibond,proj_type)+ &
                    dfeat_iat0(2,i2b,itype,iatom_nonzero_feat)*dxpp(2,idir,iatom_or_ibond,proj_type)+ &
                    dfeat_iat0(3,i2b,itype,iatom_nonzero_feat)*dxpp(3,idir,iatom_or_ibond,proj_type)
                
                ! dfeat part
                dfeat_M1(num,iatom_nonzero_feat,:,:)=0.d0
                
                do ineigh=1,num_neigh(itype,iatom_nonzero_feat)
                    jj=ind_all_neigh(ineigh,itype,iatom_nonzero_feat)
                    
                    real_bin(:)=ddfeat_iat1(:,jj,1,i2b,itype,iatom_nonzero_feat)*dxpp(1,idir,iatom_or_ibond,proj_type)+ &
                        ddfeat_iat1(:,jj,2,i2b,itype,iatom_nonzero_feat)*dxpp(2,idir,iatom_or_ibond,proj_type)+ &
                        ddfeat_iat1(:,jj,3,i2b,itype,iatom_nonzero_feat)*dxpp(3,idir,iatom_or_ibond,proj_type)
                        
                    dfeat_M1(num,iatom_nonzero_feat,jj,:)=real_bin(:)
                    dfeat_M1(num,iatom_nonzero_feat,1,:)=dfeat_M1(num,iatom_nonzero_feat,1,:)-real_bin(:)
                enddo
                
                if(idir.eq.1) then
                    real_bin(:)=dfeat_iat0(1,i2b,itype,iatom_nonzero_feat)*dxpp1_diat2(1,:,iatom_or_ibond,proj_type)/A_AU_1+ &
                        dfeat_iat0(2,i2b,itype,iatom_nonzero_feat)*dxpp1_diat2(2,:,iatom_or_ibond,proj_type)/A_AU_1+ &
                        dfeat_iat0(3,i2b,itype,iatom_nonzero_feat)*dxpp1_diat2(3,:,iatom_or_ibond,proj_type)/A_AU_1
                    
                    dfeat_M1(num,iatom_nonzero_feat,jj1,:)=dfeat_M1(num,iatom_nonzero_feat,jj1,:)-real_bin(:)
                    dfeat_M1(num,iatom_nonzero_feat,jj2,:)=dfeat_M1(num,iatom_nonzero_feat,jj2,:)+real_bin(:)
                endif
                if(idir.eq.2) then
                    real_bin(:)=dfeat_iat0(1,i2b,itype,iatom_nonzero_feat)*dxpp2_diat4(1,:,iatom_or_ibond,proj_type)/A_AU_1+ &
                        dfeat_iat0(2,i2b,itype,iatom_nonzero_feat)*dxpp2_diat4(2,:,iatom_or_ibond,proj_type)/A_AU_1+ &
                        dfeat_iat0(3,i2b,itype,iatom_nonzero_feat)*dxpp2_diat4(3,:,iatom_or_ibond,proj_type)/A_AU_1
                    
                    dfeat_M1(num,iatom_nonzero_feat,jj3,:)=dfeat_M1(num,iatom_nonzero_feat,jj3,:)-real_bin(:)
                    dfeat_M1(num,iatom_nonzero_feat,jj4,:)=dfeat_M1(num,iatom_nonzero_feat,jj4,:)+real_bin(:)
                endif
                if(idir.eq.3) then
                    ! PLEASE pay special attention to the meaning of dxpp1_diat2, dxpp2_diat4, dxpp3_diat2, dxpp3_diat4
                    ! for example
                    ! dxpp3_diat2(1,2) is dxpp(1,3)/dx22
                    ! dxpp3_diat4(2,3) is dxpp(2,3)/dx43
                    ! Misusing the first two indexes cause trouble and waste time
                    real_bin(:)=dfeat_iat0(1,i2b,itype,iatom_nonzero_feat)*dxpp3_diat2(1,:,iatom_or_ibond,proj_type)/A_AU_1+ &
                        dfeat_iat0(2,i2b,itype,iatom_nonzero_feat)*dxpp3_diat2(2,:,iatom_or_ibond,proj_type)/A_AU_1+ &
                        dfeat_iat0(3,i2b,itype,iatom_nonzero_feat)*dxpp3_diat2(3,:,iatom_or_ibond,proj_type)/A_AU_1
                    
                    dfeat_M1(num,iatom_nonzero_feat,jj1,:)=dfeat_M1(num,iatom_nonzero_feat,jj1,:)-real_bin(:)
                    dfeat_M1(num,iatom_nonzero_feat,jj2,:)=dfeat_M1(num,iatom_nonzero_feat,jj2,:)+real_bin(:)
                    
                    real_bin(:)=dfeat_iat0(1,i2b,itype,iatom_nonzero_feat)*dxpp3_diat4(1,:,iatom_or_ibond,proj_type)/A_AU_1+ &
                        dfeat_iat0(2,i2b,itype,iatom_nonzero_feat)*dxpp3_diat4(2,:,iatom_or_ibond,proj_type)/A_AU_1+ &
                        dfeat_iat0(3,i2b,itype,iatom_nonzero_feat)*dxpp3_diat4(3,:,iatom_or_ibond,proj_type)/A_AU_1
                    
                    dfeat_M1(num,iatom_nonzero_feat,jj3,:)=dfeat_M1(num,iatom_nonzero_feat,jj3,:)-real_bin(:)
                    dfeat_M1(num,iatom_nonzero_feat,jj4,:)=dfeat_M1(num,iatom_nonzero_feat,jj4,:)+real_bin(:)
                endif
                
            enddo
        enddo
        
        nfeat0M1=nfeat0m
    end subroutine calc_new_feat_and_dfeat_2b_ppp
    
end module calc_ftype1_ppp
