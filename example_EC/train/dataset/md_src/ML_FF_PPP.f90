subroutine ML_FF_PPP(itype_mol,xatom,AL,ppp,ppp1,ppp2,dppp,dppp1,dppp2)
! ------------------------------------------------------------
! added_by_fx
! ------------------------------------------------------------
! Inference code for the ppp0 prediction using PWmat MLFF.
! ------------------------------------------------------------
! I     itype_mol
! I     xatom
! I     AL
! O     ppp
! O     ppp1
! O     ppp2
! O     dppp
! O     dppp1
! O     dppp2
! ------------------------------------------------------------
    
    use mod_mpi
    use mod_data, only: natom_m,nmol,ntype_m,natom_mm
    use mod_md, only: A_AU_1
    use mod_param_densityFF, only: numb,ipol,indb
    use mod_ppp, only: &
        num_ppp_model_atom,num_ppp_model_bond, &
        ind_ppp_model_atom,ind_ppp_model_bond, &
        sign_ppp_model_atom,sign_ppp_model_bond, &
        is_fitted_atom,is_fitted_bond, &
        iatom_m_atom,iatom_m_bond, &
        load_model_type1_ppp_from_mem,load_model_type2_ppp_from_mem,load_neighbore_from_mem
    use calc_ftype1_ppp, only: feat_M1,dfeat_M1,nfeat0M1, &
        set_image_info_type1_ppp,calc_orig_dfeat_and_ddfeat_2b_ppp,calc_new_feat_and_dfeat_2b_ppp, &
        m_neigh
    use calc_ftype2_ppp, only: feat_M2,dfeat_M2,nfeat0M2, &
        set_image_info_type2_ppp,calc_orig_dfeat_and_ddfeat_3b_ppp,calc_new_feat_and_dfeat_3b_ppp
    ! use calc_ftype2_ppp, only: feat_M2,dfeat_M2,nfeat0M2,num_neigh_alltypeM2,list_neigh_alltypeM2, &
        ! set_image_info_type2_ppp,gen_feature_type2_ppp, &
        ! iat11_ftype2,iat22_ftype2,iat33_ftype2,iat44_ftype2,idir_ftype2,path_3b
    use calc_lin_ppp, only: natom_nn,Etot_pred_lin,force_pred_lin, &
        set_image_info_lin_ppp_atom_part,set_image_info_lin_ppp_bond_part, &
        cal_energy_force_lin_ppp_atom_part,cal_energy_force_lin_ppp_bond_part
    use find_neighbore_ppp, only: num_neigh,iat_neigh,list_neigh,dR_neigh,num_neigh_alltype, &
        list_neigh_alltype,calc_dR_neigh
    
    implicit double precision (a-h,o-z)
    
    logical :: iflag_timer=.false.
    real*8 t_ppp_start,t_ppp_end,dt_ppp
    real*8 t_atom_start,t_atom_end,t_atom_read_file_start,t_atom_read_file_end, &
        t_atom_load_model_start,t_atom_load_model_end,t_atom_set_info_start,t_atom_set_info_end, &
        t_atom_gen_feature_start,t_atom_gen_feature_end,t_atom_inference_start,t_atom_inference_end, &
        dt_atom,dt_atom_read_file,dt_atom_load_model,dt_atom_set_info, &
        dt_atom_gen_feature,dt_atom_inference
    real*8 t_bond_start,t_bond_end,t_bond_read_file_start,t_bond_read_file_end, &
        t_bond_load_model_start,t_bond_load_model_end,t_bond_set_info_start,t_bond_set_info_end, &
        t_bond_gen_feature_start,t_bond_gen_feature_end,t_bond_inference_start,t_bond_inference_end, &
        dt_bond,dt_bond_read_file,dt_bond_load_model,dt_bond_set_info, &
        dt_bond_gen_feature,dt_bond_inference
    
    real*8 xatom(3,natom_mm)
    real*8 AL(3,3)
    real*8 ppp(3,2*natom_mm)
    real*8 ppp1(2,200)
    real*8 ppp2(2,200)
    real*8 dppp(3,2*natom_mm,natom_mm,3)
    real*8 dppp1(3,200,2,200)
    real*8 dppp2(3,200,2,200)
    
    integer nfeat0
    real*8,allocatable,dimension (:,:) :: feat
    real*8,allocatable,dimension (:,:,:,:) :: dfeat
    integer ii,jj,iat
    logical is_reset

    integer iatom
    integer ibond
    integer ippp
    integer ind_ppp_model
    integer iatom_nonzero_feat
    
    character(len=20) ind_ppp_model_string
    character(len=20) iatom_string
    character(len=20) ibond_string
    character(len=100) first_line
    character(len=1000) path
    
    integer iatom_m(natom_mm)
    integer ntype_ML
        
    character(len=20) str
    character(len=20) str2
    
    if(iflag_timer) then
        t_ppp_start=0.d0
        t_ppp_end=0.d0
        dt_ppp=0.d0
        
        t_atom_start=0.d0
        t_atom_end=0.d0
        t_atom_read_file_start=0.d0
        t_atom_read_file_end=0.d0
        t_atom_load_model_start=0.d0
        t_atom_load_model_end=0.d0
        t_atom_set_info_start=0.d0
        t_atom_set_info_end=0.d0
        t_atom_gen_feature_start=0.d0
        t_atom_gen_feature_end=0.d0
        t_atom_inference_start=0.d0
        t_atom_inference_end=0.d0
        dt_atom=0.d0
        dt_atom_read_file=0.d0
        dt_atom_load_model=0.d0
        dt_atom_set_info=0.d0
        dt_atom_gen_feature=0.d0
        dt_atom_inference=0.d0
        
        t_bond_start=0.d0
        t_bond_end=0.d0
        t_bond_read_file_start=0.d0
        t_bond_read_file_end=0.d0
        t_bond_load_model_start=0.d0
        t_bond_load_model_end=0.d0
        t_bond_set_info_start=0.d0
        t_bond_set_info_end=0.d0
        t_bond_gen_feature_start=0.d0
        t_bond_gen_feature_end=0.d0
        t_bond_inference_start=0.d0
        t_bond_inference_end=0.d0
        dt_bond=0.d0
        dt_bond_read_file=0.d0
        dt_bond_load_model=0.d0
        dt_bond_set_info=0.d0
        dt_bond_gen_feature=0.d0
        dt_bond_inference=0.d0
    endif
    
    if(iflag_timer) t_ppp_start=mpi_wtime()
    
    ppp=0.d0
    ppp1=0.d0
    ppp2=0.d0
    dppp=0.d0
    dppp1=0.d0
    dppp2=0.d0
    
    natom_nn=natom_m(itype_mol)  ! in mod_linear
    
    ! ---------- load neighbore list and precalculate orig dfeat, ddfeat ----------
    call load_neighbore_from_mem(itype_mol)
    call calc_dR_neigh(xatom,AL)
    
    ! ---------- load 2b and 3b model from mem ----------
    call load_model_type1_ppp_from_mem(itype_mol)
    call load_model_type2_ppp_from_mem(itype_mol)
    
    ! ---------- set atom type ----------
    iatom_m(:)=iatom_m_atom(:,1,1,itype_mol)
    call set_image_info_type1_ppp(iatom_m(1),natom_m(itype_mol))
    call set_image_info_type2_ppp(iatom_m(1),natom_m(itype_mol))
    
    ! ---------- precalculate dfeat/dR_self and d(dfeat/dR_self)/dR ----------
    call calc_orig_dfeat_and_ddfeat_2b_ppp()
    call calc_orig_dfeat_and_ddfeat_3b_ppp()
    
    if(iflag_timer) t_atom_start=mpi_wtime()
    do 1001 iatom=1,natom_m(itype_mol)
        do 1000 ippp=1,6
            ind_ppp_model=ind_ppp_model_atom(ippp,iatom,itype_mol)
            write(ind_ppp_model_string,*) ind_ppp_model
            ind_ppp_model_string=adjustl(ind_ppp_model_string)
            
            write(str,*) itype_mol
            str=adjustl(str)
            
            write(str2,*) (ippp-1)/2+1
            str2=adjustl(str2)
            
            write (iatom_string,*) iatom
            iatom_string = adjustl(iatom_string)
            
            if(iflag_timer) t_atom_read_file_start=mpi_wtime()
            ! open(400,file='mol.'//trim(str)//'/'//'ppp_model'//'/'//'atom_'//trim(ind_ppp_model_string)//'/'//'fit_linearMM.input')
            ! rewind(400)
            ! read(400,*) first_line
            ! close(400)
            if(is_fitted_atom(ind_ppp_model,itype_mol).eq..false.) goto 1000
            
            ! open(400,file='mol.'//trim(str)//'/'//'ppp_model'//'/'//'atom_'//trim(ind_ppp_model_string)//'/'//'fit_linearMM.input')
            ! rewind(400)
            ! read(400,*) ntype_ML
            ! close(400)
            ! ntype_ML=ntype_ML_atom(ind_ppp_model,itype_mol)
            
            ! open(400,file='mol.'//trim(str)//'/'//'ppp_neighbore'//'/'//'atom_'//trim(iatom_string)//'_'//trim(str2)//'/'//'input'//'/'//'MOVEMENT.type')
            ! rewind(400)
            ! read(400,*) natom_movement_type
            ! if(natom_movement_type.ne.natom_m(itype_mol)) then
                ! write(6,*) "natom_movement_type.ne.natom_m(itype_mol),stop",natom_movement_type,natom_m(itype_mol)
                ! stop
            ! endif
            ! do i=1,natom_m(itype_mol)
                ! read(400,*) ielement,iatom_m(i)
            ! enddo
            ! close(400)
            iatom_m(:)=iatom_m_atom(:,(ippp-1)/2+1,iatom,itype_mol)
            if(iflag_timer) t_atom_read_file_end=mpi_wtime()
            if(iflag_timer) dt_atom_read_file=dt_atom_read_file+t_atom_read_file_end-t_atom_read_file_start
            
            ! path='mol.'//trim(str)//'/'//'ppp_model'//'/'//'atom_'//trim(ind_ppp_model_string)
            if(iflag_timer) t_atom_load_model_start=mpi_wtime()
            ! call load_model_lin(ntype_ML,path)
            ! call load_model_type1_ppp(path)
            ! call load_model_type2_ppp(path)
            ! call load_model_lin_ppp_atom_part(ippp,iatom,itype_mol)
            ! call load_model_type1_ppp_from_mem(itype_mol)
            ! call load_model_type2_ppp_from_mem(itype_mol)
            if(iflag_timer) t_atom_load_model_end=mpi_wtime()
            if(iflag_timer) dt_atom_load_model=dt_atom_load_model+t_atom_load_model_end-t_atom_load_model_start
            
            is_reset= .true.
            if(iflag_timer) t_atom_set_info_start=mpi_wtime()
            call set_image_info_lin_ppp_atom_part(iatom_m(1),is_reset,natom_m(itype_mol),ippp,iatom,itype_mol)       ! iatom_m(1,itype_mol) here is strange but correct. Is this because when passing an array, one only need to pass the reference of the first element?
            ! call set_image_info_type1_ppp(iatom_m(1),natom_m(itype_mol))     ! To be consistent with ML_FF_EF.f90, I do not make any change here.
            ! call set_image_info_type2_ppp(iatom_m(1),is_reset,natom_m(itype_mol))     ! Please be attention!
            if(iflag_timer) t_atom_set_info_end=mpi_wtime()
            if(iflag_timer) dt_atom_set_info=dt_atom_set_info+t_atom_set_info_end-t_atom_set_info_start
            
            ! path_3b='mol.'//trim(str)//'/'//'ppp_neighbore'//'/'//'atom_'//trim(iatom_string)//'_'//trim(str2)
            ! iat11_ftype1=ipol(1,iatom,itype_mol)
            ! iat11_ftype2=ipol(1,iatom,itype_mol)
            ! iat22_ftype1=ipol(2,iatom,itype_mol)
            ! iat22_ftype2=ipol(2,iatom,itype_mol)
            ! iat33_ftype1=ipol(3,iatom,itype_mol)
            ! iat33_ftype2=ipol(3,iatom,itype_mol)
            ! iat44_ftype1=ipol(4,iatom,itype_mol)
            ! iat44_ftype2=ipol(4,iatom,itype_mol)
            ! if((ippp.eq.1).or.(ippp.eq.2)) then
                ! idir_ftype1=1
                ! idir_ftype2=1
            ! endif
            ! if((ippp.eq.3).or.(ippp.eq.4)) then
                ! idir_ftype1=2
                ! idir_ftype2=2
            ! endif
            ! if((ippp.eq.5).or.(ippp.eq.6)) then
                ! idir_ftype1=3
                ! idir_ftype2=3
            ! endif
            if(iflag_timer) t_atom_gen_feature_start=mpi_wtime()
            ! ---------- project to get new feat and dfeat ----------
            call calc_new_feat_and_dfeat_2b_ppp(itype_mol,1,iatom,ippp)
            call calc_new_feat_and_dfeat_3b_ppp(itype_mol,1,iatom,ippp)
            if(iflag_timer) t_atom_gen_feature_end=mpi_wtime()
            if(iflag_timer) dt_atom_gen_feature=dt_atom_gen_feature+t_atom_gen_feature_end-t_atom_gen_feature_start
            
            nfeat0=nfeat0M1+nfeat0M2
            
            iatom_nonzero_feat=iatom
            
            allocate(feat(nfeat0,natom_nn))
            allocate(dfeat(nfeat0,natom_nn,m_neigh,3))
            
            ! do iat=1,natom_nn
                ! do ii=1,nfeat0M1
                    ! feat(ii,iat)=feat_M1(ii,iat)
                ! enddo
                ! do ii=1,nfeat0M2
                    ! feat(ii+nfeat0M1,iat)=feat_M2(ii,iat)
                ! enddo
            ! enddo
            
            feat(1:nfeat0M1,iatom_nonzero_feat)=feat_M1(1:nfeat0M1,iatom_nonzero_feat)
            feat(1+nfeat0M1:nfeat0M2+nfeat0M1,iatom_nonzero_feat)=feat_M2(1:nfeat0M2,iatom_nonzero_feat)
            
            ! do jj=1,m_neigh
                ! do iat=1,natom_nn
                    ! do ii=1,nfeat0M1
                        ! dfeat(ii,iat,jj,1)=dfeat_M1(ii,iat,jj,1)
                        ! dfeat(ii,iat,jj,2)=dfeat_M1(ii,iat,jj,2)
                        ! dfeat(ii,iat,jj,3)=dfeat_M1(ii,iat,jj,3)
                    ! enddo
                ! enddo
                ! do iat=1,natom_nn
                    ! do ii=1,nfeat0M2
                        ! dfeat(ii+nfeat0M1,iat,jj,1)=dfeat_M2(ii,iat,jj,1)
                        ! dfeat(ii+nfeat0M1,iat,jj,2)=dfeat_M2(ii,iat,jj,2)
                        ! dfeat(ii+nfeat0M1,iat,jj,3)=dfeat_M2(ii,iat,jj,3)
                    ! enddo
                ! enddo
            ! enddo
            
            dfeat(1:nfeat0M1,iatom_nonzero_feat,1:m_neigh,1:3)=dfeat_M1(1:nfeat0M1,iatom_nonzero_feat,1:m_neigh,1:3)
            dfeat(1+nfeat0M1:nfeat0M2+nfeat0M1,iatom_nonzero_feat,1:m_neigh,1:3)=dfeat_M2(1:nfeat0M2,iatom_nonzero_feat,1:m_neigh,1:3)
            
            if(iflag_timer) t_atom_inference_start=mpi_wtime()
            ! if(iflag_model_ppp.eq.1) then
                call cal_energy_force_lin_ppp_atom_part(feat,dfeat,num_neigh_alltype,  &
                    list_neigh_alltype,AL,xatom,natom_nn,nfeat0,m_neigh, &
                    ippp,iatom,itype_mol)
                ppp((ippp-1)/2+1,(iatom-1)*2+mod(ippp-1,2)+1)=Etot_pred_lin*sign_ppp_model_atom(ippp,iatom,itype_mol)
                do iat=1,natom_nn
                    dppp((ippp-1)/2+1,(iatom-1)*2+mod(ippp-1,2)+1,iat,:)=force_pred_lin(:,iat)*sign_ppp_model_atom(ippp,iatom,itype_mol)
                enddo
            ! endif
            if(iflag_timer) t_atom_inference_end=mpi_wtime()
            if(iflag_timer) dt_atom_inference=dt_atom_inference+t_atom_inference_end-t_atom_inference_start
            
            deallocate(feat)
            deallocate(dfeat)
        1000 continue
    1001 continue
    if(iflag_timer) t_atom_end=mpi_wtime()
    if(iflag_timer) dt_atom=t_atom_end-t_atom_start
    
    if(iflag_timer) t_bond_start=mpi_wtime()
    do 2001 ibond=1,numb(itype_mol)
        do 2000 ippp=1,4
            ind_ppp_model=ind_ppp_model_bond(ippp,ibond,itype_mol)
            write(ind_ppp_model_string,*) ind_ppp_model
            ind_ppp_model_string=adjustl(ind_ppp_model_string)
            
            write(str,*) itype_mol
            str=adjustl(str)
            
            write(str2,*) (ippp-1)/2+1
            str2=adjustl(str2)
            
            write (ibond_string,*) ibond
            ibond_string = adjustl(ibond_string)
            
            if(iflag_timer) t_bond_read_file_start=mpi_wtime()
            ! open(400,file='mol.'//trim(str)//'/'//'ppp_model'//'/'//'bond_'//trim(ind_ppp_model_string)//'/'//'fit_linearMM.input')
            ! rewind(400)
            ! read(400,*) first_line
            ! close(400)
            if(is_fitted_bond(ind_ppp_model,itype_mol).eq..false.) goto 2000
            
            ! open(400,file='mol.'//trim(str)//'/'//'ppp_model'//'/'//'bond_'//trim(ind_ppp_model_string)//'/'//'fit_linearMM.input')
            ! rewind(400)
            ! read(400,*) ntype_ML
            ! close(400)
            ! ntype_ML=ntype_ML_bond(ind_ppp_model,itype_mol)
            
            ! open(400,file='mol.'//trim(str)//'/'//'ppp_neighbore'//'/'//'bond_'//trim(ibond_string)//'_'//trim(str2)//'/'//'input'//'/'//'MOVEMENT.type')
            ! rewind(400)
            ! read(400,*) natom_movement_type
            ! if(natom_movement_type.ne.natom_m(itype_mol)) then
                ! write(6,*) "natom_movement_type.ne.natom_m(itype_mol),stop",natom_movement_type,natom_m(itype_mol)
                ! stop
            ! endif
            ! do i=1,natom_m(itype_mol)
                ! read(400,*) ielement,iatom_m(i)
            ! enddo
            ! close(400)
            iatom_m(:)=iatom_m_bond(:,(ippp-1)/2+1,ibond,itype_mol)
            if(iflag_timer) t_bond_read_file_end=mpi_wtime()
            if(iflag_timer) dt_bond_read_file=dt_bond_read_file+t_bond_read_file_end-t_bond_read_file_start
            
            ! path='mol.'//trim(str)//'/'//'ppp_model'//'/'//'bond_'//trim(ind_ppp_model_string)
            if(iflag_timer) t_bond_load_model_start=mpi_wtime()
            ! call load_model_lin(ntype_ML,path)
            ! call load_model_type1_ppp(path)
            ! call load_model_type2_ppp(path)
            ! call load_model_lin_ppp_bond_part(ippp,ibond,itype_mol)
            ! call load_model_type1_ppp_from_mem(itype_mol)
            ! call load_model_type2_ppp_from_mem(itype_mol)
            if(iflag_timer) t_bond_load_model_end=mpi_wtime()
            if(iflag_timer) dt_bond_load_model=dt_bond_load_model+t_bond_load_model_end-t_bond_load_model_start
            
            is_reset= .true.
            if(iflag_timer) t_bond_set_info_start=mpi_wtime()
            call set_image_info_lin_ppp_bond_part(iatom_m(1),is_reset,natom_m(itype_mol),ippp,ibond,itype_mol)       ! iatom_m(1,itype_mol) here is strange but correct. Is this because when passing an array, one only need to pass the reference of the first element?
            ! call set_image_info_type1_ppp(iatom_m(1),natom_m(itype_mol))     ! To be consistent with ML_FF_EF.f90, I do not make any change here.
            ! call set_image_info_type2_ppp(iatom_m(1),is_reset,natom_m(itype_mol))     ! Please be attention!
            if(iflag_timer) t_bond_set_info_end=mpi_wtime()
            if(iflag_timer) dt_bond_set_info=dt_bond_set_info+t_bond_set_info_end-t_bond_set_info_start
            
            ! path_3b='mol.'//trim(str)//'/'//'ppp_neighbore'//'/'//'bond_'//trim(ibond_string)//'_'//trim(str2)
            ! iat11_ftype1=indb(1,ibond,itype_mol)
            ! iat11_ftype2=indb(1,ibond,itype_mol)
            ! iat22_ftype1=indb(2,ibond,itype_mol)
            ! iat22_ftype2=indb(2,ibond,itype_mol)
            ! iat33_ftype1=0
            ! iat33_ftype2=0
            ! iat44_ftype1=0
            ! iat44_ftype2=0
            ! idir_ftype1=1
            ! idir_ftype2=1
            if(iflag_timer) t_bond_gen_feature_start=mpi_wtime()
            call calc_new_feat_and_dfeat_2b_ppp(itype_mol,2,ibond,ippp)
            call calc_new_feat_and_dfeat_3b_ppp(itype_mol,2,ibond,ippp)
            if(iflag_timer) t_bond_gen_feature_end=mpi_wtime()
            if(iflag_timer) dt_bond_gen_feature=dt_bond_gen_feature+t_bond_gen_feature_end-t_bond_gen_feature_start
            
            nfeat0=nfeat0M1+nfeat0M2
            
            if((ippp.eq.1).or.(ippp.eq.2)) iatom_nonzero_feat=indb(1,ibond,itype_mol)
            if((ippp.eq.3).or.(ippp.eq.4)) iatom_nonzero_feat=indb(2,ibond,itype_mol)
            
            allocate(feat(nfeat0,natom_nn))
            allocate(dfeat(nfeat0,natom_nn,m_neigh,3))
            
            ! do iat=1,natom_nn
                ! do ii=1,nfeat0M1
                    ! feat(ii,iat)=feat_M1(ii,iat)
                ! enddo
                ! do ii=1,nfeat0M2
                    ! feat(ii+nfeat0M1,iat)=feat_M2(ii,iat)
                ! enddo
            ! enddo
            
            feat(1:nfeat0M1,iatom_nonzero_feat)=feat_M1(1:nfeat0M1,iatom_nonzero_feat)
            feat(1+nfeat0M1:nfeat0M2+nfeat0M1,iatom_nonzero_feat)=feat_M2(1:nfeat0M2,iatom_nonzero_feat)
            
            ! do jj=1,m_neigh
                ! do iat=1,natom_nn
                    ! do ii=1,nfeat0M1
                        ! dfeat(ii,iat,jj,1)=dfeat_M1(ii,iat,jj,1)
                        ! dfeat(ii,iat,jj,2)=dfeat_M1(ii,iat,jj,2)
                        ! dfeat(ii,iat,jj,3)=dfeat_M1(ii,iat,jj,3)
                    ! enddo
                ! enddo
                ! do iat=1,natom_nn
                    ! do ii=1,nfeat0M2
                        ! dfeat(ii+nfeat0M1,iat,jj,1)=dfeat_M2(ii,iat,jj,1)
                        ! dfeat(ii+nfeat0M1,iat,jj,2)=dfeat_M2(ii,iat,jj,2)
                        ! dfeat(ii+nfeat0M1,iat,jj,3)=dfeat_M2(ii,iat,jj,3)
                    ! enddo
                ! enddo
            ! enddo
            
            dfeat(1:nfeat0M1,iatom_nonzero_feat,1:m_neigh,1:3)=dfeat_M1(1:nfeat0M1,iatom_nonzero_feat,1:m_neigh,1:3)
            dfeat(1+nfeat0M1:nfeat0M2+nfeat0M1,iatom_nonzero_feat,1:m_neigh,1:3)=dfeat_M2(1:nfeat0M2,iatom_nonzero_feat,1:m_neigh,1:3)
            
            if(iflag_timer) t_bond_inference_start=mpi_wtime()
            ! if(iflag_model_ppp.eq.1) then
                call cal_energy_force_lin_ppp_bond_part(feat,dfeat,num_neigh_alltype,  &
                    list_neigh_alltype,AL,xatom,natom_nn,nfeat0,m_neigh, &
                    ippp,ibond,itype_mol)
                if(((ippp-1)/2+1).eq.1) then
                    ppp1(mod((ippp-1),2)+1,ibond)=Etot_pred_lin*sign_ppp_model_bond(ippp,ibond,itype_mol)
                    do iat=1,natom_nn
                        dppp1(:,iat,mod((ippp-1),2)+1,ibond)=force_pred_lin(:,iat)*sign_ppp_model_bond(ippp,ibond,itype_mol)
                    enddo
                else
                    ppp2(mod((ippp-1),2)+1,ibond)=Etot_pred_lin*sign_ppp_model_bond(ippp,ibond,itype_mol)
                    do iat=1,natom_nn
                        dppp2(:,iat,mod((ippp-1),2)+1,ibond)=force_pred_lin(:,iat)*sign_ppp_model_bond(ippp,ibond,itype_mol)
                    enddo
                endif
            ! endif
            if(iflag_timer) t_bond_inference_end=mpi_wtime()
            if(iflag_timer) dt_bond_inference=dt_bond_inference+t_bond_inference_end-t_bond_inference_start
            
            deallocate(feat)
            deallocate(dfeat)
        2000 continue
    2001 continue
    if(iflag_timer) t_bond_end=mpi_wtime()
    if(iflag_timer) dt_bond=t_bond_end-t_bond_start
    
    dppp=dppp*A_AU_1
    dppp1=dppp1*A_AU_1
    dppp2=dppp2*A_AU_1
    
    if(iflag_timer) t_ppp_end=mpi_wtime()
    if(iflag_timer) dt_ppp=t_ppp_end-t_ppp_start
    
    if(iflag_timer.and.(inode_tot.eq.1)) then
        open(400,file='PROFILE.ML_FF_PPP',position='append')
        200 FORMAT('', A25, F15.10, F15.10)
        300 FORMAT('', A25, F15.10)
        write(400,*)
        write(400,*) '------------------ ppp model timer -------------------'
        write(400,300) 'ppp model',dt_ppp
        write(400,*) '------------------------------------------------------'
        write(400,200) 'atom part',dt_atom,dt_atom/dt_ppp
        write(400,*) '                      ----------                      '
        write(400,200) 'atom part, read file',dt_atom_read_file,dt_atom_read_file/dt_atom
        write(400,200) 'atom part, load model',dt_atom_load_model,dt_atom_load_model/dt_atom
        write(400,200) 'atom part, set info',dt_atom_set_info,dt_atom_set_info/dt_atom
        write(400,200) 'atom part, gen feature',dt_atom_gen_feature,dt_atom_gen_feature/dt_atom
        write(400,200) 'atom part, inference',dt_atom_inference,dt_atom_inference/dt_atom
        write(400,*) '------------------------------------------------------'
        write(400,200) 'bond part',dt_bond,dt_bond/dt_ppp
        write(400,*) '                      ----------                      '
        write(400,200) 'bond part, read file',dt_bond_read_file,dt_bond_read_file/dt_bond
        write(400,200) 'bond part, load model',dt_bond_load_model,dt_bond_load_model/dt_bond
        write(400,200) 'bond part, set info',dt_bond_set_info,dt_bond_set_info/dt_bond
        write(400,200) 'bond part, gen feature',dt_bond_gen_feature,dt_bond_gen_feature/dt_bond
        write(400,200) 'bond part, inference',dt_bond_inference,dt_bond_inference/dt_bond
        write(400,*) '------------------------------------------------------'
        write(400,*)
        close(400)
    endif
    
end subroutine ML_FF_PPP
