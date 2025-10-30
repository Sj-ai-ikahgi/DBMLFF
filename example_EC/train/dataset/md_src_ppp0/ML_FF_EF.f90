subroutine ML_FF_EF(Etot_allM,fatom,xatom,AL,natom_tmp)

        use mod_mpi
        use mod_data, only :iflag_model,iatom_m,natom_m,nmol,ntype_m,ntype_ML
        use calc_ftype1, only : feat_M1,dfeat_M1,nfeat0M1,gen_feature_type1,  &
            nfeat0M1,num_neigh_alltypeM1,list_neigh_alltypeM1,  &
            m_neigh,load_model_type1,set_image_info_type1,path_2b
        use calc_ftype2, only : feat_M2,dfeat_M2,nfeat0M2,gen_feature_type2,  &
            nfeat0M2,num_neigh_alltypeM2,list_neigh_alltypeM2,  &
            load_model_type2,set_image_info_type2,path_3b
    
        !  Note: num_neigh)alltypeM1,2; list_neigh_altypeM1,2 should be the same for 1 & 2
        !  It is important to use only, there are many variable, which has the
        !  same name, but different meaning in different subrtouine, a
        !  consequence of merging different codes
        use calc_lin, only : cal_energy_force_lin,Etot_pred_lin,force_pred_lin,&
            natom_nn,load_model_lin,set_image_info_lin
        use calc_VV, only : cal_energy_force_VV,Etot_pred_VV,force_pred_VV
        use calc_NN, only : cal_energy_force_NN,Etot_pred_NN,force_pred_NN
        implicit none
    
        integer natom_tmp  ! conflict with the one in calc_lin
        real*8 Etot
        real*8 fatom(3,natom_tmp)
        real*8 fatom_t(3,natom_tmp)
        real*8 xatom(3,natom_tmp)
        real*8 AL(3,3)
        real*8,allocatable,dimension (:,:) :: feat
        real*8,allocatable,dimension (:,:,:,:) :: dfeat
        integer nfeat0
        integer ii,jj,iat
        integer kkk1,kkk2,kk,imol,ierr
        real*8 Etot_allM,Etot_allM_t
        real*8 tt1,tt2,tt3,tt4,tt5
        real*8,allocatable,dimension(:,:) :: fatom_mt,xatom_mt
        integer iflag_mol(1000)
        integer kkk_mol,itype_mol
        logical is_reset
          
        real*8 E_eachM(1000), E_eachM_t(1000)

        character(len=1000) path
        character(len=20) str

        E_eachM = 0.0d0
        E_eachM_t = 0.0d0

        Etot_allM=0.d0
        fatom=0.d0
    
        kkk1=0
        kkk2=0
        kkk_mol=0
        iflag_mol(1:ntype_m)=0
        do 6001 itype_mol=1,ntype_m
    
            natom_nn=natom_m(itype_mol)  ! in mod_linear
            natom_n=natom_m(itype_mol)   ! in mod_mpi
    
            if(natom_n.eq.1) then       ! There was a original bug, fixed:Aug.2
                do imol=1,nmol(itype_mol)
                    kkk1=kkk1+natom_nn
                    kkk2=kkk2+natom_nn
                enddo
                goto 6001   ! special, skip this for one atom case
            endif
    
    
            allocate(fatom_mt(3,natom_nn))
            allocate(xatom_mt(3,natom_nn))
            
    
    
            do 6000 imol=1,nmol(itype_mol)
                kkk_mol=kkk_mol+1
    
                if(mod(kkk_mol-1,nnodes_tot).ne.inode_tot-1) then
                    kkk1=kkk1+natom_nn
                    kkk2=kkk2+natom_nn
                    goto 6000
                endif
    
                !cccccccccccccccccccccccccc
                write(str,*) itype_mol
                str=adjustl(str)
                path = 'mol.'//trim(str)
                if(iflag_mol(itype_mol).eq.0) then
                    iflag_mol(itype_mol)=1    ! read in the data
                    call load_model_lin(ntype_ML(itype_mol),path)       ! only works for lin
                    call load_model_type1(path)
                    call load_model_type2(path)
                    is_reset= .true.
                    call set_image_info_lin(iatom_m(1,itype_mol),is_reset,natom_m(itype_mol))
                    call set_image_info_type1(iatom_m(1,itype_mol),is_reset,natom_m(itype_mol))
                    call set_image_info_type2(iatom_m(1,itype_mol),is_reset,natom_m(itype_mol))
                endif
                !cccccccccccccccccccccccccc
    
                do ii=1,natom_nn
                    kkk1=kkk1+1
                    xatom_mt(:,ii)=xatom(:,kkk1)
                enddo
    
                !      all other information has been uploaded from load_mode_lin
                write(str,*) itype_mol
                str=adjustl(str)
                path_2b = 'mol.'//trim(str)
                path_3b = path_2b

                tt1=mpi_wtime()
                call gen_feature_type1(AL,xatom_mt)
                tt2=mpi_wtime()
                call gen_feature_type2(AL,xatom_mt)
                tt3=mpi_wtime()
    
                nfeat0=nfeat0M1+nfeat0M2
    
                !cccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !  Assemble different feature types
                allocate(feat(nfeat0,natom_nn))
                allocate(dfeat(nfeat0,natom_nn,m_neigh,3))
            
                do iat=1,natom_nn
                    do ii=1,nfeat0M1
                        feat(ii,iat)=feat_M1(ii,iat)
                    enddo
                    do ii=1,nfeat0M2
                        feat(ii+nfeat0M1,iat)=feat_M2(ii,iat)
                    enddo
                enddo
    
    
                do jj=1,m_neigh
                    do iat=1,natom_nn
                        do ii=1,nfeat0M1
                            dfeat(ii,iat,jj,1)=dfeat_M1(ii,iat,jj,1)
                            dfeat(ii,iat,jj,2)=dfeat_M1(ii,iat,jj,2)
                            dfeat(ii,iat,jj,3)=dfeat_M1(ii,iat,jj,3)
                        enddo
                    enddo
                    do iat=1,natom_nn
                        do ii=1,nfeat0M2
                            dfeat(ii+nfeat0M1,iat,jj,1)=dfeat_M2(ii,iat,jj,1)
                            dfeat(ii+nfeat0M1,iat,jj,2)=dfeat_M2(ii,iat,jj,2)
                            dfeat(ii+nfeat0M1,iat,jj,3)=dfeat_M2(ii,iat,jj,3)
                        enddo
                    enddo
                enddo
                !cccccccccccccccccccccccccccccccccccccccccccccccccccccc
                
                if(iflag_model.eq.1) then
                    call cal_energy_force_lin(feat,dfeat,num_neigh_alltypeM1,  &
                        list_neigh_alltypeM1,AL,xatom_mt,natom_nn,nfeat0,m_neigh)
                    Etot=Etot_pred_lin    ! unit?
                    fatom_mt(:,1:natom_nn)=force_pred_lin(:,1:natom_nn)   ! unit, and - sign?
                endif

                !        if(iflag_model.eq.2) then
                !        call cal_energy_force_VV(feat,dfeat,num_neigh_alltypeM1,  &
                !         list_neigh_alltypeM1,AL,xatom_mt,natom_n,nfeat0,m_neigh)
                !        Etot=Etot_pred_VV    ! unit?
                !        fatom_mt(:,1:natom_m)=force_pred_VV(:,1:natom_m)   ! unit, and - sign?
                !        endif
    
                !        if(iflag_model.eq.3) then
                !        call cal_energy_force_NN(feat,dfeat,num_neigh_alltypeM1,  &
                !         list_neigh_alltypeM1,AL,xatom_mt,natom_n,nfeat0,m_neigh)
                !        Etot=Etot_pred_NN    ! unit?
                !        fatom_mt(:,1:natom_m)=force_pred_NN(:,1:natom_m)   ! unit, and - sign?
                !        endif
    
                tt4=mpi_wtime()
    
                !        if(inode.eq.1) then
                !        write(6,*) "time gen_f1,2       ",tt2-tt1,tt3-tt2
                !        write(6,*) "time f1,2 to E_force", tt4-tt3
                !        write(6,*) "time tot_E_force    ", tt4-tt1
                !        endif
     
                Etot_allM=Etot_allM+Etot

                E_eachM(kkk_mol) = Etot

    
                do ii=1,natom_nn
                    kkk2=kkk2+1
                    fatom(:,kkk2)=fatom_mt(:,ii)
                enddo
    
                deallocate(feat)
                deallocate(dfeat)
    
    
    6000    continue
            deallocate(fatom_mt)
            deallocate(xatom_mt)
    6001 continue
    
         call mpi_allreduce(Etot_allM,Etot_allM_t,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
         Etot_allM=Etot_allM_t
         call mpi_allreduce(fatom,fatom_t,3*natom_tmp,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
         fatom=fatom_t

         call mpi_allreduce(E_eachM,E_eachM_t,1000,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
         E_eachM = E_eachM_t
         
        !  if(inode_tot.eq.1) then
        !     open(36, file='Binding_Energy_MLFF.txt', position='append')
        !     write(36, "(4(F16.6, 1X))") Etot_allM, E_eachM(1), E_eachM(2), E_eachM(3)
        !     close(36)
        !  endif 

         
    
         return
     end subroutine ML_FF_EF
            

    
