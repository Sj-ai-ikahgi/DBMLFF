module mod_ppp
    
    use mod_setting, only: max_num_feature_type
    use mod_data, only: ntype_mm,natom_m,natom_mm
    use mod_param_densityFF, only: numb,max_numb
    use calc_ftype1_ppp, only: Rc_M_1=>Rc_M,m_neigh_1=>m_neigh,ntype_1=>ntype, &
        iat_type_1=>iat_type,Rc_type_1=>Rc_type,n2b_type_1=>n2b_type,iflag_ftype_1=>iflag_ftype, &
        n2bm_1=>n2bm,nfeat0m_1=>nfeat0m,grid2_2_1=>grid2_2
    use calc_ftype2_ppp, only: Rc_M_2=>Rc_M,m_neigh_2=>m_neigh,ntype_2=>ntype, &
        iat_type_2=>iat_type,Rc_type_2=>Rc_type,Rc2_type_2=>Rc2_type,n3b1_type_2=>n3b1_type, &
        n3b2_type_2=>n3b2_type,iflag_ftype_2=>iflag_ftype,n3b1m_2=>n3b1m,n3b2m_2=>n3b2m, &
        nfeat0m_2=>nfeat0m,grid31_2_2=>grid31_2,grid32_2_2=>grid32_2
    use find_neighbore_ppp, only: natom_find_neigh=>natom,ntype_find_neigh=>ntype, &
        m_neigh_find_neigh=>m_neigh,num_neigh_find_neigh=>num_neigh, &
        num_neigh_alltype_find_neigh=>num_neigh_alltype,iat_neigh_find_neigh=>iat_neigh, &
        iat_neigh_alltype_find_neigh=>iat_neigh_alltype,list_neigh_find_neigh=>list_neigh, &
        list_neigh_alltype_find_neigh=>list_neigh_alltype,ind_all_neigh_find_neigh=>ind_all_neigh
    
    implicit none 
    
    interface
        double precision function ddot(n,dx,incx,dy,incy)
            integer :: n,incx,incy
            double precision,dimension(*) :: dx,dy
        end function ddot
    end interface
    
    ! read at main_MD.f90 from polar_param.input (molecule-wise)
    integer,allocatable,dimension(:) :: num_ppp_model_atom
    integer,allocatable,dimension(:) :: num_ppp_model_bond
    integer,allocatable,dimension(:,:,:) :: ind_ppp_model_atom
    integer,allocatable,dimension(:,:,:) :: ind_ppp_model_bond
    integer,allocatable,dimension(:,:,:) :: sign_ppp_model_atom
    integer,allocatable,dimension(:,:,:) :: sign_ppp_model_bond
    
    ! read at main_MD.f90 from MD.input
    integer iflag_save_ppp
    integer iflag_save_ppp_m
    
    ! read at main_MD.f90 (using subroutine load_ppp_model) from fit_linearMM.input (model-wise)
    logical,allocatable,dimension(:,:) :: is_fitted_atom
    logical,allocatable,dimension(:,:) :: is_fitted_bond
    integer,allocatable,dimension(:,:) :: ntype_atom
    integer,allocatable,dimension(:,:) :: ntype_bond
    integer,allocatable,dimension(:,:) :: m_neigh_atom
    integer,allocatable,dimension(:,:) :: m_neigh_bond
    integer,allocatable,dimension(:,:) :: natom_mn_atom
    integer,allocatable,dimension(:,:) :: natom_mn_bond
    integer,allocatable,dimension(:,:,:) :: itype_atom_atom
    integer,allocatable,dimension(:,:,:) :: itype_atom_bond
    
    ! read at main_MD.f90 (using subroutine load_ppp_model) from MOVEMENT.type (ppp-wise)
    integer,allocatable,dimension(:,:,:,:) :: iatom_m_atom
    integer,allocatable,dimension(:,:,:,:) :: iatom_m_bond
    
    ! read at main_MD.f90 (using subroutine load_ppp_model) from feat.info (model-wise)
    integer,allocatable,dimension(:,:,:) :: ifeat_type_atom
    integer,allocatable,dimension(:,:,:) :: ifeat_type_bond
    integer,allocatable,dimension(:,:,:) :: nfeat1_atom
    integer,allocatable,dimension(:,:,:) :: nfeat1_bond
    integer,allocatable,dimension(:,:,:) :: nfeat2_atom
    integer,allocatable,dimension(:,:,:) :: nfeat2_bond
    ! calculated based on the variables above
    integer,allocatable,dimension(:,:) :: nfeat1m_atom
    integer,allocatable,dimension(:,:) :: nfeat1m_bond
    integer,allocatable,dimension(:,:) :: nfeat2m_atom
    integer,allocatable,dimension(:,:) :: nfeat2m_bond
    integer,allocatable,dimension(:,:) :: nfeat2tot_atom
    integer,allocatable,dimension(:,:) :: nfeat2tot_bond
    integer,allocatable,dimension(:,:,:) :: nfeat2i_atom
    integer,allocatable,dimension(:,:,:) :: nfeat2i_bond
    
    ! read at main_MD.f90 (using subroutine load_ppp_model) from linear_fitB.ntype (model-wise)
    real*8,allocatable,dimension(:,:,:) :: bb_atom
    real*8,allocatable,dimension(:,:,:) :: bb_bond
    ! calculated based on the variables above
    real*8,allocatable,dimension(:,:,:,:) :: bb_type_atom
    real*8,allocatable,dimension(:,:,:,:) :: bb_type_bond
    real*8,allocatable,dimension(:,:,:,:) :: bb_type0_atom
    real*8,allocatable,dimension(:,:,:,:) :: bb_type0_bond
    
    ! read at main_MD.f90 (using subroutine load_ppp_model) from feat_PV.x (model-wise)
    real*8,allocatable,dimension(:,:,:,:,:) :: pv_atom
    real*8,allocatable,dimension(:,:,:,:,:) :: pv_bond
    real*8,allocatable,dimension(:,:,:,:) :: feat2_shift_atom
    real*8,allocatable,dimension(:,:,:,:) :: feat2_shift_bond
    real*8,allocatable,dimension(:,:,:,:) :: feat2_scale_atom
    real*8,allocatable,dimension(:,:,:,:) :: feat2_scale_bond
    ! calculated based on the variables above
    real*8,allocatable,dimension(:,:,:,:) :: pv_scale_bb_atom
    real*8,allocatable,dimension(:,:,:,:) :: pv_scale_bb_bond
    real*8,allocatable,dimension(:,:,:) :: shift_scale_bb_atom
    real*8,allocatable,dimension(:,:,:) :: shift_scale_bb_bond
    
    ! read at main_MD.f90 (using subroutine load_ppp_model) from gen_2b_feature.in (molecule-wise)
    real*8,allocatable,dimension(:) :: Rc_M_mem
    integer,allocatable,dimension(:) :: m_neigh_mem
    integer,allocatable,dimension(:) :: ntype_mem
    integer,allocatable,dimension(:,:) :: iat_type_mem
    real*8,allocatable,dimension(:,:) :: Rc_type_1_mem
    integer,allocatable,dimension(:,:) :: n2b_type_1_mem
    ! calculated based on the variables above
    integer,allocatable,dimension(:) :: n2bm_1_mem
    integer,allocatable,dimension(:) :: nfeat0m_1_mem
    ! read at main_MD.f90 (using subroutine load_ppp_model) from grid2b_type3.x (molecule-wise)
    real*8,allocatable,dimension(:,:,:,:) :: grid2_2_1_mem
    
    ! read at main_MD.f90 (using subroutine load_ppp_model) from gen_3b_feature.in (molecule-wise)
    real*8,allocatable,dimension(:,:) :: Rc_type_2_mem
    real*8,allocatable,dimension(:,:) :: Rc2_type_2_mem
    integer,allocatable,dimension(:,:) :: n3b1_type_2_mem
    integer,allocatable,dimension(:,:) :: n3b2_type_2_mem
    ! calculated based on the variables above
    integer,allocatable,dimension(:) :: n3b1m_2_mem
    integer,allocatable,dimension(:) :: n3b2m_2_mem
    integer,allocatable,dimension(:) :: nfeat0m_2_mem
    ! read at main_MD.f90 (using subroutine load_ppp_model) from grid3b_cb12_type3.x (molecule-wise)
    real*8,allocatable,dimension(:,:,:,:) :: grid31_2_2_mem
    ! read at main_MD.f90 (using subroutine load_ppp_model) from grid3b_b1b2_type3.x (molecule-wise)
    real*8,allocatable,dimension(:,:,:,:) :: grid32_2_2_mem
    
    ! read at main_MD.f90 (using subroutine load_ppp_model) from IN.NEIGHBORE (molecule-wise)
    integer,allocatable,dimension(:,:,:) :: num_neigh_mem
    integer,allocatable,dimension(:,:,:,:) :: iat_neigh_mem
    integer,allocatable,dimension(:,:,:,:) :: list_neigh_mem
    ! calculated based on the variables above
    integer,allocatable,dimension(:,:) :: num_neigh_alltype_mem
    integer,allocatable,dimension(:,:,:) :: iat_neigh_alltype_mem
    integer,allocatable,dimension(:,:,:) :: list_neigh_alltype_mem
    integer,allocatable,dimension(:,:,:,:) :: ind_all_neigh_mem
    
    character(len=20) str
    character(len=20) str2
    character(len=20) str3
    
contains
    subroutine load_ppp_model()
        implicit none
        
        integer i,j,k
        integer itype_mol
        integer imodel
        integer iatom,iat
        integer ibond
        integer itype
        integer counter,itype1,itype2,k1,k2,k12,ii_f
        integer num
        
        character(len=100) str_bin,str_bin2,str_bin3
        integer int_bin,int_bin2
        integer io_stat
        real*8 real_bin
        
        integer max_num_ppp_model_atom
        integer max_num_ppp_model_bond
        integer max_ntype_atom
        integer max_ntype_bond
        integer max_nfeat1m_atom
        integer max_nfeat1m_bond
        integer max_nfeat2m_atom
        integer max_nfeat2m_bond
        integer max_nfeat2tot_atom
        integer max_nfeat2tot_bond
        integer max_n2bm_1_mem
        integer max_ntype_mem
        integer max_m_neigh_mem
        integer max_n3b1m_2_mem
        integer max_n3b2m_2_mem
        
        integer natom_movement_type
        integer nfeat_type
        
        real*8,allocatable,dimension(:) :: scale_bb
        
        max_num_ppp_model_atom=0
        max_num_ppp_model_bond=0
        do itype_mol=1,ntype_mm
            if(num_ppp_model_atom(itype_mol).gt.max_num_ppp_model_atom) max_num_ppp_model_atom=num_ppp_model_atom(itype_mol)
            if(num_ppp_model_bond(itype_mol).gt.max_num_ppp_model_bond) max_num_ppp_model_bond=num_ppp_model_bond(itype_mol)
        enddo
        
        ! ---------- read variable is_fitted, ntype, m_neigh, natom_mn, and itype_atom ----------
        ! ---------- calc max_ntype_atom, max_ntype_bond ----------
        allocate(is_fitted_atom(max_num_ppp_model_atom,ntype_mm))
        allocate(is_fitted_bond(max_num_ppp_model_bond,ntype_mm))
        allocate(ntype_atom(max_num_ppp_model_atom,ntype_mm))
        allocate(ntype_bond(max_num_ppp_model_bond,ntype_mm))
        allocate(m_neigh_atom(max_num_ppp_model_atom,ntype_mm))
        allocate(m_neigh_bond(max_num_ppp_model_bond,ntype_mm))
        allocate(natom_mn_atom(max_num_ppp_model_atom,ntype_mm))
        allocate(natom_mn_bond(max_num_ppp_model_bond,ntype_mm))
        allocate(itype_atom_atom(natom_mm,max_num_ppp_model_atom,ntype_mm)) ! size of the first dimension should be max_ntype_atom for optimal
        allocate(itype_atom_bond(natom_mm,max_num_ppp_model_bond,ntype_mm)) ! natom_mm is actually much larger
        is_fitted_atom=.true.
        is_fitted_bond=.true.
        ntype_atom=0
        max_ntype_atom=0
        ntype_bond=0
        max_ntype_bond=0
        m_neigh_atom=0
        m_neigh_bond=0
        natom_mn_atom=0
        natom_mn_bond=0
        itype_atom_atom=0
        itype_atom_bond=0
        do itype_mol=1,ntype_mm
            write(str,*) itype_mol
            str=adjustl(str)
            
            if (natom_m(itype_mol) .lt. 2) cycle

            do imodel=1,num_ppp_model_atom(itype_mol)
                write(str2,*) imodel
                str2=adjustl(str2)
                
                open(400,file='mol.'//trim(str)//'/'//'ppp_model'//'/'//'atom_'//trim(str2)//'/'//'fit_linearMM.input')
                rewind(400)
                read(400,*) str_bin,str_bin2,str_bin3
                if(trim(str_bin).eq.'Notfittedyet!') then
                    is_fitted_atom(imodel,itype_mol)=.false.
                else
                    read(str_bin,*,iostat=io_stat) ntype_atom(imodel,itype_mol)
                    if(io_stat.ne.0) then
                        write(6,*) 'read(str_bin,*,iostat=io_stat) io_stat.ne.0,stop',io_stat
                        stop
                    endif
                    if(ntype_atom(imodel,itype_mol).gt.max_ntype_atom) max_ntype_atom=ntype_atom(imodel,itype_mol)
                    read(str_bin2,*,iostat=io_stat) m_neigh_atom(imodel,itype_mol)
                    if(io_stat.ne.0) then
                        write(6,*) 'read(str_bin2,*,iostat=io_stat) io_stat.ne.0,stop',io_stat
                        stop
                    endif
                    read(str_bin3,*,iostat=io_stat) natom_mn_atom(imodel,itype_mol)
                    if(io_stat.ne.0) then
                        write(6,*) 'read(str_bin3,*,iostat=io_stat) io_stat.ne.0,stop',io_stat
                        stop
                    endif
                    do i=1,ntype_atom(imodel,itype_mol)
                        read(400,*) itype_atom_atom(i,imodel,itype_mol)
                    enddo
                endif
                close(400)
            enddo
            
            do imodel=1,num_ppp_model_bond(itype_mol)
                write(str2,*) imodel
                str2=adjustl(str2)
                
                open(401,file='mol.'//trim(str)//'/'//'ppp_model'//'/'//'bond_'//trim(str2)//'/'//'fit_linearMM.input')
                rewind(401)
                read(401,*) str_bin,str_bin2,str_bin3
                if(trim(str_bin).eq.'Notfittedyet!') then
                    is_fitted_bond(imodel,itype_mol)=.false.
                else
                    read(str_bin,*,iostat=io_stat) ntype_bond(imodel,itype_mol)
                    if(io_stat.ne.0) then
                        write(6,*) 'read(str_bin,*,iostat=io_stat) io_stat.ne.0,stop',io_stat
                        stop
                    endif
                    if(ntype_bond(imodel,itype_mol).gt.max_ntype_bond) max_ntype_bond=ntype_bond(imodel,itype_mol)
                    read(str_bin2,*,iostat=io_stat) m_neigh_bond(imodel,itype_mol)
                    if(io_stat.ne.0) then
                        write(6,*) 'read(str_bin2,*,iostat=io_stat) io_stat.ne.0,stop',io_stat
                        stop
                    endif
                    read(str_bin3,*,iostat=io_stat) natom_mn_bond(imodel,itype_mol)
                    if(io_stat.ne.0) then
                        write(6,*) 'read(str_bin3,*,iostat=io_stat) io_stat.ne.0,stop',io_stat
                        stop
                    endif
                    do i=1,ntype_bond(imodel,itype_mol)
                        read(401,*) itype_atom_bond(i,imodel,itype_mol)
                    enddo
                endif
                close(401)
            enddo
        enddo
        
        ! ---------- read iatom_m ----------
        allocate(iatom_m_atom(natom_mm,3,natom_mm,ntype_mm))
        allocate(iatom_m_bond(natom_mm,2,max_numb,ntype_mm))
        iatom_m_atom=0
        iatom_m_bond=0
        do itype_mol=1,ntype_mm
            write(str,*) itype_mol
            str=adjustl(str)
            
            if (natom_m(itype_mol) .lt. 2) cycle

            do iatom=1,natom_m(itype_mol)
                write(str2,*) iatom
                str2=adjustl(str2)
                
                do i=1,3
                    if((is_fitted_atom(ind_ppp_model_atom((i-1)*2+1,iatom,itype_mol),itype_mol).eq..false.) &
                        .and.(is_fitted_atom(ind_ppp_model_atom((i-1)*2+2,iatom,itype_mol),itype_mol).eq..false.)) cycle
                    
                    write(str3,*) i
                    str3=adjustl(str3)
                    
                    ! open(402,file='mol.'//trim(str)//'/'//'ppp_neighbore'//'/'//'atom_'//trim(str2)//'_'//trim(str3)//'/'//'input'//'/'//'MOVEMENT.type')
                    open(402,file='mol.'//trim(str)//'/'//'ppp_neighbore/MOVEMENT.type')
                    rewind(402)
                    read(402,*) natom_movement_type
                    if(natom_movement_type.ne.natom_m(itype_mol)) then
                        write(6,*) 'natom_movement_type.ne.natom_m(itype_mol),stop',natom_movement_type,natom_m(itype_mol)
                        stop
                    endif
                    do j=1,natom_m(itype_mol)
                        read(402,*) int_bin,iatom_m_atom(j,i,iatom,itype_mol)
                    enddo
                    close(402)
                enddo
            enddo
            
            do ibond=1,numb(itype_mol)
                write(str2,*) ibond
                str2=adjustl(str2)
                
                do i=1,2
                    if((is_fitted_bond(ind_ppp_model_bond((i-1)*2+1,ibond,itype_mol),itype_mol).eq..false.) &
                        .and.(is_fitted_bond(ind_ppp_model_bond((i-1)*2+2,ibond,itype_mol),itype_mol).eq..false.)) cycle
                    
                    write(str3,*) i
                    str3=adjustl(str3)
                    
                    ! open(403,file='mol.'//trim(str)//'/'//'ppp_neighbore'//'/'//'bond_'//trim(str2)//'_'//trim(str3)//'/'//'input'//'/'//'MOVEMENT.type')
                    open(403,file='mol.'//trim(str)//'/'//'ppp_neighbore/MOVEMENT.type')
                    rewind(403)
                    read(403,*) natom_movement_type
                    if(natom_movement_type.ne.natom_m(itype_mol)) then
                        write(6,*) 'natom_movement_type.ne.natom_m(itype_mol),stop',natom_movement_type,natom_m(itype_mol)
                        stop
                    endif
                    do j=1,natom_m(itype_mol)
                        read(403,*) int_bin,iatom_m_bond(j,i,ibond,itype_mol)
                    enddo
                    close(403)
                enddo
            enddo
        enddo
        
        ! ---------- read ifeat_type, nfeat1, nfeat2 ----------
        ! ---------- calc nfeat1m, nfeat2m, nfeat2tot, nfeat2i ----------
        ! ---------- calc max_nfeat1m_atom, max_nfeat1m_bond, max_nfeat2m_atom, max_nfeat2m_bond, max_nfeat2tot_atom, max_nfeat2tot_bond ----------
        allocate(ifeat_type_atom(max_num_feature_type,max_num_ppp_model_atom,ntype_mm))
        allocate(ifeat_type_bond(max_num_feature_type,max_num_ppp_model_bond,ntype_mm))
        allocate(nfeat1_atom(max_ntype_atom,max_num_ppp_model_atom,ntype_mm))
        allocate(nfeat1_bond(max_ntype_bond,max_num_ppp_model_bond,ntype_mm))
        allocate(nfeat2_atom(max_ntype_atom,max_num_ppp_model_atom,ntype_mm))
        allocate(nfeat2_bond(max_ntype_bond,max_num_ppp_model_bond,ntype_mm))
        allocate(nfeat1m_atom(max_num_ppp_model_atom,ntype_mm))
        allocate(nfeat1m_bond(max_num_ppp_model_bond,ntype_mm))
        allocate(nfeat2m_atom(max_num_ppp_model_atom,ntype_mm))
        allocate(nfeat2m_bond(max_num_ppp_model_bond,ntype_mm))
        allocate(nfeat2tot_atom(max_num_ppp_model_atom,ntype_mm))
        allocate(nfeat2tot_bond(max_num_ppp_model_bond,ntype_mm))
        allocate(nfeat2i_atom(max_ntype_atom,max_num_ppp_model_atom,ntype_mm))
        allocate(nfeat2i_bond(max_ntype_bond,max_num_ppp_model_bond,ntype_mm))
        ifeat_type_atom=0
        ifeat_type_bond=0
        nfeat1_atom=0
        nfeat1_bond=0
        nfeat2_atom=0
        nfeat2_bond=0
        nfeat1m_atom=0
        max_nfeat1m_atom=0
        nfeat1m_bond=0
        max_nfeat1m_bond=0
        nfeat2m_atom=0
        max_nfeat2m_atom=0
        nfeat2m_bond=0
        max_nfeat2m_bond=0
        nfeat2tot_atom=0
        max_nfeat2tot_atom=0
        nfeat2tot_bond=0
        max_nfeat2tot_bond=0
        nfeat2i_atom=0
        nfeat2i_bond=0
        do itype_mol=1,ntype_mm
            write(str,*) itype_mol
            str=adjustl(str)
            
            if (natom_m(itype_mol) .lt. 2) cycle

            do imodel=1,num_ppp_model_atom(itype_mol)
                write(str2,*) imodel
                str2=adjustl(str2)
                
                if(is_fitted_atom(imodel,itype_mol).eq..false.) cycle
                
                open(404,file='mol.'//trim(str)//'/'//'ppp_model'//'/'//'atom_'//trim(str2)//'/'//'feat.info')
                rewind(404)
                read(404,*)
                read(404,*) nfeat_type
                do i=1,nfeat_type
                    read(404,*) ifeat_type_atom(i,imodel,itype_mol)
                enddo
                read(404,*) int_bin
                if(int_bin.ne.ntype_atom(imodel,itype_mol)) then
                    write(6,*) 'int_bin.ne.ntype_atom(imodel,itype_mol)',int_bin,ntype_atom(imodel,itype_mol)
                    stop
                endif
                do i=1,ntype_atom(imodel,itype_mol)
                    read(404,*) int_bin,nfeat1_atom(i,imodel,itype_mol),nfeat2_atom(i,imodel,itype_mol)
                    if(int_bin.ne.itype_atom_atom(i,imodel,itype_mol)) then
                        write(6,*) 'int_bin.ne.itype_atom_atom(i,imodel,itype_mol),stop',int_bin,itype_atom_atom(i,imodel,itype_mol)
                        stop
                    endif
                    if(nfeat1_atom(i,imodel,itype_mol).gt.nfeat1m_atom(imodel,itype_mol)) nfeat1m_atom(imodel,itype_mol)=nfeat1_atom(i,imodel,itype_mol)
                    if(nfeat2_atom(i,imodel,itype_mol).gt.nfeat2m_atom(imodel,itype_mol)) nfeat2m_atom(imodel,itype_mol)=nfeat2_atom(i,imodel,itype_mol)
                    nfeat2tot_atom(imodel,itype_mol)=nfeat2tot_atom(imodel,itype_mol)+nfeat2_atom(i,imodel,itype_mol)
                    if(i.gt.1) then
                        nfeat2i_atom(i,imodel,itype_mol)=nfeat2i_atom(i-1,imodel,itype_mol)+nfeat2_atom(i-1,imodel,itype_mol)
                    endif
                enddo
                if(nfeat1m_atom(imodel,itype_mol).gt.max_nfeat1m_atom) max_nfeat1m_atom=nfeat1m_atom(imodel,itype_mol)
                if(nfeat2m_atom(imodel,itype_mol).gt.max_nfeat2m_atom) max_nfeat2m_atom=nfeat2m_atom(imodel,itype_mol)
                if(nfeat2tot_atom(imodel,itype_mol).gt.max_nfeat2tot_atom) max_nfeat2tot_atom=nfeat2tot_atom(imodel,itype_mol)
                close(404)
            enddo
            
            do imodel=1,num_ppp_model_bond(itype_mol)
                write(str2,*) imodel
                str2=adjustl(str2)
                
                if(is_fitted_bond(imodel,itype_mol).eq..false.) cycle
                
                open(405,file='mol.'//trim(str)//'/'//'ppp_model'//'/'//'bond_'//trim(str2)//'/'//'feat.info')
                rewind(405)
                read(405,*)
                read(405,*) nfeat_type
                do i=1,nfeat_type
                    read(405,*) ifeat_type_bond(i,imodel,itype_mol)
                enddo
                read(405,*) int_bin
                if(int_bin.ne.ntype_bond(imodel,itype_mol)) then
                    write(6,*) 'int_bin.ne.ntype_bond(imodel,itype_mol)',int_bin,ntype_bond(imodel,itype_mol)
                    stop
                endif
                do i=1,ntype_bond(imodel,itype_mol)
                    read(405,*) int_bin,nfeat1_bond(i,imodel,itype_mol),nfeat2_bond(i,imodel,itype_mol)
                    if(int_bin.ne.itype_atom_atom(i,imodel,itype_mol)) then
                        write(6,*) 'int_bin.ne.itype_atom_atom(i,imodel,itype_mol),stop',int_bin,itype_atom_atom(i,imodel,itype_mol)
                        stop
                    endif
                    if(nfeat1_bond(i,imodel,itype_mol).gt.nfeat1m_bond(imodel,itype_mol)) nfeat1m_bond(imodel,itype_mol)=nfeat1_bond(i,imodel,itype_mol)
                    if(nfeat2_bond(i,imodel,itype_mol).gt.nfeat2m_bond(imodel,itype_mol)) nfeat2m_bond(imodel,itype_mol)=nfeat2_bond(i,imodel,itype_mol)
                    nfeat2tot_bond(imodel,itype_mol)=nfeat2tot_bond(imodel,itype_mol)+nfeat2_bond(i,imodel,itype_mol)
                    if(i.gt.1) then
                        nfeat2i_bond(i,imodel,itype_mol)=nfeat2i_bond(i-1,imodel,itype_mol)+nfeat2_bond(i-1,imodel,itype_mol)
                    endif
                enddo
                if(nfeat1m_bond(imodel,itype_mol).gt.max_nfeat1m_bond) max_nfeat1m_bond=nfeat1m_bond(imodel,itype_mol)
                if(nfeat2m_bond(imodel,itype_mol).gt.max_nfeat2m_bond) max_nfeat2m_bond=nfeat2m_bond(imodel,itype_mol)
                if(nfeat2tot_bond(imodel,itype_mol).gt.max_nfeat2tot_bond) max_nfeat2tot_bond=nfeat2tot_bond(imodel,itype_mol)
                close(405)
            enddo
        enddo
        
        ! ---------- read bb ----------
        ! ---------- calc bb_type, bb_type0 ----------
        allocate(bb_atom(max_nfeat2tot_atom,max_num_ppp_model_atom,ntype_mm))
        allocate(bb_bond(max_nfeat2tot_bond,max_num_ppp_model_bond,ntype_mm))
        allocate(bb_type_atom(max_nfeat2m_atom,max_ntype_atom,max_num_ppp_model_atom,ntype_mm))
        allocate(bb_type_bond(max_nfeat2m_bond,max_ntype_bond,max_num_ppp_model_bond,ntype_mm))
        allocate(bb_type0_atom(max_nfeat2m_atom,max_ntype_atom,max_num_ppp_model_atom,ntype_mm))
        allocate(bb_type0_bond(max_nfeat2m_bond,max_ntype_bond,max_num_ppp_model_bond,ntype_mm))
        bb_atom=0.d0
        bb_bond=0.d0
        bb_type_atom=0.d0
        bb_type_bond=0.d0
        bb_type0_atom=0.d0
        bb_type0_bond=0.d0
        do itype_mol=1,ntype_mm
            write(str,*) itype_mol
            str=adjustl(str)
            
            if (natom_m(itype_mol) .lt. 2) cycle

            do imodel=1,num_ppp_model_atom(itype_mol)
                write(str2,*) imodel
                str2=adjustl(str2)
                
                if(is_fitted_atom(imodel,itype_mol).eq..false.) cycle
                
                open(406,file='mol.'//trim(str)//'/'//'ppp_model'//'/'//'atom_'//trim(str2)//'/'//'linear_fitB.ntype')
                rewind(406)
                read(406,*) int_bin
                if(int_bin.ne.nfeat2tot_atom(imodel,itype_mol)) then
                    write(6,*) 'int_bin.ne.nfeat2tot_atom(imodel,itype_mol),stop',int_bin,nfeat2tot_atom(imodel,itype_mol)
                    stop
                endif
                do i=1,nfeat2tot_atom(imodel,itype_mol)
                    read(406,*) int_bin,bb_atom(i,imodel,itype_mol)
                enddo
                close(406)
                do i=1,ntype_atom(imodel,itype_mol)
                    do j=1,nfeat2_atom(i,imodel,itype_mol)
                        bb_type0_atom(j,i,imodel,itype_mol)=bb_atom(j+nfeat2i_atom(i,imodel,itype_mol),imodel,itype_mol)
                    enddo
                enddo
                close(406)
            enddo
            
            do imodel=1,num_ppp_model_bond(itype_mol)
                write(str2,*) imodel
                str2=adjustl(str2)
                
                if(is_fitted_bond(imodel,itype_mol).eq..false.) cycle
                
                open(407,file='mol.'//trim(str)//'/'//'ppp_model'//'/'//'bond_'//trim(str2)//'/'//'linear_fitB.ntype')
                rewind(407)
                read(407,*) int_bin
                if(int_bin.ne.nfeat2tot_bond(imodel,itype_mol)) then
                    write(6,*) 'int_bin.ne.nfeat2tot_bond(imodel,itype_mol),stop',int_bin,nfeat2tot_bond(imodel,itype_mol)
                    stop
                endif
                do i=1,nfeat2tot_bond(imodel,itype_mol)
                    read(407,*) int_bin,bb_bond(i,imodel,itype_mol)
                enddo
                close(407)
                do i=1,ntype_bond(imodel,itype_mol)
                    do j=1,nfeat2_bond(i,imodel,itype_mol)
                        bb_type0_bond(j,i,imodel,itype_mol)=bb_bond(j+nfeat2i_bond(i,imodel,itype_mol),imodel,itype_mol)
                    enddo
                enddo
                close(407)
            enddo
        enddo
        
        ! ---------- read pv, feat2_shift, feat2_scale ----------
        ! ---------- calc pv_scale_bb_atom, pv_scale_bb_bond, shift_scale_bb_atom, shift_scale_bb_bond ----------
        allocate(pv_atom(max_nfeat1m_atom,max_nfeat2m_atom,max_ntype_atom,max_num_ppp_model_atom,ntype_mm))
        allocate(pv_bond(max_nfeat1m_bond,max_nfeat2m_bond,max_ntype_bond,max_num_ppp_model_bond,ntype_mm))
        allocate(feat2_shift_atom(max_nfeat2m_atom,max_ntype_atom,max_num_ppp_model_atom,ntype_mm))
        allocate(feat2_shift_bond(max_nfeat2m_bond,max_ntype_bond,max_num_ppp_model_bond,ntype_mm))
        allocate(feat2_scale_atom(max_nfeat2m_atom,max_ntype_atom,max_num_ppp_model_atom,ntype_mm))
        allocate(feat2_scale_bond(max_nfeat2m_bond,max_ntype_bond,max_num_ppp_model_bond,ntype_mm))
        allocate(pv_scale_bb_atom(max_nfeat1m_atom,max_ntype_atom,max_num_ppp_model_atom,ntype_mm))
        allocate(pv_scale_bb_bond(max_nfeat1m_atom,max_ntype_bond,max_num_ppp_model_bond,ntype_mm))
        allocate(shift_scale_bb_atom(max_ntype_atom,max_num_ppp_model_atom,ntype_mm))
        allocate(shift_scale_bb_bond(max_ntype_bond,max_num_ppp_model_bond,ntype_mm))
        pv_atom=0.d0
        pv_bond=0.d0
        feat2_shift_atom=0.d0
        feat2_shift_bond=0.d0
        feat2_scale_atom=0.d0
        feat2_scale_bond=0.d0
        pv_scale_bb_atom=0.d0
        pv_scale_bb_bond=0.d0
        shift_scale_bb_atom=0.d0
        shift_scale_bb_bond=0.d0
        do itype_mol=1,ntype_mm
            write(str,*) itype_mol
            str=adjustl(str)
            
            if (natom_m(itype_mol) .lt. 2) cycle

            do imodel=1,num_ppp_model_atom(itype_mol)
                write(str2,*) imodel
                str2=adjustl(str2)
                
                if(is_fitted_atom(imodel,itype_mol).eq..false.) cycle
                
                do itype=1,ntype_atom(imodel,itype_mol)
                    write(str3,*) itype
                    str3=adjustl(str3)
                    
                    open(408,file='mol.'//trim(str)//'/'//'ppp_model'//'/'//'atom_'//trim(str2)//'/'//'feat_PV.'//trim(str3),form='unformatted')
                    rewind(408)
                    read(408) int_bin,int_bin2
                    if(int_bin2.ne.nfeat2_atom(itype,imodel,itype_mol)) then
                        write(6,*) 'int_bin2.ne.nfeat2_atom(itype,imodel,itype_mol),stop',int_bin2,nfeat2_atom(itype,imodel,itype_mol)
                        stop
                    endif
                    if(int_bin.ne.nfeat1_atom(itype,imodel,itype_mol)) then
                        write(6,*) 'int_bin.ne.nfeat1_atom(itype,imodel,itype_mol),stop',int_bin,nfeat1_atom(itype,imodel,itype_mol)
                        stop
                    endif
                    read(408) pv_atom(1:nfeat1_atom(itype,imodel,itype_mol),1:nfeat2_atom(itype,imodel,itype_mol),itype,imodel,itype_mol)
                    read(408) feat2_shift_atom(1:nfeat2_atom(itype,imodel,itype_mol),itype,imodel,itype_mol)
                    read(408) feat2_scale_atom(1:nfeat2_atom(itype,imodel,itype_mol),itype,imodel,itype_mol)
                    close(408)
                    
                    if(nfeat2_atom(itype,imodel,itype_mol).ne.0) then
                        allocate(scale_bb(nfeat2_atom(itype,imodel,itype_mol)-1)) ! nfeat2m_atom(imodel,itype_mol)-1 because the last feat2 will be set to 1.d0
                        scale_bb=0.d0
                        scale_bb(:)=feat2_scale_atom(1:nfeat2_atom(itype,imodel,itype_mol)-1,itype,imodel,itype_mol)*bb_type0_atom(1:nfeat2_atom(itype,imodel,itype_mol)-1,itype,imodel,itype_mol)
                        call dgemv('N',nfeat1_atom(itype,imodel,itype_mol),nfeat2_atom(itype,imodel,itype_mol)-1,1.d0, &
                            pv_atom(1,1,itype,imodel,itype_mol),max_nfeat1m_atom,scale_bb,1,0.d0, &
                            pv_scale_bb_atom(1,itype,imodel,itype_mol),1)
                        shift_scale_bb_atom(itype,imodel,itype_mol)=ddot(nfeat2_atom(itype,imodel,itype_mol)-1,feat2_shift_atom(1,itype,imodel,itype_mol),1,scale_bb,1)
                        deallocate(scale_bb)
                    endif
                enddo
            enddo
            
            do imodel=1,num_ppp_model_bond(itype_mol)
                write(str2,*) imodel
                str2=adjustl(str2)
                
                if(is_fitted_bond(imodel,itype_mol).eq..false.) cycle
                
                do itype=1,ntype_bond(imodel,itype_mol)
                    write(str3,*) itype
                    str3=adjustl(str3)
                    
                    open(409,file='mol.'//trim(str)//'/'//'ppp_model'//'/'//'bond_'//trim(str2)//'/'//'feat_PV.'//trim(str3),form='unformatted')
                    rewind(409)
                    read(409) int_bin,int_bin2
                    if(int_bin2.ne.nfeat2_bond(itype,imodel,itype_mol)) then
                        write(6,*) 'int_bin2.ne.nfeat2_bond(itype,imodel,itype_mol),stop',int_bin2,nfeat2_bond(itype,imodel,itype_mol)
                        stop
                    endif
                    if(int_bin.ne.nfeat1_bond(itype,imodel,itype_mol)) then
                        write(6,*) 'int_bin.ne.nfeat1_bond(itype,imodel,itype_mol),stop',int_bin,nfeat1_bond(itype,imodel,itype_mol)
                        stop
                    endif
                    read(409) pv_bond(1:nfeat1_bond(itype,imodel,itype_mol),1:nfeat2_bond(itype,imodel,itype_mol),itype,imodel,itype_mol)
                    read(409) feat2_shift_bond(1:nfeat2_bond(itype,imodel,itype_mol),itype,imodel,itype_mol)
                    read(409) feat2_scale_bond(1:nfeat2_bond(itype,imodel,itype_mol),itype,imodel,itype_mol)
                    close(409)
                    
                    if(nfeat2_bond(itype,imodel,itype_mol).ne.0) then
                        allocate(scale_bb(nfeat2_bond(itype,imodel,itype_mol)-1)) ! nfeat2_bond(itype,imodel,itype_mol)-1 because the last feat2 will be set to 1.d0
                        scale_bb=0.d0
                        scale_bb(:)=feat2_scale_bond(1:nfeat2_bond(itype,imodel,itype_mol)-1,itype,imodel,itype_mol)*bb_type0_bond(1:nfeat2_bond(itype,imodel,itype_mol)-1,itype,imodel,itype_mol)
                        call dgemv('N',nfeat1_bond(itype,imodel,itype_mol),nfeat2_bond(itype,imodel,itype_mol)-1,1.d0, &
                            pv_bond(1,1,itype,imodel,itype_mol),max_nfeat1m_bond,scale_bb,1,0.d0, &
                            pv_scale_bb_bond(1,itype,imodel,itype_mol),1)
                        shift_scale_bb_bond(itype,imodel,itype_mol)=ddot(nfeat2_bond(itype,imodel,itype_mol)-1,feat2_shift_bond(1,itype,imodel,itype_mol),1,scale_bb,1)
                        deallocate(scale_bb)
                    endif
                enddo
            enddo
        enddo
        
        ! ---------- read Rc_M_mem, m_neigh_mem, ntype_mem, iat_type_mem, Rc_type_1_mem, n2b_type_1_mem ----------
        ! ---------- calc n2bm_1_mem, nfeat0m_1_mem ----------
        ! ---------- calc max_n2bm_1_mem, max_ntype_mem, max_m_neigh_mem ----------
        allocate(Rc_M_mem(ntype_mm))
        allocate(m_neigh_mem(ntype_mm))
        allocate(ntype_mem(ntype_mm))
        allocate(iat_type_mem(natom_mm,ntype_mm))
        allocate(Rc_type_1_mem(natom_mm,ntype_mm))
        allocate(n2b_type_1_mem(natom_mm,ntype_mm))
        allocate(n2bm_1_mem(ntype_mm))
        allocate(nfeat0m_1_mem(ntype_mm))
        Rc_M_mem=0.d0
        m_neigh_mem=0
        ntype_mem=0
        iat_type_mem=0
        Rc_type_1_mem=0.d0
        n2b_type_1_mem=0
        n2bm_1_mem=0
        nfeat0m_1_mem=0
        max_n2bm_1_mem=0
        max_ntype_mem=0
        do itype_mol=1,ntype_mm
            write(str,*) itype_mol
            str=adjustl(str)
            
            if (natom_m(itype_mol) .lt. 2) cycle

            open(410,file='mol.'//trim(str)//'/'//'ppp_model'//'/'//'atom_1'//'/'//'input'//'/'//'gen_2b_feature.in')
            rewind(410)
            read(410,*) Rc_M_mem(itype_mol),m_neigh_mem(itype_mol)
            read(410,*) ntype_mem(itype_mol)
            do i=1,ntype_mem(itype_mol)
                read(410,*) iat_type_mem(i,itype_mol)
                read(410,*) Rc_type_1_mem(i,itype_mol)
                read(410,*) n2b_type_1_mem(i,itype_mol)
                
                if(Rc_type_1_mem(i,itype_mol).gt.Rc_M_mem(itype_mol)) then
                    write(6,*) 'Rc_type_1_mem(i,itype_mol).gt.Rc_M_mem(itype_mol),stop',Rc_type_1_mem(i,itype_mol),Rc_M_mem(itype_mol)
                    stop
                endif
            enddo
            close(410)
            if(m_neigh_mem(itype_mol).gt.max_m_neigh_mem) max_m_neigh_mem=m_neigh_mem(itype_mol)
            if(ntype_mem(itype_mol).gt.max_ntype_mem) max_ntype_mem=ntype_mem(itype_mol)
            
            do i=1,ntype_mem(itype_mol)
                if(n2b_type_1_mem(i,itype_mol).gt.n2bm_1_mem(itype_mol)) n2bm_1_mem(itype_mol)=n2b_type_1_mem(i,itype_mol)
            enddo
            if(n2bm_1_mem(itype_mol).gt.max_n2bm_1_mem) max_n2bm_1_mem=n2bm_1_mem(itype_mol)
            nfeat0m_1_mem(itype_mol)=ntype_mem(itype_mol)*n2bm_1_mem(itype_mol)
        enddo
        
        ! ---------- read grid2_2_1_mem ----------
        allocate(grid2_2_1_mem(2,max_n2bm_1_mem,max_ntype_mem,ntype_mm))
        grid2_2_1_mem=0.d0
        do itype_mol=1,ntype_mm
            write(str,*) itype_mol
            str=adjustl(str)
            
            if (natom_m(itype_mol) .lt. 2) cycle

            do itype=1,ntype_mem(itype_mol)
                write(str2,*) itype
                str2=adjustl(str2)
                
                open(411,file='mol.'//trim(str)//'/'//'ppp_model'//'/'//'atom_1'//'/'//'output'//'/'//'grid2b_type3.'//trim(str2))
                rewind(411)
                read(411,*) int_bin
                if(int_bin.ne.n2b_type_1_mem(itype,itype_mol)) then
                    write(6,*) 'int_bin.ne.n2b_type_1_mem(itype,itype_mol),stop',int_bin,n2b_type_1_mem(itype,itype_mol)
                    stop
                endif
                do i=1,n2b_type_1_mem(itype,itype_mol)
                    read(411,*) int_bin,grid2_2_1_mem(1,i,itype,itype_mol),grid2_2_1_mem(2,i,itype,itype_mol)
                    if(grid2_2_1_mem(2,i,itype,itype_mol).gt.Rc_type_1_mem(itype,itype_mol)) then
                        write(6,*) 'grid2_2_1_mem(2,i,itype,itype_mol).gt.Rc_type_1_mem(itype,itype_mol)',grid2_2_1_mem(2,i,itype,itype_mol),Rc_type_1_mem(itype,itype_mol)
                    endif
                enddo
                close(411)
            enddo
        enddo
        
        ! ---------- read Rc_type_2_mem, Rc2_type_2_mem, n3b1_type_2_mem, n3b2_type_2_mem ----------
        ! ---------- calc n3b1m_2_mem, n3b2m_2_mem, nfeat0m_2_mem ----------
        ! ---------- calc max_n3b1m_2_mem, max_n3b2m_2_mem ----------
        allocate(Rc_type_2_mem(natom_mm,ntype_mm))
        allocate(Rc2_type_2_mem(natom_mm,ntype_mm))
        allocate(n3b1_type_2_mem(natom_mm,ntype_mm))
        allocate(n3b2_type_2_mem(natom_mm,ntype_mm))
        allocate(n3b1m_2_mem(ntype_mm))
        allocate(n3b2m_2_mem(ntype_mm))
        allocate(nfeat0m_2_mem(ntype_mm))
        Rc_type_2_mem=0.d0
        Rc2_type_2_mem=0.d0
        n3b1_type_2_mem=0
        n3b2_type_2_mem=0
        n3b1m_2_mem=0
        n3b2m_2_mem=0
        nfeat0m_2_mem=0
        max_n3b1m_2_mem=0
        max_n3b2m_2_mem=0
        do itype_mol=1,ntype_mm
            write(str,*) itype_mol
            str=adjustl(str)
            
            if (natom_m(itype_mol) .lt. 2) cycle

            open(412,file='mol.'//trim(str)//'/'//'ppp_model'//'/'//'atom_1'//'/'//'input'//'/'//'gen_3b_feature.in')
            rewind(412)
            read(412,*) real_bin,int_bin
            read(412,*) int_bin
            do i=1,ntype_mem(itype_mol)
                read(412,*) int_bin
                read(412,*) Rc_type_2_mem(i,itype_mol),Rc2_type_2_mem(i,itype_mol)
                read(412,*) n3b1_type_2_mem(i,itype_mol),n3b2_type_2_mem(i,itype_mol)
                
                if(Rc_type_2_mem(i,itype_mol).gt.Rc_M_mem(itype_mol)) then
                    write(6,*) 'Rc_type_2_mem(i,itype_mol).gt.Rc_M_mem(itype_mol),stop',Rc_type_2_mem(i,itype_mol),Rc_M_mem(itype_mol)
                    stop
                endif
                if(Rc2_type_2_mem(i,itype_mol).gt.2*Rc_type_2_mem(i,itype_mol)) then
                    write(6,*) 'Rc2_type_2_mem(i,itype_mol).gt.2*Rc_type_2_mem(i,itype_mol),stop',Rc2_type_2_mem(i,itype_mol),Rc_type_2_mem(i,itype_mol)
                endif
            enddo
            close(412)
            
            do i=1,ntype_mem(itype_mol)
                if(n3b1_type_2_mem(i,itype_mol).gt.n3b1m_2_mem(itype_mol)) n3b1m_2_mem(itype_mol)=n3b1_type_2_mem(i,itype_mol)
                if(n3b2_type_2_mem(i,itype_mol).gt.n3b2m_2_mem(itype_mol)) n3b2m_2_mem(itype_mol)=n3b2_type_2_mem(i,itype_mol)
            enddo
            if(n3b1m_2_mem(itype_mol).gt.max_n3b1m_2_mem) max_n3b1m_2_mem=n3b1m_2_mem(itype_mol)
            if(n3b2m_2_mem(itype_mol).gt.max_n3b2m_2_mem) max_n3b2m_2_mem=n3b2m_2_mem(itype_mol)
            
            counter=0
            do itype2=1,ntype_mem(itype_mol)
                do itype1=1,itype2
                    do k1=1,n3b1m_2_mem(itype_mol)
                        do k2=1,n3b1m_2_mem(itype_mol)
                            do k12=1,n3b2m_2_mem(itype_mol)
                                ii_f=0
                                if(itype1.ne.itype2) ii_f=1
                                if((itype1.eq.itype2).and.(k1.le.k2)) ii_f=1
                                if(ii_f.gt.0) then
                                    counter=counter+1
                                endif
                            enddo
                        enddo
                    enddo
                enddo
            enddo
            nfeat0m_2_mem(itype_mol)=counter
        enddo
        
        ! ---------- read grid31_2_2_mem ----------
        allocate(grid31_2_2_mem(2,max_n3b1m_2_mem,max_ntype_mem,ntype_mm))
        grid31_2_2_mem=0
        do itype_mol=1,ntype_mm
            write(str,*) itype_mol
            str=adjustl(str)
            
            if (natom_m(itype_mol) .lt. 2) cycle

            do itype=1,ntype_mem(itype_mol)
                write(str2,*) itype
                str2=adjustl(str2)
                
                open(413,file='mol.'//trim(str)//'/'//'ppp_model'//'/'//'atom_1'//'/'//'output'//'/'//'grid3b_cb12_type3.'//trim(str2))
                rewind(413)
                read(413,*) int_bin
                if(int_bin.ne.n3b1_type_2_mem(itype,itype_mol)) then
                    write(6,*) 'int_bin.ne.n3b1_type_2_mem(itype,itype_mol),stop',int_bin,n3b1_type_2_mem(itype,itype_mol)
                    stop
                endif
                do i=1,n3b1_type_2_mem(itype,itype_mol)
                    read(413,*) int_bin,grid31_2_2_mem(1,i,itype,itype_mol),grid31_2_2_mem(2,i,itype,itype_mol)
                    if(grid31_2_2_mem(2,i,itype,itype_mol).gt.Rc_type_2_mem(itype,itype_mol)) then
                        write(6,*) 'grid31_2_2_mem(2,i,itype,itype_mol).gt.Rc_type_2_mem(itype,itype_mol)',grid31_2_2_mem(2,i,itype,itype_mol),Rc_type_2_mem(itype,itype_mol)
                    endif
                enddo
                close(413)
            enddo
        enddo
        
        ! ---------- read grid32_2_2_mem ----------
        allocate(grid32_2_2_mem(2,max_n3b2m_2_mem,max_ntype_mem,ntype_mm))
        grid32_2_2_mem=0
        do itype_mol=1,ntype_mm
            write(str,*) itype_mol
            str=adjustl(str)
            
            if (natom_m(itype_mol) .lt. 2) cycle

            do itype=1,ntype_mem(itype_mol)
                write(str2,*) itype
                str2=adjustl(str2)
                
                open(414,file='mol.'//trim(str)//'/'//'ppp_model'//'/'//'atom_1'//'/'//'output'//'/'//'grid3b_b1b2_type3.'//trim(str2))
                rewind(414)
                read(414,*) int_bin
                if(int_bin.ne.n3b2_type_2_mem(itype,itype_mol)) then
                    write(6,*) 'int_bin.ne.n3b2_type_2_mem(itype,itype_mol),stop',int_bin,n3b2_type_2_mem(itype,itype_mol)
                    stop
                endif
                do i=1,n3b2_type_2_mem(itype,itype_mol)
                    read(414,*) int_bin,grid32_2_2_mem(1,i,itype,itype_mol),grid32_2_2_mem(2,i,itype,itype_mol)
                    if(grid32_2_2_mem(2,i,itype,itype_mol).gt.Rc2_type_2_mem(itype,itype_mol)) then
                        write(6,*) 'grid32_2_2_mem(2,i,itype,itype_mol).gt.Rc2_type_2_mem(itype,itype_mol)',grid32_2_2_mem(2,i,itype,itype_mol),Rc2_type_2_mem(itype,itype_mol)
                    endif
                enddo
                close(414)
            enddo
        enddo
        
        ! ---------- read num_neigh_mem, iat_neigh_mem, list_neigh_mem ----------
        ! ---------- calc num_neigh_alltype_mem, iat_neigh_alltype_mem, list_neigh_alltype_mem, ind_all_neigh_mem ----------
        allocate(num_neigh_mem(max_ntype_mem,natom_mm,ntype_mm))
        allocate(iat_neigh_mem(max_m_neigh_mem,max_ntype_mem,natom_mm,ntype_mm))
        allocate(list_neigh_mem(max_m_neigh_mem,max_ntype_mem,natom_mm,ntype_mm))
        allocate(num_neigh_alltype_mem(natom_mm,ntype_mm))
        allocate(iat_neigh_alltype_mem(max_m_neigh_mem,natom_mm,ntype_mm))
        allocate(list_neigh_alltype_mem(max_m_neigh_mem,natom_mm,ntype_mm))
        allocate(ind_all_neigh_mem(max_m_neigh_mem,max_ntype_mem,natom_mm,ntype_mm))
        num_neigh_mem=0
        iat_neigh_mem=0
        list_neigh_mem=0
        num_neigh_alltype_mem=0
        iat_neigh_alltype_mem=0
        list_neigh_alltype_mem=0
        ind_all_neigh_mem=0
        do itype_mol=1,ntype_mm
            write(str,*) itype_mol
            str=adjustl(str)
            
            if (natom_m(itype_mol) .lt. 2) cycle

            open(415,file='mol.'//trim(str)//'/'//'ppp_neighbore'//'/'//'IN.NEIGHBORE')
            rewind(415)
            read(415,*) int_bin,int_bin2
            if(int_bin.ne.natom_m(itype_mol)) then
                write(6,*) "int_bin.ne.natom_m(itype_mol),stop",int_bin,natom_m(itype_mol)
                stop
            endif
            if(int_bin2.ne.ntype_mem(itype_mol)) then
                write(6,*) "int_bin2.ne.ntype_mem(itype_mol),stop",int_bin2,ntype_mem(itype_mol)
                stop
            endif

            do i=1,natom_m(itype_mol)
                read(415,*) int_bin,(num_neigh_mem(j,i,itype_mol),j=1,ntype_mem(itype_mol))
            enddo
            do i=1,natom_m(itype_mol)
                do j=1,ntype_mem(itype_mol)
                    do k=1,num_neigh_mem(j,i,itype_mol)
                        read(415,*) int_bin,int_bin2,iat_neigh_mem(k,j,i,itype_mol),list_neigh_mem(k,j,i,itype_mol)
                    enddo
                enddo
            enddo
            close(415)
            
            do iat=1,natom_m(itype_mol)
                num=1
                list_neigh_alltype_mem(1,iat,itype_mol)=iat    ! the first neighbore is itself
                
                do itype=1,ntype_mem(itype_mol)
                    do j=1,num_neigh_mem(itype,iat,itype_mol)
                        num=num+1
                        if(num.gt.m_neigh_mem(itype_mol)) then
                            write(6,*) "total num_neigh.gt.m_neigh_mem(itype_mol),stop",m_neigh_mem(itype_mol)
                            stop
                        endif
                        ind_all_neigh_mem(j,itype,iat,itype_mol)=num
                        iat_neigh_alltype_mem(num,iat,itype_mol)=iat_neigh_mem(j,itype,iat,itype_mol)    ! iat_neigh_alltype(1,iat) is meaningless, do not know the type of the center
                        list_neigh_alltype_mem(num,iat,itype_mol)=list_neigh_mem(j,itype,iat,itype_mol)
                    enddo
                enddo
                num_neigh_alltype_mem(iat,itype_mol)=num
            enddo
        enddo
    end subroutine load_ppp_model
    
    
    subroutine unload_ppp_model()
        implicit none
        
        deallocate(is_fitted_atom)
        deallocate(is_fitted_bond)
        deallocate(ntype_atom)
        deallocate(ntype_bond)
        deallocate(m_neigh_atom)
        deallocate(m_neigh_bond)
        deallocate(natom_mn_atom)
        deallocate(natom_mn_bond)
        deallocate(itype_atom_atom)
        deallocate(itype_atom_bond)
        
        deallocate(iatom_m_atom)
        deallocate(iatom_m_bond)
        
        deallocate(ifeat_type_atom)
        deallocate(ifeat_type_bond)
        deallocate(nfeat1_atom)
        deallocate(nfeat1_bond)
        deallocate(nfeat2_atom)
        deallocate(nfeat2_bond)
        deallocate(nfeat1m_atom)
        deallocate(nfeat1m_bond)
        deallocate(nfeat2m_atom)
        deallocate(nfeat2m_bond)
        deallocate(nfeat2tot_atom)
        deallocate(nfeat2tot_bond)
        deallocate(nfeat2i_atom)
        deallocate(nfeat2i_bond)
        
        deallocate(bb_atom)
        deallocate(bb_bond)
        deallocate(bb_type_atom)
        deallocate(bb_type_bond)
        deallocate(bb_type0_atom)
        deallocate(bb_type0_bond)
        
        deallocate(pv_atom)
        deallocate(pv_bond)
        deallocate(feat2_shift_atom)
        deallocate(feat2_shift_bond)
        deallocate(feat2_scale_atom)
        deallocate(feat2_scale_bond)
        deallocate(pv_scale_bb_atom)
        deallocate(pv_scale_bb_bond)
        deallocate(shift_scale_bb_atom)
        deallocate(shift_scale_bb_bond)
        
        deallocate(Rc_M_mem)
        deallocate(m_neigh_mem)
        deallocate(ntype_mem)
        deallocate(iat_type_mem)
        deallocate(Rc_type_1_mem)
        deallocate(n2b_type_1_mem)
        deallocate(n2bm_1_mem)
        deallocate(nfeat0m_1_mem)
        
        deallocate(grid2_2_1_mem)
        
        deallocate(Rc_type_2_mem)
        deallocate(Rc2_type_2_mem)
        deallocate(n3b1_type_2_mem)
        deallocate(n3b2_type_2_mem)
        deallocate(n3b1m_2_mem)
        deallocate(n3b2m_2_mem)
        deallocate(nfeat0m_2_mem)
        
        deallocate(grid31_2_2_mem)
        
        deallocate(grid32_2_2_mem)
        
        deallocate(num_neigh_mem)
        deallocate(iat_neigh_mem)
        deallocate(list_neigh_mem)
        deallocate(num_neigh_alltype_mem)
        deallocate(iat_neigh_alltype_mem)
        deallocate(list_neigh_alltype_mem)
        deallocate(ind_all_neigh_mem)
    end subroutine
    
    
    subroutine load_model_type1_ppp_from_mem(itype_mol)
        ! considered variables: Rc_M(Rc_M_1),m_neigh(m_neigh_1),ntype(ntype_1),iat_type(iat_type_1),
        ! Rc_type(Rc_type_1),n2b_type(n2b_type_1),iflag_ftype(iflag_ftype_1),n2bm(n2bm_1),
        ! nfeat0m(nfeat0m_1),grid2_2(grid2_2_1)
        
        implicit none
        
        integer itype_mol
        integer i,j,k
        
        Rc_M_1=Rc_M_mem(itype_mol)
        m_neigh_1=m_neigh_mem(itype_mol)
        ntype_1=ntype_mem(itype_mol)
        iat_type_1=0
        iat_type_1(1:ntype_1)=iat_type_mem(1:ntype_1,itype_mol)
        Rc_type_1=0.d0
        Rc_type_1(1:ntype_1)=Rc_type_1_mem(1:ntype_1,itype_mol)
        n2b_type_1=0
        n2b_type_1(1:ntype_1)=n2b_type_1_mem(1:ntype_1,itype_mol)
        iflag_ftype_1=4
        n2bm_1=n2bm_1_mem(itype_mol)
        nfeat0m_1=nfeat0m_1_mem(itype_mol)
        if(allocated(grid2_2_1)) then
            deallocate(grid2_2_1)
        endif
        ! must be consistent with gen_2b_feature.f90, otherwise the size of grid2_2 in 
        ! gen_2b_feature.f90 will be different with grid2_2 in find_feature_2b_type3.f90, 
        ! and error will occur when passing grid2_2
        allocate(grid2_2_1(2,n2bm_1+1,ntype_1))
        grid2_2_1=0.d0
        grid2_2_1(:,1:n2bm_1,1:ntype_1)=grid2_2_1_mem(:,1:n2bm_1,1:ntype_1,itype_mol)
    end subroutine load_model_type1_ppp_from_mem
    
    
    subroutine load_model_type2_ppp_from_mem(itype_mol)
        ! considered variables: Rc_M(Rc_M_2),m_neigh(m_neigh_2),ntype(ntype_2),iat_type(iat_type_2),
        ! Rc_type(Rc_type_2),Rc2_type(Rc2_type_2),n3b1_type(n3b1_type_2),n3b2_type(n3b2_type_2),
        ! iflag_ftype(iflag_ftype_2),n3b1m(n3b1m_2),n3b2m(n3b2m_2),nfeat0m(nfeat0m_2)
        
        implicit none
        
        integer itype_mol
        integer i,j,k
        
        Rc_M_2=Rc_M_mem(itype_mol)
        m_neigh_2=m_neigh_mem(itype_mol)
        ntype_2=ntype_mem(itype_mol)
        iat_type_2=0
        iat_type_2(1:ntype_2)=iat_type_mem(1:ntype_2,itype_mol)
        Rc_type_2=0.d0
        Rc_type_2(1:ntype_2)=Rc_type_2_mem(1:ntype_2,itype_mol)
        Rc2_type_2=0.d0
        Rc2_type_2(1:ntype_2)=Rc2_type_2_mem(1:ntype_2,itype_mol)
        n3b1_type_2=0
        n3b1_type_2(1:ntype_2)=n3b1_type_2_mem(1:ntype_2,itype_mol)
        n3b2_type_2=0
        n3b2_type_2(1:ntype_2)=n3b2_type_2_mem(1:ntype_2,itype_mol)
        iflag_ftype_2=4
        n3b1m_2=n3b1m_2_mem(itype_mol)
        n3b2m_2=n3b2m_2_mem(itype_mol)
        nfeat0m_2=nfeat0m_2_mem(itype_mol)
        if(allocated(grid31_2_2)) then
            deallocate (grid31_2_2)
            deallocate (grid32_2_2)
        endif
        allocate(grid31_2_2(2,n3b1m_2,ntype_2))
        allocate(grid32_2_2(2,n3b2m_2,ntype_2))
        grid31_2_2=0.d0
        grid31_2_2(:,1:n3b1m_2,1:ntype_2)=grid31_2_2_mem(:,1:n3b1m_2,1:ntype_2,itype_mol)
        grid32_2_2=0.d0
        grid32_2_2(:,1:n3b2m_2,1:ntype_2)=grid32_2_2_mem(:,1:n3b2m_2,1:ntype_2,itype_mol)
    end subroutine load_model_type2_ppp_from_mem
    
    
    subroutine load_neighbore_from_mem(itype_mol)
        implicit none
        
        integer,intent(in) :: itype_mol
        
        if(allocated(num_neigh_find_neigh)) then
            deallocate(num_neigh_find_neigh)
            deallocate(num_neigh_alltype_find_neigh)
            deallocate(iat_neigh_find_neigh)
            deallocate(iat_neigh_alltype_find_neigh)
            deallocate(list_neigh_find_neigh)
            deallocate(list_neigh_alltype_find_neigh)
            deallocate(ind_all_neigh_find_neigh)
        endif
        allocate(num_neigh_find_neigh(ntype_mem(itype_mol),natom_m(itype_mol)))
        allocate(num_neigh_alltype_find_neigh(natom_m(itype_mol)))
        allocate(iat_neigh_find_neigh(m_neigh_mem(itype_mol),ntype_mem(itype_mol),natom_m(itype_mol)))
        allocate(iat_neigh_alltype_find_neigh(m_neigh_mem(itype_mol),natom_m(itype_mol)))
        allocate(list_neigh_find_neigh(m_neigh_mem(itype_mol),ntype_mem(itype_mol),natom_m(itype_mol)))
        allocate(list_neigh_alltype_find_neigh(m_neigh_mem(itype_mol),natom_m(itype_mol)))
        allocate(ind_all_neigh_find_neigh(m_neigh_mem(itype_mol),ntype_mem(itype_mol),natom_m(itype_mol)))

        natom_find_neigh=natom_m(itype_mol)
        ntype_find_neigh=ntype_mem(itype_mol)
        m_neigh_find_neigh=m_neigh_mem(itype_mol)
        num_neigh_find_neigh(1:ntype_find_neigh,1:natom_find_neigh)=num_neigh_mem(1:ntype_find_neigh,1:natom_find_neigh,itype_mol)
        num_neigh_alltype_find_neigh(1:natom_find_neigh)=num_neigh_alltype_mem(1:natom_find_neigh,itype_mol)
        iat_neigh_find_neigh(1:m_neigh_find_neigh,1:ntype_find_neigh,1:natom_find_neigh)=iat_neigh_mem(1:m_neigh_find_neigh,1:ntype_find_neigh,1:natom_find_neigh,itype_mol)
        iat_neigh_alltype_find_neigh(1:m_neigh_find_neigh,1:natom_find_neigh)=iat_neigh_alltype_mem(1:m_neigh_find_neigh,1:natom_find_neigh,itype_mol)
        list_neigh_find_neigh(1:m_neigh_find_neigh,1:ntype_find_neigh,1:natom_find_neigh)=list_neigh_mem(1:m_neigh_find_neigh,1:ntype_find_neigh,1:natom_find_neigh,itype_mol)
        list_neigh_alltype_find_neigh(1:m_neigh_find_neigh,1:natom_find_neigh)=list_neigh_alltype_mem(1:m_neigh_find_neigh,1:natom_find_neigh,itype_mol)
        ind_all_neigh_find_neigh(1:m_neigh_find_neigh,1:ntype_find_neigh,1:natom_find_neigh)=ind_all_neigh_mem(1:m_neigh_find_neigh,1:ntype_find_neigh,1:natom_find_neigh,itype_mol)
    end subroutine load_neighbore_from_mem
    
end module mod_ppp
