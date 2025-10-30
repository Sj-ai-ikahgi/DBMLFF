program main_MD

    use mod_mpi
    use mod_data
    use mod_control
    use mod_md
    use calc_ftype1, only: load_model_type1, set_image_info_type1
    use calc_ftype2, only: load_model_type2, set_image_info_type2
    use calc_lin, only: load_model_lin, set_image_info_lin
    use calc_VV, only: set_paths_VV, load_model_VV, set_image_info_VV
    use calc_NN, only: set_paths_NN, load_model_NN, set_image_info_NN
    use mod_param_densityFF, only: Ecut2, fact_kin2, Rbox, n1, n2, n3, &
        nm1_all, nm2_all, nm3_all, ncent, ntype_cent, itype_cent, icent, nat_cent, z_cent, w_cent, &
        read_param_densityFF, ion_type, z_ion, imax_ion, r_ion, rho_ion, &
        ion_type_pol, ion_type_atomp, ion_type_cent, num_term_pol, ind, BB_poldir, &
        ntype_ion, ntype_ion_pol, ipol, atom_charge_param, iflag_polar, &
        itype_cent_d, nat_cent_d, icent_d, z_cent_d, w_cent_d, r_ion_d, atom_charge_param_d, &
        rho_ion_d, imax_ion_d, z_ion_d, ion_type_cent_d, ion_type_atomp_d, &
        c_devptr_itype_cent, c_devptr_nat_cent, c_devptr_icent, c_devptr_z_cent, &
        c_devptr_w_cent, c_devptr_iatom_m, c_devptr_mass_mol, c_devptr_r_ion, &
        c_devptr_rho_ion, c_devptr_imax_ion, c_devptr_z_ion, c_devptr_ion_type_cent, &
        c_devptr_ion_type_atomp, c_devptr_atom_charge_param, E_polar_max, nm1_max, nm2_max, &
        nm3_max, indb, numb, max_numb, BB_direct
    use mod_cuinterface
    use mod_mem_preallocated

    use cudafor

    use mod_ppp, only: num_ppp_model_atom,num_ppp_model_bond, &
        ind_ppp_model_atom,ind_ppp_model_bond, &
        sign_ppp_model_atom,sign_ppp_model_bond, &
        load_ppp_model,unload_ppp_model

    implicit none

    integer ierr
    CHARACTER(LEN=200) :: right
    integer iMD, MDstep
    real*8 dtMD, Temperature1, Temperature2
    logical right_logical
    logical is_reset
    integer nstep_temp_VVMD
    integer iscale_temp_VVMD
    integer i, iat1, natom_t
    integer ii, jj, num
    integer icolor
    integer itype
    integer iter
    integer(kind=8) :: h_bytes, d_bytes
    character*50 file_ionrhoR(20)

    integer iatom_cent(1000, 20)
    integer ip, i1, i2, itype_mol, j
    integer, allocatable, dimension(:, :) :: iatom_type
    integer nparam, indp(200), indpp(3, 3, 200)
    integer numpp, ip1, j1, j2, kk, kkk

    real*8 aL1, aL2, aL3

    integer smooth_ion

    integer(kind=8) :: bytes2GB, free, total

    character(len=20) str
 
    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, inode_tot, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, nnodes_tot, ierr)
    inode_tot = inode_tot + 1

    icolor = inode_tot - 1

    call mpi_comm_split(MPI_COMM_WORLD, icolor, inode, MPI_COMM_MOL, ierr)
    inode = inode + 1
    nnodes = 1

    ! write(6,*) " @@@ icolor ", icolor, inode, nnodes, inode_tot

    if (inode_tot .eq. 1) write (6, *) "Program start !"

    ! Note, in this code, (inode,nnodes) is within each molecule, here
    ! inode=1,nnodes=1, size(MPI_COMM_MOL)=1

    iflag_debug = .true.

    open (9, file="MD.input")
    rewind (9)
    read (9, *) f_xatom
    read (9, *) iMD, MDstep, dtMD, Temperature1, Temperature2
    read (9, *) right_logical  ! IN.MDOPT:q
    read (9, *) MCTRL_stress_step
    read (9, *) iflag_mlff, iflag_model, iflag_charge, iflag_polar, iflag_training, iflag_debug  ! 1: lineear; 2: VV; 3: NN
    read (9, *) Ecut2
    read (9, *) fact_kin2
    read (9, *) Rbox
    Rbox = Rbox/A_AU_1
    Ecut2 = Ecut2/2 ! Ryd to Hartree
    read (9, *) MCTRL_output_nstep
    read (9, *) n1, n2, n3
    read (9, *) ntype_mm   ! number of type of molecule
    read (9, *) (E_polar_max(i), i=1, ntype_mm)
    ! Note, ntype_mm is just a parameter holder, the actual type of molecule
    ! is ntype_m read-in from the xatom.config file, and ntype_m can be
    ! smaller than ntype_mm
    do itype = 1, ntype_mm
        read (9, *) nm1_all(itype), nm2_all(itype), nm3_all(itype)
    end do
    read (9, *) ntype_ion, smooth_ion
    do itype = 1, ntype_ion
        read (9, *) ion_type(itype), z_ion(itype), file_ionrhoR(itype)
    end do
    read(9, *) iflag_dftd3
    read(9, *) iflag_pc
    read(9, *) polar_factor
    close (9)

    if (inode_tot .eq. 1) then
        write(6,*) "@@@ iflag_training, iflag_debug ", iflag_training, iflag_debug, polar_factor
    endif

    do itype = 1, ntype_ion
        open (9, file=file_ionrhoR(itype))
        rewind (9)
        read (9, *)
        read (9, *)
        read (9, *)

        if (smooth_ion .eq. 1) then
            do i = 1, 4000
                read (9, *, iostat=ierr) r_ion(i, itype), rho_ion(i, itype)
                if (ierr .ne. 0) exit
            end do
            imax_ion(itype) = i - 1
        else
            do i = 2, 3000
                read (9, *, iostat=ierr) r_ion(i, itype), rho_ion(i, itype)
                if (ierr .ne. 0) exit
            end do
            imax_ion(itype) = i - 1
            r_ion(1, itype) = 0.d0
            rho_ion(1, itype) = rho_ion(2, itype)
        end if

        close (9)
    end do

    do itype_mol = 1, ntype_mm
        open (9, file="mol."//char(itype_mol + 48)//"/mol.info")
        rewind (9)
        read (9, *) natom_m(itype_mol), ntype_ML(itype_mol)
        close (9)
    end do

    natom_mm = 0
    do itype_mol = 1, ntype_mm
        if (natom_m(itype_mol) .gt. natom_mm) natom_mm = natom_m(itype_mol)
    end do
    ! ntype_ML is the ML_linear model atom type
    allocate (iatom_m(natom_mm, ntype_mm))
    allocate (iatom_type(natom_mm, ntype_mm))
    allocate (mass_mol(natom_mm, ntype_mm))
    allocate (itype_cent(1000, ntype_mm))
    allocate (icent(2, 1000, ntype_mm))
    allocate (nat_cent(1000, ntype_mm))
    allocate (z_cent(1000, ntype_mm))
    allocate (w_cent(2, 1000, ntype_mm))

    do itype_mol = 1, ntype_mm
        open (9, file="mol."//char(itype_mol + 48)//"/mol.info")
        rewind (9)
        read (9, *) natom_m(itype_mol), ntype_ML(itype_mol)
        do i = 1, natom_m(itype_mol)
            read (9, *) i1, iatom_m(i, itype_mol), iatom_type(i, itype_mol), mass_mol(i, itype_mol)
           ! iatom_m is the type used in ML_FF
           ! iatom_type is the type (real type) for ion
        end do
        read (9, *)
        read (9, *) ncent(itype_mol), ntype_cent(itype_mol)
        ! ntype_cent is the spherical atom charge type, including bond centers.
        do ii = 1, ncent(itype_mol)
            read (9, *) itype_cent(ii, itype_mol), z_cent(ii, itype_mol), num, (icent(jj, ii, itype_mol), jj=1, num), &
                (w_cent(jj, ii, itype_mol), jj=1, num), iatom_cent(ii, itype_mol)
            nat_cent(ii, itype_mol) = num
        end do
        close (9)
    end do

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    allocate (ind(2, 4*natom_mm*(natom_mm + 1)/2, ntype_mm))
    allocate (BB_poldir(9*4*natom_mm*(natom_mm + 1)/2, ntype_mm))

    if (iflag_polar) then

        do itype_mol = 1, ntype_mm

            write(str,*) itype_mol
            str=adjustl(str)

            open (25, file="mol."//trim(str)//"/polar_param.input")
            rewind (25)
            read (25, *) natom

            if (natom .ne. natom_m(itype_mol)) then
                write (6, *) "natom in fit_BB_poldir not correct", itype_mol
                stop
            end if
            
            do ii = 1, natom
                read(25,*)
                read(25,*)
                read(25,*)
            enddo

            read (25,*) numb(itype_mol)
            close(25)
        enddo

        max_numb=0
        do i=1,ntype_mm
            if(numb(i).gt.max_numb) max_numb=numb(i)
        enddo

        allocate(num_ppp_model_atom(ntype_mm))
        allocate(num_ppp_model_bond(ntype_mm))
        allocate(ind_ppp_model_atom(6,natom_mm,ntype_mm))
        allocate(ind_ppp_model_bond(4,max_numb,ntype_mm))
        allocate(sign_ppp_model_atom(6,natom_mm,ntype_mm))
        allocate(sign_ppp_model_bond(4,max_numb,ntype_mm))

        do itype_mol = 1, ntype_mm

            write(str,*) itype_mol
            str=adjustl(str)

            open (25, file="mol."//trim(str)//"/polar_param.input")
            rewind (25)
            read (25, *) natom

            if (natom .ne. natom_m(itype_mol)) then
                write (6, *) "natom in fit_BB_poldir not correct", itype_mol
                stop
            end if

            nparam = 0
            num_ppp_model_atom(itype_mol)=0
            num_ppp_model_bond(itype_mol)=0
            do ii = 1, natom
                read (25, *) i1, ipol(1, ii, itype_mol), ipol(2, ii, itype_mol), ipol(3, ii, itype_mol), ipol(4, ii, itype_mol)
                read (25, *) (indp(i), i=1, 6)
                read(25,*) (ind_ppp_model_atom(i,ii,itype_mol),i=1,6),(sign_ppp_model_atom(i,ii,itype_mol),i=1,6)
                do i = 1, 6
                    if (indp(i) .gt. nparam) nparam = indp(i)
                    if(ind_ppp_model_atom(i,ii,itype_mol).gt.num_ppp_model_atom(itype_mol)) num_ppp_model_atom(itype_mol)=ind_ppp_model_atom(i,ii,itype_mol)
                end do
                indpp(1, 1, ii) = indp(1)
                indpp(2, 2, ii) = indp(2)
                indpp(3, 3, ii) = indp(3)
                indpp(1, 2, ii) = indp(4)
                indpp(2, 1, ii) = indp(4)
                indpp(1, 3, ii) = indp(5)
                indpp(3, 1, ii) = indp(5)
                indpp(2, 3, ii) = indp(6)
                indpp(3, 2, ii) = indp(6)
            end do
            read (25, *) numb(itype_mol)
            do ii=1,numb(itype_mol)
                read(25,*) indb(1,ii,itype_mol),indb(2,ii,itype_mol),indb(3,ii,itype_mol), &
                    (ind_ppp_model_bond(i,ii,itype_mol),i=1,4),(sign_ppp_model_bond(i,ii,itype_mol),i=1,4)
                if(indb(3,ii,itype_mol).gt.nparam) nparam=indb(3,ii,itype_mol)
                do i=1,4
                    if(ind_ppp_model_bond(i,ii,itype_mol).gt.num_ppp_model_bond(itype_mol)) num_ppp_model_bond(itype_mol)=ind_ppp_model_bond(i,ii,itype_mol)
                enddo
            enddo

            read (25, *) ntype_ion_pol(itype_mol)
            do ii = 1, ntype_ion_pol(itype_mol)
                read (25, *) &
                    ion_type_pol(ii, itype_mol), (atom_charge_param(j, ii, itype_mol), j=1, 4)
            end do
            close (25)

            if (natom .eq. 1) then
                num_term_pol(itype_mol) = 0
                goto 330
            end if

            numpp = 4*nparam

            open (25, file="mol."//char(itype_mol + 48)//"/fit.BB")
            rewind (25)
            do ip = 1, numpp
                read (25, *) ip1, BB_direct(ip, itype_mol)
            end do
            close (25)
            ip = 0
            do ii = 1, natom
                do i1 = 1, 2
                    do i2 = 1, 2
                        ip = ip + 1
                        ind(1, ip, itype_mol) = (ii - 1)*2 + i1
                        ind(2, ip, itype_mol) = (ii - 1)*2 + i2
                    end do
                end do
            end do
            num_term_pol(itype_mol) = ip
            !    nump = ip

            !1111111111This is really wasteful, since many are zero
            !  We should change this, simplify it to the simple model form
            BB_poldir(:, itype_mol) = 0.d0
            ip = 0
            do jj = 1, natom
                do j1 = 1, 2
                    do j2 = 1, 2
                        do i1 = 1, 3
                            do i2 = 1, 3
                                ip = ip + 1
                                kk = indpp(i1, i2, jj)
                                if (kk .ne. 0) then
                                    kkk = (kk - 1)*4 + j2 + (j1 - 1)*2
                                    BB_poldir(ip, itype_mol) = BB_direct(kkk, itype_mol)
                                end if
                            end do
                        end do
                    end do
                end do
            end do
    330 continue
        end do

        call load_ppp_model()
    endif
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccc

    do itype_mol = 1, ntype_mm
        do ii = 1, ncent(itype_mol)
            ion_type_cent(ii, itype_mol) = 0
            do jj = 1, ntype_ion
                if (iatom_cent(ii, itype_mol) .eq. ion_type(jj)) &
                    ion_type_cent(ii, itype_mol) = jj
            end do
        end do

        if (iflag_polar) then
            do ii = 1, natom_m(itype_mol)
                ion_type_atomp(ii, itype_mol) = 0
                do jj = 1, ntype_ion_pol(itype_mol)
                    if (iatom_type(ii, itype_mol) .eq. ion_type_pol(jj, itype_mol)) &
                        ion_type_atomp(ii, itype_mol) = jj
                end do
                if (ion_type_atomp(ii, itype_mol) .eq. 0) then
                    write (6, *) "atom type not found in ion_type", iatom_type(ii, itype_mol)
                    write (6, *) "ion_type_pol", ion_type_pol(1:ntype_ion_pol(itype_mol), itype_mol)
                    stop
                end if
            end do
        endif
        
    end do

    deallocate (iatom_type)

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccc

    nstep_temp_VVMD = 100
    iscale_temp_VVMD = 0

    call readf_xatom_new(iMD)

    call get_ALI(AL, ALI)

    if (right_logical) then
        call read_mdopt(dtMD, Temperature1, temperature2, imov_at, natom, AL)
    else
        call default_mdopt(dtMD, Temperature1, temperature2, imov_at, natom, AL)
    end if

    !ccccccccccccccccccccccccccccccccccccccc
    if (iMD == 4 .or. iMD == 44) MCTRL_stress = 1
    if (iMD == 5 .or. iMD == 55) MCTRL_stress = 1
    if (iMD == 7 .or. iMD == 77) MCTRL_stress = 1
    if (iMD == 8 .or. iMD == 88) MCTRL_stress = 1
    if (iMD == 8 .or. iMD == 88) then
        MCTRL_is_MSST = .true.
    else
        MCTRL_is_MSST = .false.
    end if

    !ccccccccccccccccccccccccccccccccccccccc
    !init resources and bind process with a gpu
    h_bytes = (3*3*2 + 4)*8
    d_bytes = 3 * 3 * 8 + int(n1, 8) * n2 * n3 * 16 + int(n1, 8) * n2 * n3 * 8 + &
              int(n1, 8) * n2 * n3 * 16 + int(n1, 8) * n2 * n3 * 8 + &
              int(n1, 8) * n2 * n3 * 16 + 4 * 8 + ceiling((int(n1, 8) * n2 * n3) / (1024.d0 * 8)) * 4 * 8
  
    call init_resources(h_bytes, d_bytes, n1, n2, n3, nm1_all, nm2_all, nm3_all, ntype_m)
    ! calculate imax_ncent_itype_mol
    imax_ncent_itype_mol = 0
    do itype_mol = 1, ntype_m
        imax_ncent_itype_mol = max(imax_ncent_itype_mol, ncent(itype_mol))
    end do
    ! calculate nm1_max, nm2_max, nm3_max
    nm1_max = 0
    nm2_max = 0
    nm3_max = 0
    do iter = 1, ntype_m
        nm1_max = max(nm1_max, nm1_all(iter))
        nm2_max = max(nm2_max, nm2_all(iter))
        nm3_max = max(nm3_max, nm3_all(iter))
    end do
    ! memory prealloated

    aL1 = dsqrt(AL(1, 1)**2 + AL(2, 1)**2 + AL(3, 1)**2)
    aL2 = dsqrt(AL(1, 2)**2 + AL(2, 2)**2 + AL(3, 2)**2)
    aL3 = dsqrt(AL(1, 3)**2 + AL(2, 3)**2 + AL(3, 3)**2)
    id1_max = (Rbox/aL1)*n1
    id2_max = (Rbox/aL2)*n2
    id3_max = (Rbox/aL3)*n3

    if (inode_tot .eq. 1) then
        write (6, *) "@@@ id1_max ", id1_max, id2_max, id3_max
        write (6, *) "@@@ imax_ncent_itype_mol ", imax_ncent_itype_mol
    endif

    !id1_max = 60
    !id2_max = 60
    !id3_max = 60
    ierr = cudaMalloc(c_devptr_dxyz_box, 3*(2*id1_max + 1)*(2*id2_max + 1)*(2*id3_max + 1)*8)
    ierr = cudaMalloc(c_devptr_rho, n1*n2*n3*8)
    ierr = cudaMalloc(c_devptr_rho_z, n1*n2*n3*8)
    ierr = cudaMalloc(c_devptr_vxc, n1*n2*n3*8)
    ierr = cudaMalloc(c_devptr_vxc2, n1*n2*n3*8)
    ierr = cudaMalloc(c_devptr_vcoul, n1*n2*n3*8)
    ierr = cudaMalloc(c_devptr_xatom_cent, 3*10*natom_mm*8)
    ierr = cudaMalloc(c_devptr_AL, 3*3*8)
    ierr = cudaMalloc(c_devptr_ALm, 3*3*8)
    ierr = cudaMalloc(c_devptr_xatom_m, 3*natom_mm*nmolm*ntype_m*8)
    ierr = cudaMalloc(c_devptr_force_m, 3*natom_mm*nmolm*ntype_m*8)
    ierr = cudaMalloc(c_devptr_icorner, 3*4)
    ierr = cudaMalloc(c_devptr_force_cent, 3*10*natom_mm*8)

    ierr = cudaMalloc(c_devptr_box, (2*id1_max + 1)*(2*id2_max + 1)*(2*id3_max + 1)*imax_ncent_itype_mol*8)
    ierr = cudaMalloc(c_devptr_box2, (2*id1_max + 1)*(2*id2_max + 1)*(2*id3_max + 1)*imax_ncent_itype_mol*8)
    ierr = cudaMalloc(c_devptr_box3, (2*id1_max + 1)*(2*id2_max + 1)*(2*id3_max + 1)*imax_ncent_itype_mol*8)

    ierr = cudaMalloc(c_devptr_dbox, (2*id1_max + 1)*(2*id2_max + 1)*(2*id3_max + 1)*3*imax_ncent_itype_mol*8)
    ierr = cudaMalloc(c_devptr_dbox2, (2*id1_max + 1)*(2*id2_max + 1)*(2*id3_max + 1)*3*imax_ncent_itype_mol*8)
    ierr = cudaMalloc(c_devptr_dbox3, (2*id1_max + 1)*(2*id2_max + 1)*(2*id3_max + 1)*3*imax_ncent_itype_mol*8)

    if (iflag_polar .eq. 1) then
        ierr = cudaMalloc(c_devptr_dxyzp_store, 3*3*natom_mm*8)
        ierr = cudaMalloc(c_devptr_force_polar_m, 3*natom_mm*nmolm*ntype_mm*8)
        ierr = cudaMalloc(c_devptr_ppp, 3*2*natom_mm*8)
        ierr = cudaMalloc(c_devptr_dppp, 3*2*natom_mm*natom_mm*3*8)
        ierr = cudaMalloc(c_devptr_pxyz, 3*2*natom_mm*8)
        ierr = cudaMalloc(c_devptr_dpxyz, 3*3*2*natom_mm*8)
    end if
    ierr = cudaMalloc(c_devptr_rho_m, nm1_max*nm2_max*nm3_max*8)
    ierr = cudaMalloc(c_devptr_rho_mz, nm1_max*nm2_max*nm3_max*8)
    ierr = cudaMalloc(c_devptr_vxc_m, nm1_max*nm2_max*nm3_max*8)
    ierr = cudaMalloc(c_devptr_vxc2_m, nm1_max*nm2_max*nm3_max*8)
    ierr = cudaMalloc(c_devptr_vcoul_m, nm1_max*nm2_max*nm3_max*8)
    
    ierr = cudaMalloc(c_devptr_dbox_c, (2*id1_max + 1)*(2*id2_max + 1)*(2*id3_max + 1)*3*imax_ncent_itype_mol*8)
    ierr = cudaMalloc(c_devptr_dbox3_c, (2*id1_max + 1)*(2*id2_max + 1)*(2*id3_max + 1)*3*imax_ncent_itype_mol*8)
    ! calc_polar_mem
    ierr = cudaMalloc(c_devptr_fact2_store, 1000*8)
    ierr = cudaMalloc(c_devptr_fact_store, 1000*8)
    ierr = cudaMalloc(c_devptr_CC_pol, 6*natom_mm*nmolm*ntype_mm*8)
    ierr = cudaMalloc(c_devptr_xc_cent, 3*8)
    ierr = cudaMalloc(c_devptr_rhop_tot, n1*n2*n3*8)
    ierr = cudaMalloc(c_devptr_vrhop_tot_coul, n1*n2*n3*8)
    ierr = cudaMalloc(c_devptr_dvxc_tot, n1*n2*n3*8)
    ierr = cudaMalloc(c_devptr_dvxc_m, nm1_max*nm2_max*nm3_max*8)
    ierr = cudaMalloc(c_devptr_rhop_m, nm1_max*nm2_max*nm3_max*8)
    ierr = cudaMalloc(c_devptr_vrhop_m_coul, nm1_max*nm2_max*nm3_max*8)
    !cuda malloc
    ierr = cudaMalloc(c_devptr_itype_cent, 1000*ntype_mm*4)
    ierr = cudaMalloc(c_devptr_nat_cent, 1000*ntype_mm*4)
    ierr = cudaMalloc(c_devptr_icent, 2*1000*ntype_mm*4)
    ierr = cudaMalloc(c_devptr_z_cent, 1000*ntype_mm*8)
    ierr = cudaMalloc(c_devptr_w_cent, 2*1000*ntype_mm*8)
    ierr = cudaMalloc(c_devptr_iatom_m, natom_mm*ntype_mm*4)
    ierr = cudaMalloc(c_devptr_mass_mol, natom_mm*ntype_mm*8)
    ierr = cudaMalloc(c_devptr_r_ion, 5000*20*8)
    ierr = cudaMalloc(c_devptr_rho_ion, 5000*20*8)
    ierr = cudaMalloc(c_devptr_imax_ion, 20*4)
    ierr = cudaMalloc(c_devptr_z_ion, 20*8)
    ierr = cudaMalloc(c_devptr_ion_type_cent, 1000*50*4)
    ierr = cudaMalloc(c_devptr_ion_type_atomp, 1000*50*4)
    ierr = cudaMalloc(c_devptr_atom_charge_param, 4*20*50*8)

    ierr = cudaMemGetInfo(free,total)
    bytes2GB = 1024*1024*1024

    if (inode_tot .eq. 1) then
        write(6,*) "@@@ Cuda Mem(free/Total GB): ", inode_tot, 1.0*free/bytes2GB,"GB/",1.0*total/bytes2GB,"GB"
    endif
    
    call c_f_pointer(c_devptr_itype_cent, itype_cent_d, (/1000, ntype_mm/))
    call c_f_pointer(c_devptr_nat_cent, nat_cent_d, (/1000, ntype_mm/))
    call c_f_pointer(c_devptr_icent, icent_d, (/2, 1000, ntype_mm/))
    call c_f_pointer(c_devptr_z_cent, z_cent_d, (/1000, ntype_mm/))
    call c_f_pointer(c_devptr_w_cent, w_cent_d, (/2, 1000, ntype_mm/))
    call c_f_pointer(c_devptr_iatom_m, iatom_m_d, (/natom_mm, ntype_mm/))
    call c_f_pointer(c_devptr_mass_mol, mass_mol_d, (/natom_mm, ntype_mm/))
    call c_f_pointer(c_devptr_r_ion, r_ion_d, (/5000, 20/))
    call c_f_pointer(c_devptr_rho_ion, rho_ion_d, (/5000, 20/))
    call c_f_pointer(c_devptr_imax_ion, imax_ion_d, (/20/))
    call c_f_pointer(c_devptr_z_ion, z_ion_d, (/20/))
    call c_f_pointer(c_devptr_ion_type_cent, ion_type_cent_d, (/1000, 50/))
    call c_f_pointer(c_devptr_ion_type_atomp, ion_type_atomp_d, (/1000, 50/))
    call c_f_pointer(c_devptr_atom_charge_param, atom_charge_param_d, (/4, 20, 50/))

    iatom_m_d = iatom_m
    mass_mol_d = mass_mol
    itype_cent_d = itype_cent
    icent_d = icent
    nat_cent_d = nat_cent
    z_cent_d = z_cent
    w_cent_d = w_cent
    r_ion_d = r_ion
    rho_ion_d = rho_ion
    imax_ion_d = imax_ion
    z_ion_d = z_ion
    ion_type_cent_d = ion_type_cent
    ion_type_atomp_d = ion_type_atomp
    atom_charge_param_d = atom_charge_param

    call read_param_densityFF(ntype_m)

    ! cpu memmory preallocated
    if (iflag_polar .eq. 1) then
        allocate (CC_pol(6, natom_mm, nmolm, ntype_mm))
        allocate (dxyzp_store(3, 3, natom_mm))
        allocate (force_polar_m(3, natom_mm, nmolm, ntype_mm))
        allocate (fact_polar(nmolm, ntype_mm))

        allocate (ppp(3, 2*natom_mm))
        allocate (dppp(3, 2*natom_mm, natom_mm, 3))
        allocate (pxyz(3, 2*natom_mm))
        allocate (dpxyz(3, 3, 2*natom_mm))
    end if

    call molecular_dynamics_new(dtMD,iMD,MDstep,xatom,iMDatom,&
           temperature1,temperature2,iscale_temp_VVMD,nstep_temp_VVMD,&
           imov_at,natom,AL,ALI,iatom,f_xatom)

    call mpi_finalize(ierr)
    ierr = cudaFree(c_devptr_itype_cent)
    ierr = cudaFree(c_devptr_nat_cent)
    ierr = cudaFree(c_devptr_icent)
    ierr = cudaFree(c_devptr_z_cent)
    ierr = cudaFree(c_devptr_w_cent)
    ierr = cudaFree(c_devptr_iatom_m)
    ierr = cudaFree(c_devptr_mass_mol)
    ierr = cudaFree(c_devptr_r_ion)
    ierr = cudaFree(c_devptr_rho_ion)
    ierr = cudaFree(c_devptr_imax_ion)
    ierr = cudaFree(c_devptr_z_ion)
    ierr = cudaFree(c_devptr_ion_type_cent)
    ierr = cudaFree(c_devptr_ion_type_atomp)
    ierr = cudaFree(c_devptr_atom_charge_param)
    call destroy_resources()
    ierr = cudaFree(c_devptr_box)
    ierr = cudaFree(c_devptr_box2)
    ierr = cudaFree(c_devptr_box3)
    ierr = cudaFree(c_devptr_dbox)
    ierr = cudaFree(c_devptr_dbox2)
    ierr = cudaFree(c_devptr_dbox3)
    ierr = cudaFree(c_devptr_vxc)
    ierr = cudaFree(c_devptr_vxc2)
    ierr = cudaFree(c_devptr_vcoul)
    ierr = cudaFree(c_devptr_dxyz_box)
    ierr = cudaFree(c_devptr_AL)
    ierr = cudaFree(c_devptr_ALm)
    ierr = cudaFree(c_devptr_xatom_cent)
    ierr = cudaFree(c_devptr_force_m)
    ierr = cudaFree(c_devptr_icorner)
    ierr = cudaFree(c_devptr_force_cent)
    ierr = cudaFree(c_devptr_rho)
    ierr = cudaFree(c_devptr_rho_z)
    ierr = cudaFree(c_devptr_xatom_m)
    if (iflag_polar .eq. 1) then
        ierr = cudaFree(c_devptr_dxyzp_store)
        ierr = cudaFree(c_devptr_ppp)
        ierr = cudaFree(c_devptr_dppp)
        ierr = cudaFree(c_devptr_pxyz)
        ierr = cudaFree(c_devptr_dpxyz)
        ierr = cudaFree(c_devptr_force_polar_m)
    end if
    ierr = cudaFree(c_devptr_rho_m)
    ierr = cudaFree(c_devptr_rho_mz)
    ierr = cudaFree(c_devptr_vxc_m)
    ierr = cudaFree(c_devptr_vxc2_m)
    ierr = cudaFree(c_devptr_vcoul_m)
    ierr = cudaFree(c_devptr_dbox_c)
    ierr = cudaFree(c_devptr_dbox3_c)
    ! calc_polar mem free
    ierr = cudaFree(c_devptr_fact2_store)
    ierr = cudaFree(c_devptr_fact_store)
    ierr = cudaFree(c_devptr_CC_pol)
    ierr = cudaFree(c_devptr_xc_cent)
    ierr = cudaFree(c_devptr_rhop_tot)
    ierr = cudaFree(c_devptr_vrhop_tot_coul)
    ierr = cudaFree(c_devptr_dvxc_tot)
    ierr = cudaFree(c_devptr_dvxc_m)
    ierr = cudaFree(c_devptr_rhop_m)
    ierr = cudaFree(c_devptr_vrhop_m_coul)
    ! cpu memory free
    if (iflag_polar .eq. 1) then
        deallocate (dxyzp_store)
        deallocate (ppp)
        deallocate (dppp)
        deallocate (pxyz)
        deallocate (dpxyz)
        deallocate (CC_pol)
        deallocate (force_polar_m)
        deallocate (fact_polar)

        deallocate(num_ppp_model_atom)
        deallocate(num_ppp_model_bond)
        deallocate(ind_ppp_model_atom)
        deallocate(ind_ppp_model_bond)
        deallocate(sign_ppp_model_atom)
        deallocate(sign_ppp_model_bond)
        call unload_ppp_model()
    end if

    if (inode_tot .eq. 1) write (6, *) "Program end !"
end

