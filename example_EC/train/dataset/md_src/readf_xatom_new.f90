subroutine readf_xatom_new(iMD)

    use mod_data
    use mod_mpi
    use mod_control
    use mod_md, only: A_AU_1, Hartree_eV

    implicit none
    ! real*8, parameter :: A_AU_1=0.52917721067d0
    ! real*8, parameter :: Hartree_eV=27.21138602d0 !(eV)
    ! real*8, parameter :: A_AU_1 = 0.529177d0
    ! real*8, parameter :: Hartree_eV = 27.211396d0 !(eV)

    real*8, allocatable, dimension(:, :) :: xatom_tmp
    integer, allocatable, dimension(:, :) :: imov_latt_vectmp
    integer, allocatable, dimension(:)   ::  iatom_tmp, iatom_tmp2
    real*8, allocatable, dimension(:) :: VX_1_tmp, VY_1_tmp, VZ_1_tmp
    real*8, allocatable, dimension(:) :: weight_mag_tmp
    real*8, allocatable, dimension(:, :) :: weight_mag_tmpxyz
    real*8, allocatable, dimension(:) ::langevin_factT_tmp, langevin_factG_tmp
    real*8, allocatable, dimension(:) :: constraint_mag_atom_tmp
    real*8, allocatable, dimension(:) :: constraint_mag_alpha_tmp
    real*8, allocatable, dimension(:, :) :: LDAU_lambda_tmp
    integer, allocatable, dimension(:)   :: ind_order_old2new_tmp, ind_order_new2old_tmp
    real*8, allocatable, dimension(:) :: weight_atom_tmp
    integer natom_m_t(1000)
    character*200 message

    integer ierr, i, ia, jjj, imol, i1, jj, kk
    integer iMD
    integer itype_tmp, imov_sum, ii, ncount
    integer num_test
    integer itype

    logical :: scanit
    intrinsic :: sum

    open (10, file=f_xatom, status='old', action='read', iostat=ierr)

    if (ierr .ne. 0) then
        if (inode_tot .eq. 1) write (message, *) "IN.ATOM", f_xatom, "not exist, stop"
        write (*, *) message
        stop
    end if

    rewind (10)
    !
    read (10, *) natom, ntype_m, (natom_m_t(ii), ii=1, ntype_m), (nmol(ii), ii=1, ntype_m)
    if (ntype_m .gt. ntype_mm) then
        write (6, *) "ntype_m.gt.ntype_mm,stop", ntype_m, ntype_mm
        stop
    end if

    do ii = 1, ntype_m
        if (natom_m_t(ii) .ne. natom_m(ii)) then   ! natom_m from mod_data
            write (6, *) "natom_m from xatom.config not same as in MD.input", natom_m_t, natom_m
            stop
        end if
    end do

    num_test = 0
    do ii = 1, ntype_m
        num_test = num_test + natom_m(ii)*nmol(ii)
    end do

    if (natom .ne. num_test) then
        write (6, *) "natom.ne.nmol*natom_m,stop", natom, num_test
        stop
    end if

    if (natom .gt. matom_1) then
        if (inode_tot .eq. 1) then
            write (message, *) "natom.gt.matom_1, increase matom_1 in data.f,stop ", f_xatom, natom, matom_1
            write (*, *) message
            stop
        end if
    end if

    nmolm = 0
    do itype = 1, ntype_m
        if (nmol(itype) .gt. nmolm) nmolm = nmol(itype)
    end do

    allocate (xatom_m(3, natom_mm, nmolm, ntype_m))

    allocate (xatom_tmp(3, natom))
    allocate (imov_latt_vectmp(3, natom))
    allocate (iatom_tmp(natom))
    allocate (iatom_tmp2(natom))
    allocate (VX_1_tmp(natom))
    allocate (VY_1_tmp(natom))
    allocate (VZ_1_tmp(natom))
    allocate (weight_mag_tmp(natom))
    allocate (weight_mag_tmpxyz(natom, 3))
    allocate (langevin_factT_tmp(natom))
    allocate (langevin_factG_tmp(natom))
    allocate (constraint_mag_atom_tmp(natom))
    allocate (constraint_mag_alpha_tmp(natom))
    allocate (ind_order_old2new_tmp(natom))
    allocate (ind_order_new2old_tmp(natom))
    allocate (LDAU_lambda_tmp(natom, 2))
    allocate (weight_atom_tmp(natom))

    call scan_key_words(10, "LATTICE", len("LATTICE"), scanit)
    if (scanit) then
        read (10, *) (AL(i, 1), i=1, 3)
        read (10, *) (AL(i, 2), i=1, 3)
        read (10, *) (AL(i, 3), i=1, 3)
    else
        write (message, *) "Must provide LATTICE in IN.ATOM file ", ADJUSTL(trim(f_xatom))
        write (*, *) message
        stop
    end if

    AL = AL/A_AU_1 !zhilin

    ! print *,"AL raw",AL
    ! print *,"A_AU_1 raw",A_AU_1

    call scan_key_words(10, "STRESS_MASK", len("STRESS_MASK"), scanit)
    if (.not. scanit) then
        stress_mask = 1
    else
        read (10, *) (stress_mask(i, 1), i=1, 3)
        read (10, *) (stress_mask(i, 2), i=1, 3)
        read (10, *) (stress_mask(i, 3), i=1, 3)
    end if
    call scan_key_words(10, "STRESS_EXTERNAL", len("STRESS_EXTERNAL"), scanit)
    if (.not. scanit) then
        exist_stress_ext = .false.
    else
        exist_stress_ext = .true.
        read (10, *) (stress_ext(i, 1), i=1, 3)
        read (10, *) (stress_ext(i, 2), i=1, 3)
        read (10, *) (stress_ext(i, 3), i=1, 3)
        stress_ext = stress_ext/Hartree_ev*natom
    end if

    call scan_key_words(10, "POSITION", len("POSITION"), scanit)
    if (.not. scanit) then
        write (*, *) "keyword 'position' is needed at", f_xatom
        if (inode_tot .eq. 1) then
            call mpi_abort(mpi_comm_world, ierr)
        end if
    else

        i = 0
        do itype = 1, ntype_m
            do imol = 1, nmol(itype)
                do i1 = 1, natom_m(itype)
                    i = i + 1
                    ! read(10,*) iatom_tmp(i),xatom_tmp(1,i),xatom_tmp(2,i),xatom_tmp(3,i),imov_latt_vectmp(1,i),imov_latt_vectmp(2,i),imov_latt_vectmp(3,i), jj
                    ! if(jj.ne.imol) then
                    !     write(6,*) "jj.me.imol in xatom.config, safety check"
                    !     stop
                    ! endif
                    read(10,*) iatom_tmp(i),xatom_tmp(1,i),xatom_tmp(2,i),xatom_tmp(3,i),imov_latt_vectmp(1,i),imov_latt_vectmp(2,i),imov_latt_vectmp(3,i)

                    xatom_m(1, i1, imol, itype) = xatom_tmp(1, i)
                    xatom_m(2, i1, imol, itype) = xatom_tmp(2, i)
                    xatom_m(3, i1, imol, itype) = xatom_tmp(3, i)
                    !   The order of atoms will not be changed in this version.

                end do
            end do
        end do
        !!!!!  check imov avoiding imov==0 in POSITION-RELAX MD NEB TDDFT NAMD
        imov_sum = sum(imov_latt_vectmp)
        !!!!!
    end if

    call scan_key_words(10, "VELOCITY", len("VELOCITY"), scanit)
    if (.not. scanit) then
        !write(*,*) "keyword 'velocity' is needed at",f_xatom
        !stop
        if ((iMD .eq. 11) .or. (iMD .eq. 22) .or. (iMD .eq. 33)) then
            write (*, *) "keyword 'velocity' is needed at", f_xatom
            if (inode_tot .eq. 1) call mpi_abort(mpi_comm_world, ierr)
        end if
        exist_velocity = .false.
    end if
    if (scanit) then
        if (inode_tot .eq. 1) then
            write (*, *) "velocity read in from ", f_xatom
        end if
        do i = 1, natom
            read (10, *) iatom_tmp2(i), VX_1_tmp(i), VY_1_tmp(i), VZ_1_tmp(i)
            if (iatom_tmp2(i) .ne. iatom_tmp(i)) then
                write (6, *) "order of iatom in position/velocity not the same", i
                if (inode_tot .eq. 1) then
                    call mpi_abort(mpi_comm_world, ierr)
                end if
            end if
        end do
        exist_velocity = .true.
    end if

    call scan_key_words(10, "LANGEVIN_ATOMFACT_TG", len("LANGEVIN_ATOMFACT_TG"), scanit)
    if (.not. scanit) then

        do i = 1, natom
            langevin_factT_tmp(i) = 1.d0
            langevin_factG_tmp(i) = 1.d0
        end do
    else
        do i = 1, natom
            read (10, *) iatom_tmp2(i), langevin_factT_tmp(i), langevin_factG_tmp(i)
        end do
    end if

    !********************************************************
    !********************************************************
    call scan_key_words(10, "Weight_atom", len("Weight_atom"), scanit)
    if (.not. scanit) then
        do i = 1, natom
            weight_atom_tmp(i) = 1.d0
        end do
    else
        do i = 1, natom
            read (10, *) iatom_tmp2(i), weight_atom_tmp(i)
        end do
    end if
    !********************************************************

    close (10)

    do i = 1, natom
        ii = i
        iatom(ii) = iatom_tmp(i)
        xatom(:, ii) = xatom_tmp(:, i)
        imov_at(:, ii) = imov_latt_vectmp(:, i)
        VX_1(ii) = VX_1_tmp(i)
        VY_1(ii) = VY_1_tmp(i)
        VZ_1(ii) = VZ_1_tmp(i)
        langevin_factT(ii) = langevin_factT_tmp(i)
        langevin_factG(ii) = langevin_factG_tmp(i)
        ind_order_new2old_tmp(ii) = i
        ind_order_old2new_tmp(i) = ii
        weight_atom_1(ii) = weight_atom_tmp(i)

    end do

    iMDatom(1:natom) = 1.d0
    kk = 0
    do itype = 1, ntype_m
        do imol = 1, nmol(itype)
            do ii = 1, natom_m(itype)
                kk = kk + 1
                iMDatom(kk) = mass_mol(ii, itype)
            end do
        end do
    end do

    if (inode_tot .eq. 1) then
        open (unit=2200, file='final.config')
        rewind (2200)
        write (2200, *) natom
        write (2200, *) "Lattice vector (Angstrom), stress(eV/natom)"
        write (2200, "(3(E19.10,1x))") A_AU_1*AL(1, 1), A_AU_1*AL(2, 1), A_AU_1*AL(3, 1)
        write (2200, "(3(E19.10,1x))") A_AU_1*AL(1, 2), A_AU_1*AL(2, 2), A_AU_1*AL(3, 2)
        write (2200, "(3(E19.10,1x))") A_AU_1*AL(1, 3), A_AU_1*AL(2, 3), A_AU_1*AL(3, 3)
        write (2200, *) "Position, move_x, move_y, move_z"
        do ia = 1, natom
            write (2200, 1113) iatom(ia), xatom(1, ia), xatom(2, ia), xatom(3, ia), imov_at(1, ia), imov_at(2, ia), imov_at(3, ia)
        end do
        close (2200)
1113    format(i4, 1x, 3(f14.9, 1x), 4x, 3(i1, 2x))
    end if

    deallocate (xatom_tmp)
    deallocate (imov_latt_vectmp)
    deallocate (iatom_tmp)
    deallocate (iatom_tmp2)
    deallocate (VX_1_tmp)
    deallocate (VY_1_tmp)
    deallocate (VZ_1_tmp)
    deallocate (weight_mag_tmp)
    deallocate (weight_mag_tmpxyz)
    deallocate (langevin_factT_tmp)
    deallocate (langevin_factG_tmp)
    deallocate (constraint_mag_atom_tmp)
    deallocate (constraint_mag_alpha_tmp)
    deallocate (ind_order_old2new_tmp)
    deallocate (ind_order_new2old_tmp)
    deallocate (LDAU_lambda_tmp)
    deallocate (weight_atom_tmp)

    return
end subroutine readf_xatom_new

