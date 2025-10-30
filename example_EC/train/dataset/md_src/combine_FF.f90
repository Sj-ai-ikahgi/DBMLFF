subroutine combine_FF(Etot, fatom, xatom, iatom, AL, natom_tot, natom_m, nmol)

    use mod_data, only: natom_mm, nmolm, ntype_m, iflag_mlff, iflag_charge, iflag_pc, &
        iflag_debug
    use mod_profile
    use mod_mpi
    use mod_md, only: A_AU_1, Hartree_eV

    implicit none

    integer natom_tot, natom_m(ntype_m), nmol(ntype_m)

    real*8 Etot, Etot1, Etot2
    real*8 fatom(3, natom_tot)
    real*8 xatom(3, natom_tot)
    real*8 AL(3, 3), AL_tmp(3, 3)
    real*8 xatom_m(3, natom_mm, nmolm, ntype_m)
    real*8 force_m(3, natom_mm, nmolm, ntype_m)
    integer iatom_m(natom_mm, nmolm, ntype_m)
    integer iatom(natom_tot)
    integer kkk1, imol, ii
    integer itype_mol
    real*8 Etot_molecule_pair_correction
    integer :: i

    integer ierr

    if (switch) combine_FF_start = mpi_wtime()

    kkk1 = 0
    do itype_mol = 1, ntype_m
        do imol = 1, nmol(itype_mol)
            do ii = 1, natom_m(itype_mol)
                kkk1 = kkk1 + 1
                xatom_m(:, ii, imol, itype_mol) = xatom(:, kkk1)
                iatom_m(ii, imol, itype_mol) = iatom(kkk1)
            end do
        end do
    end do

    AL_tmp = AL*A_AU_1

    Etot = 0.d0
    fatom = 0.d0
    
    if (iflag_mlff) then
        if (switch) ML_FF_EF_start = mpi_wtime()

        call ML_FF_EF(Etot1, fatom, xatom, AL_tmp, natom_tot)

        Etot = Etot1/Hartree_eV  ! convert back to Hartree
        fatom(:, 1:natom_tot) = fatom(:, 1:natom_tot)*A_AU_1/Hartree_eV  ! convert to Hatree/Bohr

        if (switch) ML_FF_EF_end = mpi_wtime()
    end if


    if ((inode_tot .eq. 1) .and. iflag_debug .and. iflag_mlff) then
        open (36, file="OUT.ML_FF", position="append")
        write (36, *) Etot1
        do i = 1, natom_tot
           write (36, "(i4, 1x, 3(E25.17,1X))") i, fatom(:, i)*Hartree_eV/A_AU_1
        end do
        close (36)
    end if

    if (iflag_charge .eq. 1) then

        if (switch) nonbond_force_start = mpi_wtime()

        force_m = 0.d0
        call nonbond_force(natom_m, nmol, xatom_m, AL, force_m, Etot2)

        Etot = Etot + Etot2
        kkk1 = 0
        do itype_mol = 1, ntype_m
            do imol = 1, nmol(itype_mol)
                do ii = 1, natom_m(itype_mol)
                    kkk1 = kkk1 + 1
                    fatom(:, kkk1) = fatom(:, kkk1) + force_m(:, ii, imol, itype_mol)
                end do
            end do
        end do

        if (switch) nonbond_force_end = mpi_wtime()


        if ((inode_tot .eq. 1) .and. iflag_debug) then
            open (36, file="OUT.nonbond_force", position="append")
            write (36, *) Etot2*Hartree_eV
            i = 0
            do itype_mol = 1, ntype_m
                do imol = 1, nmol(itype_mol)
                    do ii = 1, natom_m(itype_mol)
                        i = i + 1
                        write (36, "(i4, 1x, 3(E25.17,1X))") i, force_m(:, ii, imol, itype_mol)*Hartree_eV/A_AU_1
                    end do
                end do
            end do
            close (36)
        end if


    end if


!    !======================================================================
!    ! 20230224 molecule pair correction
!    !======================================================================
!     if (iflag_pc .eq. 1) then

!         force_m = 0.d0
!         ! AL have been scaled with A_AU_1 in readf_xatom_new.f90
!         call ML_FF_EF_PC(natom_m, nmol, xatom_m, AL_tmp, force_m, Etot_molecule_pair_correction)
  
!         Etot = Etot + Etot_molecule_pair_correction/Hartree_eV
!         kkk1 = 0
!         do itype_mol = 1, ntype_m
!            do imol = 1, nmol(itype_mol)
!               do ii = 1, natom_m(itype_mol)
!                  kkk1 = kkk1 + 1
!                  fatom(:, kkk1) = fatom(:, kkk1) + force_m(:, ii, imol, itype_mol)*A_AU_1/Hartree_eV
!               end do
!            end do
!         end do


!         if ((inode_tot .eq. 1) .and. iflag_debug) then
!             open (36, file="OUT.PC", position="append")
!             write (36, *) Etot_molecule_pair_correction
!             i = 0
!             do itype_mol = 1, ntype_m
!                 do imol = 1, nmol(itype_mol)
!                     do ii = 1, natom_m(itype_mol)
!                         i = i + 1
!                         write (36, "(i4, 1x, 3(E25.17,1X))") i, force_m(:, ii, imol, itype_mol)
!                     end do
!                 end do
!             end do
!             close (36)
!         end if


!     endif


    if ((inode_tot .eq. 1) .and. iflag_debug) then
        open (36, file="OUT.energy_force", position="append")
        write (36, *) Etot*Hartree_eV
        do i = 1, natom_tot
           write (36, "(i4, 1x, 3(E25.17,1X))") i, fatom(:, i)*Hartree_eV/A_AU_1
        end do
        close (36)
     end if


    if (switch) combine_FF_end = mpi_wtime()

    return
end subroutine combine_FF

