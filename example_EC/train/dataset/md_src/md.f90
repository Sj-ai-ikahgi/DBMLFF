subroutine molecular_dynamics_kernel(Etot, fatom, e_stress)
    use mod_md
    use mod_mpi
    use mod_control, only: MCTRL_iMD, MCTRL_AL, MCTRL_output_nstep, &
        MCTRL_stress, MCTRL_stress_step
    use mod_data
    use mod_profile
    ! This is a bit crazy: mod_md also used mod_control, in that way, through
    ! mod_md, we have all the contains of mod_control, so not really needed the
    ! mod_control here.
    ! The MCTRL_xatom, and MCTRL_iatom are accessed from mod_md

    implicit none
    integer precision_flag
    real(8) fatom_old(3, matom_1)
    real(8) fatom(3, matom_1)
    real(8) e_stress(3, 3), e_stress_old(3, 3)
    real*8 Etot
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !mymd
    type(type_md) :: md
    integer istep
    character(len=200) message
    real*8 tt
    integer ierr

    !data for interpolation
    real*8 :: stress_plumed(3, 3)
    real*8 AL_tmp(3, 3), AL_old(3, 3)

    real*8 Etot_1, Etot_2, Etot_3

    real*8 tt1, tt2

    
    call init_type_md(md)

    MCTRL_ido_stop = 0
    MCTRL_ido_ns = 1

    if (MCTRL_stress .eq. 1) then
        AL_old = MCTRL_AL
        AL(1, :) = AL_old(1, :)*1.01   ! numerical stress calculation
        call combine_FF(Etot_1, fatom, MCTRL_xatom, MCTRL_iatom, AL, &
            MCTRL_natom, natom_m, nmol)
        AL(2, :) = AL_old(2, :)*1.01   ! numerical stress calculation
        call combine_FF(Etot_2, fatom, MCTRL_xatom, MCTRL_iatom, AL, &
            MCTRL_natom, natom_m, nmol)
        AL(3, :) = AL_old(3, :)*1.01   ! numerical stress calculation
        call combine_FF(Etot_3, fatom, MCTRL_xatom, MCTRL_iatom, AL, &
            MCTRL_natom, natom_m, nmol)
        AL = AL_old
    end if

    call combine_FF(Etot, fatom, MCTRL_xatom, MCTRL_iatom, MCTRL_AL, &
        MCTRL_natom, natom_m, nmol)

    if (MCTRL_stress .eq. 1) then
        e_stress = 0.d0
        e_stress(1, 1) = (Etot_1 - Etot)/0.01
        e_stress(2, 2) = (Etot_2 - Etot)/0.01
        e_stress(3, 3) = (Etot_3 - Etot)/0.01
        e_stress_old = e_stress
    else
        e_stress = 0.d0
    end if

    call exchange_data_scf2md(md, fatom, Etot, e_stress)
    !
    if (inode_tot .eq. 1) then
        print *, '***********************************************'
        print *, '*                                             *'
        print *, '*               MD starts Now!                *'
        print *, '*                                             *'
        print *, '***********************************************'
    end if
    
    call md_init(md)

    !ccccccccccc  Need to do one fatom before MD step, to get the fatom
    ! It will be used in the first step to update velocity
    !
    !TODO: initialize plumed
    !
    !--------------------------------
    ! MD LOOP
    !--------------------------------
    do istep = 1, MCTRL_MDstep
        if (inode_tot .eq. 1) then
            write (*, *) "MD step=", istep
        end if

        tt1 = mpi_wtime()

        call update_time(md)
 
        tt = mpi_wtime()
 
        call update_T(md, istep)

        !--------------------------------
        ! move atom
        ! Note, the velocity is also updated in update_r, it update the
        ! second half with a give mymd.f
        !--------------------------------
        call update_r(md, istep)
        !
        !--------------------------------
        ! rho,ug interpolation
        !--------------------------------
        !  if(MCTRL_iMD.ne.101) then
        ! call interpolation()
        ! endif
        !
        !--------------------------------
        ! new force and Etot
        !--------------------------------
        call exchange_data_md2scf(md)

        call MPI_Bcast(MCTRL_xatom, MCTRL_natom*3, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(md%v, MCTRL_natom*3, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

        if (switch) then
            combine_FF_total_start = mpi_wtime()
        end if

        if (MCTRL_stress .eq. 1 .and. mod(istep - 1, MCTRL_stress_step) .eq. 0) then
            AL_old = MCTRL_AL
            AL(1, :) = AL_old(1, :)*1.01   ! numerical stress calculation
            call combine_FF(Etot_1, fatom, MCTRL_xatom, MCTRL_iatom, AL, &
                MCTRL_natom, natom_m, nmol)
            AL(2, :) = AL_old(2, :)*1.01   ! numerical stress calculation
            call combine_FF(Etot_2, fatom, MCTRL_xatom, MCTRL_iatom, AL, &
                MCTRL_natom, natom_m, nmol)
            AL(3, :) = AL_old(3, :)*1.01   ! numerical stress calculation
            call combine_FF(Etot_3, fatom, MCTRL_xatom, MCTRL_iatom, AL, &
                MCTRL_natom, natom_m, nmol)
            AL = AL_old
        end if

        call combine_FF(Etot, fatom, MCTRL_xatom, MCTRL_iatom, MCTRL_AL, &
            MCTRL_natom, natom_m, nmol)

        if (switch) then
            combine_FF_total_end = mpi_wtime()
        end if

        if (MCTRL_stress .eq. 1) then
            if (mod(istep - 1, MCTRL_stress_step) .eq. 0) then
                e_stress = 0.d0
                e_stress(1, 1) = (Etot_1 - Etot)/0.01
                e_stress(2, 2) = (Etot_2 - Etot)/0.01
                e_stress(3, 3) = (Etot_3 - Etot)/0.01
                e_stress_old = e_stress
            else
                e_stress = e_stress_old
            end if

        else
            e_stress = 0.d0
        end if

        call exchange_data_scf2md(md, fatom, Etot, e_stress)

        !--------------------------------
        ! new velocity
        ! Note, the velocity is also updated in update_r, it update half,
        ! and here, it update another half
        !--------------------------------
        call update_v(md)
        !--------------------------------
        ! post processing
        !--------------------------------
        call get_energy_kinetic(md)
        call get_temperature(md)
        ! the scaling must follow get_energy_kinetic & get_temperature
        call energy_scaling(md, istep)
        
        ! if ( mod(istep,1000) .eq. 0) then
        call zero_momentum(md, istep)
        ! endif

        !call get_diffusion_coeff(md)
        !call get_cell_info(md)
        call get_average_temperature(md)
        call get_average_pressure(md)
        call post_check(md)
        !--------------------------------
        ! output
        !--------------------------------
        call write_MDSTEPS(md)
        if (mod(istep - 1, MCTRL_output_nstep) .eq. 0) then
            call write_MOVEMENT(md, istep, MCTRL_MDstep)
        end if
        call write_finalconfig(md)
        !if(.true.) then
        !    call write_diffusion_coeff(md)
        !    call write_average_temperature(md)
        !endif
        !
        tt = mpi_wtime() - tt

        tt2 = mpi_wtime()
!         if (inode_tot .eq. 1) then
!             write (6, *) "time MDstep", tt2 - tt1
!             if (switch) then
!                 write (*, "(A60)") "---function name-----------spend time(s)------ratio------------"
!                 write(*, 200) "   combine_FF_total        ", combine_FF_total_end - combine_FF_total_start, (combine_FF_total_end - combine_FF_total_start) / (tt2 - tt1)
!                 write (*, "(A60)") "——---------------------subroutine:combine_FF----------------------------"
!                 write (*, 300) "   combine_FF time         ", combine_FF_end - combine_FF_start
!                 write (*, "(A60)") "---function name-----------spend time(s)------ratio------------"
!                 write(*, 200) "   ML_FF_EF                ", ML_FF_EF_end - ML_FF_EF_start, (ML_FF_EF_end - ML_FF_EF_start) / (combine_FF_end - combine_FF_start)
!                 write(*, 200) "   nonbond_force           ", nonbond_force_end - nonbond_force_start, (nonbond_force_end - nonbond_force_start) / (combine_FF_end - combine_FF_start)
!                 write (*, "(A60)") "----------------------subroutine:nonbond_force--------------------------"
!                 write (*, 300) "   nonbond_force time", subroutine_nonbond_force_end - subroutine_nonbond_force_start
!                 write (*, "(A60)") "---function name-----------spend time(s)------ratio------------"
!                 write(*, 200) "   loop_5001               ", loop_5001_end - loop_5001_start, (loop_5001_end - loop_5001_start) / (subroutine_nonbond_force_end - subroutine_nonbond_force_start)
!                 write(*, 200) "   loop_6001               ", loop_6001_end - loop_6001_start, (loop_6001_end - loop_6001_start) / (subroutine_nonbond_force_end - subroutine_nonbond_force_start)
!                 write(*, 200) "   loop_400                ", loop_400_end - loop_400_start, (loop_400_end - loop_400_start) / (subroutine_nonbond_force_end - subroutine_nonbond_force_start)
!                 write(*, 200) "   loop_500                ", loop_500_end - loop_500_start, (loop_500_end - loop_500_start) / (subroutine_nonbond_force_end - subroutine_nonbond_force_start)
!                 write(*, 200) "   loop_polarization       ", loop_polarization_end - loop_polarization_start, (loop_polarization_end - loop_polarization_start) / (subroutine_nonbond_force_end - subroutine_nonbond_force_start)
!                 write(*, 200) "   calc_polar_F            ", calc_polar_F_end - calc_polar_F_start, (calc_polar_F_end - calc_polar_F_start) / (subroutine_nonbond_force_end - subroutine_nonbond_force_start)
!                 write(*, 200) "   calc_coul_xc_kin_V      ", calc_coul_xc_kin_V_total_time, (calc_coul_xc_kin_V_total_time) / (subroutine_nonbond_force_end - subroutine_nonbond_force_start)
!                 write(*, 200) "   ML_FF_PPP               ", ML_FF_PPP_total_time, (ML_FF_PPP_total_time) / (subroutine_nonbond_force_end - subroutine_nonbond_force_start)
! 200             FORMAT('', A25, F15.10, F15.10)
! 300             FORMAT('', A25, F15.10)
!             end if
!         end if


        if (inode_tot .eq. 1) then
            write (6, *) "time MDstep", tt2 - tt1
            if (switch) then
                write (*, "(A60)") "——---------------------subroutine:combine_FF----------------------------"
                write (*, 300) "   combine_FF time         ", combine_FF_end - combine_FF_start
                write (*, "(A60)") "---function name-----------spend time(s)------ratio------------"
                write(*, 200) "   ML_FF_EF                ", ML_FF_EF_end - ML_FF_EF_start, (ML_FF_EF_end - ML_FF_EF_start) / (combine_FF_end - combine_FF_start)
                write(*, 200) "   nonbond_force           ", nonbond_force_end - nonbond_force_start, (nonbond_force_end - nonbond_force_start) / (combine_FF_end - combine_FF_start)
                write (*, "(A60)") "---function name-----------spend time(s)------ratio------------"
                write(*, 200) "   nonbond_polar           ", loop_polarization_total_time, loop_polarization_total_time / (subroutine_nonbond_force_end - subroutine_nonbond_force_start)
                write(*, 200) "   calc_polar              ", calc_polar_F_end - calc_polar_F_start, (calc_polar_F_end - calc_polar_F_start) / (subroutine_nonbond_force_end - subroutine_nonbond_force_start)
                write(*, 200) "   calc_coul_xc_kin_V      ", calc_coul_xc_kin_V_total_time, (calc_coul_xc_kin_V_total_time) / (subroutine_nonbond_force_end - subroutine_nonbond_force_start)
                write(*, 200) "   ML_FF_PPP               ", ML_FF_PPP_total_time, (ML_FF_PPP_total_time) / (subroutine_nonbond_force_end - subroutine_nonbond_force_start)
200             FORMAT('', A25, F15.10, F15.10)
300             FORMAT('', A25, F15.10)
            end if
        end if


    end do

    return
end subroutine molecular_dynamics_kernel

