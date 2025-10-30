module mod_profile
implicit none 
    logical :: switch = .true.
    ! variables for md's profile
    real*8 combine_FF_start, combine_FF_end
    ! variables for combine_FF's profile
    real*8 combine_FF_total_start, combine_FF_total_end, &
    ML_FF_EF_start, ML_FF_EF_end, &
    nonbond_force_start, nonbond_force_end
    ! variables for nonbond_force's profile
    real*8 subroutine_nonbond_force_start, subroutine_nonbond_force_end, &
    calc_coul_xc_kin_V_00_start, calc_coul_xc_kin_V_00_end, &
    calc_coul_xc_kin_V_01_start, calc_coul_xc_kin_V_01_end, &
    calc_coul_xc_kin_V_total_time, &
    loop_5001_start, loop_5001_end, &
    loop_6001_start, loop_6001_end, &
    loop_400_start, loop_400_end, &
    loop_500_start, loop_500_end, &
    loop_510_start, loop_510_end, &
    loop_polarization_start, loop_polarization_end, &
    calc_polar_F_start, calc_polar_F_end
    ! variables for calc_polar's profile
    real*8 subroutine_calc_polar_loop500_start, subroutine_calc_polar_loop500_end, &
    subroutine_calc_polar_loop600_start, subroutine_calc_polar_loop600_end, &
    subroutine_calc_polar_loop700_start, subroutine_calc_polar_loop700_end, &
    subroutine_calc_polar_loop800_start, subroutine_calc_polar_loop800_end
    real*8 ML_FF_PPP_start, ML_FF_PPP_end, ML_FF_PPP_total_time
    real*8 nonbond_others_start, nonbond_others_end, nonbond_others_total_time
    real*8 loop_polarization_total_time


contains
    subroutine init_var_nonbond()
        calc_coul_xc_kin_V_total_time = 0
        ML_FF_PPP_total_time = 0
        nonbond_others_total_time = 0
        loop_polarization_total_time = 0
    end subroutine init_var_nonbond
end module mod_profile
