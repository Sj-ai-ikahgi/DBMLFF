module mod_cuinterface
    use iso_c_binding
    interface

    subroutine cudakernel_nonbond_loop_82to96(dxyz_box, AL, id3, id2, id1, n3, n2, n1) bind(c, name = 'fcc_nonbond_loop_82to96')
        use iso_c_binding 
        use cudafor
        type(c_devptr), value :: dxyz_box, AL
        integer(c_int), value :: id3, id2, id1, n3, n2, n1
    end subroutine

    subroutine cudakernel_nonbond_loop_134to149(icent, w_cent, nat_cent, xatom_m, itype_mol, imol, xatom_cent, ncent_itype_mol, natom_mm, nmolm) bind(c, name = 'fcc_nonbond_loop_134to149')
        use iso_c_binding   
        use cudafor
        type(c_devptr), value :: icent, w_cent, nat_cent, xatom_m, xatom_cent
        integer(c_int), value :: itype_mol, imol, ncent_itype_mol, natom_mm, nmolm
    end subroutine

    subroutine cudakernel_nonbond_loop_149to238(xatom_cent, z_cent, itype_cent, funcr2, AL, ncent_itype_mol, nr, Rm2, dxyz_box, &
                                        rho, rho_z, id3, id2, id1, n3, n2, n1, Rbox2, itype_mol, pi, ion_type_cent, imax_ion, &
                                        r_ion, rho_ion, vol, box, box2, box3, Q_type, z_ion, imax_nr, imax_ntype_cent) bind(c, name = 'fcc_nonbond_loop_149to238')
        use iso_c_binding
        use cudafor
        type(c_devptr), value :: xatom_cent, z_cent, itype_cent, funcr2, AL, dxyz_box, rho, rho_z, r_ion, rho_ion, box, box2, box3, Q_type, z_ion, imax_ion, ion_type_cent
        integer(c_int), value :: id3, id2, id1, n3, n2, n1, nr, itype_mol, ncent_itype_mol, imax_nr, imax_ntype_cent
        real(c_double), value :: Rm2, Rbox2, pi, vol
    end subroutine

    subroutine cudakernel_nonbond_loop_400(itype_cent, xatom_cent, z_cent, icorner, dbox_c, dbox3_c, dxyz_box, AL, rho_m, rho_mz, funcr2, nr, Rm2, &
                                    id3, id2, id1, nm3, nm2, nm1, n3, n2, n1, Rbox2, itype_mol, ncent_itype_mol, pi, ion_type_cent, imax_ion, &
                                    r_ion, rho_ion, box, box2, box3, vol, Q_type, z_ion, dbox, dbox2, dbox3, imax_nr, imax_ntype_cent) bind(c, name = 'fcc_nonbond_loop_400')
        use iso_c_binding
        use cudafor
        type(c_devptr), value :: itype_cent, xatom_cent, z_cent, icorner, dxyz_box, AL, funcr2, dbox_c, dbox3_c, rho_m, rho_mz, box, box2, box3, Q_type, z_ion, &
                                imax_ion, ion_type_cent, r_ion, rho_ion, dbox, dbox2, dbox3
        integer(c_int), value :: nr, id3, id2, id1, n3, n2, n1, nm3, nm2, nm1, itype_mol, ncent_itype_mol, imax_nr, imax_ntype_cent
        real(c_double), value :: Rm2, Rbox2, pi, vol
    end subroutine

    subroutine cudakernel_nonbond_loop_500(itype_cent, xatom_cent, force_cent, icorner, vxc, vxcm, vcoul, vcoulm, dbox_c, dbox3_c, vol_n, &
                                        id3, id2, id1, n3, n2, n1, nm3, nm2, nm1, itype_mol, ncent_itype_mol) bind(c, name = 'fcc_nonbond_loop_500')
        use iso_c_binding
        use cudafor
        type(c_devptr), value :: itype_cent, xatom_cent, vxc, vcoul, vxcm, vcoulm, icorner, dbox_c, dbox3_c, force_cent
        integer(c_int), value :: id3, id2, id1, n3, n2, n1, nm3, nm2, nm1, itype_mol, ncent_itype_mol
        real(c_double), value :: vol_n
    end subroutine

    subroutine cudakernel_nonbond_loop_510(natom_m_itype_mol, itype_mol, ion_type_atomp, atom_charge_param, xatom_m, n3, n2, n1, nm3, nm2, nm1, icorner, AL, &
                                      vol_n, id3, id2, id1, dxyz_box, Rbox2, pxyz, dpxyz, vxc2, vxc2_m, vcoul, vcoul_m, imol, natom_mm, nmolm) bind(c, name = 'fcc_nonbond_loop_510')
        use iso_c_binding
        use cudafor
        type(c_devptr), value :: ion_type_atomp, atom_charge_param, xatom_m, icorner, AL, dxyz_box, pxyz, dpxyz, vxc2, vxc2_m, vcoul, vcoul_m
        integer(c_int), value :: natom_m_itype_mol, itype_mol, n3, n2, n1, id3, id2, id1, imol, natom_mm, nmolm, nm3, nm2, nm1
        real(c_double), value :: vol_n, Rbox2
    end subroutine

    subroutine cudakernel_calc_coul_xc_kin_V(rho, rho_z, vxc, vxc2, vcoul, n1, n2, n3, AL, E_coul,  &
                                            E_xc, E_kin1, E_kin2, fact_kin2, pi) bind(c, name = 'fcc_calc_coul_xc_kin_V')           
        use iso_c_binding
        use cudafor
        type(c_devptr), value :: rho, rho_z, vxc, vxc2, vcoul, AL, AL_d
        integer(c_int), value :: n1, n2, n3
        real(c_double) :: E_coul, E_xc, E_kin1, E_kin2
        real(c_double), value :: pi, fact_kin2
    end subroutine

    subroutine cudakernel_polar_loop_113to129(dxyz_box, AL, id3, id2, id1, n3, n2, n1) bind (c, name = 'fcc_polar_loop_113to129')
        use iso_c_binding
        use cudafor
        type(c_devptr), value :: dxyz_box, AL
        integer(c_int), value :: id3, id2, id1, n3, n2, n1
    end subroutine

    subroutine cudakernel_polar_loop_500(natom_m_itype_mol, ion_type_atomp, atom_charge_param, xatom_m, imol, itype_mol, n3, n2, n1, id3, id2, id1, natom_mm, nmolm, ntype_mm, AL, dxyz_box, CC_pol, rhop_tot, Rbox2) bind (c, name = 'fcc_polar_loop_500')
        use iso_c_binding
        use cudafor
        type(c_devptr), value :: ion_type_atomp, atom_charge_param, xatom_m, AL, dxyz_box, CC_pol, rhop_tot
        integer(c_int), value :: natom_m_itype_mol, imol, n3 ,n2, n1, id3, id2, id1, natom_mm, nmolm, ntype_mm, itype_mol
        real(c_double), value :: Rbox2
    end subroutine

    subroutine cudakernel_calc_coul(rho, vcoul, n1, n2, n3, AL, pi) bind (c, name = 'fcc_calc_coul')
        use iso_c_binding
        use cudafor
        type(c_devptr), value :: rho, vcoul, AL
        integer(c_int), value :: n1, n2, n3
        real(c_double), value :: pi
    end subroutine

    subroutine cudakernel_calc_dvxc(rho, rho_p, dvxc, n1, n2, n3, pi) bind (c, name = 'fcc_calc_dvxc')
        use iso_c_binding
        use cudafor
        type(c_devptr), value :: rho, rho_p, dvxc
        integer(c_int), value :: n1, n2, n3
        real(c_double), value :: pi
    end subroutine

    subroutine cudakernel_polar_loop_312to339(ncent_itype_mol, icent, imol, itype_mol, w_cent, xc_cent, xatom_cent, icorner, n3, n2, n1, &
                                            nm3, nm2, nm1, xatom_m, natom_mm, nmolm, nat_cent) bind (c, name = 'fcc_polar_loop_312to339')
        use iso_c_binding
        use cudafor
        type(c_devptr), value :: icent, w_cent, xc_cent, xatom_cent, icorner, xatom_m, nat_cent
        integer(c_int), value :: ncent_itype_mol, imol, itype_mol, n3, n2, n1, nm3, nm2, nm1, natom_mm, nmolm
    end subroutine

    subroutine cudakernel_polar_loop_600(ncent_itype_mol, itype_cent, itype_mol, xatom_cent, icorner, AL, ion_type_cent, imax_ion, r_ion, id3, &
                                        id2, id1, box, box2, box3, dxyz_box, funcr2, fact_store, fact2_store, n3, n2, n1, nr, Rm2, Rbox2, &
                                        rho_ion, vol, nm3, nm2, nm1, rho_m, Q_type, z_ion, imax_nr, imax_ntype_cent) bind (c, name = 'fcc_polar_loop_600')
        use iso_c_binding
        use cudafor
        type(c_devptr), value :: itype_cent, xatom_cent, icorner, AL, ion_type_cent, imax_ion, r_ion, box, box2, box3, dxyz_box, funcr2, fact_store, &
                                fact2_store, rho_ion, rho_m, z_ion, Q_type
        integer(c_int), value :: ncent_itype_mol, itype_mol, id3, id2, id1, n3, n2, n1, nr, nm3, nm2, nm1, imax_nr, imax_ntype_cent
        real(c_double), value :: Rm2, Rbox2, vol
    end subroutine

    subroutine cudakernel_polar_loop_700(natom_m_itype_mol, ion_type_atomp, imol, itype_mol, atom_charge_param, xatom_m, &
                                        icorner, n3, n2, n1, id3, id2, id1, CC_pol, rhop_m, dxyz_box, Rbox2, natom_mm, nmolm, &
                                        AL, nm3, nm2, nm1) bind (c, name = 'fcc_polar_loop700')
        use iso_c_binding
        use cudafor        
        type(c_devptr), value :: ion_type_atomp, atom_charge_param, xatom_m, icorner, CC_pol, rhop_m, dxyz_box, AL
        integer(c_int), value :: natom_m_itype_mol, imol, itype_mol, n3, n2, n1, id3, id2, id1, natom_mm, nmolm, nm3, nm2, nm1
        real(c_double), value :: Rbox2
    end subroutine

    subroutine cudakernel_polar_loop_800(ncent_itype_mol, itype_cent, xatom_cent, icorner, AL, n3, n2, n1, ion_type_cent, imax_ion, &
                                    r_ion, dxyz_box, funcr2, box, box2, box3, dbox, dbox2, dbox3, fact_store, fact2_store, Rbox2, nr, &
                                    nm3, nm2, nm1, vrhop_tot_coul, vrhop_m_coul, dvxc_tot, dvxc_m, vol_n, id3, id2, id1, force_cent, Rm2, &
                                    itype_mol, rho_ion, imax_nr, imax_ntype_cent) bind (c, name = 'fcc_polar_loop800')
        use iso_c_binding
        use cudafor        
        type(c_devptr), value :: itype_cent, xatom_cent, icorner, AL, ion_type_cent, imax_ion, r_ion, dxyz_box, funcr2, box, box2, box3, &
                                dbox, dbox2, dbox3, fact_store,fact2_store, vrhop_tot_coul, vrhop_m_coul, dvxc_tot, dvxc_m, rho_ion, force_cent
        integer(c_int), value :: ncent_itype_mol, n3, n2, n1, nr, nm3, nm2, nm1, id3, id2, id1, itype_mol, imax_nr, imax_ntype_cent
        real(c_double), value :: Rbox2, vol_n, Rm2
    end subroutine

    subroutine init_resources(h_bytes, d_bytes, n1, n2, n3, nm1_all, nm2_all, nm3_all, ntype_m) bind (c, name = 'init_resources')
        use iso_c_binding
        integer(c_int), value :: n1, n2, n3, ntype_m
        integer(c_int64_t), value :: h_bytes, d_bytes
        integer, dimension(1 : 1000) :: nm1_all, nm2_all, nm3_all
    end subroutine
    
    subroutine destroy_resources() bind (c, name = 'destroy_resources')
        use iso_c_binding
    end subroutine

    subroutine zero_gpu_mem(ptr_device, len) bind (c, name = 'zero_gpu_mem')
        use iso_c_binding
        use cudafor
        integer(c_int), value :: len
        type(c_devptr), value :: ptr_device
    end subroutine
    subroutine print_gpu_array(ptr_device, len) bind (c, name = 'print_gpu_array')
        use iso_c_binding
        use cudafor
        integer(c_int), value :: len
        type(c_devptr), value :: ptr_device
    end subroutine

    subroutine nccl_mpi_allreduce(sendbuff, recvbuff, count, datatype, op) bind (c, name = 'nccl_mpi_allreduce')
        use iso_c_binding
        use cudafor
        type(c_devptr), value :: sendbuff, recvbuff
        integer(c_int), value :: count, datatype, op
    end subroutine

    subroutine cudakernel_vec_add(src1, src2, dst, count) bind (c, name = 'fcc_vec_add')
        use iso_c_binding
        use cudafor
        type(c_devptr), value :: src1, src2, dst
        integer(c_int), value :: count
    end subroutine

    end interface
end module mod_cuinterface
