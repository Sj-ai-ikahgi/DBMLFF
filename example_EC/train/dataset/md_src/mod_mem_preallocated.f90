module mod_mem_preallocated
use cudafor
! memory in GPU to be preallocated
integer imax_ncent_itype_mol, id1_max, id2_max, id3_max
type(c_devptr)  c_devptr_dxyz_box, c_devptr_rho, c_devptr_rho_z, c_devptr_vxc, c_devptr_vxc2, c_devptr_vcoul, c_devptr_xatom_cent, c_devptr_AL, &
                c_devptr_ALm, c_devptr_xatom_m, c_devptr_force_m, c_devptr_icorner, c_devptr_force_cent, c_devptr_box, c_devptr_box2, c_devptr_box3, &
                c_devptr_dbox, c_devptr_dbox2, c_devptr_dbox3, c_devptr_dxyzp_store, c_devptr_force_polar_m, c_devptr_ppp, c_devptr_dppp, c_devptr_pxyz, &
                c_devptr_dpxyz, c_devptr_rho_m, c_devptr_rho_mz, c_devptr_vxc_m, c_devptr_vxc2_m, c_devptr_vcoul_m, c_devptr_dbox_c, c_devptr_dbox3_c
type(c_devptr)  c_devptr_fact2_store, c_devptr_fact_store, c_devptr_CC_pol, c_devptr_xc_cent, c_devptr_rhop_tot, c_devptr_vrhop_tot_coul, c_devptr_dvxc_tot, &
                c_devptr_dvxc_m, c_devptr_rhop_m, c_devptr_vrhop_m_coul

real*8, allocatable,dimension(:,:,:,:) :: CC_pol
real*8, allocatable,dimension(:,:,:) :: dxyzp_store
real*8, allocatable,dimension(:,:,:,:) :: force_polar_m
real*8, allocatable,dimension(:,:) :: fact_polar
real*8, allocatable,dimension(:,:) ::  pxyz,ppp
real*8, allocatable,dimension(:,:,:,:) ::  dppp
real*8, allocatable,dimension(:,:,:) ::  dpxyz

end module