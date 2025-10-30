subroutine nonbond_force(natom_m, nmol, xatom_m, AL, force_m, Etot)

    !  first assume all the molecules are in the same type
    !  There is no molecule type right now
    !  And all the molecules have the same order and type etc
    use mod_mpi
    use mod_param_densityFF
    use mod_data, only: &
        natom_mm, nmolm, ntype_m, ntype_mm, Etot_nonbond, E_polar, E_kin2_nonbond, &
        iflag_training, iflag_debug, polar_factor
    use mod_profile
    use mod_cuinterface
    use cudafor
    use mod_mem_preallocated
    use mod_md, only: A_AU_1, Hartree_eV
    use mod_calc_pol_dir, only: dxpp,dxpp1diat2=>dxpp1_diat2,dxpp2diat4=>dxpp2_diat4, &
        dxpp3diat2=>dxpp3_diat2,dxpp3diat4=>dxpp3_diat4,calc_pol_dir

    implicit double precision(a - h, o - z)

    real*8 AL(3, 3), vol, vol_n
    real*8, device, allocatable :: AL_d(:, :)
    real*8 ALm(3, 3)
    real*8, device, allocatable :: ALm_d(:, :)
    real*8 xc_cent(3)
    integer icorner(3)
    integer, device, allocatable :: icorner_d(:)
    real*8 Etot
    real*8 pi

    integer natom_m(ntype_m), nmol(ntype_m)
    real*8 xatom_m(3, natom_mm, nmolm, ntype_m)      ! NEED TO ALLOCATE
    real*8, device, allocatable :: xatom_m_d(:, :, :, :)
    real*8 force_m(3, natom_mm, nmolm, ntype_m)
    real*8, device, allocatable :: force_m_d(:, :, :, :)
    integer iatom_mm(natom_mm, nmolm, ntype_m)

    integer itype_ion

    real*8 force_m_t(3, natom_mm, nmolm, ntype_m)

    real*8 xatom_cent(3, 10*natom_mm)
    real*8, device, allocatable :: xatom_cent_d(:, :)
    real*8 force_cent(3, 10*natom_mm)
    real*8, device, allocatable :: force_cent_d(:, :)

    real*8, allocatable, dimension(:, :, :, :) :: dxyz_box
    real*8, allocatable, device, dimension(:, :, :, :) :: dxyz_box_d
    real*8, allocatable, dimension(:, :, :, :) :: dbox, dbox2, dbox3
    real*8, allocatable, dimension(:, :, :, :, :) :: dbox_c, dbox3_c
    real*8, allocatable, device, dimension(:, :, :, :, :) :: dbox_c_d, dbox3_c_d
    real*8, allocatable, dimension(:, :, :, :) :: box_c, box2_c, box3_c
    real*8, allocatable, device, dimension(:, :, :, :) :: box_c_d, box2_c_d, box3_c_d
    real*8, allocatable, dimension(:, :, :) :: rho, rho_z, rho_t
    real*8, allocatable, device, dimension(:, :, :) :: rho_d, rho_z_d
    real*8, allocatable, dimension(:, :, :) :: rho_m, rho_mz
    real*8, allocatable, device, dimension(:, :, :) :: rho_m_d, rho_mz_d
    real*8, allocatable, dimension(:, :, :) :: vxc, vcoul, vxc2
    real*8, allocatable, device, dimension(:, :, :) :: vxc_d, vcoul_d, vxc2_d
    real*8, allocatable, dimension(:, :, :) :: vxc_m, vcoul_m, vxc2_m
    real*8, allocatable, device, dimension(:, :, :) :: vxc_m_d, vcoul_m_d, vxc2_m_d

    real*8, allocatable, device, dimension(:, :, :, :) :: CC_pol_d
    real*8, allocatable, device, dimension(:, :, :) :: dxyzp_store_d
    real*8, allocatable, device, dimension(:, :, :, :) :: force_polar_m_d

    real*8, allocatable, device, dimension(:, :) :: pxyz_d, ppp_d
    real*8, allocatable, device, dimension(:, :, :, :) :: dppp_d
    real*8, allocatable, device, dimension(:, :, :) :: dpxyz_d

    real*8 dxp1(3), dxp2(3), dxp3(3)
    real*8 dxpp1(3), dxpp2(3), dxpp3(3)
    real*8 dxpp3_diat2(3, 3), dxpp3_diat4(3, 3)
    real*8 dxpp1_diat2(3, 3), dxpp2_diat4(3, 3)
    real*8 dxt(3), dwx1(3), dwx2(3), px1(3), px2(3)
    real*8 dpx1(3, 3), dpx2(3, 3)

    real*8 ppp1(2, 200), ppp2(2, 200)
    real*8 dppp1(3, 2, 200), dppp2(3, 2, 200)
    real*8 Epp(3, 200), dEpp(3, 3, 200)
    real*8 force_Epp(3, 200)

    integer iflag_mol(1000)
    integer ierr

    real*8 AL_tmp(3,3)
    real*8 ppp1_m(2,200),ppp2_m(2,200)
    real*8,allocatable,dimension(:,:) :: ppp_m
    real*8,allocatable,dimension(:,:,:,:) :: dppp_m
    real*8 dppp1_m(3,200,2,200),dppp2_m(3,200,2,200)

    integer nzdp1, nzdp2, nzdp

    call init_var_nonbond()

    if (switch) subroutine_nonbond_force_start = mpi_wtime()

    AL_tmp = AL*A_AU_1

    vol = al(3, 1)*(al(1, 2)*al(2, 3) - al(1, 3)*al(2, 2)) &
        + al(3, 2)*(al(1, 3)*al(2, 1) - al(1, 1)*al(2, 3)) &
        + al(3, 3)*(al(1, 1)*al(2, 2) - al(1, 2)*al(2, 1))

    vol = abs(vol)
    vol_n = vol/(n1*n2*n3)

    pi = 4*datan(1.d0)

    aL1 = dsqrt(AL(1, 1)**2 + AL(2, 1)**2 + AL(3, 1)**2)
    aL2 = dsqrt(AL(1, 2)**2 + AL(2, 2)**2 + AL(3, 2)**2)
    aL3 = dsqrt(AL(1, 3)**2 + AL(2, 3)**2 + AL(3, 3)**2)
    id1 = (Rbox/aL1)*n1
    id2 = (Rbox/aL2)*n2
    id3 = (Rbox/aL3)*n3
    Rbox2 = Rbox**2

    call c_f_pointer(c_devptr_dxyz_box, dxyz_box_d, (/3, -id1:id1, -id2:id2, -id3:id3/))
    call c_f_pointer(c_devptr_rho, rho_d, (/n1, n2, n3/))
    call c_f_pointer(c_devptr_rho_z, rho_z_d, (/n1, n2, n3/))
    call c_f_pointer(c_devptr_vxc, vxc_d, (/n1, n2, n3/))
    call c_f_pointer(c_devptr_vxc2, vxc2_d, (/n1, n2, n3/))
    call c_f_pointer(c_devptr_vcoul, vcoul_d, (/n1, n2, n3/))
    call c_f_pointer(c_devptr_xatom_cent, xatom_cent_d, (/3, 10*natom_mm/))
    call c_f_pointer(c_devptr_AL, AL_d, (/3, 3/))
    call c_f_pointer(c_devptr_ALm, ALm_d, (/3, 3/))
    call c_f_pointer(c_devptr_icorner, icorner_d, (/3/))
    call c_f_pointer(c_devptr_xatom_m, xatom_m_d, (/3, natom_mm, nmolm, ntype_m/))
    call c_f_pointer(c_devptr_force_m, force_m_d, (/3, natom_mm, nmolm, ntype_m/))
    call c_f_pointer(c_devptr_force_cent, force_cent_d, (/3, 10*natom_mm/))
    ! call c_f_pointer(c_devptr_box3, box3_c_d, (/-id1:id1, -id2:id2, -id3:id3, imax_ncent_itype_mol/))

    AL_d = AL
    xatom_m_d = xatom_m

    call cudakernel_nonbond_loop_82to96(c_devptr_dxyz_box, c_devptr_AL, id3, id2, id1, n3, n2, n1)

    !cccccccccccccccc  first, get the whole charge density
    E_polar = 0.d0

    call zero_gpu_mem(c_devptr_rho, n1*n2*n3*8)
    call zero_gpu_mem(c_devptr_rho_z, n1*n2*n3*8)

    kkk_mol = 0
    iflag_mol(1:ntype_m) = 0

    if (switch) loop_5001_start = mpi_wtime()

    do 5001 itype_mol = 1, ntype_m

        do 5000 imol = 1, nmol(itype_mol)
            kkk_mol = kkk_mol + 1

            if (mod(kkk_mol - 1, nnodes_tot) .ne. inode_tot - 1) goto 5000

            nm1 = nm1_all(itype_mol)
            nm2 = nm2_all(itype_mol)
            nm3 = nm3_all(itype_mol)

            call cudakernel_nonbond_loop_134to149(c_devptr_icent, c_devptr_w_cent, c_devptr_nat_cent, c_devptr_xatom_m, itype_mol, imol, c_devptr_xatom_cent, &
                ncent(itype_mol), natom_mm, nmolm)

            call cudakernel_nonbond_loop_149to238(c_devptr_xatom_cent, c_devptr_z_cent, c_devptr_itype_cent, c_devptr_funcr2, c_devptr_AL, ncent(itype_mol), nr, Rm2,  &
                c_devptr_dxyz_box, c_devptr_rho, c_devptr_rho_z, id3, id2, id1, n3, n2, n1, Rbox2, itype_mol, pi, c_devptr_ion_type_cent, &
                c_devptr_imax_ion, c_devptr_r_ion, c_devptr_rho_ion, vol, c_devptr_box, c_devptr_box2, c_devptr_box3, c_devptr_Q_type, &
                c_devptr_z_ion, imax_nr, imax_ntype_cent)

5000    continue  ! imol=1,nmol
5001 continue

     if (switch) loop_5001_end = mpi_wtime()


     if (switch) nonbond_others_start = mpi_wtime()

     call nccl_mpi_allreduce(c_devptr_rho, c_devptr_rho, n1*n2*n3, 10, 0)
     call nccl_mpi_allreduce(c_devptr_rho_z, c_devptr_rho_z, n1*n2*n3, 10, 0)
     call cudakernel_vec_add(c_devptr_rho_z, c_devptr_rho, c_devptr_rho_z, n1*n2*n3)


     if (switch) then
        nonbond_others_end = mpi_wtime()
        nonbond_others_total_time = nonbond_others_total_time + nonbond_others_end - nonbond_others_start
     endif


     ! Right now, all processors are calling this. In the future,
     ! this should be done with parallel FFT, and this box can be
     ! big. One can even use the small box FFT

     if (switch) calc_coul_xc_kin_V_00_start = mpi_wtime()

     allocate (rho(n1, n2, n3))
     allocate (rho_z(n1, n2, n3))

     call cudakernel_calc_coul_xc_kin_V(c_devptr_rho, c_devptr_rho_z, c_devptr_vxc, c_devptr_vxc2, c_devptr_vcoul, n1, n2, n3, c_devptr_AL, &
         E_coul, E_xc, E_kin1, E_kin2, fact_kin2, pi)

     if (switch) then
         calc_coul_xc_kin_V_00_end = mpi_wtime()
         calc_coul_xc_kin_V_total_time = calc_coul_xc_kin_V_total_time + calc_coul_xc_kin_V_00_end - calc_coul_xc_kin_V_00_start
     end if

     Etot = E_coul + E_xc + E_kin1 + E_kin2*fact_kin2
    
     E_kin2_nonbond = E_kin2
     E_kin2_M = 0.d0

     !ccccccccccccccccccccccccccccccccccccccccccccccccccccc
     !ccccccccccccccccccccccccccccccccccccccccccccccccccccc
     !ccccccccccccccccccccccccccccccccccccccccccccccccccccc

     if (iflag_polar .eq. 1) then
         call c_f_pointer(c_devptr_pxyz, pxyz_d, (/3, 2*natom_mm/))
         call c_f_pointer(c_devptr_dpxyz, dpxyz_d, (/3, 3, 2*natom_mm/))

         allocate(ppp_m(3, 2*natom_mm))
         allocate(dppp_m(3,2*natom_mm,natom_mm,3))

         CC_pol = 0.d0
         force_polar_m = 0.d0
         fact_polar = 1.d0
     end if

     force_m = 0.d0
     EtotM = 0.d0

     if (switch) loop_6001_start = mpi_wtime()

     kkk_mol = 0
     iflag_mol(1:ntype_m) = 0

     do 6001 itype_mol = 1, ntype_m

         nm1 = nm1_all(itype_mol)
         nm2 = nm2_all(itype_mol)
         nm3 = nm3_all(itype_mol)

         ALm(:, 1) = AL(:, 1)*nm1*1.d0/n1
         ALm(:, 2) = AL(:, 2)*nm2*1.d0/n2
         ALm(:, 3) = AL(:, 3)*nm3*1.d0/n3
         ierr = cudaMemcpy(ALm_d, ALm, 9, cudaMemcpyHostToDevice)

         if (nm1 .lt. 2*id1 + 1 .or. nm2 .lt. 2*id2 + 1 .or. nm3 .lt. 2*id3 + 1) then
             write (6, *) "nm1,nm2,nm3 too small", nm1, nm2, nm3
             write (6, *) "id1,id2,id3", id1, id2, id3
             stop
         end if

         do 6000 imol = 1, nmol(itype_mol)
             kkk_mol = kkk_mol + 1

             if (mod(kkk_mol - 1, nnodes_tot) .ne. inode_tot - 1) goto 6000

             xc_cent = 0.d0
             do ii = 1, ncent(itype_mol)
                 do ixyz = 1, 3
                     x1 = xatom_m(ixyz, icent(1, ii, itype_mol), imol, itype_mol)
                     xc = x1*w_cent(1, ii, itype_mol)
                     w_sum = w_cent(1, ii, itype_mol)
                     do jj = 2, nat_cent(ii, itype_mol)
                         x2 = xatom_m(ixyz, icent(jj, ii, itype_mol), imol, itype_mol)
                         if (abs(x2 + 1 - x1) .lt. abs(x2 - x1)) x2 = x2 + 1
                         if (abs(x2 - 1 - x1) .lt. abs(x2 - x1)) x2 = x2 - 1
                         xc = xc + x2*w_cent(jj, ii, itype_mol)
                         w_sum = w_sum + w_cent(jj, ii, itype_mol)
                     end do
                     xc = xc/w_sum
                     xatom_cent(ixyz, ii) = xc
                     dx = xatom_cent(ixyz, ii) - xatom_cent(ixyz, 1)
                     if (abs(dx + 1) .lt. abs(dx)) dx = dx + 1    ! only correct when molecule is small
                     if (abs(dx - 1) .lt. abs(dx)) dx = dx - 1

                     xc_cent(ixyz) = xc_cent(ixyz) + (dx + xatom_cent(ixyz, 1))
                 end do
             end do

             xc_cent = xc_cent/ncent(itype_mol)    ! center of the molecule
             icorner(1) = mod(xc_cent(1) + 1.d0, 1.d0)*n1 - int(0.6 + nm1/2.d0) + 1   ! the corner of the rho in (n1,n2,n3)
             icorner(2) = mod(xc_cent(2) + 1.d0, 1.d0)*n2 - int(0.6 + nm2/2.d0) + 1
             icorner(3) = mod(xc_cent(3) + 1.d0, 1.d0)*n3 - int(0.6 + nm3/2.d0) + 1
             !  Note, the corner can be negative. (1,1,1) is the original corner

             ierr = cudaMemcpy(icorner_d, icorner, 3, cudaMemcpyHostToDevice); 
             ierr = cudaMemcpy(xatom_cent_d, xatom_cent, 3*10*natom_mm, cudaMemcpyHostToDevice); 
            !  call zero_gpu_mem(c_devptr_rho_m, nm1_max*nm2_max*nm3_max*8)
            !  call zero_gpu_mem(c_devptr_rho_mz, nm1_max*nm2_max*nm3_max*8)

             if (switch) loop_400_start = mpi_wtime()

             call cudakernel_nonbond_loop_400(c_devptr_itype_cent, c_devptr_xatom_cent, c_devptr_z_cent, c_devptr_icorner, c_devptr_dbox_c, c_devptr_dbox3_c, c_devptr_dxyz_box,  &
                 c_devptr_AL, c_devptr_rho_m, c_devptr_rho_mz, c_devptr_funcr2, nr, Rm2, id3, id2, id1, nm3, nm2, nm1, n3, n2, n1, Rbox2, itype_mol, ncent(itype_mol), pi, &
                 c_devptr_ion_type_cent, c_devptr_imax_ion, c_devptr_r_ion, c_devptr_rho_ion, c_devptr_box, c_devptr_box2, c_devptr_box3, vol, c_devptr_Q_type, c_devptr_z_ion, &
                 c_devptr_dbox, c_devptr_dbox2, c_devptr_dbox3, imax_nr, imax_ntype_cent)

             if (switch) loop_400_end = mpi_wtime()

            !  if (switch) calc_coul_xc_kin_V_01_start = mpi_wtime()

            !  call cudakernel_calc_coul_xc_kin_V(c_devptr_rho_m, c_devptr_rho_mz, c_devptr_vxc_m, c_devptr_vxc2_m, c_devptr_vcoul_m, nm1, nm2, nm3, c_devptr_ALm, &
            !      E_coul, E_xc, E_kin1, E_kin2, fact_kin2, pi)

            !  if (switch) then
            !      calc_coul_xc_kin_V_01_end = mpi_wtime()
            !      calc_coul_xc_kin_V_total_time = calc_coul_xc_kin_V_total_time + calc_coul_xc_kin_V_01_end - calc_coul_xc_kin_V_01_start
            !  end if

            !  EtotM = EtotM + E_coul + E_xc + (E_kin1 + E_kin2*fact_kin2)
            !  E_kin2_M = E_kin2_M + E_kin2

             !ccccccccccccccccccccccccccccccccccccccccccccccccccc

             if (switch) loop_500_start = mpi_wtime()

             call cudakernel_nonbond_loop_500(c_devptr_itype_cent, c_devptr_xatom_cent, c_devptr_force_cent, c_devptr_icorner, c_devptr_vxc, c_devptr_vxc_m, c_devptr_vcoul,  &
                 c_devptr_vcoul_m, c_devptr_dbox_c, c_devptr_dbox3_c, vol_n, id3, id2, id1, n3, n2, n1, nm3, nm2, nm1, itype_mol, ncent(itype_mol))

             ierr = cudaMemcpy(force_cent, force_cent_d, 3*10*natom_mm, cudaMemcpyDeviceToHost)


             if (switch) loop_500_end = mpi_wtime()


             do ii = 1, ncent(itype_mol)
                 do jj = 1, nat_cent(ii, itype_mol)
                     do ixyz = 1, 3

                         force_m(ixyz, icent(jj, ii, itype_mol), imol, itype_mol) = force_m(ixyz, icent(jj, ii, itype_mol), imol, itype_mol) &
                             + w_cent(jj, ii, itype_mol)*force_cent(ixyz, ii)
                     end do
                 end do
             end do


             !!!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
             !!!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
             !!!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
             !!!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
             if ((iflag_polar .ne. 1) .or. (num_term_pol(itype_mol) .eq. 0)) goto 6010

             if (switch) loop_polarization_start = mpi_wtime()

             if (switch) loop_510_start = mpi_wtime()

             call cudakernel_nonbond_loop_510(natom_m(itype_mol), itype_mol, c_devptr_ion_type_atomp, c_devptr_atom_charge_param, c_devptr_xatom_m, n3, n2, n1, nm3, nm2, nm1, c_devptr_icorner, c_devptr_AL, &
                 vol_n, id3, id2, id1, c_devptr_dxyz_box, Rbox2, c_devptr_pxyz, c_devptr_dpxyz, c_devptr_vxc2, c_devptr_vxc2_m, c_devptr_vcoul, c_devptr_vcoul_m, imol, natom_mm, &
                 nmolm)
             ierr = cudaMemcpy(pxyz, pxyz_d, 3*2*natom_mm, cudaMemcpyDeviceToHost)
             ierr = cudaMemcpy(dpxyz, dpxyz_d, 3*3*2*natom_mm, cudaMemcpyDeviceToHost)

             if (switch) loop_510_end = mpi_wtime()

             !ccccccccccc  Now, convert pxyz to p_dir (the molecule direction defined by ipol

            ! ---------- calc dir and ddir ----------
             call calc_pol_dir(itype_mol,xatom_m(1:3,1:natom_m(itype_mol),imol,itype_mol),AL)

             ppp_m=0.d0
             ppp1_m=0.d0
             ppp2_m=0.d0
             dppp_m=0.d0
             dppp1_m=0.d0
             dppp2_m=0.d0
            !  write(6,*) " @@@ start ML_FF_PPP"

             if (switch) ML_FF_PPP_start = mpi_wtime()

            if (iflag_training .eq. 0) then
                call ML_FF_PPP(itype_mol,xatom_m(:,:,imol,itype_mol),AL_tmp,ppp_m,ppp1_m,ppp2_m,dppp_m,dppp1_m,dppp2_m)
            endif

             if (switch) then 
                ML_FF_PPP_end = mpi_wtime()
                ML_FF_PPP_total_time = ML_FF_PPP_total_time + ML_FF_PPP_end - ML_FF_PPP_start
             endif
            !  write(6,*) " @@@ end ML_FF_PPP"


             ! ccccccccc NEW, BOND
             do jj = 1, numb(itype_mol)
                !  jj1 = indb(1, jj, itype_mol)
                !  jj2 = indb(2, jj, itype_mol)
                !  dx1 = xatom_m(1, jj1, imol, itype_mol) - xatom_m(1, jj2, imol, itype_mol)
                !  dx2 = xatom_m(2, jj1, imol, itype_mol) - xatom_m(2, jj2, imol, itype_mol)
                !  dx3 = xatom_m(3, jj1, imol, itype_mol) - xatom_m(3, jj2, imol, itype_mol)
                !  if (abs(dx1 + 1) .lt. abs(dx1)) dx1 = dx1 + 1
                !  if (abs(dx1 - 1) .lt. abs(dx1)) dx1 = dx1 - 1
                !  if (abs(dx2 + 1) .lt. abs(dx2)) dx2 = dx2 + 1
                !  if (abs(dx2 - 1) .lt. abs(dx2)) dx2 = dx2 - 1
                !  if (abs(dx3 + 1) .lt. abs(dx3)) dx3 = dx3 + 1
                !  if (abs(dx3 - 1) .lt. abs(dx3)) dx3 = dx3 - 1
                !  dx = AL(1, 1)*dx1 + AL(1, 2)*dx2 + AL(1, 3)*dx3
                !  dy = AL(2, 1)*dx1 + AL(2, 2)*dx2 + AL(2, 3)*dx3
                !  dz = AL(3, 1)*dx1 + AL(3, 2)*dx2 + AL(3, 3)*dx3
                !  d = dsqrt(dx**2 + dy**2 + dz**2)
                !  Epp(1, jj) = dx/d
                !  Epp(2, jj) = dy/d
                !  Epp(3, jj) = dz/d

                !  dEpp(1, 1, jj) = 1/d - dx**2/d**3
                !  dEpp(2, 1, jj) = -dx*dy/d**3
                !  dEpp(3, 1, jj) = -dx*dz/d**3
                !  dEpp(1, 2, jj) = -dy*dx/d**3
                !  dEpp(2, 2, jj) = 1/d - dy**2/d**3
                !  dEpp(3, 2, jj) = -dy*dz/d**3
                !  dEpp(1, 3, jj) = -dz*dx/d**3
                !  dEpp(2, 3, jj) = -dz*dy/d**3
                !  dEpp(3, 3, jj) = 1/d - dz**2/d**3
                Epp(:,jj)=-dxpp(:,1,jj,2)
                dEpp(:,:,jj)=dxpp1diat2(:,:,jj,2)
             end do

             ppp1 = 0.d0
             ppp2 = 0.d0
             dppp1 = 0.d0
             dppp2 = 0.d0
             Ep_bond = 0.d0
             force_Epp = 0.d0
             do jj = 1, numb(itype_mol)
                 jj1 = indb(1, jj, itype_mol)
                 jj2 = indb(2, jj, itype_mol)
                 jj11 = (jj1 - 1)*2
                 jj22 = (jj2 - 1)*2

                 do ii1 = 1, 2
                     ij1 = ii1 + jj11
                     ij2 = ii1 + jj22
                     ppp1(ii1, jj) = pxyz(1, ij1)*Epp(1, jj) + pxyz(2, ij1)*Epp(2, jj) + pxyz(3, ij1)*Epp(3, jj)
                     ppp2(ii1, jj) = pxyz(1, ij2)*Epp(1, jj) + pxyz(2, ij2)*Epp(2, jj) + pxyz(3, ij2)*Epp(3, jj)
                     dppp1(:, ii1, jj) = pxyz(1, ij1)*dEpp(:, 1, jj) + pxyz(2, ij1)*dEpp(:, 2, jj) + pxyz(3, ij1)*dEpp(:, 3, jj)
                     dppp2(:, ii1, jj) = pxyz(1, ij2)*dEpp(:, 1, jj) + pxyz(2, ij2)*dEpp(:, 2, jj) + pxyz(3, ij2)*dEpp(:, 3, jj)
                 end do

                 ip = (indb(3, jj, itype_mol) - 1)*4
                 do ii1 = 1, 2
                     do ii2 = 1, 2
                         ij1 = ii1 + jj11
                         ij2 = ii2 + jj22
                         ip = ip + 1

                         Ep_bond=Ep_bond+(ppp1(ii1,jj)-ppp1_m(ii1,jj))*(ppp2(ii2,jj)-ppp2_m(ii2,jj))*BB_direct(ip,itype_mol)
                         force_Epp(:,jj1)=force_Epp(:,jj1)+(dppp1(:,ii1,jj)*(ppp2(ii2,jj)-ppp2_m(ii2,jj))+(ppp1(ii1,jj)-ppp1_m(ii1,jj))*dppp2(:,ii2,jj))*BB_direct(ip,itype_mol)
                         force_Epp(:,jj2)=force_Epp(:,jj2)-(dppp1(:,ii1,jj)*(ppp2(ii2,jj)-ppp2_m(ii2,jj))+(ppp1(ii1,jj)-ppp1_m(ii1,jj))*dppp2(:,ii2,jj))*BB_direct(ip,itype_mol)

                         !  The force is defined as dE/dR
                         force_Epp(:,jj1)=force_Epp(:,jj1)+(dpxyz(1,:,ij1)*Epp(1,jj)+dpxyz(2,:,ij1)*Epp(2,jj)+dpxyz(3,:,ij1)*Epp(3,jj)) &
                             *(ppp2(ii2,jj)-ppp2_m(ii2,jj))*BB_direct(ip,itype_mol)
                         force_Epp(:,jj2)=force_Epp(:,jj2)+(ppp1(ii1,jj)-ppp1_m(ii1,jj))* &
                             (dpxyz(1,:,ij2)*Epp(1,jj)+dpxyz(2,:,ij2)*Epp(2,jj)+dpxyz(3,:,ij2)*Epp(3,jj))*BB_direct(ip,itype_mol)
                         do iat=1,natom_m(itype_mol)
                             force_polar_m(:,iat,imol,itype_mol)=force_polar_m(:,iat,imol,itype_mol)+ &
                                 ((-dppp1_m(:,iat,ii1,jj))*(ppp2(ii2,jj)-ppp2_m(ii2,jj))+(ppp1(ii1,jj)-ppp1_m(ii1,jj))*(-dppp2_m(:,iat,ii2,jj)))*BB_direct(ip,itype_mol)
                         enddo
                        !  Ep_bond = Ep_bond + ppp1(ii1, jj)*ppp2(ii2, jj)*BB_direct(ip, itype_mol)
                        !  force_Epp(:, jj1) = force_Epp(:, jj1) + (dppp1(:, ii1, jj)*ppp2(ii2, jj) + ppp1(ii1, jj)*dppp2(:, ii2, jj))*BB_direct(ip, itype_mol)
                        !  force_Epp(:, jj2) = force_Epp(:, jj2) - (dppp1(:, ii1, jj)*ppp2(ii2, jj) + ppp1(ii1, jj)*dppp2(:, ii2, jj))*BB_direct(ip, itype_mol)
                        !  ! The force is defined as dE/dR
                        !  force_Epp(:, jj1) = force_Epp(:, jj1) + (dpxyz(1, :, ij1)*Epp(1, jj) + dpxyz(2, :, ij1)*Epp(2, jj) + dpxyz(3, :, ij1)*Epp(3, jj)) &
                        !      *ppp2(ii2, jj)*BB_direct(ip, itype_mol)
                        !  force_Epp(:, jj2) = force_Epp(:, jj2) + ppp1(ii1, jj)* &
                        !      (dpxyz(1, :, ij2)*Epp(1, jj) + dpxyz(2, :, ij2)*Epp(2, jj) + dpxyz(3, :, ij2)*Epp(3, jj))*BB_direct(ip, itype_mol)

                     end do
                 end do
             end do
             !cccccccccccccccccccccccccccccccccccccccccccccccccccccc


            !  ! ---------------------------------------------------------------
            !  ! ---------------------------------------------------------------
            !  ! ---------------------------------------------------------------

            !  ! output dppp -- bond part
            !  if ((inode_tot.eq.1) .and. (iflag_training .eq. 1)) then

            !     nzdp1 = 0
            !     nzdp2 = 0
            !     nzdp = 0

            !     iat = 2

            !     open(400,file='OUT.dppp_bond',position='append')

            !     do jj = 1, numb(itype_mol)
            !         do ii1 = 1, 2
            !             if ( (abs(dppp1_m(1,iat,ii1,jj)) .gt. 1.0E-9) .or. &
            !                  (abs(dppp1_m(2,iat,ii1,jj)) .gt. 1.0E-9) .or. &
            !                  (abs(dppp1_m(3,iat,ii1,jj)) .gt. 1.0E-9) ) then
            !                 nzdp1 = nzdp1 + 1
            !             endif
            !         enddo
            !     enddo
            !     do jj = 1, numb(itype_mol)
            !         do ii1 = 1, 2
            !             if ( (abs(dppp2_m(1,iat,ii1,jj)) .gt. 1.0E-9) .or. &
            !                  (abs(dppp2_m(2,iat,ii1,jj)) .gt. 1.0E-9) .or. &
            !                  (abs(dppp2_m(3,iat,ii1,jj)) .gt. 1.0E-9) ) then
            !                 nzdp2 = nzdp2 + 1
            !             endif
            !         enddo
            !     enddo

            !     write(400,*) nzdp1, nzdp2

            !     do jj = 1, numb(itype_mol)
            !         do ii1 = 1, 2
            !             if ( (abs(dppp1_m(1,iat,ii1,jj)) .gt. 1.0E-9) .or. &
            !                  (abs(dppp1_m(2,iat,ii1,jj)) .gt. 1.0E-9) .or. &
            !                  (abs(dppp1_m(3,iat,ii1,jj)) .gt. 1.0E-9) ) then
            !                 write(400,401) jj, ii1, dppp1_m(:,iat,ii1,jj)
            !             endif
            !         enddo
            !     enddo
            !     do jj = 1, numb(itype_mol)
            !         do ii1 = 1, 2
            !             if ( (abs(dppp2_m(1,iat,ii1,jj)) .gt. 1.0E-9) .or. &
            !                  (abs(dppp2_m(2,iat,ii1,jj)) .gt. 1.0E-9) .or. &
            !                  (abs(dppp2_m(3,iat,ii1,jj)) .gt. 1.0E-9) ) then
            !                 write(400,401) jj, ii1, dppp2_m(:,iat,ii1,jj)
            !             endif
            !         enddo
            !     enddo
            !     close(400)

            !     ! output dppp -- atomic part
            !     open(400,file='OUT.dppp_atom',position='append')

            !     do jj = 1, 2*natom_m(itype_mol)
            !         do i1 = 1, 3
            !             if ( (abs(dppp_m(i1,jj,iat,1)) .gt. 1.0E-9) .or. &
            !                  (abs(dppp_m(i1,jj,iat,2)) .gt. 1.0E-9) .or. &
            !                  (abs(dppp_m(i1,jj,iat,3)) .gt. 1.0E-9) ) then
            !                 nzdp = nzdp + 1
            !             endif
            !         enddo
            !     enddo

            !     write(400,*) nzdp

            !     do jj = 1, 2*natom_m(itype_mol)
            !         do i1 = 1, 3
            !             if ( (abs(dppp_m(i1,jj,iat,1)) .gt. 1.0E-9) .or. &
            !                  (abs(dppp_m(i1,jj,iat,2)) .gt. 1.0E-9) .or. &
            !                  (abs(dppp_m(i1,jj,iat,3)) .gt. 1.0E-9) ) then
            !                 write(400,401) i1, jj, dppp_m(i1,jj,iat,:)
            !             endif
            !         enddo
            !     enddo
            !     close(400)

            !  endif

! 401 format(2(I7,1X),3(F20.9,1X)) 

             ppp = 0.d0
             dppp = 0.d0

             iat_count = 0
             do iat = 1, natom_m(itype_mol)

                 iat1 = ipol(1, iat, itype_mol)
                 iat2 = ipol(2, iat, itype_mol)
                 iat3 = ipol(3, iat, itype_mol)
                 iat4 = ipol(4, iat, itype_mol)

                !  do i = 1, 3
                !      dxp1(i) = xatom_m(i, iat2, imol, itype_mol) - xatom_m(i, iat1, imol, itype_mol)
                !      dxp2(i) = xatom_m(i, iat4, imol, itype_mol) - xatom_m(i, iat3, imol, itype_mol)
                !      if (abs(dxp1(i) + 1) .lt. abs(dxp1(i))) dxp1(i) = dxp1(i) + 1
                !      if (abs(dxp1(i) - 1) .lt. abs(dxp1(i))) dxp1(i) = dxp1(i) - 1
                !      if (abs(dxp2(i) + 1) .lt. abs(dxp2(i))) dxp2(i) = dxp2(i) + 1
                !      if (abs(dxp2(i) - 1) .lt. abs(dxp2(i))) dxp2(i) = dxp2(i) - 1
                !  end do

                !  dxpp1(:) = AL(:, 1)*dxp1(1) + AL(:, 2)*dxp1(2) + AL(:, 3)*dxp1(3)
                !  dxpp2(:) = AL(:, 1)*dxp2(1) + AL(:, 2)*dxp2(2) + AL(:, 3)*dxp2(3)

                !  dxpp3(1) = dxpp1(2)*dxpp2(3) - dxpp1(3)*dxpp2(2)   ! dpp3=dpp1 x dpp2
                !  dxpp3(2) = dxpp1(3)*dxpp2(1) - dxpp1(1)*dxpp2(3)
                !  dxpp3(3) = dxpp1(1)*dxpp2(2) - dxpp1(2)*dxpp2(1)

                !  d1 = dsqrt(dxpp1(1)**2 + dxpp1(2)**2 + dxpp1(3)**2)

                !  do i1 = 1, 3
                !      do i2 = 1, 3
                !          dxpp1_diat2(i1, i2) = -dxpp1(i1)*dxpp1(i2)/d1**3
                !      end do
                !      dxpp1_diat2(i1, i1) = dxpp1_diat2(i1, i1) + 1/d1
                !  end do

                !  d2 = dsqrt(dxpp2(1)**2 + dxpp2(2)**2 + dxpp2(3)**2)

                !  do i1 = 1, 3
                !      do i2 = 1, 3
                !          dxpp2_diat4(i1, i2) = -dxpp2(i1)*dxpp2(i2)/d2**3
                !      end do
                !      dxpp2_diat4(i1, i1) = dxpp2_diat4(i1, i1) + 1/d2
                !  end do

                !  d3 = dsqrt(dxpp3(1)**2 + dxpp3(2)**2 + dxpp3(3)**2)

                !  do i1 = 1, 3
                !      !cccc derivative for the 1/d part
                !      dxpp3_diat2(i1, 1) = -dxpp3(i1)/d3**3*(-dxpp3(2)*dxpp2(3) + dxpp3(3)*dxpp2(2))
                !      dxpp3_diat2(i1, 2) = -dxpp3(i1)/d3**3*(dxpp3(1)*dxpp2(3) - dxpp3(3)*dxpp2(1))
                !      dxpp3_diat2(i1, 3) = -dxpp3(i1)/d3**3*(-dxpp3(1)*dxpp2(2) + dxpp3(2)*dxpp2(1))

                !      dxpp3_diat4(i1, 1) = -dxpp3(i1)/d3**3*(dxpp3(2)*dxpp1(3) - dxpp3(3)*dxpp1(2))
                !      dxpp3_diat4(i1, 2) = -dxpp3(i1)/d3**3*(-dxpp3(1)*dxpp1(3) + dxpp3(3)*dxpp1(1))
                !      dxpp3_diat4(i1, 3) = -dxpp3(i1)/d3**3*(dxpp3(1)*dxpp1(2) - dxpp3(2)*dxpp1(1))
                !  end do

                !  dxpp3_diat2(1, 2) = dxpp3_diat2(1, 2) + dxpp2(3)/d3
                !  dxpp3_diat2(1, 3) = dxpp3_diat2(1, 3) - dxpp2(2)/d3
                !  dxpp3_diat2(2, 1) = dxpp3_diat2(2, 1) - dxpp2(3)/d3
                !  dxpp3_diat2(2, 3) = dxpp3_diat2(2, 3) + dxpp2(1)/d3
                !  dxpp3_diat2(3, 1) = dxpp3_diat2(3, 1) + dxpp2(2)/d3
                !  dxpp3_diat2(3, 2) = dxpp3_diat2(3, 2) - dxpp2(1)/d3

                !  dxpp3_diat4(1, 2) = dxpp3_diat4(1, 2) - dxpp1(3)/d3
                !  dxpp3_diat4(1, 3) = dxpp3_diat4(1, 3) + dxpp1(2)/d3
                !  dxpp3_diat4(2, 1) = dxpp3_diat4(2, 1) + dxpp1(3)/d3
                !  dxpp3_diat4(2, 3) = dxpp3_diat4(2, 3) - dxpp1(1)/d3
                !  dxpp3_diat4(3, 1) = dxpp3_diat4(3, 1) - dxpp1(2)/d3
                !  dxpp3_diat4(3, 2) = dxpp3_diat4(3, 2) + dxpp1(1)/d3

                !  dxpp1 = dxpp1/d1
                !  dxpp2 = dxpp2/d2
                !  dxpp3 = dxpp3/d3

                 ! use precalculated data
                 dxpp1(:)=dxpp(:,1,iat,1)
                 dxpp2(:)=dxpp(:,2,iat,1)
                 dxpp3(:)=dxpp(:,3,iat,1)
                 dxpp1_diat2(:,:)=dxpp1diat2(:,:,iat,1)
                 dxpp2_diat4(:,:)=dxpp2diat4(:,:,iat,1)
                 dxpp3_diat2(:,:)=dxpp3diat2(:,:,iat,1)
                 dxpp3_diat4(:,:)=dxpp3diat4(:,:,iat,1)

                 dxyzp_store(:, 1, iat) = dxpp1(:)
                 dxyzp_store(:, 2, iat) = dxpp2(:)
                 dxyzp_store(:, 3, iat) = dxpp3(:)

                 ! dpp1=(dxpp1,dypp1,dzpp1),dpp2,dpp3 are the three polarization
                 ! directions, to be applied with the polarization  model with
                 !------------------------------------------------

                 do kkk = 1, 2
                     iat_count = iat_count + 1

                     ppp(1, iat_count) = dxpp1(1)*pxyz(1, iat_count) + dxpp1(2)*pxyz(2, iat_count) + dxpp1(3)*pxyz(3, iat_count)
                     ppp(2, iat_count) = dxpp2(1)*pxyz(1, iat_count) + dxpp2(2)*pxyz(2, iat_count) + dxpp2(3)*pxyz(3, iat_count)
                     ppp(3, iat_count) = dxpp3(1)*pxyz(1, iat_count) + dxpp3(2)*pxyz(2, iat_count) + dxpp3(3)*pxyz(3, iat_count)

                     dppp(1, iat_count, iat, :) = dppp(1, iat_count, iat, :) + &
                         dxpp1(1)*dpxyz(1, :, iat_count) + dxpp1(2)*dpxyz(2, :, iat_count) + dxpp1(3)*dpxyz(3, :, iat_count)
                     dppp(2, iat_count, iat, :) = dppp(2, iat_count, iat, :) + &
                         dxpp2(1)*dpxyz(1, :, iat_count) + dxpp2(2)*dpxyz(2, :, iat_count) + dxpp2(3)*dpxyz(3, :, iat_count)
                     dppp(3, iat_count, iat, :) = dppp(3, iat_count, iat, :) + &
                         dxpp3(1)*dpxyz(1, :, iat_count) + dxpp3(2)*dpxyz(2, :, iat_count) + dxpp3(3)*dpxyz(3, :, iat_count)

                     dppp(1, iat_count, iat2, :) = dppp(1, iat_count, iat2, :) + &
                         dxpp1_diat2(1, :)*pxyz(1, iat_count) + dxpp1_diat2(2, :)*pxyz(2, iat_count) + &
                         dxpp1_diat2(3, :)*pxyz(3, iat_count)

                     dppp(1, iat_count, iat1, :) = dppp(1, iat_count, iat1, :) - &
                         dxpp1_diat2(1, :)*pxyz(1, iat_count) - dxpp1_diat2(2, :)*pxyz(2, iat_count) - &
                         dxpp1_diat2(3, :)*pxyz(3, iat_count)

                     dppp(2, iat_count, iat4, :) = dppp(2, iat_count, iat4, :) + &
                         dxpp2_diat4(1, :)*pxyz(1, iat_count) + dxpp2_diat4(2, :)*pxyz(2, iat_count) + &
                         dxpp2_diat4(3, :)*pxyz(3, iat_count)

                     dppp(2, iat_count, iat3, :) = dppp(2, iat_count, iat3, :) - &
                         dxpp2_diat4(1, :)*pxyz(1, iat_count) - dxpp2_diat4(2, :)*pxyz(2, iat_count) - &
                         dxpp2_diat4(3, :)*pxyz(3, iat_count)

                     dppp(3, iat_count, iat2, :) = dppp(3, iat_count, iat2, :) + &
                         dxpp3_diat2(1, :)*pxyz(1, iat_count) + dxpp3_diat2(2, :)*pxyz(2, iat_count) + &
                         dxpp3_diat2(3, :)*pxyz(3, iat_count)

                     dppp(3, iat_count, iat1, :) = dppp(3, iat_count, iat1, :) - &
                         dxpp3_diat2(1, :)*pxyz(1, iat_count) - dxpp3_diat2(2, :)*pxyz(2, iat_count) - &
                         dxpp3_diat2(3, :)*pxyz(3, iat_count)

                     dppp(3, iat_count, iat4, :) = dppp(3, iat_count, iat4, :) + &
                         dxpp3_diat4(1, :)*pxyz(1, iat_count) + dxpp3_diat4(2, :)*pxyz(2, iat_count) + &
                         dxpp3_diat4(3, :)*pxyz(3, iat_count)

                     dppp(3, iat_count, iat3, :) = dppp(3, iat_count, iat3, :) - &
                         dxpp3_diat4(1, :)*pxyz(1, iat_count) - dxpp3_diat4(2, :)*pxyz(2, iat_count) - &
                         dxpp3_diat4(3, :)*pxyz(3, iat_count)

                 end do  ! kkk

             end do


             !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

             ip = 0
             E_polar_tmp2 = 0.d0
             do jj = 1, num_term_pol(itype_mol)
                 do i1 = 1, 3
                     do i2 = 1, 3
                         ip = ip + 1

                         ! E_polar=E_polar+ppp(i1,ind(1,jj,itype_mol))*ppp(i2,ind(2,jj,itype_mol))*BB_poldir(ip,itype_mol)
                        !  E_polar_tmp2 = E_polar_tmp2 + ppp(i1, ind(1, jj, itype_mol))*ppp(i2, ind(2, jj, itype_mol))*BB_poldir(ip, itype_mol)
                         E_polar_tmp2=E_polar_tmp2+(ppp(i1,ind(1,jj,itype_mol))-ppp_m(i1,ind(1,jj,itype_mol)))*(ppp(i2,ind(2,jj,itype_mol))-ppp_m(i2,ind(2,jj,itype_mol)))*BB_poldir(ip,itype_mol)

                         do iat = 1, natom_m(itype_mol)
                            !  force_polar_m(:, iat, imol, itype_mol) = force_polar_m(:, iat, imol, itype_mol) + &
                            !      (dppp(i1, ind(1, jj, itype_mol), iat, :)*ppp(i2, ind(2, jj, itype_mol)) + &
                            !      ppp(i1, ind(1, jj, itype_mol))*dppp(i2, ind(2, jj, itype_mol), iat, :))*BB_poldir(ip, itype_mol)
                            force_polar_m(:,iat,imol,itype_mol)=force_polar_m(:,iat,imol,itype_mol)+ &
                                ((dppp(i1,ind(1,jj,itype_mol),iat,:)-dppp_m(i1,ind(1,jj,itype_mol),iat,:))*(ppp(i2,ind(2,jj,itype_mol))-ppp_m(i2,ind(2,jj,itype_mol)))+ &
                                (ppp(i1,ind(1,jj,itype_mol))-ppp_m(i1,ind(1,jj,itype_mol)))*(dppp(i2,ind(2,jj,itype_mol),iat,:)-dppp_m(i2,ind(2,jj,itype_mol),iat,:)))*BB_poldir(ip,itype_mol)

                                end do
                     end do
                 end do
             end do


            ! ------------------------------------------------------------
            ! If iflag_save_ppp.eq.1, save ppp, ppp1, ppp2 in OUT.PPP
            ! ------------------------------------------------------------
             if ((inode_tot.eq.1) .and. (iflag_training .eq. 1)) then
                open(400,file='OUT.PPP',position='append')
                ! rewind(400)
                
                write(400,*) natom_m(itype_mol)
                write(400,*) numb(itype_mol)
                write(400,*) natom_m(itype_mol)*6+numb(itype_mol)*4
                
                iat_count=0
                do iat=1,natom_m(itype_mol)
                    iat_count=iat_count+1
                    write(400,*) ppp(1,iat_count)
                    write(400,*) ppp(1,iat_count+1)
                    write(400,*) ppp(2,iat_count)
                    write(400,*) ppp(2,iat_count+1)
                    write(400,*) ppp(3,iat_count)
                    write(400,*) ppp(3,iat_count+1)
                    iat_count=iat_count+1
                enddo
                
                do jj=1,numb(itype_mol)
                    write(400,*) ppp1(1,jj)
                    write(400,*) ppp1(2,jj)
                    write(400,*) ppp2(1,jj)
                    write(400,*) ppp2(2,jj)
                enddo
                
                close(400)
             endif


            !  if ((inode_tot.eq.1) .and. (iflag_training .eq. 0)) then
            !     open(400,file='OUT.PPP',position='append')
            !     ! rewind(400)
                
            !     write(400,*) natom_m(itype_mol)
            !     write(400,*) numb(itype_mol)
            !     write(400,*) natom_m(itype_mol)*6+numb(itype_mol)*4
                
            !     iat_count=0
            !     do iat=1,natom_m(itype_mol)
            !         iat_count=iat_count+1
            !         write(400,*) ppp(1,iat_count),   ppp_m(1,iat_count)
            !         write(400,*) ppp(1,iat_count+1), ppp_m(1,iat_count+1)
            !         write(400,*) ppp(2,iat_count),   ppp_m(2,iat_count)
            !         write(400,*) ppp(2,iat_count+1), ppp_m(2,iat_count+1)
            !         write(400,*) ppp(3,iat_count),   ppp_m(3,iat_count)
            !         write(400,*) ppp(3,iat_count+1), ppp_m(3,iat_count+1)
            !         iat_count=iat_count+1
            !     enddo
                
            !     do jj=1,numb(itype_mol)
            !         write(400,*) ppp1(1,jj), ppp1_m(1,jj)
            !         write(400,*) ppp1(2,jj), ppp1_m(2,jj)
            !         write(400,*) ppp2(1,jj), ppp2_m(1,jj)
            !         write(400,*) ppp2(2,jj), ppp2_m(2,jj)
            !     enddo
                
            !     close(400)
            !  endif


             E_polar_tmp2 = E_polar_tmp2 + Ep_bond

             !occcccccccccccccccccccccccccccccccc
             do iat = 1, natom_m(itype_mol)
                 force_polar_m(:, iat, imol, itype_mol) = force_polar_m(:, iat, imol, itype_mol) + force_Epp(:, iat)
             end do

             ! different molecule can have different E_polar_max, here we just used a
             ! unified one
             if (E_polar_tmp2 .lt. 0.d0) then
                 yy = -E_polar_tmp2/E_polar_max(itype_mol)
                 E_polar_tmp3 = -E_polar_max(itype_mol)*(1 - (1 + yy)*exp(-2*yy))
                 fact_polar(imol, itype_mol) = (1 + 2*yy)*exp(-2*yy)
             else
                 E_polar_tmp3 = E_polar_tmp2
                 fact_polar(imol, itype_mol) = 1.d0
             end if

             E_polar = E_polar + E_polar_tmp3

             do iat = 1, natom_m(itype_mol)
                 force_polar_m(:, iat, imol, itype_mol) = force_polar_m(:, iat, imol, itype_mol)*fact_polar(imol, itype_mol)
             end do

             !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
             !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

             !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
             !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

             do jj = 1, numb(itype_mol)
                 jj1 = indb(1, jj, itype_mol)
                 jj2 = indb(2, jj, itype_mol)
                 ip = (indb(3, jj, itype_mol) - 1)*4
                 do ii1 = 1, 2
                     do ii2 = 1, 2
                         ip = ip + 1
                        !  pp_tmp = ppp2(ii2, jj)*BB_direct(ip, itype_mol)*fact_polar(imol, itype_mol)
                         pp_tmp=(ppp2(ii2,jj)-ppp2_m(ii2,jj))*BB_direct(ip,itype_mol)*fact_polar(imol,itype_mol)
                         CC_pol(1 + (ii1 - 1)*3, jj1, imol, itype_mol) = CC_pol(1 + (ii1 - 1)*3, jj1, imol, itype_mol) + Epp(1, jj)*pp_tmp
                         CC_pol(2 + (ii1 - 1)*3, jj1, imol, itype_mol) = CC_pol(2 + (ii1 - 1)*3, jj1, imol, itype_mol) + Epp(2, jj)*pp_tmp
                         CC_pol(3 + (ii1 - 1)*3, jj1, imol, itype_mol) = CC_pol(3 + (ii1 - 1)*3, jj1, imol, itype_mol) + Epp(3, jj)*pp_tmp

                        !  pp_tmp = ppp1(ii1, jj)*BB_direct(ip, itype_mol)*fact_polar(imol, itype_mol)
                         pp_tmp=(ppp1(ii1,jj)-ppp1_m(ii1,jj))*BB_direct(ip,itype_mol)*fact_polar(imol,itype_mol)
                         CC_pol(1 + (ii2 - 1)*3, jj2, imol, itype_mol) = CC_pol(1 + (ii2 - 1)*3, jj2, imol, itype_mol) + Epp(1, jj)*pp_tmp
                         CC_pol(2 + (ii2 - 1)*3, jj2, imol, itype_mol) = CC_pol(2 + (ii2 - 1)*3, jj2, imol, itype_mol) + Epp(2, jj)*pp_tmp
                         CC_pol(3 + (ii2 - 1)*3, jj2, imol, itype_mol) = CC_pol(3 + (ii2 - 1)*3, jj2, imol, itype_mol) + Epp(3, jj)*pp_tmp
                     end do
                 end do
             end do
             !cccccccccccccccccccccccccccccccccccccccccccccccccccccc
             !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

             ip = 0
             do jj = 1, num_term_pol(itype_mol)
                 do i1 = 1, 3
                     do i2 = 1, 3
                         ip = ip + 1
                         ! Note, ind(1,jj), ind(2,jj) are not equivalent

                         iat_tmp = (ind(1, jj, itype_mol) - 1)/2 + 1
                         ii1 = ind(1, jj, itype_mol) - 1 - (iat_tmp - 1)*2
                         ! jj = 0 or 1
                        !  pp_tmp = ppp(i2, ind(2, jj, itype_mol))*BB_poldir(ip, itype_mol)
                         pp_tmp=(ppp(i2,ind(2,jj,itype_mol))-ppp_m(i2,ind(2,jj,itype_mol)))*BB_poldir(ip,itype_mol)

                         pp_tmp = pp_tmp*fact_polar(imol, itype_mol)   ! this is for polarization force on other molecules
                         if (i1 .eq. 1) then
                             CC_pol(1 + ii1*3, iat_tmp, imol, itype_mol) = CC_pol(1 + ii1*3, iat_tmp, imol, itype_mol) + &
                                 dxyzp_store(1, 1, iat_tmp)*pp_tmp
                             CC_pol(2 + ii1*3, iat_tmp, imol, itype_mol) = CC_pol(2 + ii1*3, iat_tmp, imol, itype_mol) + &
                                 dxyzp_store(2, 1, iat_tmp)*pp_tmp
                             CC_pol(3 + ii1*3, iat_tmp, imol, itype_mol) = CC_pol(3 + ii1*3, iat_tmp, imol, itype_mol) + &
                                 dxyzp_store(3, 1, iat_tmp)*pp_tmp
                         elseif (i1 .eq. 2) then
                             CC_pol(1 + ii1*3, iat_tmp, imol, itype_mol) = CC_pol(1 + ii1*3, iat_tmp, imol, itype_mol) + &
                                 dxyzp_store(1, 2, iat_tmp)*pp_tmp
                             CC_pol(2 + ii1*3, iat_tmp, imol, itype_mol) = CC_pol(2 + ii1*3, iat_tmp, imol, itype_mol) + &
                                 dxyzp_store(2, 2, iat_tmp)*pp_tmp
                             CC_pol(3 + ii1*3, iat_tmp, imol, itype_mol) = CC_pol(3 + ii1*3, iat_tmp, imol, itype_mol) + &
                                 dxyzp_store(3, 2, iat_tmp)*pp_tmp
                         elseif (i1 .eq. 3) then
                             CC_pol(1 + ii1*3, iat_tmp, imol, itype_mol) = CC_pol(1 + ii1*3, iat_tmp, imol, itype_mol) + &
                                 dxyzp_store(1, 3, iat_tmp)*pp_tmp
                             CC_pol(2 + ii1*3, iat_tmp, imol, itype_mol) = CC_pol(2 + ii1*3, iat_tmp, imol, itype_mol) + &
                                 dxyzp_store(2, 3, iat_tmp)*pp_tmp
                             CC_pol(3 + ii1*3, iat_tmp, imol, itype_mol) = CC_pol(3 + ii1*3, iat_tmp, imol, itype_mol) + &
                                 dxyzp_store(3, 3, iat_tmp)*pp_tmp
                         end if

                         iat_tmp = (ind(2, jj, itype_mol) - 1)/2 + 1
                         ii2 = ind(2, jj, itype_mol) - 1 - (iat_tmp - 1)*2
                         ! jj=0 or 1
                        !  pp_tmp = ppp(i1, ind(1, jj, itype_mol))*BB_poldir(ip, itype_mol)
                         pp_tmp=(ppp(i1,ind(1,jj,itype_mol))-ppp_m(i1,ind(1,jj,itype_mol)))*BB_poldir(ip,itype_mol)

                         pp_tmp = pp_tmp*fact_polar(imol, itype_mol)   ! this is for polarization force on other molecules
                         if (i2 .eq. 1) then
                             CC_pol(1 + ii2*3, iat_tmp, imol, itype_mol) = CC_pol(1 + ii2*3, iat_tmp, imol, itype_mol) + &
                                 dxyzp_store(1, 1, iat_tmp)*pp_tmp
                             CC_pol(2 + ii2*3, iat_tmp, imol, itype_mol) = CC_pol(2 + ii2*3, iat_tmp, imol, itype_mol) + &
                                 dxyzp_store(2, 1, iat_tmp)*pp_tmp
                             CC_pol(3 + ii2*3, iat_tmp, imol, itype_mol) = CC_pol(3 + ii2*3, iat_tmp, imol, itype_mol) + &
                                 dxyzp_store(3, 1, iat_tmp)*pp_tmp
                         elseif (i2 .eq. 2) then
                             CC_pol(1 + ii2*3, iat_tmp, imol, itype_mol) = CC_pol(1 + ii2*3, iat_tmp, imol, itype_mol) + &
                                 dxyzp_store(1, 2, iat_tmp)*pp_tmp
                             CC_pol(2 + ii2*3, iat_tmp, imol, itype_mol) = CC_pol(2 + ii2*3, iat_tmp, imol, itype_mol) + &
                                 dxyzp_store(2, 2, iat_tmp)*pp_tmp
                             CC_pol(3 + ii2*3, iat_tmp, imol, itype_mol) = CC_pol(3 + ii2*3, iat_tmp, imol, itype_mol) + &
                                 dxyzp_store(3, 2, iat_tmp)*pp_tmp
                         elseif (i2 .eq. 3) then
                             CC_pol(1 + ii2*3, iat_tmp, imol, itype_mol) = CC_pol(1 + ii2*3, iat_tmp, imol, itype_mol) + &
                                 dxyzp_store(1, 3, iat_tmp)*pp_tmp
                             CC_pol(2 + ii2*3, iat_tmp, imol, itype_mol) = CC_pol(2 + ii2*3, iat_tmp, imol, itype_mol) + &
                                 dxyzp_store(2, 3, iat_tmp)*pp_tmp
                             CC_pol(3 + ii2*3, iat_tmp, imol, itype_mol) = CC_pol(3 + ii2*3, iat_tmp, imol, itype_mol) + &
                                 dxyzp_store(3, 3, iat_tmp)*pp_tmp
                         end if

                     end do
                 end do
             end do

             if (switch) then
                 loop_polarization_end = mpi_wtime()
                 loop_polarization_total_time = loop_polarization_total_time + loop_polarization_end - loop_polarization_start
             endif

6010     continue

6000 continue   ! imolecule
6001 continue

     if (switch) loop_6001_end = mpi_wtime()

     call mpi_allreduce(force_m, force_m_t, 3*natom_mm*nmolm*ntype_m, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
     force_m = force_m_t

    !  call mpi_allreduce(EtotM, EtotM_t, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    !  EtotM = EtotM_t

     call mpi_allreduce(E_kin2_M, E_kin2_M_t, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

     E_kin2_M = E_kin2_M_t
     E_kin2_nonbond = E_kin2_nonbond - E_kin2_M


     !$ for new strategy
     ! Etot = Etot - EtotM

     if ((inode_tot .eq. 1) .and. iflag_debug) then
        open (36, file='OUT.DB', position='append')
        write (36, "(3(E25.17,1X))") Etot*Hartree_eV, EtotM*Hartree_eV, (Etot - EtotM)*Hartree_eV
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


     ! Note, CC_Pol is not full, but we don;t need it to be full,
     ! sinc the job is divided in the same way in the calc_polar_F
     if ((iflag_polar .eq. 1) .and. (iflag_training .eq. 0)) then

         call mpi_allreduce(E_polar, E_polar_t, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

         E_polar = E_polar_t


         E_polar = E_polar * polar_factor


         Etot = Etot + E_polar/Hartree_eV

         Etot_nonbond = Etot*Hartree_eV

         if (switch) calc_polar_F_start = mpi_wtime()

         call calc_polar_F(natom_mm, nmolm, ntype_mm, &
             natom_m, nmol, xatom_m, AL, force_polar_m, rho_z, rho, CC_pol, c_devptr_xatom_m, c_devptr_rho_z, c_devptr_rho)

         if (switch) calc_polar_F_end = mpi_wtime()

         ! The unit for force_polar_m is eV/Bohr
         call mpi_allreduce(force_polar_m, force_m_t, 3*natom_mm*nmolm*ntype_m, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

         force_polar_m = force_m_t

         
         force_polar_m = force_polar_m*polar_factor


         force_m = force_m + force_polar_m/Hartree_eV   ! convert to Hartree/Bohr


         if ((inode_tot .eq. 1) .and. iflag_debug) then
             open (36, file='OUT.POL', position='append')
             write (36, "(1(E25.17,1X))") E_polar
             i = 0
             do itype_mol = 1, ntype_m
                 do imol = 1, nmol(itype_mol)
                     do ii = 1, natom_m(itype_mol)
                         i = i + 1
                         write (36, "(i4, 1x, 3(E25.17,1X))") i, force_polar_m(:, ii, imol, itype_mol)/A_AU_1
                     end do
                 end do
             end do
             close (36)

             write (6, "('Etot,E_polar,DB_binding = ', 3(F25.6,1X))") Etot*Hartree_eV, E_polar, (Etot - EtotM)*Hartree_eV
         end if

         
     end if

     if (switch) subroutine_nonbond_force_end = mpi_wtime()

     deallocate (rho)
     deallocate (rho_z)

     if (iflag_polar .eq. 1) then
        deallocate(ppp_m)
        deallocate(dppp_m)
     endif

     return
 end subroutine nonbond_force

