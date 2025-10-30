subroutine calc_polar_F(natom_mm, nmolm, ntype_mm, &
      natom_m, nmol, xatom_m, AL, force_m, rho_totz, rho_tot, &
      CC_pol, c_devptr_xatom_m, c_devptr_rho_totz, c_devptr_rho_tot)
      !         call calc_polar_F(natom_mm,nmolm,ntype_mm, &
      !  natom_m,nmol,xatom_m,AL,force_polar_m,rho_z,rho,CC_pol,c_devptr_xatom_m)
      ! This subroutine take the charge density of the previous nonbond
      ! calculation, but only calculate the force of the polarization term
      use mod_mpi
      use mod_profile
      use mod_param_densityFF
      use cudafor
      use mod_cuinterface
      use mod_mem_preallocated
  
      implicit double precision(a - h, o - z)
  
      real*8 rho_tot(n1, n2, n3), rho_totz(n1, n2, n3)
      real*8 AL(3, 3), vol, vol_n
      real*8, allocatable, device, dimension(:, :) :: AL_d
      real*8 ALm(3, 3)
      real*8, allocatable, device, dimension(:, :) :: ALm_d
      real*8 ALtt(3, 3)
      real*8 xc_cent(3)
      real*8, allocatable, device, dimension(:) :: xc_cent_d
      integer icorner(3)
      integer, allocatable, device, dimension(:) :: icorner_d
      real*8 Etot
      real*8 pi
  
      integer natom_m_t(1000)
  
      integer iatom(200)
      real*8 xatom(3, 200)
      real*8 pp(2000)
      real*8 pp2(2000)
      character*40 file_ionrhoR(20)
      real*8 CC_pol(6, natom_mm, nmolm, ntype_mm)
      real*8, allocatable, device, dimension(:, :, :, :) :: CC_pol_d
      real*8 fact_store(1000), fact2_store(1000)
      real*8, allocatable, device, dimension(:) :: fact_store_d, fact2_store_d
      type(c_devptr) :: c_devptr_xatom_m, c_devptr_rho_totz, c_devptr_rho_tot
  
      real*8, allocatable, dimension(:) :: vr_tmp
      real*8, allocatable, dimension(:, :, :) :: vr_tmp2
  
      integer natom_m(ntype_mm), nmol(ntype_mm)
      real*8 xatom_m(3, natom_mm, nmolm, ntype_mm)
      real*8 force_m(3, natom_mm, nmolm, ntype_mm)
  
      real*8, allocatable, dimension(:, :) :: xatom_cent
      real*8, allocatable, device, dimension(:, :) :: xatom_cent_d
  
      real*8, allocatable, dimension(:, :, :) :: box, box2, box3
      real*8, allocatable, dimension(:, :, :, :) :: dxyz_box
      real*8, allocatable, device, dimension(:, :, :, :) :: dxyz_box_d
      real*8, allocatable, dimension(:, :, :, :) :: dbox, dbox2, dbox3
      real*8, allocatable, dimension(:, :, :, :) :: box_c, box2_c, box3_c
      real*8, allocatable, dimension(:, :, :) :: rhop_tot, vrhop_tot_coul
      real*8, allocatable, device, dimension(:, :, :) :: rhop_tot_d, vrhop_tot_coul_d
      real*8, allocatable, dimension(:, :, :) :: rho_t
      real*8, allocatable, dimension(:, :, :) :: dvxc_tot
      real*8, allocatable, device, dimension(:, :, :) :: dvxc_tot_d
      real*8, allocatable, dimension(:, :, :) :: dvxc_m
      real*8, allocatable, device, dimension(:, :, :) :: dvxc_m_d
      real*8, allocatable, dimension(:, :, :) :: rhop_m, vrhop_m_coul
      real*8, allocatable, device, dimension(:, :, :) :: rhop_m_d, vrhop_m_coul_d
      real*8, allocatable, dimension(:, :, :) :: rho_m
      real*8, allocatable, device, dimension(:, :, :) :: rho_m_d
  
      ! real*8 pxyz(3,1000),ppp(3,1000)
      ! real*8 dxp1(3),dxp2(3),dxp3(3)
  
      integer iflag_mol(1000)
      integer ierr
      character*40 f_xatom
  
      real*8 force_cent(3, 1000)
      real*8, allocatable, device, dimension(:, :) :: force_cent_d
  

      real*8 cp_total_time_begin, cp_total_time_end
      real*8 cp_copy_1_begin, cp_copy_1_end
      real*8 cp_113to129_begin, cp_113to129_end
      real*8 cp_500_begin, cp_500_end, cp_500_total
      real*8 cp_calc_begin, cp_calc_end
      real*8 cp_312to339_begin, cp_312to339_end, cp_312to339_total
      real*8 cp_800_begin, cp_800_end, cp_800_total
      real*8 cp_copy_2_begin, cp_copy_2_end, cp_copy_2_total


      cp_500_total = 0.d0
      cp_312to339_total = 0.d0
      cp_800_total = 0.d0
      cp_copy_2_total = 0.d0

      if (switch) cp_total_time_begin = mpi_wtime()





      ntype_m = ntype_mm
  
      imax_ncent_itype_mol = 0
      do itype_mol = 1, ntype_m
          imax_ncent_itype_mol = max(imax_ncent_itype_mol, ncent(itype_mol))
      end do
  
      call c_f_pointer(c_devptr_xatom_cent, xatom_cent_d, (/3, 10*natom_mm/))
      call c_f_pointer(c_devptr_icorner, icorner_d, (/3/))
      call c_f_pointer(c_devptr_fact2_store, fact2_store_d, (/1000/))
      call c_f_pointer(c_devptr_fact_store, fact_store_d, (/1000/))
      call c_f_pointer(c_devptr_force_cent, force_cent_d, (/3, 10*natom_mm/))
  
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

  
      call c_f_pointer(c_devptr_AL, AL_d, (/3, 3/))
      call c_f_pointer(c_devptr_ALm, ALm_d, (/3, 3/))
      call c_f_pointer(c_devptr_dxyz_box, dxyz_box_d, (/3, -id1:id1, -id2:id2, -id3:id3/))
      call c_f_pointer(c_devptr_CC_pol, CC_pol_d, (/6, natom_mm, nmolm, ntype_mm/))
      call c_f_pointer(c_devptr_xc_cent, xc_cent_d, (/3/))
  

      if (switch) cp_copy_1_begin = mpi_wtime()

      ierr = cudaMemcpy(AL_d, AL, 3*3, cudaMemcpyHostToDevice)
      ierr = cudaMemcpy(CC_pol_d, CC_pol, 6*natom_mm*nmolm*ntype_mm, cudaMemcpyHostToDevice)
  
      if (switch) cp_copy_1_end = mpi_wtime()

      if (switch) cp_113to129_begin = mpi_wtime()

      call cudakernel_polar_loop_113to129(c_devptr_dxyz_box, c_devptr_AL, id3, id2, id1, n3, n2, n1)
  
      if (switch) cp_113to129_end = mpi_wtime()


      call c_f_pointer(c_devptr_rhop_tot, rhop_tot_d, (/n1, n2, n3/))
      call c_f_pointer(c_devptr_vrhop_tot_coul, vrhop_tot_coul_d, (/n1, n2, n3/))
      call c_f_pointer(c_devptr_dvxc_tot, dvxc_tot_d, (/n1, n2, n3/))
      !cccccccccccccccc  first, get the whole charge density
  
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !  First do loop, just to get rhop_tot
  
      ! rhop_tot=0.d0
      call zero_gpu_mem(c_devptr_rhop_tot, n1*n2*n3*8)
  
      kkk_mol = 0
      iflag_mol(1:ntype_m) = 0
      do 5001 itype_mol = 1, ntype_m
  
          do 5000 imol = 1, nmol(itype_mol)
  
              kkk_mol = kkk_mol + 1
              if (mod(kkk_mol - 1, nnodes_tot) .ne. inode_tot - 1) goto 5000
  
              if (switch) cp_500_begin = mpi_wtime()

              call cudakernel_polar_loop_500(natom_m(itype_mol), c_devptr_ion_type_atomp, c_devptr_atom_charge_param, c_devptr_xatom_m, imol, itype_mol, n3, n2, n1, id3, id2, id1, natom_mm, nmolm, ntype_mm, c_devptr_AL, c_devptr_dxyz_box, c_devptr_CC_pol, c_devptr_rhop_tot, Rbox2)
  
              if (switch) then 
                cp_500_end = mpi_wtime()
                cp_500_total = cp_500_total + cp_500_end - cp_500_begin
              endif

  5000    continue   ! imolecule
  5001 continue
       !ccccccccccccccccccccccccccccccccccccccccccccccccccccc
       ! ierr = cudaMemcpy(rhop_tot, rhop_tot_d, n1 * n2 * n3, cudaMemcpyDeviceToHost)
  

       if (switch) cp_calc_begin = mpi_wtime()

       call nccl_mpi_allreduce(c_devptr_rhop_tot, c_devptr_rhop_tot, n1*n2*n3, 10, 0)

       call cudakernel_calc_coul(c_devptr_rhop_tot, c_devptr_vrhop_tot_coul, n1, n2, n3, c_devptr_AL, pi)

       call cudakernel_calc_dvxc(c_devptr_rho_tot, c_devptr_rhop_tot, c_devptr_dvxc_tot, n1, n2, n3, pi)

       if (switch) cp_calc_end = mpi_wtime()

  
       kkk_mol = 0
       iflag_mol(1:ntype_m) = 0

       do 6001 itype_mol = 1, ntype_m
  
        nm1 = nm1_all(itype_mol)
        nm2 = nm2_all(itype_mol)
        nm3 = nm3_all(itype_mol)
  
           do 6000 imol = 1, nmol(itype_mol)
               kkk_mol = kkk_mol + 1
  
                if (mod(kkk_mol - 1, nnodes_tot) .ne. inode_tot - 1) goto 6000
    
                call zero_gpu_mem(c_devptr_xc_cent, 3*8)

                if(switch) cp_312to339_begin = mpi_wtime()

                call cudakernel_polar_loop_312to339(ncent(itype_mol), c_devptr_icent, imol, itype_mol, c_devptr_w_cent, c_devptr_xc_cent, c_devptr_xatom_cent, c_devptr_icorner, n3, n2, n1, nm3, nm2, nm1, c_devptr_xatom_m, natom_mm, nmolm, c_devptr_nat_cent)
    
                if (switch) then
                    cp_312to339_end = mpi_wtime()
                    cp_312to339_total = cp_312to339_total + cp_312to339_end - cp_312to339_begin
                endif

                if (switch) cp_800_begin = mpi_wtime()

                call cudakernel_polar_loop_800(ncent(itype_mol), c_devptr_itype_cent, c_devptr_xatom_cent, c_devptr_icorner, c_devptr_AL, n3, n2, n1, c_devptr_ion_type_cent, c_devptr_imax_ion, &
                    c_devptr_r_ion, c_devptr_dxyz_box, c_devptr_funcr2, c_devptr_box, c_devptr_box2, c_devptr_box3, c_devptr_dbox, c_devptr_dbox2, c_devptr_dbox3, &
                    c_devptr_fact_store, c_devptr_fact2_store, Rbox2, nr, nm3, nm2, nm1, c_devptr_vrhop_tot_coul, c_devptr_vrhop_m_coul, c_devptr_dvxc_tot, c_devptr_dvxc_m, &
                    vol_n, id3, id2, id1, c_devptr_force_cent, Rm2, itype_mol, c_devptr_rho_ion, imax_nr, imax_ntype_cent)

                if (switch) then
                    cp_800_end = mpi_wtime()
                    cp_800_total = cp_800_total + cp_800_end - cp_800_begin
                endif


                if (switch) cp_copy_2_begin = mpi_wtime()

                ierr = cudaMemcpy(force_cent, force_cent_d, 3*10*natom_mm, cudaMemcpyDeviceToHost)

                if (switch) then
                    cp_copy_2_end = mpi_wtime()
                    cp_copy_2_total = cp_copy_2_total + cp_copy_2_end - cp_copy_2_begin
                endif

                do ii = 1, ncent(itype_mol)
                    do jj = 1, nat_cent(ii, itype_mol)
                        do ixyz = 1, 3
                            force_m(ixyz, icent(jj, ii, itype_mol), imol, itype_mol) = force_m(ixyz, icent(jj, ii, itype_mol), imol, itype_mol) &
                                + w_cent(jj, ii, itype_mol)*force_cent(ixyz, ii)
                        end do
                    end do
                end do
  
  6000     continue   ! imolecule
  
  6001 continue
  
       if (switch) cp_total_time_end = mpi_wtime()


    !    write(6,*) " @@@ total   ", cp_total_time_end - cp_total_time_begin
    !    write(6,*) " @@@ copy 1  ", cp_copy_1_end - cp_copy_1_begin
    !    write(6,*) " @@@ 113     ", cp_113to129_end - cp_113to129_begin
    !    write(6,*) " @@@ 500     ", cp_500_total
    !    write(6,*) " @@@ calc    ", cp_calc_end - cp_calc_begin
    !    write(6,*) " @@@ 312     ", cp_312to339_total
    !    write(6,*) " @@@ 800     ", cp_800_total
    !    write(6,*) " @@@ copy 2  ", cp_copy_2_total

       return
   end subroutine calc_polar_F
  
  