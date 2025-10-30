subroutine ML_FF_EF_PC(matom_m, nmol, xatom, AL, fatom, Etot)

   use mod_mpi
   use mod_data, only: natom_mm, ntype_m, nmolm

   use calc_ftype2b_3b, only : load_model_type2, gen_feature_type2, &
   & m_neigh, nfeat0m2, nfeat0m3, feat2, feat3, dfeat2, dfeat3, neigh_numR, neigh_listR
   use calc_lin_PC, only : set_paths_lin, load_model_lin, set_image_info_lin, cal_energy_force_lin, &
   & energy_pred_lin, etot_pred_lin, force_pred_lin, natom_nn

   implicit none

   integer, intent(in) :: matom_m(ntype_m)
   integer, intent(in) :: nmol(ntype_m)
   real*8, intent(in)  :: xatom(3, natom_mm, nmolm, ntype_m)      ! NEED TO ALLOCATE
   real*8,  intent(in) :: AL(3, 3)
   real*8, intent(out) :: fatom(3, natom_mm, nmolm, ntype_m)
   real*8, intent(out) :: Etot

   real*8 fatom_t(3, natom_mm, nmolm, ntype_m)
   real*8 Etot_t

   integer npair, ipair, tnum
   integer iatom_tmp_list(100,20)
   logical is_reset
   integer kkk_mol

   integer itype, jtype, imol, jmol
   integer mol1_type, mol2_type, mol1_Id, mol2_Id
   integer mol1_atom_num
   integer switch_pair

   integer, allocatable, dimension(:) :: iatom_tmp
   real*8, allocatable, dimension(:,:) :: xatom_tmp
   integer, allocatable, dimension(:) :: atom_mol_type, atom_ind_mol

   integer i, ii, jj, iat
   integer nfeat0
   real*8 Rc, Rc_list(5)

   real*8, allocatable, dimension(:,:) :: feat_all
   real*8, allocatable, dimension(:,:,:,:) :: dfeat_all

   integer ierr

   logical need_calc
   real*8 mol_min_dist, mol_max_dist

   character(len=2) :: element(50,5)
   real*8 rx, ry, rz

   element(1:11,1) = (/'O','O','O','C','C','C','H','H','H','H','Li'/)
   element(1:17,2) = (/'O','O','O','C','C','C','H','H','H','H','P','F','F','F','F','F','F'/)
   element(1:8,3)  = (/'P','F','F','F','F','F','F','Li'/)
   element(1:19,4)  = (/'O','O','O','C','C',  'C','C','C','H','H',  'H','H','H','H','H', 'H','H','H','Li'/)
   element(1:25,5)  = (/'O','O','O','C','C',  'C','C','C','H','H',  'H','H','H','H','H', 'H','H','H','P','F',  'F','F','F','F','F'/)


   open(112, file="pair_iatom_list.dat")
   read(112,*) npair
   read(112,*) Rc_list(1:npair)
   do i = 1, npair
      read(112,*) ipair, tnum, iatom_tmp_list(1:tnum, ipair)
   enddo
   close(112)

   is_reset = .true.
   Etot = 0.d0
   fatom = 0.d0
   kkk_mol = 0

   do 1001 itype = 1, ntype_m

      do 1002 imol = 1, nmol(itype)

         kkk_mol = kkk_mol + 1
         if (mod(kkk_mol-1, nnodes_tot) .ne. inode_tot - 1) then
            cycle
         endif

         do 1003 jtype = itype+1, ntype_m

            natom_nn = matom_m(itype) + matom_m(jtype)

            do 1004 jmol = 1, nmol(jtype)

               if ((itype .eq. 1) .and. (jtype .eq. 4)) cycle

               ! write(6,*) "@@@ start PC " , itype, imol, jtype, jmol


               etot_pred_lin = 0.d0
               force_pred_lin = 0.d0

               switch_pair = 0
               if ( (itype .eq. 1) .and. (jtype .eq. 2) ) then
                  ipair = 1
               elseif ( (itype .eq. 1) .and. (jtype .eq. 3) ) then
                  ipair = 2
               elseif ( (itype .eq. 2) .and. (jtype .eq. 3) ) then
                  ipair = 3
                  switch_pair = 1
               elseif ( (itype .eq. 2) .and. (jtype .eq. 4) ) then
                  ipair = 4
                  switch_pair = 1
               elseif ( (itype .eq. 3) .and. (jtype .eq. 4) ) then
                  ipair = 5
                  switch_pair = 1
               endif


               ! set atom type and coordiantion
               ! ensure Li in the last, and PF6 in the last when there is no Li.

               ! ------------------------------
               ! set temporary configuration
               allocate(iatom_tmp(natom_nn))
               allocate(xatom_tmp(3, natom_nn))
               allocate(atom_mol_type(natom_nn))
               allocate(atom_ind_mol(natom_nn))


               Rc = Rc_list(ipair)
               do i = 1, natom_nn
                  iatom_tmp(i) = iatom_tmp_list(i,ipair)
               enddo

               if (switch_pair) then
                  mol1_type = jtype
                  mol2_type = itype
                  mol1_Id = jmol
                  mol2_Id = imol
               else
                  mol1_type = itype
                  mol2_type = jtype
                  mol1_Id = imol
                  mol2_Id = jmol
               endif

               i = 0
               mol1_atom_num = matom_m(mol1_type)
               do ii = 1, matom_m(mol1_type)
                  i = i+1
                  atom_mol_type(i) = mol1_type
                  atom_ind_mol(i) = ii
                  xatom_tmp(1:3,i) = xatom(1:3,ii,mol1_Id,mol1_type)
               enddo
               do ii = 1, matom_m(mol2_type)
                  i = i+1
                  atom_mol_type(i) = mol2_type
                  atom_ind_mol(i) = ii
                  xatom_tmp(1:3,i) = xatom(1:3,ii,mol2_Id,mol2_type)
               enddo


               call move_two_mol_together(matom_m(mol1_type), matom_m(mol2_type), natom_nn, xatom_tmp)


               ! -------------------------------
               ! generature feature
               call load_model_type2(ipair)
               call gen_feature_type2(AL, natom_nn, iatom_tmp, xatom_tmp, &
               & Rc, mol1_type, atom_mol_type, atom_ind_mol, mol1_atom_num)
 


               ! -------------------------------
               need_calc = .False.
               if (ipair .eq. 1) then
                  if (neigh_numR(11) .gt. 1) need_calc = .True.
               elseif (ipair .eq. 2) then
                  do i = 11, natom_nn
                     if (neigh_numR(i) .gt. 7) need_calc = .True.
                  enddo
               elseif (ipair .eq. 3) then
                  if (neigh_numR(8) .gt. 1) need_calc = .True.
               elseif (ipair .eq. 4) then
                  if (neigh_numR(19) .gt. 1) need_calc = .True.
               elseif (ipair .eq. 5) then
                  do i = 19, natom_nn
                     if (neigh_numR(i) .gt. 7) need_calc = .True.
                  enddo   
               endif

               ! write(6,*)
               ! do i = 1, natom_nn
               !    write(6,"(I5,3(1X,F12.9))") iatom_tmp_list(i,ipair), xatom_tmp(1:3,i)
               ! enddo
               ! write(6,*)

               ! open(123, file='pc_pair_all.xyz', position="append")
               ! write(123,"(I5)") natom_nn
               ! write(123,"('Iteration 1')")
               ! do i = 1, natom_nn
               !    rx=AL(1,1)*xatom_tmp(1,i)+AL(1,2)*xatom_tmp(2,i)+AL(1,3)*xatom_tmp(3,i)
               !    ry=AL(2,1)*xatom_tmp(1,i)+AL(2,2)*xatom_tmp(2,i)+AL(2,3)*xatom_tmp(3,i)
               !    rz=AL(3,1)*xatom_tmp(1,i)+AL(3,2)*xatom_tmp(2,i)+AL(3,3)*xatom_tmp(3,i)
               !    write(123,"(A2,3(1X,F12.9))") element(i,ipair), rx, ry, rz
               ! enddo

               ! call mol_atom_distance(natom_nn, xatom_tmp, AL, mol1_atom_num, mol_min_dist, mol_max_dist)
               ! write(6,*) "mol atom distance ", mol_min_dist, mol_max_dist

               if (need_calc) then

                  ! -------------------------------
                  ! combine different feature
                  ! nfeat0m2, nfeat0m3 set in load_model_type2
                  nfeat0 = nfeat0m2 + nfeat0m3

                  allocate(feat_all(nfeat0, natom_nn))
                  allocate(dfeat_all(nfeat0, natom_nn, m_neigh, 3))

                  do iat = 1, natom_nn
                     do ii = 1, nfeat0m2
                        feat_all(ii, iat) = feat2(ii, iat)
                     enddo
                     do ii = 1, nfeat0m3
                        feat_all(ii+nfeat0m2, iat) = feat3(ii, iat)
                     enddo
                  enddo

                  do jj = 1, m_neigh
                     do iat = 1, natom_nn
                        do ii = 1, nfeat0m2
                           dfeat_all(ii, iat, jj, :) = dfeat2(ii,iat,jj,:)
                        enddo
                     enddo
                     do iat = 1, natom_nn
                        do ii = 1, nfeat0m3
                           dfeat_all(ii+nfeat0m2, iat, jj, :) = dfeat3(ii,iat,jj,:)
                        enddo
                     enddo
                  enddo

                  ! -----------------------------
                  ! calculate energy and force
                  call set_paths_lin(ipair)
                  call load_model_lin()
                  call set_image_info_lin(iatom_tmp, is_reset, natom_nn)
                  call cal_energy_force_lin(feat_all, dfeat_all, neigh_numR, neigh_listR, &
                    AL, xatom_tmp, natom_nn, nfeat0, m_neigh, atom_mol_type, Rc)
                  
                  deallocate(feat_all)
                  deallocate(dfeat_all)


                  ! -----------------------------
                  ! sum energy and force
                  Etot = Etot + etot_pred_lin

                  i = 0
                  do ii = 1, matom_m(mol1_type)
                     i = i+1
                     fatom(1:3,ii,mol1_Id,mol1_type) = fatom(1:3,ii,mol1_Id,mol1_type) + force_pred_lin(1:3,i)
                  enddo
                  do ii = 1, matom_m(mol2_type)
                     i = i+1
                     fatom(1:3,ii,mol2_Id,mol2_type) = fatom(1:3,ii,mol2_Id,mol2_type) + force_pred_lin(1:3,i)
                  enddo

                  ! write(6,"('@@@ Etot_pred_lin ', 4(I5,1X), F12.6)") , itype, imol, jtype, jmol, Etot_pred_lin

                  ! write(6,*) "Position, move_x, move_y, move_z"
                  ! do i = 1, natom_nn
                  !    write(6,"(I5,3(1X,F12.9))") iatom_tmp(i), xatom_tmp(1:3,i)
                  ! enddo
                  ! write(6,*) "force x y z"
                  ! do i = 1, natom_nn
                  !    write(6,"(I5,3(1X,F12.9))") iatom_tmp(i), force_pred_lin(1:3,i)
                  ! enddo

                  ! if (abs(Etot_pred_lin) .gt. 0.8) then
                  !    iatom_tmp(1:11) = (/8,8,8,6,6,6,1,1,1,1,3/)
                  !    write(6,*) "11  4 10 1 7 18  1 1 0 0"
                  !    write(6,*) "Lattice vector"
                  !    write(6,*) AL(1,1), "      0.0000000000      0.0000000000"
                  !    write(6,*) "0.0000000000      ", AL(2,2), "     0.0000000000"
                  !    write(6,*) "0.0000000000      0.0000000000      ", AL(3,3)
                  !    write(6,*) "Position, move_x, move_y, move_z"
                  !    do i = 1, natom_nn
                  !       write(6,"(I5,3(1X,F12.9))") iatom_tmp(i), xatom_tmp(1:3,i)
                  !    enddo
                  ! endif

                  open(123, file='pc_pair.xyz', position="append")
                  write(123,"(I5)") natom_nn
                  write(123,"('Iteration 1')")
                  do i = 1, natom_nn
                     rx=AL(1,1)*xatom_tmp(1,i)+AL(1,2)*xatom_tmp(2,i)+AL(1,3)*xatom_tmp(3,i)
                     ry=AL(2,1)*xatom_tmp(1,i)+AL(2,2)*xatom_tmp(2,i)+AL(2,3)*xatom_tmp(3,i)
                     rz=AL(3,1)*xatom_tmp(1,i)+AL(3,2)*xatom_tmp(2,i)+AL(3,3)*xatom_tmp(3,i)
                     write(123,"(A2,3(1X,F12.9))") element(i,ipair), rx, ry, rz
                  enddo

               endif

               deallocate(iatom_tmp, xatom_tmp)
               deallocate(atom_mol_type, atom_ind_mol)

1004        enddo ! itype
1003     enddo ! imol
1002  enddo ! jtype
1001 enddo ! jmol


   call mpi_allreduce(Etot, Etot_t, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
   Etot = Etot_t
   call mpi_allreduce(fatom, fatom_t, 3*natom_mm*nmolm*ntype_m, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
   fatom = fatom_t

end subroutine ML_FF_EF_PC

