!// forquill v1.01 beta www.fcode.cn
module calc_lin_ppp
   use mod_mpi
   use mod_ppp, only: &
       ind_ppp_model_atom,ind_ppp_model_bond, &
       ntype_atom,ntype_bond, &
       m_neigh_atom,m_neigh_bond, &
       itype_atom_atom,itype_atom_bond, &
       nfeat1_atom,nfeat1_bond, &
       nfeat2_atom,nfeat2_bond, &
       nfeat1m_atom,nfeat1m_bond, &
       nfeat2m_atom,nfeat2m_bond, &
       bb_type0_atom,bb_type0_bond, &
       pv_atom,pv_bond, &
       feat2_shift_atom,feat2_shift_bond, &
       feat2_scale_atom,feat2_scale_bond, &
       pv_scale_bb_atom,pv_scale_bb_bond, &
       shift_scale_bb_atom,shift_scale_bb_bond
   use mod_param_densityFF, only: indb
       
   implicit none
   
   interface
       double precision function ddot(n,dx,incx,dy,incy)
           integer :: n,incx,incy
           double precision,dimension(*) :: dx,dy
       end function ddot
   end interface

   !!!!!!!!!!!!!          以下为  module variables     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   character(80), parameter :: fit_input_path = "fit_linearMM.input"
   character(80), parameter :: feat_info_path = "feat.info"
   character(80), parameter :: model_coefficients_path = "linear_fitB.ntype"
   character(80), parameter :: model_bond_path = "bond.molecule"
   character(80), parameter :: weight_feat_path_header = "weight_feat."
   character(80), parameter :: feat_pv_path_header = "feat_PV."
   character(80), parameter :: vdw_path = "vdw_fitB.ntype"

   integer(4) :: ntype                                      !模型所有涉及的原子种类
   integer(4) :: m_neigh                                    !模型所使用的最大近邻数(考虑这个数是否可以不用)
   integer(4) :: nfeat1m                                    !不同种原子的原始feature数目中最大者(目前似无意义)
   integer(4) :: nfeat2m                                    !不同种原子的PCA之后feature数目中最大者
   integer(4) :: nfeat2tot                                  !PCA之后各种原子的feature数目之和
   integer(4), allocatable, dimension(:) :: nfeat1          !各种原子的原始feature数目
   integer(4), allocatable, dimension(:) :: nfeat2          !各种原子PCA之后的feature数目
   integer(4), allocatable, dimension(:) :: nfeat2i         !用来区分计算时各段各属于哪种原子的分段端点序号

   real(8), allocatable, dimension(:) :: bb                 !计算erergy和force时与new feature相乘的系数向量w
   real(8), allocatable, dimension(:, :) :: bb_type         !不明白有何作用,似乎应该是之前用的变量
   real(8), allocatable, dimension(:, :) :: bb_type0        !将bb分别归类到不同种类的原子中，第二维才是代表原子种类

   real(8), allocatable, dimension(:, :, :) :: pv             !PCA所用的转换矩阵
   real(8), allocatable, dimension(:, :) :: feat2_shift     !PCA之后用于标准化feat2的平移矩阵
   real(8), allocatable, dimension(:, :) :: feat2_scale     !PCA之后用于标准化feat2的伸缩系数矩阵

   integer(4) :: natom_nn                                   !image的原子个数
   integer(4) :: natom_mn                                 ! the number of atom in the current molecule
   integer(4) :: nexp_vdw                                 ! the number of atom in the current molecule
   integer(4), allocatable, dimension(:) :: num             !属于每种原子的原子个数，但似乎在calc_linear中无用
   integer(4), allocatable, dimension(:) :: num_atomtype    !属于每种原子的原子个数，似是目前所用的
   integer(4), allocatable, dimension(:) :: itype_atom      !每一种原子的原子属于第几种原子
   integer(4), allocatable, dimension(:) :: iatom           !每种原子的原子序数列表，即atomTypeList
   integer(4), allocatable, dimension(:) :: iatom_type      !每种原子的种类，即序数在种类列表中的序数

   real(8), allocatable, dimension(:) :: energy_pred_lin       !每个原子的能量预测值
   real(8), allocatable, dimension(:) :: energy_pred_tmp        !每个原子的能量预测值
   real(8), allocatable, dimension(:, :) :: force_pred_lin       !每个原子的受力预测值
   real(8), allocatable, dimension(:, :) :: force_pred_tmp       !每个原子的受力预测值
   real(8) :: etot_pred_lin
   character(200) :: error_msg
   integer(4) :: istat
   real(8), allocatable, dimension(:) ::  const_f
   integer(4), allocatable, dimension(:) :: direction, add_force_atom
   integer(4) :: add_force_num, power, axis

   real*8, allocatable, dimension(:) :: rad_atom, E_ave_vdw
   real*8, allocatable, dimension(:, :, :) :: wp_atom
   integer(4) :: nfeat1tm(100), ifeat_type(100), nfeat1t(100)

   real*8, allocatable, dimension(:) :: bond_length, bond_alpha
   real*8, allocatable, dimension(:, :) :: dist_pair1, dist_pair2
   integer, allocatable, dimension(:, :) :: iat_bond
   integer nbond
   real*8 Epair0
   integer iflag_bvp(3)    ! 1: bond, 2: vdw, 3: pair
   
   character(len=20) str

!!!!!!!!!!!!!          以上为  module variables     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains
   subroutine set_image_info_lin_ppp_atom_part(iatom_tmp, is_reset, natom_tmp, &
       ippp,imodel,itype_mol)
       
       integer(4) :: i, j, itype, iitype
       integer iatom_tmp(natom_tmp)
       logical, intent(in) :: is_reset
       integer(4) :: image_size
       integer :: natom_tmp
       integer ippp,imodel,itype_mol
       
       integer ind_ppp_model
       ind_ppp_model=ind_ppp_model_atom(ippp,imodel,itype_mol)
       
       image_size = natom_tmp
       if (is_reset .or. (.not. allocated(iatom)) .or. image_size /= natom_nn) then

           if (allocated(iatom)) then
               if (image_size == natom_nn .and. maxval(abs(iatom_tmp - iatom)) == 0) then
                   return
               end if
               deallocate (iatom)
               deallocate (iatom_type)
               deallocate (energy_pred_lin)
               deallocate (energy_pred_tmp)
               deallocate (force_pred_lin)
               deallocate (force_pred_tmp)
           end if

           natom_nn = image_size
           allocate (iatom(natom_nn))
           allocate (iatom_type(natom_nn))
           allocate (energy_pred_lin(natom_nn))
           allocate (energy_pred_tmp(natom_nn))
           allocate (force_pred_lin(3, natom_nn))
           allocate (force_pred_tmp(3, natom_nn))

           iatom = iatom_tmp

           do i = 1, natom_nn
               iitype = 0
               do itype = 1, ntype_atom(ind_ppp_model,itype_mol)
                   if (itype_atom_atom(itype,ind_ppp_model,itype_mol) == iatom(i)) then
                       iitype = itype
                   end if
               end do
               if (iitype == 0) then
                   write (6, *) 'this type not found', iatom(i)
               end if
               iatom_type(i) = iitype
           end do
       end if

   end subroutine set_image_info_lin_ppp_atom_part


   subroutine set_image_info_lin_ppp_bond_part(iatom_tmp, is_reset, natom_tmp, &
       ippp,imodel,itype_mol)
       
       integer(4) :: i, j, itype, iitype
       integer iatom_tmp(natom_tmp)
       logical, intent(in) :: is_reset
       integer(4) :: image_size
       integer :: natom_tmp
       integer ippp,imodel,itype_mol
       
       integer ind_ppp_model
       ind_ppp_model=ind_ppp_model_bond(ippp,imodel,itype_mol)
       
       image_size = natom_tmp
       if (is_reset .or. (.not. allocated(iatom)) .or. image_size /= natom_nn) then

           if (allocated(iatom)) then
               if (image_size == natom_nn .and. maxval(abs(iatom_tmp - iatom)) == 0) then
                   return
               end if
               deallocate (iatom)
               deallocate (iatom_type)
               deallocate (energy_pred_lin)
               deallocate (energy_pred_tmp)
               deallocate (force_pred_lin)
               deallocate (force_pred_tmp)
           end if

           natom_nn = image_size
           allocate (iatom(natom_nn))
           allocate (iatom_type(natom_nn))
           allocate (energy_pred_lin(natom_nn))
           allocate (energy_pred_tmp(natom_nn))
           allocate (force_pred_lin(3, natom_nn))
           allocate (force_pred_tmp(3, natom_nn))

           iatom = iatom_tmp

           do i = 1, natom_nn
               iitype = 0
               do itype = 1, ntype_bond(ind_ppp_model,itype_mol)
                   if (itype_atom_bond(itype,ind_ppp_model,itype_mol) == iatom(i)) then
                       iitype = itype
                   end if
               end do
               if (iitype == 0) then
                   write (6, *) 'this type not found', iatom(i)
               end if
               iatom_type(i) = iitype
           end do
       end if

   end subroutine set_image_info_lin_ppp_bond_part


   subroutine cal_energy_force_lin_ppp_atom_part(feat, dfeat, num_neigh, list_neigh, AL, xatom, natom_tmp, nfeat0_tmp, m_neigh_tmp, &
       ippp,imodel,itype_mol)
       
       integer(4)  :: itype, ixyz, i, j, jj
       integer natom_tmp, nfeat0_tmp, m_neigh_tmp
       real(8) :: sum
       real(8), intent(in) :: feat(nfeat0_tmp, natom_nn)
       real*8, intent(in) :: dfeat(nfeat0_tmp, natom_nn, m_neigh_tmp, 3)
       integer(4), intent(in) :: num_neigh(natom_tmp)
       integer(4), intent(in) :: list_neigh(m_neigh_tmp, natom_tmp)
       real(8), intent(in) :: AL(3, 3)
       real(8), dimension(:, :), intent(in) :: xatom
       integer ippp,imodel,itype_mol

       real(8), allocatable, dimension(:, :) :: feat2
       real(8), allocatable, dimension(:, :, :) :: feat_type
       real(8), allocatable, dimension(:, :, :) :: feat2_type
       integer(4), allocatable, dimension(:, :) :: ind_type
       real(8), allocatable, dimension(:, :, :) :: dfeat_type
       real(8), allocatable, dimension(:, :, :) :: dfeat2_type
       real(8), allocatable, dimension(:, :, :, :) :: dfeat2

       real*8 pi, dE, dFx, dFy, dFz, yy, Ebond
       real*8 rad1, rad2, rad, dx1, dx2, dx3, dx, dy, dz, dd, w22, dEdd, d, w22_1, w22_2, w22F_1, w22F_2
       integer iat1, iat2, ierr, ibond

       real*8 Epair
       
       integer iatom_nonzero_feat
       integer itype_nonzero_feat
       integer counter
       integer,allocatable,dimension(:) :: natom_per_type
       real*8,allocatable,dimension(:) :: force_pred_type
       
       integer ind_ppp_model
       ind_ppp_model=ind_ppp_model_atom(ippp,imodel,itype_mol)
       
       pi = 4*datan(1.d0)

       ! allocate (num(ntype_atom(ind_ppp_model,itype_mol)))
       ! allocate (feat2(nfeat2m_atom(ind_ppp_model,itype_mol), natom_nn))
       ! allocate (feat_type(nfeat1m_atom(ind_ppp_model,itype_mol), natom_nn, ntype_atom(ind_ppp_model,itype_mol)))
       ! allocate (feat2_type(nfeat2m_atom(ind_ppp_model,itype_mol), natom_nn, ntype_atom(ind_ppp_model,itype_mol)))
       ! allocate (ind_type(natom_nn, ntype_atom(ind_ppp_model,itype_mol)))
       allocate (dfeat_type(nfeat1m_atom(ind_ppp_model,itype_mol), natom_nn*m_neigh_atom(ind_ppp_model,itype_mol)*3, ntype_atom(ind_ppp_model,itype_mol)))
       ! allocate (dfeat2_type(nfeat2m_atom(ind_ppp_model,itype_mol), natom_nn*m_neigh_atom(ind_ppp_model,itype_mol)*3, ntype_atom(ind_ppp_model,itype_mol)))
       ! allocate (dfeat2(nfeat2m_atom(ind_ppp_model,itype_mol), natom_nn, m_neigh_atom(ind_ppp_model,itype_mol), 3))
       allocate (force_pred_type(natom_nn*3))

       istat = 0
       error_msg = ''

       if (nfeat0_tmp /= nfeat1m_atom(ind_ppp_model,itype_mol) .or. natom_tmp /= natom_nn .or. m_neigh_tmp /= m_neigh_atom(ind_ppp_model,itype_mol)) then
           write (*, *) "Shape of input arrays don't match the model!"
           stop
       end if
       
       ! ---------- energy part (old approach) ----------
       ! num = 0
       ! iat1 = 0
       ! do i = 1, natom_nn
           ! if (mod(i - 1, nnodes) .eq. inode - 1) then
               ! iat1 = iat1 + 1
               ! itype = iatom_type(i)
               ! num(itype) = num(itype) + 1
               ! ind_type(num(itype), itype) = iat1
               ! feat_type(:, num(itype), itype) = feat(:, iat1)
           ! end if
       ! end do
       
       ! do itype = 1, ntype_atom(ind_ppp_model,itype_mol)
           ! call dgemm('T','N',nfeat2_atom(itype,ind_ppp_model,itype_mol),num(itype),nfeat1_atom(itype,ind_ppp_model,itype_mol),1.d0, &
               ! pv_atom(1,1,itype,ind_ppp_model,itype_mol),nfeat1m_atom(ind_ppp_model,itype_mol),feat_type(1,1,itype),nfeat1m_atom(ind_ppp_model,itype_mol),0.d0, &
               ! feat2_type(1,1,itype),nfeat2m_atom(ind_ppp_model,itype_mol))
       ! end do

       ! do itype = 1, ntype_atom(ind_ppp_model,itype_mol)
           ! do i = 1, num(itype)
               ! do j = 1, nfeat2_atom(itype,ind_ppp_model,itype_mol) - 1
                   ! feat2_type(j, i, itype) = (feat2_type(j, i, itype) - feat2_shift_atom(j,itype,ind_ppp_model,itype_mol))*feat2_scale_atom(j,itype,ind_ppp_model,itype_mol)
               ! end do
               ! if(nfeat2_atom(itype,ind_ppp_model,itype_mol).ne.0) feat2_type(nfeat2_atom(itype,ind_ppp_model,itype_mol), i, itype) = 1.d0
           ! end do
       ! end do

       ! num = 0
       ! iat1 = 0
       ! do i = 1, natom_nn
           ! if (mod(i - 1, nnodes) .eq. inode - 1) then
               ! iat1 = iat1 + 1
               ! itype = iatom_type(i)
               ! num(itype) = num(itype) + 1
               ! feat2(:, iat1) = feat2_type(:, num(itype), itype)
           ! end if
       ! end do

       ! energy_pred_tmp = 0.d0

       ! iat1 = 0
       ! do i = 1, natom_nn
           ! if (mod(i - 1, nnodes) .eq. inode - 1) then
               ! iat1 = iat1 + 1
               ! itype = iatom_type(i)
               ! sum = 0.d0
               ! do j = 1, nfeat2_atom(itype,ind_ppp_model,itype_mol)
                   ! sum = sum + feat2(j, iat1)*bb_type0_atom(j,itype,ind_ppp_model,itype_mol)
               ! end do
               ! energy_pred_tmp(i) = sum
           ! end if
       ! end do
       
       ! ---------- energy part (new approach) ----------
       allocate(natom_per_type(ntype_atom(ind_ppp_model,itype_mol)))
       natom_per_type=0
       do i=1,natom_nn
           itype=iatom_type(i)
           natom_per_type(itype)=natom_per_type(itype)+1
       enddo
       
       do i=1,ntype_atom(ind_ppp_model,itype_mol)
           if(nfeat2_atom(i,ind_ppp_model,itype_mol).ne.0) itype_nonzero_feat=i
       enddo
       
       iatom_nonzero_feat=imodel
       
       if(iatom_type(iatom_nonzero_feat).ne.itype_nonzero_feat) then
           write(6,*) 'iatom_type(iatom_nonzero_feat).ne.itype_nonzero_feat, something goes wrong! stop',iatom_type(iatom_nonzero_feat),itype_nonzero_feat
           stop
       endif
       
       etot_pred_lin=ddot(nfeat1_atom(itype_nonzero_feat,ind_ppp_model,itype_mol),feat(1,iatom_nonzero_feat),1,pv_scale_bb_atom(1,itype_nonzero_feat,ind_ppp_model,itype_mol),1)
       etot_pred_lin=etot_pred_lin+(-shift_scale_bb_atom(itype_nonzero_feat,ind_ppp_model,itype_mol)+bb_type0_atom(nfeat2_atom(itype_nonzero_feat,ind_ppp_model,itype_mol),itype_nonzero_feat,ind_ppp_model,itype_mol))*natom_per_type(itype_nonzero_feat)
       deallocate(natom_per_type)
       
       ! ---------- force part (old approach) ----------
       ! num = 0
       ! iat1 = 0
       ! do i = 1, natom_nn
           ! if (mod(i - 1, nnodes) .eq. inode - 1) then
               ! iat1 = iat1 + 1
               ! itype = iatom_type(i)
               ! do jj = 1, num_neigh(i)
                   ! num(itype) = num(itype) + 1
                   ! dfeat_type(:, num(itype), itype) = dfeat(:, iat1, jj, 1)
                   ! num(itype) = num(itype) + 1
                   ! dfeat_type(:, num(itype), itype) = dfeat(:, iat1, jj, 2)
                   ! num(itype) = num(itype) + 1
                   ! dfeat_type(:, num(itype), itype) = dfeat(:, iat1, jj, 3)
               ! end do
           ! end if
       ! end do
       ! !cccccccc note: num(itype) is rather large, in the scane of natom*num_neigh

       ! do itype = 1, ntype_atom(ind_ppp_model,itype_mol)
           ! call dgemm('T', 'N', nfeat2_atom(itype,ind_ppp_model,itype_mol), num(itype), nfeat1_atom(itype,ind_ppp_model,itype_mol), 1.d0, &
               ! pv_atom(1,1,itype,ind_ppp_model,itype_mol), nfeat1m_atom(ind_ppp_model,itype_mol), dfeat_type(1,1,itype), nfeat1m_atom(ind_ppp_model,itype_mol), 0.d0, &
               ! dfeat2_type(1,1,itype), nfeat2m_atom(ind_ppp_model,itype_mol))
       ! end do

       ! num = 0
       ! iat1 = 0
       ! do i = 1, natom_nn
           ! if (mod(i - 1, nnodes) .eq. inode - 1) then
               ! iat1 = iat1 + 1
               ! itype = iatom_type(i)
               ! do jj = 1, num_neigh(i)
                   ! num(itype) = num(itype) + 1
                   ! do j = 1, nfeat2_atom(itype,ind_ppp_model,itype_mol) - 1
                       ! dfeat2(j, iat1, jj, 1) = dfeat2_type(j, num(itype), itype)*feat2_scale_atom(j,itype,ind_ppp_model,itype_mol)
                   ! end do
                   ! if(nfeat2_atom(itype,ind_ppp_model,itype_mol).ne.0) dfeat2(nfeat2_atom(itype,ind_ppp_model,itype_mol), iat1, jj, 1) = 0.d0
                   ! num(itype) = num(itype) + 1
                   ! do j = 1, nfeat2_atom(itype,ind_ppp_model,itype_mol) - 1
                       ! dfeat2(j, iat1, jj, 2) = dfeat2_type(j, num(itype), itype)*feat2_scale_atom(j,itype,ind_ppp_model,itype_mol)
                   ! end do
                   ! if(nfeat2_atom(itype,ind_ppp_model,itype_mol).ne.0) dfeat2(nfeat2_atom(itype,ind_ppp_model,itype_mol), iat1, jj, 2) = 0.d0
                   ! num(itype) = num(itype) + 1
                   ! do j = 1, nfeat2_atom(itype,ind_ppp_model,itype_mol) - 1
                       ! dfeat2(j, iat1, jj, 3) = dfeat2_type(j, num(itype), itype)*feat2_scale_atom(j,itype,ind_ppp_model,itype_mol)
                   ! end do
                   ! if(nfeat2_atom(itype,ind_ppp_model,itype_mol).ne.0) dfeat2(nfeat2_atom(itype,ind_ppp_model,itype_mol), iat1, jj, 3) = 0.d0
               ! end do
           ! end if
       ! end do

       ! !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       ! !cc  the new dfeat2 is:
       ! !cc dfeat2(nfeat2,natom,j_neigh,3): dfeat2(j,i,jj,3)= d/dr(jj_neigh)(feat2(j,i))
       ! !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       ! !cccc now, we have the new features, we need to calculate the distance to reference state

       ! force_pred_tmp = 0.d0

       ! iat1 = 0
       ! do i = 1, natom_nn
           ! if (mod(i - 1, nnodes) .eq. inode - 1) then
               ! iat1 = iat1 + 1
               ! itype = iatom_type(i)
               ! do jj = 1, num_neigh(i)
                   ! iat2 = list_neigh(jj, i)

                   ! do j = 1, nfeat2_atom(itype,ind_ppp_model,itype_mol)
                       ! force_pred_tmp(1, iat2) = force_pred_tmp(1, iat2) + dfeat2(j, iat1, jj, 1)*bb_type0_atom(j,itype,ind_ppp_model,itype_mol)
                       ! force_pred_tmp(2, iat2) = force_pred_tmp(2, iat2) + dfeat2(j, iat1, jj, 2)*bb_type0_atom(j,itype,ind_ppp_model,itype_mol)
                       ! force_pred_tmp(3, iat2) = force_pred_tmp(3, iat2) + dfeat2(j, iat1, jj, 3)*bb_type0_atom(j,itype,ind_ppp_model,itype_mol)
                   ! end do
               ! end do
           ! end if
       ! end do
       
       ! ---------- force part (new approach) ----------
       force_pred_type=0.d0
       force_pred_lin=0.d0
       
       counter=0
       do jj=1,num_neigh(iatom_nonzero_feat)
           counter=counter+1
           dfeat_type(:,counter,itype_nonzero_feat)=dfeat(:,iatom_nonzero_feat,jj,1)
           counter=counter+1
           dfeat_type(:,counter,itype_nonzero_feat)=dfeat(:,iatom_nonzero_feat,jj,2)
           counter=counter+1
           dfeat_type(:,counter,itype_nonzero_feat)=dfeat(:,iatom_nonzero_feat,jj,3)
       enddo
       
       call dgemv('T',nfeat1_atom(itype_nonzero_feat,ind_ppp_model,itype_mol),counter,1.d0, &
           dfeat_type(1,1,itype_nonzero_feat),nfeat1m_atom(ind_ppp_model,itype_mol),pv_scale_bb_atom(1,itype_nonzero_feat,ind_ppp_model,itype_mol),1,0.d0, &
           force_pred_type(1),1)
       
       counter=0
       do jj=1,num_neigh(iatom_nonzero_feat)
           counter=counter+1
           force_pred_lin(1,list_neigh(jj,iatom_nonzero_feat))=force_pred_type(counter)
           counter=counter+1
           force_pred_lin(2,list_neigh(jj,iatom_nonzero_feat))=force_pred_type(counter)
           counter=counter+1
           force_pred_lin(3,list_neigh(jj,iatom_nonzero_feat))=force_pred_type(counter)
       enddo
       
       !ccccccccccccccccccccccccccccccccccccccccccc

       ! call mpi_allreduce(energy_pred_tmp, energy_pred_lin, natom_nn, MPI_REAL8, MPI_SUM, MPI_COMM_MOL, ierr)
       ! call mpi_allreduce(force_pred_tmp, force_pred_lin, 3*natom_nn, MPI_REAL8, MPI_SUM, MPI_COMM_MOL, ierr)

       !ccccccccccccccccccccccccccccccccccccccccccc

       ! etot_pred_lin = 0.d0
       ! do i = 1, natom_nn
           ! etot_pred_lin = etot_pred_lin + energy_pred_lin(i)
       ! end do

       ! deallocate (num)
       ! deallocate (feat2)
       ! deallocate (feat_type)
       ! deallocate (feat2_type)
       ! deallocate (ind_type)
       deallocate (dfeat_type)
       ! deallocate (dfeat2_type)
       ! deallocate (dfeat2)
       deallocate (force_pred_type)
   end subroutine cal_energy_force_lin_ppp_atom_part


   subroutine cal_energy_force_lin_ppp_bond_part(feat, dfeat, num_neigh, list_neigh, AL, xatom, natom_tmp, nfeat0_tmp, m_neigh_tmp, &
       ippp,imodel,itype_mol)
       
       integer(4)  :: itype, ixyz, i, j, jj
       integer natom_tmp, nfeat0_tmp, m_neigh_tmp
       real(8) :: sum
       real(8), intent(in) :: feat(nfeat0_tmp, natom_nn)
       real*8, intent(in) :: dfeat(nfeat0_tmp, natom_nn, m_neigh_tmp, 3)
       integer(4), intent(in) :: num_neigh(natom_tmp)
       integer(4), intent(in) :: list_neigh(m_neigh_tmp, natom_tmp)
       real(8), intent(in) :: AL(3, 3)
       real(8), dimension(:, :), intent(in) :: xatom
       integer ippp,imodel,itype_mol

       real(8), allocatable, dimension(:, :) :: feat2
       real(8), allocatable, dimension(:, :, :) :: feat_type
       real(8), allocatable, dimension(:, :, :) :: feat2_type
       integer(4), allocatable, dimension(:, :) :: ind_type
       real(8), allocatable, dimension(:, :, :) :: dfeat_type
       real(8), allocatable, dimension(:, :, :) :: dfeat2_type
       real(8), allocatable, dimension(:, :, :, :) :: dfeat2

       real*8 pi, dE, dFx, dFy, dFz, yy, Ebond
       real*8 rad1, rad2, rad, dx1, dx2, dx3, dx, dy, dz, dd, w22, dEdd, d, w22_1, w22_2, w22F_1, w22F_2
       integer iat1, iat2, ierr, ibond

       real*8 Epair
       
       integer iatom_nonzero_feat
       integer itype_nonzero_feat
       integer counter
       integer,allocatable,dimension(:) :: natom_per_type
       real*8,allocatable,dimension(:) :: force_pred_type
       
       integer ind_ppp_model
       ind_ppp_model=ind_ppp_model_bond(ippp,imodel,itype_mol)
       
       pi = 4*datan(1.d0)

       ! allocate (num(ntype_bond(ind_ppp_model,itype_mol)))
       ! allocate (feat2(nfeat2m_bond(ind_ppp_model,itype_mol), natom_nn))
       ! allocate (feat_type(nfeat1m_bond(ind_ppp_model,itype_mol), natom_nn, ntype_bond(ind_ppp_model,itype_mol)))
       ! allocate (feat2_type(nfeat2m_bond(ind_ppp_model,itype_mol), natom_nn, ntype_bond(ind_ppp_model,itype_mol)))
       ! allocate (ind_type(natom_nn, ntype_bond(ind_ppp_model,itype_mol)))
       allocate (dfeat_type(nfeat1m_bond(ind_ppp_model,itype_mol), natom_nn*m_neigh_bond(ind_ppp_model,itype_mol)*3, ntype_bond(ind_ppp_model,itype_mol)))
       ! allocate (dfeat2_type(nfeat2m_bond(ind_ppp_model,itype_mol), natom_nn*m_neigh_bond(ind_ppp_model,itype_mol)*3, ntype_bond(ind_ppp_model,itype_mol)))
       ! allocate (dfeat2(nfeat2m_bond(ind_ppp_model,itype_mol), natom_nn, m_neigh_bond(ind_ppp_model,itype_mol), 3))
       allocate (force_pred_type(natom_nn*3))

       istat = 0
       error_msg = ''

       if (nfeat0_tmp /= nfeat1m_bond(ind_ppp_model,itype_mol) .or. natom_tmp /= natom_nn .or. m_neigh_tmp /= m_neigh_bond(ind_ppp_model,itype_mol)) then
           write (*, *) "Shape of input arrays don't match the model!"
           stop
       end if
       
       ! ---------- energy part (old approach) ----------
       ! num = 0
       ! iat1 = 0
       ! do i = 1, natom_nn
           ! if (mod(i - 1, nnodes) .eq. inode - 1) then
               ! iat1 = iat1 + 1
               ! itype = iatom_type(i)
               ! num(itype) = num(itype) + 1
               ! ind_type(num(itype), itype) = iat1
               ! feat_type(:, num(itype), itype) = feat(:, iat1)
           ! end if
       ! end do
       
       ! do itype = 1, ntype_bond(ind_ppp_model,itype_mol)
           ! call dgemm('T', 'N', nfeat2_bond(itype,ind_ppp_model,itype_mol), num(itype), nfeat1_bond(itype,ind_ppp_model,itype_mol), 1.d0, pv_bond(1,1,itype,ind_ppp_model,itype_mol), nfeat1m_bond(ind_ppp_model,itype_mol), feat_type(1,1,itype), nfeat1m_bond(ind_ppp_model,itype_mol), 0.d0,feat2_type(1,1,itype), nfeat2m_bond(ind_ppp_model,itype_mol))
       ! end do

       ! do itype = 1, ntype_bond(ind_ppp_model,itype_mol)
           ! do i = 1, num(itype)
               ! do j = 1, nfeat2_bond(itype,ind_ppp_model,itype_mol) - 1
                   ! feat2_type(j, i, itype) = (feat2_type(j, i, itype) - feat2_shift_bond(j,itype,ind_ppp_model,itype_mol))*feat2_scale_bond(j,itype,ind_ppp_model,itype_mol)
               ! end do
               ! if(nfeat2_bond(itype,ind_ppp_model,itype_mol).ne.0) feat2_type(nfeat2_bond(itype,ind_ppp_model,itype_mol), i, itype) = 1.d0
           ! end do
       ! end do

       ! num = 0
       ! iat1 = 0
       ! do i = 1, natom_nn
           ! if (mod(i - 1, nnodes) .eq. inode - 1) then
               ! iat1 = iat1 + 1
               ! itype = iatom_type(i)
               ! num(itype) = num(itype) + 1
               ! feat2(:, iat1) = feat2_type(:, num(itype), itype)
           ! end if
       ! end do

       ! energy_pred_tmp = 0.d0

       ! iat1 = 0
       ! do i = 1, natom_nn
           ! if (mod(i - 1, nnodes) .eq. inode - 1) then
               ! iat1 = iat1 + 1
               ! itype = iatom_type(i)
               ! sum = 0.d0
               ! do j = 1, nfeat2_bond(itype,ind_ppp_model,itype_mol)
                   ! sum = sum + feat2(j, iat1)*bb_type0_bond(j,itype,ind_ppp_model,itype_mol)
               ! end do
               ! energy_pred_tmp(i) = sum
           ! end if
       ! end do
       
       ! ---------- energy part (new approach) ----------
       allocate(natom_per_type(ntype_bond(ind_ppp_model,itype_mol)))
       natom_per_type=0
       do i=1,natom_nn
           itype=iatom_type(i)
           natom_per_type(itype)=natom_per_type(itype)+1
       enddo
       
       do i=1,ntype_bond(ind_ppp_model,itype_mol)
           if(nfeat2_bond(i,ind_ppp_model,itype_mol).ne.0) itype_nonzero_feat=i
       enddo
       
       if((ippp.eq.1).or.(ippp.eq.2)) iatom_nonzero_feat=indb(1,imodel,itype_mol)
       if((ippp.eq.3).or.(ippp.eq.4)) iatom_nonzero_feat=indb(2,imodel,itype_mol)
       
       if(iatom_type(iatom_nonzero_feat).ne.itype_nonzero_feat) then
           write(6,*) 'iatom_type(iatom_nonzero_feat).ne.itype_nonzero_feat, something goes wrong! stop',iatom_type(iatom_nonzero_feat),itype_nonzero_feat
           stop
       endif
       
       etot_pred_lin=ddot(nfeat1_bond(itype_nonzero_feat,ind_ppp_model,itype_mol),feat(1,iatom_nonzero_feat),1,pv_scale_bb_bond(1,itype_nonzero_feat,ind_ppp_model,itype_mol),1)
       etot_pred_lin=etot_pred_lin+(-shift_scale_bb_bond(itype_nonzero_feat,ind_ppp_model,itype_mol)+bb_type0_bond(nfeat2_bond(itype_nonzero_feat,ind_ppp_model,itype_mol),itype_nonzero_feat,ind_ppp_model,itype_mol))*natom_per_type(itype_nonzero_feat)
       deallocate(natom_per_type)
       
       ! ---------- force part (old approach) ----------
       ! num = 0
       ! iat1 = 0
       ! do i = 1, natom_nn
           ! if (mod(i - 1, nnodes) .eq. inode - 1) then
               ! iat1 = iat1 + 1
               ! itype = iatom_type(i)
               ! do jj = 1, num_neigh(i)
                   ! num(itype) = num(itype) + 1
                   ! dfeat_type(:, num(itype), itype) = dfeat(:, iat1, jj, 1)
                   ! num(itype) = num(itype) + 1
                   ! dfeat_type(:, num(itype), itype) = dfeat(:, iat1, jj, 2)
                   ! num(itype) = num(itype) + 1
                   ! dfeat_type(:, num(itype), itype) = dfeat(:, iat1, jj, 3)
               ! end do
           ! end if
       ! end do
       ! !cccccccc note: num(itype) is rather large, in the scane of natom*num_neigh

       ! do itype = 1, ntype_bond(ind_ppp_model,itype_mol)
           ! call dgemm('T', 'N', nfeat2_bond(itype,ind_ppp_model,itype_mol), num(itype), nfeat1_bond(itype,ind_ppp_model,itype_mol), 1.d0, pv_bond(1,1,itype,ind_ppp_model,itype_mol), nfeat1m_bond(ind_ppp_model,itype_mol), dfeat_type(1,1,itype), nfeat1m_bond(ind_ppp_model,itype_mol), 0.d0, dfeat2_type(1,1,itype), nfeat2m_bond(ind_ppp_model,itype_mol))
       ! end do

       ! num = 0
       ! iat1 = 0
       ! do i = 1, natom_nn
           ! if (mod(i - 1, nnodes) .eq. inode - 1) then
               ! iat1 = iat1 + 1
               ! itype = iatom_type(i)
               ! do jj = 1, num_neigh(i)
                   ! num(itype) = num(itype) + 1
                   ! do j = 1, nfeat2_bond(itype,ind_ppp_model,itype_mol) - 1
                       ! dfeat2(j, iat1, jj, 1) = dfeat2_type(j, num(itype), itype)*feat2_scale_bond(j,itype,ind_ppp_model,itype_mol)
                   ! end do
                   ! if(nfeat2_bond(itype,ind_ppp_model,itype_mol).ne.0) dfeat2(nfeat2_bond(itype,ind_ppp_model,itype_mol), iat1, jj, 1) = 0.d0
                   ! num(itype) = num(itype) + 1
                   ! do j = 1, nfeat2_bond(itype,ind_ppp_model,itype_mol) - 1
                       ! dfeat2(j, iat1, jj, 2) = dfeat2_type(j, num(itype), itype)*feat2_scale_bond(j,itype,ind_ppp_model,itype_mol)
                   ! end do
                   ! if(nfeat2_bond(itype,ind_ppp_model,itype_mol).ne.0) dfeat2(nfeat2_bond(itype,ind_ppp_model,itype_mol), iat1, jj, 2) = 0.d0
                   ! num(itype) = num(itype) + 1
                   ! do j = 1, nfeat2_bond(itype,ind_ppp_model,itype_mol) - 1
                       ! dfeat2(j, iat1, jj, 3) = dfeat2_type(j, num(itype), itype)*feat2_scale_bond(j,itype,ind_ppp_model,itype_mol)
                   ! end do
                   ! if(nfeat2_bond(itype,ind_ppp_model,itype_mol).ne.0) dfeat2(nfeat2_bond(itype,ind_ppp_model,itype_mol), iat1, jj, 3) = 0.d0
               ! end do
           ! end if
       ! end do

       ! !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       ! !cc  the new dfeat2 is:
       ! !cc dfeat2(nfeat2,natom,j_neigh,3): dfeat2(j,i,jj,3)= d/dr(jj_neigh)(feat2(j,i))
       ! !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       ! !cccc now, we have the new features, we need to calculate the distance to reference state

       ! force_pred_tmp = 0.d0

       ! iat1 = 0
       ! do i = 1, natom_nn
           ! if (mod(i - 1, nnodes) .eq. inode - 1) then
               ! iat1 = iat1 + 1
               ! itype = iatom_type(i)
               ! do jj = 1, num_neigh(i)
                   ! iat2 = list_neigh(jj, i)

                   ! do j = 1, nfeat2_bond(itype,ind_ppp_model,itype_mol)
                       ! force_pred_tmp(1, iat2) = force_pred_tmp(1, iat2) + dfeat2(j, iat1, jj, 1)*bb_type0_bond(j,itype,ind_ppp_model,itype_mol)
                       ! force_pred_tmp(2, iat2) = force_pred_tmp(2, iat2) + dfeat2(j, iat1, jj, 2)*bb_type0_bond(j,itype,ind_ppp_model,itype_mol)
                       ! force_pred_tmp(3, iat2) = force_pred_tmp(3, iat2) + dfeat2(j, iat1, jj, 3)*bb_type0_bond(j,itype,ind_ppp_model,itype_mol)
                   ! end do
               ! end do
           ! end if
       ! end do
       
       ! ---------- force part (new approach) ----------
       force_pred_type=0.d0
       force_pred_lin=0.d0
       
       counter=0
       do jj=1,num_neigh(iatom_nonzero_feat)
           counter=counter+1
           dfeat_type(:,counter,itype_nonzero_feat)=dfeat(:,iatom_nonzero_feat,jj,1)
           counter=counter+1
           dfeat_type(:,counter,itype_nonzero_feat)=dfeat(:,iatom_nonzero_feat,jj,2)
           counter=counter+1
           dfeat_type(:,counter,itype_nonzero_feat)=dfeat(:,iatom_nonzero_feat,jj,3)
       enddo
       
       call dgemv('T',nfeat1_bond(itype_nonzero_feat,ind_ppp_model,itype_mol),counter,1.d0, &
           dfeat_type(1,1,itype_nonzero_feat),nfeat1m_bond(ind_ppp_model,itype_mol),pv_scale_bb_bond(1,itype_nonzero_feat,ind_ppp_model,itype_mol),1,0.d0, &
           force_pred_type(1),1)
       
       counter=0
       do jj=1,num_neigh(iatom_nonzero_feat)
           counter=counter+1
           force_pred_lin(1,list_neigh(jj,iatom_nonzero_feat))=force_pred_type(counter)
           counter=counter+1
           force_pred_lin(2,list_neigh(jj,iatom_nonzero_feat))=force_pred_type(counter)
           counter=counter+1
           force_pred_lin(3,list_neigh(jj,iatom_nonzero_feat))=force_pred_type(counter)
       enddo
       
       !ccccccccccccccccccccccccccccccccccccccccccc

       ! call mpi_allreduce(energy_pred_tmp, energy_pred_lin, natom_nn, MPI_REAL8, MPI_SUM, MPI_COMM_MOL, ierr)
       ! call mpi_allreduce(force_pred_tmp, force_pred_lin, 3*natom_nn, MPI_REAL8, MPI_SUM, MPI_COMM_MOL, ierr)

       !ccccccccccccccccccccccccccccccccccccccccccc

       ! etot_pred_lin = 0.d0
       ! do i = 1, natom_nn
           ! etot_pred_lin = etot_pred_lin + energy_pred_lin(i)
       ! end do

       ! deallocate (num)
       ! deallocate (feat2)
       ! deallocate (feat_type)
       ! deallocate (feat2_type)
       ! deallocate (ind_type)
       deallocate (dfeat_type)
       ! deallocate (dfeat2_type)
       ! deallocate (dfeat2)
       deallocate (force_pred_type)
   end subroutine cal_energy_force_lin_ppp_bond_part

end module calc_lin_ppp
