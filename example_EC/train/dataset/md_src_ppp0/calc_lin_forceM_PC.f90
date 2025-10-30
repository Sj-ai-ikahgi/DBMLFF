!// forquill v1.01 beta www.fcode.cn
module calc_lin_PC
   implicit none

   !!!!!!!!!!!!!          以下为  module variables     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   character(80), parameter :: fit_input_path0 = "fread_dfeat/fit_linearMM.input"
   character(80), parameter :: feat_info_path0 = "fread_dfeat/feat.info"
   character(80), parameter :: model_coefficients_path0 = "fread_dfeat/linear_fitB.ntype"
   character(80), parameter :: weight_feat_path_header0 = "fread_dfeat/weight_feat."
   character(80), parameter :: feat_pv_path_header0 = "fread_dfeat/feat_PV."
   character(80), parameter :: vdw_path0="fread_dfeat/vdw_fitB.ntype"


   character(200) :: fit_input_path = trim(fit_input_path0)
   character(200) :: feat_info_path = trim(feat_info_path0)
   character(200) :: model_coefficients_path = trim(model_coefficients_path0)
   character(200) :: weight_feat_path_header = trim(weight_feat_path_header0)
   character(200) :: feat_pv_path_header = trim(feat_pv_path_header0)
   character(200) :: vdw_path=trim(vdw_path0)

   integer(4) :: ntype                                    !模型所有涉及的原子种类
   integer(4) :: m_neigh                                  !模型所使用的最大近邻数(考虑这个数是否可以不用)
   integer(4) :: nfeat1m                                  !不同种原子的原始feature数目中最大者(目前似无意义)
   integer(4) :: nfeat2m                                  !不同种原子的PCA之后feature数目中最大者
   integer(4) :: nfeat2tot                                !PCA之后各种原子的feature数目之和
   integer(4), allocatable, dimension(:) :: nfeat1          !各种原子的原始feature数目
   integer(4), allocatable, dimension(:) :: nfeat2          !各种原子PCA之后的feature数目
   integer(4), allocatable, dimension(:) :: nfeat2i         !用来区分计算时各段各属于哪种原子的分段端点序号

   real(8), allocatable, dimension(:) :: bb                 !计算erergy和force时与new feature相乘的系数向量w
   real(8), allocatable, dimension(:, :) :: bb_type0         !将bb分别归类到不同种类的原子中，第二维才是代表原子种类

   real(8), allocatable, dimension(:, :, :) :: pv             !PCA所用的转换矩阵
   real(8), allocatable, dimension(:, :) :: feat2_shift     !PCA之后用于标准化feat2的平移矩阵
   real(8), allocatable, dimension(:, :) :: feat2_scale     !PCA之后用于标准化feat2的伸缩系数矩阵

   real*8,allocatable,dimension(:) :: rad_atom
   real*8,allocatable,dimension(:,:,:) :: wp_atom

   integer(4) :: natom_nn                                    !image的原子个数
   integer(4), allocatable, dimension(:) :: itype_atom      !每一种原子的原子属于第几种原子
   integer(4), allocatable, dimension(:) :: iatom           !每种原子的原子序数列表，即atomTypeList
   integer(4), allocatable, dimension(:) :: iatom_type      !每种原子的种类，即序数在种类列表中的序数

   integer(4) :: ifeat_type_l(100)
   integer(4) :: nfeat_type_l

   real(8), allocatable, dimension(:) :: energy_pred_lin       !每个原子的能量预测值
   real(8), allocatable, dimension(:, :) :: force_pred_lin       !每个原子的受力预测值
   real(8) :: etot_pred_lin

contains

   subroutine set_paths_lin(ipair_type)
      integer ipair_type

      fit_input_path = 'pair.' // char(48 + ipair_type)//'/'//trim(fit_input_path0)
      feat_info_path = 'pair.' // char(48 + ipair_type)//'/'//trim(feat_info_path0)
      model_coefficients_path = 'pair.' // char(48 + ipair_type)//'/'//trim(model_coefficients_path0)
      weight_feat_path_header = 'pair.' // char(48 + ipair_type)//'/'//trim(weight_feat_path_header0)
      feat_pv_path_header = 'pair.' // char(48 + ipair_type)//'/'//trim(feat_pv_path_header0)
      vdw_path='pair.' // char(48 + ipair_type)//'/'//trim(vdw_path0)
   end subroutine set_paths_lin

   subroutine load_model_lin()

      integer(4) :: nfeat1_tmp, nfeat2_tmp, itype, i, k, ntmp, itmp
      integer(4) :: iflag_PCA, kkk, ntype_tmp, iatom_tmp
      integer itype_t, ntype_t, nterm, j1, itype1
      real*8 E_ave_vdw

      ! **************** read fit_linearMM.input ********************
      open (10, file=trim(fit_input_path))
      rewind (10)
      read (10, *) ntype, m_neigh

      if (allocated(itype_atom)) deallocate (itype_atom)
      if (allocated(nfeat1)) deallocate (nfeat1)
      if (allocated(nfeat2)) deallocate (nfeat2)
      if (allocated(nfeat2i)) deallocate (nfeat2i)
      if (allocated(bb)) deallocate (bb)
      if (allocated(bb_type0)) deallocate (bb_type0)
      if (allocated(pv)) deallocate (pv)
      if (allocated(feat2_shift)) deallocate (feat2_shift)
      if (allocated(feat2_scale)) deallocate (feat2_scale)

      allocate (itype_atom(ntype))
      allocate (nfeat1(ntype))
      allocate (nfeat2(ntype))
      allocate (nfeat2i(ntype))

      do i = 1, ntype
         read (10, *) itype_atom(i)!,rad_atom(i),wp_atom(i)
      end do
      ! read(10,*) weight_E,weight_E0,weight_F
      close (10)


      if (allocated(rad_atom)) deallocate (rad_atom)
      if (allocated(wp_atom)) deallocate (wp_atom)
      allocate(rad_atom(ntype))
      allocate(wp_atom(ntype,ntype,2))
      wp_atom=0.d0

      ! --- read vdw
      open(10,file=trim(vdw_path))
      rewind(10)
      read(10,*) ntype_t,nterm
      if(nterm.gt.2) then
         write(6,*) "nterm.gt.2,stop"
         stop
      endif
      if(ntype_t.ne.ntype) then
         write(6,*) "ntype not same in vwd_fitB.ntype,something wrong"
         stop
      endif
      do itype1=1,ntype
         read(10,*) itype_t,rad_atom(itype1),E_ave_vdw,((wp_atom(i,itype1,j1),i=1,ntype),j1=1,nterm)
      enddo
      close(10)


      ! **************** read feat.info ********************
      open (10, file=trim(feat_info_path))
      rewind (10)
      read (10, *) iflag_PCA   ! this can be used to turn off degmm part
      read (10, *) nfeat_type_l
      do kkk = 1, nfeat_type_l
         read (10, *) ifeat_type_l(kkk)   ! the index (1,2,3) of the feature type
      end do
      read (10, *) ntype_tmp
      if (ntype_tmp .ne. ntype) then
         write (*, *) "ntype of atom not same, fit_linearMM.input, feat.info, stop"
         write (*, *) ntype, ntype_tmp
         stop
      end if

      do i = 1, ntype
         read (10, *) iatom_tmp, nfeat1(i), nfeat2(i)   ! these nfeat1,nfeat2 include all ftype
         if (iatom_tmp .ne. itype_atom(i)) then
            write (*, *) "iatom not same, fit_linearMM.input, feat.info"
            write (*, *) iatom_tmp, itype_atom(i)
            stop
         end if
      end do
      close (10)

      ! cccccccc Right now, nfeat1,nfeat2,for different types
      ! cccccccc must be the same. We will change that later, allow them
      ! cccccccc to be different
      nfeat1m = 0   ! the original feature
      nfeat2m = 0   ! the new PCA, PV feature
      nfeat2tot = 0 ! tht total feature of diff atom type
      nfeat2i = 0   ! the starting point
      nfeat2i(1) = 0
      do i = 1, ntype
         if (nfeat1(i) .gt. nfeat1m) nfeat1m = nfeat1(i)
         if (nfeat2(i) .gt. nfeat2m) nfeat2m = nfeat2(i)
         nfeat2tot = nfeat2tot + nfeat2(i)
         if (i .gt. 1) then
            nfeat2i(i) = nfeat2i(i - 1) + nfeat2(i - 1)
         end if

      end do

      allocate (bb(nfeat2tot))
      allocate (bb_type0(nfeat2m, ntype))

      open (12, file=trim(model_coefficients_path))
      rewind (12)
      read (12, *) ntmp
      if (ntmp /= nfeat2tot) then
         write (6, *) 'ntmp.not.right,linear_fitb.ntype', ntmp, nfeat2tot
         stop
      end if
      do i = 1, nfeat2tot
         read (12, *) itmp, bb(i)
      end do
      close (12)
      do itype = 1, ntype
         do k = 1, nfeat2(itype)
            bb_type0(k, itype) = bb(k + nfeat2i(itype))
         end do
      end do

      allocate (pv(nfeat1m, nfeat2m, ntype))
      allocate (feat2_shift(nfeat2m, ntype))
      allocate (feat2_scale(nfeat2m, ntype))

      do itype = 1, ntype

         if (nfeat2(itype) .eq. 0) cycle

         open (11, file=trim(feat_pv_path_header)//char(itype + 48), form='unformatted')
         rewind (11)
         read (11) nfeat1_tmp, nfeat2_tmp
         if (nfeat2_tmp /= nfeat2(itype)) then
            write (6, *) 'nfeat2.not.same,feat2_ref', itype, nfeat2_tmp, nfeat2(itype)
            stop
         end if
         if (nfeat1_tmp /= nfeat1(itype)) then
            write (6, *) 'nfeat1.not.same,feat2_ref', itype, nfeat1_tmp, nfeat1(itype)
            stop
         end if
         read (11) pv(1:nfeat1(itype), 1:nfeat2(itype), itype)
         read (11) feat2_shift(1:nfeat2(itype), itype)
         read (11) feat2_scale(1:nfeat2(itype), itype)
         close (11)
      end do
      !********************add_force****************

   end subroutine load_model_lin

   subroutine set_image_info_lin(iatom_tmp, is_reset, natom_tmp)
      integer(4) :: i, j, itype, iitype
      integer iatom_tmp(natom_tmp)
      logical, intent(in) :: is_reset
      integer(4) :: image_size
      integer :: natom_tmp

      image_size = natom_tmp
      if (is_reset .or. (.not. allocated(iatom)) .or. image_size /= natom_nn) then

         if (allocated(iatom)) then
            if (image_size == natom_nn .and. maxval(abs(iatom_tmp - iatom)) == 0) then
               return
            end if
            deallocate (iatom)
            deallocate (iatom_type)
            deallocate (energy_pred_lin)
            deallocate (force_pred_lin)
         end if

         natom_nn = image_size
         allocate (iatom(natom_nn))
         allocate (iatom_type(natom_nn))
         allocate (energy_pred_lin(natom_nn))
         allocate (force_pred_lin(3, natom_nn))
         !allocate()

         iatom = iatom_tmp

         do i = 1, natom_nn
            iitype = 0
            do itype = 1, ntype
               if (itype_atom(itype) == iatom(i)) then
                  iitype = itype
               end if
            end do
            if (iitype == 0) then
               write (6, *) 'this type not found', iatom(i)
            end if
            iatom_type(i) = iitype
         end do

      end if

   end subroutine set_image_info_lin


   subroutine cal_energy_force_lin(feat, dfeat, num_neigh, list_neigh, &
        AL, xatom, natom_tmp, nfeat0_tmp, m_neigh_tmp, atom_mol_type, neigh_cutoff)

      real(8), intent(in) :: feat(nfeat0_tmp, natom_tmp)
      real*8, intent(in) :: dfeat(nfeat0_tmp, natom_tmp, m_neigh_tmp, 3)
      integer(4), intent(in) :: num_neigh(natom_tmp)
      integer(4), intent(in) :: list_neigh(m_neigh_tmp, natom_tmp)
      integer(4), intent(in) :: natom_tmp, nfeat0_tmp, m_neigh_tmp
      real(8), intent(in) :: AL(3,3)
      real(8), intent(in) :: xatom(3, natom_tmp)
      integer, intent(in) :: atom_mol_type(natom_tmp)
      real(8), intent(in) :: neigh_cutoff

      integer(4) i, j, jj
      integer(4)  :: itype, ixyz
      real(8) :: sum

      integer(4), allocatable, dimension(:) :: num
      real(8), allocatable, dimension(:, :) :: feat2
      real(8), allocatable, dimension(:, :, :) :: feat_type
      real(8), allocatable, dimension(:, :, :) :: feat2_type
      integer(4), allocatable, dimension(:, :) :: ind_type
      real(8), allocatable, dimension(:, :, :) :: dfeat_type
      real(8), allocatable, dimension(:, :, :) :: dfeat2_type
      real(8), allocatable, dimension(:, :, :, :) :: dfeat2
      real(8), allocatable, dimension(:, :, :, :) :: SS

      real*8 rad, rad1, rad2, dE, dEdd, dFx, dFy, dFz
      real*8 dx1, dx2, dx3, dx, dy, dz, dd
      real*8 w22_1, w22_2, w22F_1, w22F_2, yy, pi


      if (nfeat0_tmp /= nfeat1m .or. natom_tmp /= natom_nn .or. m_neigh_tmp /= m_neigh) then
         write (*, *) "Shape of input arrays don't match the model!"
         write(6,*) "nfeat0_tmp, nfeat1m, natom_tmp, natom_nn, m_neigh_tmp, m_neigh"
         write(6,*) nfeat0_tmp, nfeat1m, natom_tmp, natom_nn, m_neigh_tmp, m_neigh
         stop
      end if

      allocate (num(ntype))
      allocate (feat_type(nfeat1m, natom_nn, ntype))
      allocate (feat2_type(nfeat2m, natom_nn, ntype))
      allocate (feat2(nfeat2m, natom_nn))

      allocate (dfeat_type(nfeat1m, natom_nn*m_neigh*3, ntype))
      allocate (dfeat2_type(nfeat2m, natom_nn*m_neigh*3, ntype))
      allocate (dfeat2(nfeat2m, natom_nn, m_neigh, 3))
      allocate (SS(nfeat2m, natom_nn, 3, ntype))


      num = 0
      do i = 1, natom_nn
         itype = iatom_type(i)
         num(itype) = num(itype) + 1
         feat_type(:, num(itype), itype) = feat(:, i)
      end do

      do itype = 1, ntype
         if (nfeat2(itype) .eq. 0) cycle
         call dgemm('T', 'N', nfeat2(itype), num(itype), nfeat1(itype), 1.d0, pv(1,1,itype), &
         & nfeat1m, feat_type(1,1,itype), nfeat1m, 0.d0,feat2_type(1,1,itype), nfeat2m)
      end do

      do itype = 1, ntype
         do i = 1, num(itype)
            if (nfeat2(itype) .eq. 0) cycle
            do j = 1, nfeat2(itype) - 1
               feat2_type(j, i, itype) = (feat2_type(j, i, itype) - feat2_shift(j, itype))*feat2_scale(j, itype)
            end do
            feat2_type(nfeat2(itype), i, itype) = 1.d0
         end do
      end do

      num = 0
      do i = 1, natom_nn
         itype = iatom_type(i)
         num(itype) = num(itype) + 1
         feat2(:, i) = feat2_type(:, num(itype), itype)
      end do

      ! ------------------------------------------
      energy_pred_lin = 0.d0
      do i = 1, natom_nn
         itype = iatom_type(i)
         sum = 0.d0
         do j = 1, nfeat2(itype)
            sum = sum + feat2(j, i)*bb_type0(j, itype)
         end do
         energy_pred_lin(i) = sum
      end do
      ! ------------------------------------------

      num = 0
      do i = 1, natom_nn
         do jj = 1, num_neigh(i)
            itype = iatom_type(list_neigh(jj,i))
            num(itype) = num(itype) + 1
            dfeat_type(:, num(itype), itype) = dfeat(:, i, jj, 1)
            num(itype) = num(itype) + 1
            dfeat_type(:, num(itype), itype) = dfeat(:, i, jj, 2)
            num(itype) = num(itype) + 1
            dfeat_type(:, num(itype), itype) = dfeat(:, i, jj, 3)
         end do
      end do
      !cccccccc note: num(itype) is rather large, in the scane of natom*num_neigh

      do itype = 1, ntype
         call dgemm('T', 'N', nfeat2(itype), num(itype), nfeat1(itype), &
         & 1.d0, pv(1,1,itype), nfeat1m, dfeat_type(1,1,itype), nfeat1m, 0.d0, &
         & dfeat2_type(1,1,itype), nfeat2m)
      end do

      num = 0
      do i = 1, natom_nn
         do jj = 1, num_neigh(i)
            itype = iatom_type(list_neigh(jj,i))

            num(itype) = num(itype) + 1
            do j = 1, nfeat2(itype) - 1
               dfeat2(j, i, jj, 1) = dfeat2_type(j, num(itype), itype)*feat2_scale(j, itype)
            end do
            if (nfeat2(itype) .gt. 0) dfeat2(nfeat2(itype), i, jj, 1) = 0.d0

            num(itype) = num(itype) + 1
            do j = 1, nfeat2(itype) - 1
               dfeat2(j, i, jj, 2) = dfeat2_type(j, num(itype), itype)*feat2_scale(j, itype)
            end do
            if (nfeat2(itype) .gt. 0) dfeat2(nfeat2(itype), i, jj, 2) = 0.d0

            num(itype) = num(itype) + 1
            do j = 1, nfeat2(itype) - 1
               dfeat2(j, i, jj, 3) = dfeat2_type(j, num(itype), itype)*feat2_scale(j, itype)
            end do
            if (nfeat2(itype) .gt. 0) dfeat2(nfeat2(itype), i, jj, 3) = 0.d0
         end do
      end do

      SS = 0.d0
      do i = 1, natom_nn
         do jj = 1, num_neigh(i)
            itype = iatom_type(list_neigh(jj, i))  ! this is this neighbor's type
            do j = 1, nfeat2(itype)
               SS(j, i, 1, itype) = SS(j, i, 1, itype) + dfeat2(j, i, jj, 1)
               SS(j, i, 2, itype) = SS(j, i, 2, itype) + dfeat2(j, i, jj, 2)
               SS(j, i, 3, itype) = SS(j, i, 3, itype) + dfeat2(j, i, jj, 3)
            end do
         end do
      end do

      ! ---------------------------------------------

      force_pred_lin = 0.d0
      do i = 1, natom_nn
         do ixyz = 1, 3
            sum = 0.d0
            do itype = 1, ntype
               do j = 1, nfeat2(itype)
                  sum = sum + SS(j, i, ixyz, itype)*BB_type0(j, itype)
               end do
            end do
            force_pred_lin(ixyz, i) = sum
         end do
      end do

      ! ---------------------------------------------



      ! get back the vdw part begin

      pi=4*datan(1.d0)

      do i=1,natom_nn

         rad1=rad_atom(iatom_type(i))

         dE=0.d0
         dFx=0.d0
         dFy=0.d0
         dFz=0.d0

         do jj=1,num_neigh(i)
            j=list_neigh(jj,i)

            ! if((i .ne. j) .and. (atom_mol_type(i) .eq. atom_mol_type(j)) ) then
            !    write(6,*) "@@@ same mol ", i, j, atom_mol_type(i), atom_mol_type(j)
            ! endif

            if(atom_mol_type(i) .ne. atom_mol_type(j)) then
            ! if(i.ne.j) then
               rad2=rad_atom(iatom_type(j))
               rad=rad1+rad2
               dx1=mod(xatom(1,j)-xatom(1,i)+100.d0,1.d0)
               if(abs(dx1-1).lt.abs(dx1)) dx1=dx1-1
               dx2=mod(xatom(2,j)-xatom(2,i)+100.d0,1.d0)
               if(abs(dx2-1).lt.abs(dx2)) dx2=dx2-1
               dx3=mod(xatom(3,j)-xatom(3,i)+100.d0,1.d0)
               if(abs(dx3-1).lt.abs(dx3)) dx3=dx3-1
               dx=AL(1,1)*dx1+AL(1,2)*dx2+AL(1,3)*dx3
               dy=AL(2,1)*dx1+AL(2,2)*dx2+AL(2,3)*dx3
               dz=AL(3,1)*dx1+AL(3,2)*dx2+AL(3,3)*dx3
               dd=dsqrt(dx**2+dy**2+dz**2)

               if (dd .lt. neigh_cutoff) then
               ! if(dd.lt.2*rad) then
                  !        w22=dsqrt(wp_atom(iatom_type(i))*wp_atom(iatom_type(j)))
                  !        yy=pi*dd/(4*rad)
                  ! !       dE=dE+0.5*w22*exp((1-dd/rad)*4.0)*cos(yy)**2
                  ! !       dEdd=w22*exp((1-dd/rad)*4.d0)*((-4/rad)*cos(yy)**2
                  ! !     &   -(pi/(2*rad))*cos(yy)*sin(yy))
                  !        dE=dE+0.5*4*w22*(rad/dd)**12*cos(yy)**2
                  !        dEdd=4*w22*(-12*(rad/dd)**12/dd*cos(yy)**2  &
                  !         -(pi/(2*rad))*cos(yy)*sin(yy)*(rad/dd)**12)
                  w22_1=wp_atom(iatom_type(j),iatom_type(i),1)
                  w22_2=wp_atom(iatom_type(j),iatom_type(i),2)
                  w22F_1=(wp_atom(iatom_type(j),iatom_type(i),1)+wp_atom(iatom_type(i),iatom_type(j),1))/2     ! take the average for force calc.
                  w22F_2=(wp_atom(iatom_type(j),iatom_type(i),2)+wp_atom(iatom_type(i),iatom_type(j),2))/2     ! take the average for force calc.

                  yy=pi*dd/(4*rad)
                  ! c       dE=dE+0.5*w22*exp((1-dd/rad)*4.0)*cos(yy)**2
                  ! c       dEdd=w22*exp((1-dd/rad)*4.d0)*((-4/rad)*cos(yy)**2
                  ! c     &   -(pi/(2*rad))*cos(yy)*sin(yy))
                  dE=dE+0.5*4*(w22_1*(rad/dd)**12*cos(yy)**2+w22_2*(rad/dd)**6*cos(yy)**2)
                  dEdd=4*(w22F_1*(-12*(rad/dd)**12/dd*cos(yy)**2-(pi/(2*rad))*cos(yy)*sin(yy)*(rad/dd)**12)   &
                     +W22F_2*(-6*(rad/dd)**6/dd*cos(yy)**2-(pi/(2*rad))*cos(yy)*sin(yy)*(rad/dd)**6))

                  dFx=dFx-dEdd*dx/dd       ! note, -sign, because dx=d(j)-x(i)
                  dFy=dFy-dEdd*dy/dd
                  dFz=dFz-dEdd*dz/dd
               endif
            endif
         enddo

         energy_pred_lin(i)=energy_pred_lin(i)+dE
         force_pred_lin(1,i)=force_pred_lin(1,i)+dFx   ! Note, assume force=dE/dx, no minus sign
         force_pred_lin(2,i)=force_pred_lin(2,i)+dFy
         force_pred_lin(3,i)=force_pred_lin(3,i)+dFz

      enddo
      ! get back the vdw part end


      ! ---------------------------------------------

      etot_pred_lin = 0.d0
      do i = 1, natom_nn
         etot_pred_lin = etot_pred_lin + energy_pred_lin(i)
      end do

      ! ---------------------------------------------

      deallocate (num)
      deallocate (feat_type, feat2_type, feat2)
      deallocate (dfeat_type, dfeat2_type, dfeat2, SS)

   end subroutine cal_energy_force_lin

end module calc_lin_PC


