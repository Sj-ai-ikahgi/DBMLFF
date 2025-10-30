module calc_ftype2b_3b

   IMPLICIT NONE

   integer m_neigh, ntype
   integer n2b_type(100), n2bm
   real*8,allocatable,dimension (:,:,:) :: grid2_2
   integer n3b1_type(100),n3b2_type(100), n3b1m, n3b2m
   real*8,allocatable,dimension (:,:,:) :: grid31_2, grid32_2
   integer nfeat0m2, nfeat0m3

   integer, allocatable, dimension(:) :: itype_atom
   integer, allocatable, dimension(:) :: nfeat2_atom
   integer, allocatable, dimension(:) :: nfeat3_atom
   real*8, allocatable, dimension(:, :) :: feat2
   real*8, allocatable, dimension(:, :) :: feat3
   real*8, allocatable, dimension(:, :, :, :) :: dfeat2
   real*8, allocatable, dimension(:, :, :, :) :: dfeat3

   integer, allocatable, dimension(:) :: neigh_numR
   integer, allocatable, dimension(:,:) :: neigh_listR

   integer iat_type(100)

contains

   subroutine load_model_type2(ipair_type)

      integer, intent(in) :: ipair_type

      real*8 Rc_M
      real*8 Rc_type(100), Rc2_type(100), Rm_type(100),fact_grid_type(100),dR_grid1_type(100),dR_grid2_type(100)
      real*8 E_tolerance
      integer iflag_ftype, recalc_grid
      integer iflag_grid_type(100), iflag_grid

      integer n2b, n2b_t, n3b1, n3b2, n3b1_t, n3b2_t
      integer i, itype1, itype2, k1, k2, k12, kkk, ii_f, num, it

      character(80), parameter :: feat2_input_path0 = "input/gen_2b_feature.in"
      character(80), parameter :: grid2_input_header0 = "output/grid2b_type3."
      character(200) :: feat2_input_path = trim(feat2_input_path0)
      character(200) :: grid2_input_header = trim(grid2_input_header0)

      character(80), parameter :: feat3_input_path0 = "input/gen_3b_feature.in"
      character(80), parameter :: grid3_input_header01 = "output/grid3b_cb12_type3."
      character(80), parameter :: grid3_input_header02 = "output/grid3b_b1b2_type3."
      character(200) :: feat3_input_path = trim(feat3_input_path0)
      character(200) :: grid3_input_header1 = trim(grid3_input_header01)
      character(200) :: grid3_input_header2 = trim(grid3_input_header02)

      feat2_input_path = 'pair.' // char(48 + ipair_type)//'/'//trim(feat2_input_path0)
      grid2_input_header = 'pair.' // char(48 + ipair_type)//'/'//trim(grid2_input_header0)

      feat3_input_path = 'pair.' // char(48 + ipair_type)//'/'//trim(feat3_input_path0)
      grid3_input_header1 = 'pair.' // char(48 + ipair_type)//'/'//trim(grid3_input_header01)
      grid3_input_header2 = 'pair.' // char(48 + ipair_type)//'/'//trim(grid3_input_header02)

      ! get n2b_type, n2bm, grid2_2, nfeat0m2
      open(10,file=trim(feat2_input_path),status="old",action="read")
      rewind(10)
      read(10,*) Rc_M,m_neigh
      read(10,*) ntype

      do i=1,ntype
         read(10,*) iat_type(i)
         read(10,*) Rc_type(i),Rm_type(i),iflag_grid_type(i),fact_grid_type(i),dR_grid1_type(i)
         read(10,*) n2b_type(i)

         if(Rc_type(i).gt.Rc_M) then
            write(6,*) "Rc_type must be smaller than Rc_M, gen_3b_feature.in",i,Rc_type(i),Rc_M
            stop
         endif
      enddo

      read(10,*) E_tolerance
      read(10,*) iflag_ftype
      read(10,*) recalc_grid
      close(10)

      do i=1,ntype
         if(iflag_ftype.eq.3.and.iflag_grid_type(i).ne.3) then
            write(6,*) "if iflag_ftype.eq.3, iflag_grid must equal 3, stop"
            stop
         endif
      enddo

      n2bm=0
      do i=1,ntype
         if(n2b_type(i).gt.n2bm) n2bm=n2b_type(i)
      enddo

      nfeat0m2=ntype*n2bm

      if (allocated(grid2_2))  deallocate(grid2_2)
      allocate(grid2_2(2,n2bm+1,ntype)) 

      do kkk=1,ntype    ! center atom

         iflag_grid=iflag_grid_type(kkk)
         n2b=n2b_type(kkk)

         if(iflag_grid.eq.3) then
            ! for iflag_grid.eq.3, the graid is just read in.
            ! Its format is different from above grid31, grid32.
            ! For each point, it just have two numbers, r1,r2, indicating the region of the sin peak function.

            open(13,file=trim(grid2_input_header)//char(kkk+48))
            rewind(13)
            read(13,*) n2b_t

            if(n2b_t.ne.n2b) then
               write(6,*) "n2b_t not equivalent to n2b in grid2b_type3", n2b_t,n2b
               stop
            endif

            do i=1,n2b
               read(13,*) it,grid2_2(1,i,kkk),grid2_2(2,i,kkk)
               if(grid2_2(2,i,kkk).gt.Rc_type(kkk)) write(6,*) "grid2_2 greater than Rc",grid2_2(2,i,kkk),Rc_type(kkk)
            enddo

            close(13)
         endif

      enddo     ! kkk=1,ntype



      open(10,file=trim(feat3_input_path),status="old",action="read")
      rewind(10)
      read(10,*) Rc_M,m_neigh

      read(10,*) ntype
      do i=1,ntype
         read(10,*) iat_type(i)
         read(10,*) Rc_type(i),Rc2_type(i),Rm_type(i),iflag_grid_type(i),fact_grid_type(i),dR_grid1_type(i),dR_grid2_type(i)
         read(10,*) n3b1_type(i),n3b2_type(i)

         if(Rc_type(i).gt.Rc_M) then
            write(6,*) "Rc_type must be smaller than Rc_M, gen_3b_feature.in",i,Rc_type(i),Rc_M
            stop
         endif
         if(Rc2_type(i).gt.2*Rc_type(i)) then
            write(6,*) "Rc2_type must be smaller than 2*Rc_type, gen_3b_feature.in",i,Rc_type(i),Rc2_type(i)
            stop
         endif

      enddo
      read(10,*) E_tolerance
      read(10,*) iflag_ftype
      read(10,*) recalc_grid
      close(10)


      do i=1,ntype
         if(iflag_ftype.eq.3.and.iflag_grid_type(i).ne.3) then
            write(6,*) "if iflag_ftype.eq.3, iflag_grid must equal 3, stop"
            stop
         endif
      enddo

      n3b1m=0
      n3b2m=0
      do i=1,ntype
         if(n3b1_type(i).gt.n3b1m) n3b1m=n3b1_type(i)
         if(n3b2_type(i).gt.n3b2m) n3b2m=n3b2_type(i)
      enddo

      !cccccccccccccccccccccccccccccccccccccccccccccccc
      num=0
      do itype2=1,ntype
         do itype1=1,itype2
            do k1=1,n3b1m
               do k2=1,n3b1m
                  do k12=1,n3b2m
                     ii_f=0
                     if(itype1.ne.itype2) ii_f=1
                     if(itype1.eq.itype2.and.k1.le.k2) ii_f=1
                     if(ii_f.gt.0) then
                        num=num+1
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo

      nfeat0m3 = num

      !cccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (allocated(grid31_2)) deallocate(grid31_2)
      if (allocated(grid32_2)) deallocate(grid32_2)
      allocate(grid31_2(2,n3b1m,ntype))
      allocate(grid32_2(2,n3b2m,ntype))

      !cccccccccccccccccccccccccccccccccccccccccccccccccccc

      do kkk=1,ntype    ! center atom

         iflag_grid=iflag_grid_type(kkk)
         n3b1=n3b1_type(kkk)
         n3b2=n3b2_type(kkk)

         !cccccccccccccccccccccccccccccccccccccccc
         !cccccccccccccccccccccccccccccccccccccccccccc
         if(iflag_grid.eq.3) then
            ! for iflag_grid.eq.3, the graid is just read in.
            ! Its format is different from above grid31, grid32.
            ! For each point, it just have two numbers, r1,r2, indicating the region of the sin peak function.

            open(13,file=trim(grid3_input_header1)//char(kkk+48))
            rewind(13)
            read(13,*) n3b1_t
            if(n3b1_t.ne.n3b1) then
               write(6,*) "n3b1_t.ne.n3b1,in grid31_type3", n3b1_t,n3b1
               stop
            endif
            do i=1,n3b1
               read(13,*) it,grid31_2(1,i,kkk),grid31_2(2,i,kkk)
               if(grid31_2(2,i,kkk).gt.Rc_type(kkk)) write(6,*) "grid31_2.gt.Rc",grid31_2(2,i,kkk),Rc_type(kkk)
            enddo
            close(13)

            open(13,file=trim(grid3_input_header2)//char(kkk+48))
            rewind(13)
            read(13,*) n3b2_t
            if(n3b2_t.ne.n3b2) then
               write(6,*) "n3b2_t.ne.n3b2,in grid32_type3", n3b2_t,n3b2
               stop
            endif
            do i=1,n3b2
               read(13,*) it,grid32_2(1,i,kkk),grid32_2(2,i,kkk)
               if(grid32_2(2,i,kkk).gt.Rc2_type(kkk)) write(6,*) "grid32_2.gt.Rc",grid32_2(2,i,kkk),Rc2_type(kkk)
            enddo
            close(13)
         endif

      enddo     ! kkk=1,ntype

   end subroutine load_model_type2


   subroutine gen_feature_type2(AL, natom, iatom, xatom, &
    & Rc, mol1_type, atom_mol_type, atom_ind_mol, mol1_atom_num)

      real(8), intent(in) :: AL(3,3)
      integer(4), intent(in) :: natom
      integer(4), intent(in) :: iatom(natom)
      real(8), intent(in) :: xatom(3,natom)
      real(8), intent(in) :: Rc
      integer, intent(in) :: mol1_type
      integer, intent(in) :: atom_mol_type(natom)
      integer, intent(in) :: atom_ind_mol(natom)
      integer, intent(in) :: mol1_atom_num

      integer(4)  :: i, j
      integer tnum

      if (allocated(itype_atom))  deallocate (itype_atom)
      if (allocated(nfeat2_atom)) deallocate (nfeat2_atom)
      if (allocated(nfeat3_atom))   deallocate (nfeat3_atom)
      if (allocated(feat2))   deallocate (feat2)
      if (allocated(feat3))   deallocate (feat3)
      if (allocated(dfeat2))   deallocate (dfeat2)
      if (allocated(dfeat3))   deallocate (dfeat3)
      if (allocated(neigh_numR))   deallocate (neigh_numR)
      if (allocated(neigh_listR))   deallocate (neigh_listR)

      allocate (itype_atom(natom))
      allocate (nfeat2_atom(natom))
      allocate (nfeat3_atom(natom))
      allocate (feat2(nfeat0m2, natom))
      allocate (feat3(nfeat0m3, natom))
      allocate (dfeat2(nfeat0m2, natom, m_neigh, 3))
      allocate (dfeat3(nfeat0m3, natom, m_neigh, 3))
      allocate (neigh_numR(natom))
      allocate (neigh_listR(m_neigh, natom))

      itype_atom = 0
      do i = 1, natom
         do j = 1, ntype
            if (iatom(i) .eq. iat_type(j)) then
               itype_atom(i) = j
            end if
         end do
         if (itype_atom(i) .eq. 0) then
            write (6, *) "this atom type didn't found", itype_atom(i)
            stop
         end if
      end do

      call find_feature_2b_3b_type3(natom, xatom, AL, ntype, m_neigh, itype_atom, &
         n2b_type, n2bm, grid2_2, &
         n3b1_type, n3b2_type, n3b1m, n3b2m, grid31_2, grid32_2, &
         nfeat0m2, nfeat0m3, &
         Rc, mol1_type, atom_mol_type, atom_ind_mol, mol1_atom_num, &
         nfeat2_atom, feat2, dfeat2, &
         nfeat3_atom, feat3, dfeat3, &
         neigh_numR, neigh_listR)


      ! do i = 1, natom
      !    tnum = neigh_numR(i)
      !    write(6,*) i, tnum, neigh_listR(1:tnum,i)
      ! enddo


   end subroutine gen_feature_type2

end module calc_ftype2b_3b
