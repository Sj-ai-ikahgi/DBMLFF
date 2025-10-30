program convert_rho

!ccccc This program convert a unformatted dens (or vr), graph.R, or
!charge_out file
!ccccc into a formatted xxx.plt_f file for 3D isosurface plot
!ccccc using gOpenMol. This formatted .plt_f file needs to be
!ccccc transformed into a unformatted .plt file using pltfile
!ccccc within gopenmol (RUN), before be read in by the
!gopenmol-plot-contour
!ccccc option.

   implicit double precision (a-h,o-z)

   real*8 AL(3,3),ALI(3,3),AL0(3,3)

   character*40 filename

   real*8, allocatable, dimension (:,:,:) :: vr
   real*8, allocatable, dimension (:,:,:) :: vr_xyz
   real*8, allocatable, dimension(:) :: vr_tmp
   integer  stat, sum_atom
   character(len=256) :: right, atomfile
   logical :: readit, alive
   real*8, allocatable :: atomic_position_frac(:,:)
   integer, allocatable :: atomic_number(:)
   real*8, allocatable :: atomic_position_cart_ang(:,:)
   real*8 :: lattice_ang(3,3)


   call getarg(1,filename)

   write(6,*) filename

   call input_vr()

   open (33, file = "atom.config")
   read (33, *) sum_atom
   allocate(atomic_number(sum_atom),atomic_position_frac(3,sum_atom))
   allocate(atomic_position_cart_ang(3,sum_atom))
   read (33, *)
   do i = 1, 3
      read (33, *) lattice_ang(1:3,i)
   end do
   read (33, *)
   do i = 1, sum_atom
      read (33, *) atomic_number(i), atomic_position_frac(1:3,i)
   end do
   do i  = 1, sum_atom
      do j = 1, 3
         atomic_position_cart_ang(j,i) =           &
            atomic_position_frac(1,i)*lattice_ang(j,1) +  &
            atomic_position_frac(2,i)*lattice_ang(j,2) +  &
            atomic_position_frac(3,i)*lattice_ang(j,3)
      end do
   end do
   !

   open(12,file='RHO.xsf')
   rewind(12)
   write(12,'(x,a)') 'CRYSTAL'
   write(12,'(x,a)') 'PRIMVEC'
   do i = 1, 3
      write (12, '(3f18.10)') lattice_ang(1:3,i)
   end do
   !write(12,'(x,a)') 'CONVVEC'
   !do i = 1, 3
   !    write (12, '(3f14.10)') lattice_ang(1:3,i)
   !end do
   write(12,'(x,a)') 'PRIMCOORD'
   write(12, *) sum_atom, 1
   do i = 1, sum_atom
      write(12, '(x,i3,3f18.10)') atomic_number(i),atomic_position_cart_ang(1:3,i)
   end do
   write(12,1001)"BEGIN_BLOCK_DATAGRID_3D"
   write(12,1001)"XSF_FILE"
   write(12,1001)" BEGIN_DATAGRID_3D_XSF_FILE"
   write(12,98) n1+1,n2+1,n3+1
   write(12,200) 0.0000,0.00000,0.00000
!c       write(12,200) ALI(1,1),ALI(2,1),ALI(3,1)
!c       write(12,200) ALI(1,2),ALI(2,2),ALI(3,2)
!c       write(12,200) ALI(1,3),ALI(2,3),ALI(3,3)
   write(6,*) AL(1,1),AL(2,1),AL(3,1)
   write(6,*) AL(1,2),AL(2,2),AL(3,2)
   write(6,*) AL(1,3),AL(2,3),AL(3,3)

   write(12,200) AL(1,1),AL(2,1),AL(3,1)
   write(12,200) AL(1,2),AL(2,2),AL(3,2)
   write(12,200) AL(1,3),AL(2,3),AL(3,3)
   write(12,*)
   write(12,100) (((vr(mod(i-1,n1)+1,mod(j-1,n2)+1,mod(k-1,n3)+1),i=1,n1+1),j=1,n2+1),k=1,n3+1)
   write(12,1001)"END_DATAGRID_3D"
   write(12,1001)"END_BLOCK_DATAGRID_3D"
   close(12)

99 format(2(i2,2x))
98 format(3(i5,2x))
200 format(3(1x,E13.5,1x))
!100    format(6(F16.9,1x))
100 format(6(F12.6))
1001 format(a28)

   deallocate(vr)

   stop

contains
!************************************************

   subroutine input_vr()
      implicit double precision (a-h,o-z)
      integer ierr
!*************************************************************
      open(11,file=filename,form="unformatted")
      rewind(11)
      read(11,IOSTAT=ierr) n1,n2,n3,nnodes,nstate

      if(ierr.ne.0) then
         rewind(11)
         read(11,IOSTAT=ierr) n1,n2,n3,nnodes
         nstate=1
      endif
      read(11) AL

      nr=n1*n2*n3
      nr_n=nr/nnodes
      allocate(vr_tmp(nr_n))
      allocate(vr(n1,n2,n3))
      if(nstate.ne.1) then
         write(6,*) "there are ", nstate, " states"
         write(6,*) "input the state index to be plotted"
         read(6,*) istate
      else
         istate=1
      endif


      do ist=1,istate
         do iread=1,nnodes
            read(11) vr_tmp

            do ii=1,nr_n

               jj=ii+(iread-1)*nr_n

               i=(jj-1)/(n2*n3)+1
               j=(jj-1-(i-1)*n2*n3)/n3+1
               k=jj-(i-1)*n2*n3-(j-1)*n3

               vr(i,j,k)=vr_tmp(ii)
            enddo
         enddo
      enddo
      close(11)

!************************************************************
      return
   end subroutine input_vr

end




