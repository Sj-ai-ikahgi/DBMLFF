module mod_param_densityFF
      use cudafor
      implicit none
  
      !  50: the max number of molecule type
      !  20: the max number of ion type
      real*8 Ecut2, fact_kin2, Rbox
      integer n1, n2, n3       ! box for the whole AL
      integer ntype_ion, ntype_ion_pol(50)
      integer nm1_all(1000), nm2_all(1000), nm3_all(1000)    ! box for each molecule
      integer ncent(1000), ntype_cent(1000)
      !      integer itype_cent(1000),icent(2,1000),nat_cent(1000)
      !      real*8   z_cent(1000),w_cent(2,1000)
      integer, allocatable, dimension(:, :) ::  itype_cent, nat_cent
      integer, allocatable, device, dimension(:, :) ::  itype_cent_d, nat_cent_d
      integer, allocatable, dimension(:, :, :) ::  icent
      integer, allocatable, device, dimension(:, :, :) ::  icent_d
      real*8, allocatable, dimension(:, :) ::  z_cent
      real*8, allocatable, device, dimension(:, :) :: z_cent_d
      real*8, allocatable, dimension(:, :, :) ::  w_cent
      real*8, allocatable, device, dimension(:, :, :) ::  w_cent_d
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8, allocatable, dimension(:, :) ::  BB_poldir
      integer, allocatable, dimension(:, :, :) :: ind
      integer num_term_pol(20)
      integer ipol(4, 1000, 20)
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer nr
      real*8 Rm2
      real*8, allocatable, dimension(:, :) :: Q_type
      real*8, allocatable, device, dimension(:, :) :: Q_type_d
      real*8, allocatable, dimension(:, :, :) :: funcr2, dfuncr2
      real*8, allocatable, device, dimension(:, :, :) :: funcr2_d, dfuncr2_d
  
      real*8 r_ion(5000, 20), rho_ion(5000, 20)
      real*8, allocatable, device, dimension(:, :) ::  r_ion_d, rho_ion_d
      integer imax_ion(20)
      integer, allocatable, device, dimension(:) :: imax_ion_d
      integer ion_type(20), ion_type_pol(20, 50)
      real*8 z_ion(20)
      real*8, allocatable, device, dimension(:) :: z_ion_d
  
      integer ion_type_cent(1000, 50), ion_type_atomp(1000, 50)
      integer, allocatable, device, dimension(:, :) :: ion_type_cent_d, ion_type_atomp_d
      real*8 atom_charge_param(4, 20, 50)
      real*8, allocatable, device, dimension(:, :, :) :: atom_charge_param_d
      integer iflag_polar
      real*8 E_polar_max(50)
      integer nm1_max, nm2_max, nm3_max
      integer imax_nr, imax_ntype_cent
  
      integer indb(3, 1000, 20), numb(20), max_numb
      real*8 BB_direct(2000, 20)
  
      type(c_devptr)  c_devptr_funcr2, c_devptr_dfuncr2, c_devptr_itype_cent,c_devptr_nat_cent,c_devptr_icent,c_devptr_z_cent,c_devptr_w_cent,  &
          c_devptr_iatom_m, c_devptr_mass_mol, c_devptr_r_ion, c_devptr_rho_ion, c_devptr_imax_ion, c_devptr_z_ion, c_devptr_ion_type_cent, &
          c_devptr_ion_type_atomp, c_devptr_Q_type, c_devptr_atom_charge_param
  
  contains
  
      subroutine read_param_densityFF(ntype_m)
          implicit none
          real*8 dr
          integer itype, i, ntype_cent_tmp
          real*8 sum, sum1, r
          integer ntype_m
          integer istat
          integer itype_mol
          integer nr_all(ntype_m)
          ! In the future, provide a iflag, when iflag.eq.1, read
          !  When iflag.eq.2, copy the parameter from memory storage
  
          !ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          do itype_mol = 1, ntype_m
              open (13, file='mol.'//char(itype_mol + 48)//'/'//'funcr_atom.fit.bin', form="unformatted")
              rewind (13)
              read (13) nr, Rm2, ntype_cent_tmp
              if (ntype_cent_tmp .ne. ntype_cent(itype_mol)) then
                  write (6, *) "ntype_cent in funcr_atom.fit.bin not right", &
                      ntype_cent_tmp, ntype_cent
                  stop
              end if
              nr_all(itype_mol) = nr
              close (13)
          end do
          imax_nr = 0
          imax_ntype_cent = 0
          do itype_mol = 1, ntype_m
              imax_nr = max(nr_all(itype_mol), imax_nr)
              imax_ntype_cent = max(ntype_cent(itype_mol), imax_ntype_cent)
          end do
          allocate (Q_type(imax_ntype_cent, ntype_m))
          allocate (funcr2(imax_nr, imax_ntype_cent, ntype_m))
          allocate (dfuncr2(imax_nr, imax_ntype_cent, ntype_m))
  
          istat = cudaMalloc(c_devptr_funcr2, imax_nr*imax_ntype_cent*ntype_m*8)
          istat = cudaMalloc(c_devptr_dfuncr2, imax_nr*imax_ntype_cent*ntype_m*8)
          istat = cudaMalloc(c_devptr_Q_type, imax_ntype_cent*ntype_m*8)
  
          call c_f_pointer(c_devptr_funcr2, funcr2_d, (/imax_nr, imax_ntype_cent, ntype_m/))
          call c_f_pointer(c_devptr_dfuncr2, dfuncr2_d, (/imax_nr, imax_ntype_cent, ntype_m/))
          call c_f_pointer(c_devptr_Q_type, Q_type_d, (/imax_ntype_cent, ntype_m/))
  
          do itype_mol = 1, ntype_m
              open (13, file='mol.'//char(itype_mol + 48)//'/'//'funcr_atom.fit.bin', form="unformatted")
              rewind (13)
              read (13) nr, Rm2, ntype_cent_tmp
              if (ntype_cent_tmp .ne. ntype_cent(itype_mol)) then
                  write (6, *) "ntype_cent in funcr_atom.fit.bin not right", &
                      ntype_cent_tmp, ntype_cent
                  stop
              end if
              read (13) Q_type(1:ntype_cent(itype_mol), itype_mol)
              read (13) funcr2(1:nr, 1:ntype_cent(itype_mol), itype_mol)
              close (13)
              dr = 2*Rm2/nr
              do itype = 1, ntype_cent(itype_mol)
                  do i = 2, nr - 1
                      dfuncr2(i, itype, itype_mol) = (funcr2(i + 1, itype, itype_mol) - funcr2(i - 1, itype, itype_mol))/dr
                  end do
              end do
              dfuncr2(1, :, itype_mol) = 0.d0
              dfuncr2(nr, :, itype_mol) = 0.d0
          end do
          funcr2_d = funcr2
          dfuncr2_d = dfuncr2
          Q_type_d = Q_type
          return
      end subroutine read_param_densityFF
  
  end module mod_param_densityFF
  