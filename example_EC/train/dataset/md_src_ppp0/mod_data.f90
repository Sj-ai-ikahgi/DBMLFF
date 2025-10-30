module mod_data
    use mod_parameter
    real*8 Vx_1(matom_1),Vy_1(matom_1),Vz_1(matom_1)
    real*8 xatom(3,matom_1)
    integer imov_at(3,matom_1)
    integer iatom(matom_1)
    real*8 iMDatom(matom_1)
    real*8 weight_atom_1(matom_1)
    real*8 langevin_factT(matom_1),langevin_factG(matom_1)
!    integer ntype
    real*8 AL(3,3),ALI(3,3)
    real*8 stress_mask(3,3)
    real*8 stress_ext(3,3)
    character*20 f_xatom
    integer iflag_mlff, iflag_model,iflag_charge,iflag_training,iflag_dftd3, iflag_pc
    real*8 polar_factor
    logical iflag_debug
    real*8,allocatable,dimension(:,:,:,:) :: xatom_m
    real*8,allocatable,dimension(:,:) :: mass_mol
    integer,allocatable,dimension(:,:) :: iatom_m
    integer natom,natom_m(1000),nmol(1000)
    integer ntype_ML(1000)     ! number of type in ML model, it is the same as the
                         !  ntype read from fit_input_path, stored in calc_lin
    integer ntype_m     ! the type of molecule, read in from xatom
    integer ntype_mm    ! max number of molecule, read in from MD.input
    integer natom_mm    ! the max of natom_m
    integer nmolm       ! the maximum of nmol
    real*8, allocatable, device, dimension(:,:,:,:) :: xatom_m_d
    integer, allocatable, device, dimension(:, :) :: iatom_m_d
    real*8, allocatable, device, dimension(:, :) :: mass_mol_d
    real*8 Etot_nonbond,E_polar,E_kin2_nonbond
end module mod_data 

