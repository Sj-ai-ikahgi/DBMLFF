subroutine find_pair(num_mol1,num_mol2,num_atom_mol1,num_atom_mol2,&
  ind_ref_atom_mol1,ind_ref_atom_mol2,lattice_constant,radius,&
  x_frac_mol1,x_frac_mol2,num_pair,pair,image,distance)
!======================================================================
! added_by_fx
!
! Input the fractional coords of all mol1 and mol2, and a radius for 
! judging neighbor. Output all possible pairs and images of mol2
! PLEASE BE ATTENTION!!! Assume the lattice is orthogonal
! PLEASE BE ATTENTION!!! x_frac_mol1 and x_frac_mol2 are both IO! All 
! mol1 and mol2 will be clustered and moved into the unitcell.
!
! I     num_mol1
! I     num_mol2
! I     num_atom_mol1
! I     num_atom_mol2
! I     ind_ref_atom_mol1
! I     ind_ref_atom_mol2
! I     lattice_constant: lattice constant in A
! I     radius: radius in A
! IO    x_frac_mol1: fractional coords
! IO    x_frac_mol2: fractional coords
! O     num_pair
! O     pair
! O     image
! O     distance
!======================================================================
  implicit double precision (a-h,o-z)

  integer :: num_mol1, num_mol2, num_atom_mol1, num_atom_mol2
  integer :: ind_ref_atom_mol1, ind_ref_atom_mol2
  real*8 lattice_constant(3)
  real*8 radius

  real*8 x_frac_mol1(3,100,100)
  real*8 x_frac_mol2(3,100,100)

  integer num_pair
  
  integer pair(2,100)
  integer image(3,100)
  real*8 distance(100)
  real*8 min_d
  
  real*8,allocatable,dimension(:,:) :: molecule1_center
  real*8,allocatable,dimension(:,:) :: molecule2_center
  integer,allocatable,dimension(:,:) :: supercell_range
  
  num_pair = 0
  pair = 0
  image = 0
  distance = 0.d0

  ! ----------Calculate centers----------
  allocate(molecule1_center(3,num_mol1))
  allocate(molecule2_center(3,num_mol2))
  
  ! PLEASE BE ATTENTION!!! x_frac_mol1 and x_frac_mol2 may be modified here (All mol1 and mol2 are clustered)!
  do i=1,num_mol1
      call calc_molecule_center(num_atom_mol1,ind_ref_atom_mol1,x_frac_mol1(:,:,i),molecule1_center(:,i))
  enddo

  do i=1,num_mol2
      call calc_molecule_center(num_atom_mol2,ind_ref_atom_mol2,x_frac_mol2(:,:,i),molecule2_center(:,i))
  enddo
  
  ! ----------Move all molecules into the unitcell----------
  ! PLEASE BE ATTENTION!!! x_frac_mol1 and x_frac_mol2 may be modified here (All mol1 and mol2 are moved into the unitcell)!
  do i=1,num_mol1
      do ixyz=1,3
          do while (molecule1_center(ixyz,i).gt.1)
              molecule1_center(ixyz,i)=molecule1_center(ixyz,i)-1
              
              do j=1,num_atom_mol1
                  x_frac_mol1(ixyz,j,i)=x_frac_mol1(ixyz,j,i)-1
              enddo
          enddo
          
          do while (molecule1_center(ixyz,i).lt.0)
              molecule1_center(ixyz,i)=molecule1_center(ixyz,i)+1
              
              do j=1,num_atom_mol1
                  x_frac_mol1(ixyz,j,i)=x_frac_mol1(ixyz,j,i)+1
              enddo
          enddo
      enddo
  enddo
  
  do i=1,num_mol2
      do ixyz=1,3
          do while (molecule2_center(ixyz,i).gt.1)
              molecule2_center(ixyz,i)=molecule2_center(ixyz, i)-1
              
              do j=1,num_atom_mol2
                  x_frac_mol2(ixyz,j,i)=x_frac_mol2(ixyz,j,i)-1
              enddo
          enddo
          
          do while (molecule2_center(ixyz,i).lt.0)
              molecule2_center(ixyz,i)=molecule2_center(ixyz,i)+1
              
              do j=1,num_atom_mol2
                  x_frac_mol2(ixyz,j,i)=x_frac_mol2(ixyz,j,i)+1
              enddo
          enddo
      enddo
  enddo
  
  ! ----------Create a supercell----------
  ! TODO
  ! We need to narrow the range for searching mol2 further
  ! We can refer to the practice of pymatgen
  allocate(supercell_range(2,3))
  
  do i=1,3
      supercell_range(1,i)=-ceiling(radius/lattice_constant(i))
      supercell_range(2,i)=ceiling(radius/lattice_constant(i))
  enddo
  
  num_pair=0
  min_d = 1000.0
  do i=supercell_range(1,1),supercell_range(2,1)
      do j=supercell_range(1,2),supercell_range(2,2)
          do k=supercell_range(1,3),supercell_range(2,3)
              do imol1=1,num_mol1
                  do imol2=1,num_mol2
                      dx=abs(molecule1_center(1,imol1)-(molecule2_center(1,imol2)+i))*lattice_constant(1)
                      dy=abs(molecule1_center(2,imol1)-(molecule2_center(2,imol2)+j))*lattice_constant(2)
                      dz=abs(molecule1_center(3,imol1)-(molecule2_center(3,imol2)+k))*lattice_constant(3)
                      d=dsqrt(dx**2+dy**2+dz**2)

                      if (d < min_d) min_d = d

                      if (d<radius) then
                          num_pair=num_pair+1
                          
                          pair(1,num_pair)=imol1
                          pair(2,num_pair)=imol2
                          
                          image(1,num_pair)=i
                          image(2,num_pair)=j
                          image(3,num_pair)=k
                          
                          distance(num_pair)=d
                      endif
                  enddo
              enddo
          enddo
      enddo
  enddo

!   write(6, "('@@@ min_d is ', F12.6f)") min_d
  
  ! ----------Deallocate----------
  deallocate(molecule1_center)
  deallocate(molecule2_center)
  deallocate(supercell_range)

end subroutine find_pair


subroutine calc_molecule_center(num_atom,ind_ref_atom,x_frac,center_frac)
!======================================================================
! added_by_fx
! 
! Input the fractional coords of a molecule, and output the 
! fractional coord of the molecule center. 
! Due to PBC, a reference atom is needed to cluster the whole molecule 
! together first. Pleace choose a reference atom located near the 
! center of the molecule.
! PLEASE BE ATTENTION!!! x_frac is an IO! The output x_frac are the 
! coords being clustered, and may be out of the range of [0, 1).
! 
! I     num_atom
! I     ind_ref_atom
! IO    x_frac: fractional atomic positions
! O     center_frac: fractional center of the molecule
!======================================================================
  implicit none

  integer num_atom
  real*8 x_frac(3,100)
  integer ind_ref_atom
  real*8 center_frac(3)
  
  integer i,ixyz
  real*8 x1,x2
  
  center_frac=0.d0
  
  do i=1,num_atom
      do ixyz=1,3
          x1=x_frac(ixyz,i)
          x2=x_frac(ixyz,ind_ref_atom)
          
          do while (abs(x1-1-x2).lt.abs(x1-x2))
              x1=x1-1
          enddo
          do while (abs(x1+1-x2).lt.abs(x1-x2))
              x1=x1+1
          enddo
          
          x_frac(ixyz,i)=x1
          
          center_frac(ixyz)=center_frac(ixyz)+x1
      enddo
  enddo
  
  center_frac=center_frac/num_atom

end subroutine calc_molecule_center
