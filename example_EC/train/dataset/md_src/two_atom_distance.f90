! subroutine two_atom_distance(natom, xatom, AL, iat, jat, dr, d)
!    implicit none

!    integer natom
!    real*8 xatom(3,natom)
!    real*8 AL(3,3)
!    integer iat, jat
!    real*8 dr(3)
!    real*8 d

!    real*8 dx1, dx2, dx3
!    real*8 dd

!    dr = 0.d0
!    d = 0.d0

!    dx1=xatom(1,jat)-xatom(1,iat)
!    if(abs(dx1+1).lt.abs(dx1)) dx1=dx1+1
!    if(abs(dx1-1).lt.abs(dx1)) dx1=dx1-1

!    dx2=xatom(2,jat)-xatom(2,iat)
!    if(abs(dx2+1).lt.abs(dx2)) dx2=dx2+1
!    if(abs(dx2-1).lt.abs(dx2)) dx2=dx2-1

!    dx3=xatom(3,jat)-xatom(3,iat)
!    if(abs(dx3+1).lt.abs(dx3)) dx3=dx3+1
!    if(abs(dx3-1).lt.abs(dx3)) dx3=dx3-1

!    dr(1)=AL(1,1)*dx1+AL(1,2)*dx2+AL(1,3)*dx3
!    dr(2)=AL(2,1)*dx1+AL(2,2)*dx2+AL(2,3)*dx3
!    dr(3)=AL(3,1)*dx1+AL(3,2)*dx2+AL(3,3)*dx3

!    dd=dr(1)**2+dr(2)**2+dr(3)**2

!    d = dsqrt(dd)

! end subroutine

subroutine two_coord_distance(coord1, coord2, dd)
   implicit none 

   real*8 coord1(3), coord2(3), dd
   real*8 dr(3)
   integer i

   dr(1:3) = coord2(1:3) - coord1(1:3)

   do i = 1, 3
      if (abs(dr(i)+1.0) .lt. abs(dr(i))) dr(i) = dr(i) + 1.0
      if (abs(dr(i)-1.0) .lt. abs(dr(i))) dr(i) = dr(i) - 1.0
   enddo

   dd=dr(1)**2+dr(2)**2+dr(3)**2

end subroutine two_coord_distance


subroutine two_atom_distance(natom, xatom, AL, iat, jat, dr, d)
   implicit none

   integer natom
   real*8 xatom(3,natom)
   real*8 AL(3,3)
   integer iat, jat
   real*8 dr(3)
   real*8 d

   real*8 dx1, dx2, dx3
   real*8 dd

   dr = 0.d0
   d = 0.d0

   dx1=xatom(1,jat)-xatom(1,iat)
   dx2=xatom(2,jat)-xatom(2,iat)
   dx3=xatom(3,jat)-xatom(3,iat)

   dr(1)=AL(1,1)*dx1+AL(1,2)*dx2+AL(1,3)*dx3
   dr(2)=AL(2,1)*dx1+AL(2,2)*dx2+AL(2,3)*dx3
   dr(3)=AL(3,1)*dx1+AL(3,2)*dx2+AL(3,3)*dx3

   dd=dr(1)**2+dr(2)**2+dr(3)**2

   d = dsqrt(dd)

end subroutine




subroutine Add_Neigh(m_neigh, natom, neigh_numR, neigh_listR, iat, jat, ind)

   implicit none

   integer m_neigh, natom
   integer neigh_numR(natom)
   integer neigh_listR(m_neigh, natom)
   integer iat, jat
   integer ind

   integer j, tjat

   ind = -1
   do j = 1, neigh_numR(iat)
      tjat = neigh_listR(j, iat)
      if (tjat .eq. jat) then
         ind = j
         return
      endif
   enddo

   ! if (.not. exist) then

   neigh_numR(iat) = neigh_numR(iat) + 1
   neigh_listR(neigh_numR(iat), iat) = jat
   ind = neigh_numR(iat)

end subroutine Add_Neigh



subroutine Get_Neigh_Ind(m_neigh, natom, neigh_numR, neigh_listR, iat, jat, ind)

   implicit none

   integer m_neigh, natom
   integer neigh_numR(natom)
   integer neigh_listR(m_neigh, natom)
   integer iat, jat
   integer ind

   integer j, tjat

   ind = -1
   do j = 1, neigh_numR(iat)
      tjat = neigh_listR(j, iat)
      if (tjat .eq. jat) then
         ind = j
         return
      endif
   enddo


end subroutine Get_Neigh_Ind


subroutine mol_atom_distance(natom, xatom, AL, mol1_atom_num, min_dist, max_dist)
   implicit none

   integer natom
   real*8 xatom(3,natom)
   real*8 AL(3,3)
   integer mol1_atom_num
   real*8 min_dist, max_dist

   integer i, j
   real*8 dr(3), d

   min_dist = 1000.d0
   max_dist = 0.d0

   do i = 1, mol1_atom_num
      do j = mol1_atom_num + 1, natom

         call two_atom_distance(natom, xatom, AL, i, j, dr, d)

         if (d < min_dist) min_dist = d
         if (d > max_dist) max_dist = d
      enddo
   enddo

end subroutine mol_atom_distance


subroutine move_atom_together(natoms, coords, caId)
   implicit none

   integer natoms
   real*8 coords(3,natoms)
   integer caId

   integer i, j
   real*8 dr

   do i = 1, natoms
      do j = 1, 3
         dr = coords(j,i) - coords(j,caId)
         if(abs(dr+1).lt.abs(dr)) coords(j,i) = coords(j,i) + 1
         if(abs(dr-1).lt.abs(dr)) coords(j,i) = coords(j,i) - 1
      enddo
   enddo

end subroutine move_atom_together



subroutine two_mol_min_dist(mol1_atom_num, mol2_atom_num, coord1, coord2, caId1, caId2)
   implicit none

   integer mol1_atom_num, mol2_atom_num
   real*8 coord1(3,mol1_atom_num)
   real*8 coord2(3,mol2_atom_num)
   integer caId1, caId2

   real*8 min_d, dd, dr(3)
   integer i, j

   min_d = 1000.0

   do i = 1, mol1_atom_num
      do j = 1, mol2_atom_num

         call two_coord_distance(coord1(1:3,i), coord2(1:3,j), dd)
      
         if (dd < min_d) then
            min_d = dd
            caId1 = i 
            caId2 = j
         endif
      enddo
   enddo

end subroutine two_mol_min_dist



! subroutine move_two_mol_together(mol1_atom_num, mol2_atom_num, natom, xatom)
!    implicit none

!    integer mol1_atom_num, mol2_atom_num, natom
!    real*8 AL(3,3)
!    real*8 xatom(3,natom)

!    real*8 coord1(3,mol1_atom_num), coord2(3,mol2_atom_num)
!    integer caId1, caId2
!    integer i, ii
!    real*8 dx

!    ! get two mol coordinate
!    do i = 1, mol1_atom_num
!       coord1(:,i) = xatom(:,i)
!    enddo

!    do i = 1, mol2_atom_num
!       ii = i+mol1_atom_num
!       coord2(:,i) = xatom(:,ii)
!    enddo

!    ! get closest atom index
!    call two_mol_min_dist(mol1_atom_num, mol2_atom_num, coord1, coord2, caId1, caId2)

!    ! move the two cloeset atom together
!    do i = 1, 3
!       dx = coord2(i,caId2) - coord1(i,caId1)
!       if(abs(dx+1).lt.abs(dx)) coord2(i,caId2) = coord2(i,caId2) + 1
!       if(abs(dx-1).lt.abs(dx)) coord2(i,caId2) = coord2(i,caId2) - 1
!    enddo

!    ! move atom in each mol together
!    call move_atom_together(mol1_atom_num, coord1, caId1)
!    call move_atom_together(mol2_atom_num, coord2, caId2)

!    ! save mol position
!    do i = 1, mol1_atom_num
!       xatom(:,i) = coord1(:,i)
!    enddo

!    do i = 1, mol2_atom_num
!       ii = i + mol1_atom_num
!       xatom(:,ii) = coord2(:,i)
!    enddo

! end subroutine move_two_mol_together


subroutine mol_center(mol_atom_num, coord, center_coord)
   implicit none 

   integer mol_atom_num
   real*8 coord(3,mol_atom_num), center_coord(3)
   integer i,j

   center_coord(1:3) = 0.d0

   do i = 1, mol_atom_num
      center_coord(1:3) = center_coord(1:3) + coord(1:3,i)
   enddo

   center_coord(1:3) = center_coord(1:3)/mol_atom_num

end subroutine mol_center


subroutine move_two_mol_together(mol1_atom_num, mol2_atom_num, natom, xatom)
   implicit none

   integer mol1_atom_num, mol2_atom_num, natom
   real*8 AL(3,3)
   real*8 xatom(3,natom)

   real*8 coord1(3,mol1_atom_num), coord2(3,mol2_atom_num)
   integer caId1, caId2
   integer i, j, ii
   real*8 dx, center_coord1(3), center_coord2(3), flip(3)

   ! get two mol coordinate
   do i = 1, mol1_atom_num
      coord1(:,i) = xatom(:,i)
   enddo

   do i = 1, mol2_atom_num
      ii = i+mol1_atom_num
      coord2(:,i) = xatom(:,ii)
   enddo

   caId1 = 1
   caId2 = 1

   ! move atom in each mol together
   call move_atom_together(mol1_atom_num, coord1, caId1)
   call move_atom_together(mol2_atom_num, coord2, caId2)

   call mol_center(mol1_atom_num, coord1, center_coord1)
   call mol_center(mol2_atom_num, coord2, center_coord2)

   flip(1:3) = 0.0
   ! move the two cloeset atom together
   do i = 1, 3
      dx = center_coord2(i) - center_coord1(i)
      if(abs(dx+1).lt.abs(dx)) flip(i) = 1.d0
      if(abs(dx-1).lt.abs(dx)) flip(i) = - 1.d0
   enddo

   do i = 1, mol2_atom_num
      do j = 1, 3
         coord2(j,i) = coord2(j,i) + flip(j)
      enddo
   enddo

   ! save mol position
   do i = 1, mol1_atom_num
      xatom(:,i) = coord1(:,i)
   enddo

   do i = 1, mol2_atom_num
      ii = i + mol1_atom_num
      xatom(:,ii) = coord2(:,i)
   enddo

end subroutine move_two_mol_together



