subroutine find_feature_2b_3b_type3(natom, xatom, AL, ntype, m_neigh, itype_atom, &
   n2b, n2bm, grid2_2, &
   n3b1, n3b2, n3b1m, n3b2m, grid31_2, grid32_2, &
   nfeat0m2, nfeat0m3,&
   neigh_cutoff, mol1_type, atom_mol_type, atom_ind_mol, mol1_atom_num, &
   nfeat2_atom, feat2_all, dfeat2_allR, &
   nfeat3_atom, feat3_all, dfeat3_allR, &
   neigh_numR, neigh_listR)

   implicit none

   integer natom
   real*8 xatom(3,natom)
   real*8 AL(3,3)
   integer ntype
   integer m_neigh
   integer itype_atom(natom)
   integer n2b(ntype)
   integer n2bm
   real*8 grid2_2(2, n2bm + 1, ntype)
   integer n3b1(ntype), n3b2(ntype)
   integer n3b1m, n3b2m
   real*8 grid31_2(2, n3b1m, ntype), grid32_2(2, n3b2m, ntype)
   integer nfeat0m2, nfeat0m3
   integer nfeat2_atom(natom)
   integer nfeat3_atom(natom)
   real*8 feat2_all(nfeat0m2, natom), dfeat2_allR(nfeat0m2, natom, m_neigh, 3)
   real*8 feat3_all(nfeat0m3, natom), dfeat3_allR(nfeat0m3, natom, m_neigh, 3)
   integer neigh_numR(natom)
   integer neigh_listR(m_neigh,natom)

   real*8 feat2(n2bm, ntype, natom)
   real*8 dfeat2(n2bm, ntype, natom, m_neigh, 3)
   real*8 feat3(n3b1m*n3b1m*n3b2m, ntype*(ntype + 1)/2, natom)
   real*8 dfeat3(n3b1m*n3b1m*n3b2m, ntype*(ntype + 1)/2, natom, m_neigh, 3)


   real*8 pi, pi2
   integer i, j, k, ii, jj, ii1, ii2
   integer iat, jat, jat1, jat2, tind, jj1, jj2
   integer itype0, jtype, jtype1, jtype2, jtype12
   real*8 dr12(3), dr13(3), dr23(3), d
   real*8 x, f1, y, y2
   integer num, num_k12, num_k13, num_k23
   real*8 ind_f12(n3b1m), ind_f13(n3b1m), ind_f32(n3b2m)
   real*8 f12(n3b1m), f13(n3b1m), f32(n3b2m)
   real*8 df12(n3b1m,3), df13(n3b1m,3), df32(n3b2m,2,3)
   integer i1, i2, j12, k1, k2, k12, ii_f
   integer nneigh
   real*8 dfeat2_all(nfeat0m2, natom, m_neigh, 3)
   real*8 dfeat3_all(nfeat0m3, natom, m_neigh, 3)

   integer mol1_type, mol2_type
   integer mol1_atom_num, mol2_atom_num
   real*8  neigh_cutoff
   integer ipair
   integer atom_mol_type(natom)
   integer atom_ind_mol(natom)
   integer bond_neigh_num(natom,4)
   integer bond_neigh_list(m_neigh,natom,4)

   integer mol_type_jat1, mol_atom_ind_ja1, mol_atom_ind_ja2, bid



   pi = 4*datan(1.d0)
   pi2 = 2*pi

   nfeat2_atom = 0
   nfeat3_atom = 0
   feat2_all = 0.d0
   feat3_all = 0.d0
   dfeat2_allR = 0.d0
   dfeat3_allR = 0.d0

   neigh_numR = 0
   neigh_listR = 0

   feat2 = 0.d0
   feat3 = 0.d0
   dfeat2 = 0.d0
   dfeat3 = 0.d0

   ! write(6,*) "@@@ mol1_type ", mol1_type
   open(13, file="bond_neigh_list."//char(mol1_type + 48), status="old", action="read")
   rewind(13)
   bond_neigh_num = 0
   bond_neigh_list = 0
   do i = 1, mol1_atom_num
      read(13,*) ii, num, bond_neigh_list(1:num, i, mol1_type)
      bond_neigh_num(i, mol1_type) = num
   enddo
   close(13)

   ! initial neigh list
   do iat = 1, natom
      call Add_Neigh(m_neigh, natom, neigh_numR, neigh_listR, iat, iat, tind)
   enddo

   ! start calculate 2b & 3b feature
   do 3000 iat = mol1_atom_num + 1, natom

      itype0 = itype_atom(iat)

      do 2000 jat1 = 1, natom

         if (iat .eq. jat1) cycle

         jtype1 = itype_atom(jat1)

         ! dr12 = r(jat1) - r(iat)
         call two_atom_distance(natom, xatom, AL, iat, jat1, dr12, d)
         if (d .gt. neigh_cutoff) cycle

         call Add_Neigh(m_neigh, natom, neigh_numR, neigh_listR, iat, jat1, jj1)
         call Add_Neigh(m_neigh, natom, neigh_numR, neigh_listR, jat1, iat, ii1)

         ! 2b feature
         do k = 1, n2b(itype0)
            if (d .ge. grid2_2(1, k, itype0) .and. d .lt. grid2_2(2, k, itype0)) then
               x = (d - grid2_2(1, k, itype0))/(grid2_2(2, k, itype0) - grid2_2(1, k, itype0))
               y = (x - 0.5d0)*pi2
               f1 = 0.5d0*(cos(y) + 1)
               feat2(k, jtype1, iat) = feat2(k, jtype1, iat) + f1
               y2 = -pi*sin(y)/(d*(grid2_2(2, k, itype0) - grid2_2(1, k, itype0)))
               dfeat2(k, jtype1, iat, jj1, :) = dfeat2(k, jtype1, iat, jj1, :) + y2*dr12(:)
               dfeat2(k, jtype1, iat, 1, :) = dfeat2(k, jtype1, iat, 1, :) - y2*dr12(:)
            end if
         end do   ! k=1,n2b


         ! following is for 3b featrure
         num = 0
         do k = 1, n3b1(itype0)
            if (d .ge. grid31_2(1, k, itype0) .and. d .lt. grid31_2(2, k, itype0)) then
               num = num + 1
               x = (d - grid31_2(1, k, itype0))/(grid31_2(2, k, itype0) - grid31_2(1, k, itype0))
               y = (x - 0.5d0)*pi2
               f1 = 0.5d0*(cos(y) + 1)
               f12(num) = f1
               ind_f12(num) = k
               y2 = -pi*sin(y)/(d*(grid31_2(2, k, itype0) - grid31_2(1, k, itype0)))
               df12(num, :) = y2*dr12(:)
            end if
         end do
         num_k12 = num


         if (atom_mol_type(iat) .ne. atom_mol_type(jat1)) then
            ! iat and jat1 not in the same mol (jat1 in the neigh mol)

            mol_type_jat1 = atom_mol_type(jat1)
            mol_atom_ind_ja1 = atom_ind_mol(jat1)

            do  bid = 1, bond_neigh_num(mol_atom_ind_ja1, mol_type_jat1)

               ! because the neigh mol is the first mol
               ! 1 <= jat2 <= mol1_atom_num
               mol_atom_ind_ja2 = bond_neigh_list(bid, mol_atom_ind_ja1, mol_type_jat1)
               jat2 = mol_atom_ind_ja2

               call Add_Neigh(m_neigh, natom, neigh_numR, neigh_listR, iat, jat2, jj2)
               call Add_Neigh(m_neigh, natom, neigh_numR, neigh_listR, jat2, iat, ii2)

               jtype2 = itype_atom(jat2)

               call two_atom_distance(natom, xatom, AL, iat, jat2, dr13, d)
               num = 0
               do k = 1, n3b1(itype0)
                  if (d .ge. grid31_2(1, k, itype0) .and. d .lt. grid31_2(2, k, itype0)) then
                     num = num + 1
                     x = (d - grid31_2(1, k, itype0))/(grid31_2(2, k, itype0) - grid31_2(1, k, itype0))
                     y = (x - 0.5d0)*pi2
                     f1 = 0.5d0*(cos(y) + 1)
                     f13(num) = f1
                     ind_f13(num) = k
                     y2 = -pi*sin(y)/(d*(grid31_2(2, k, itype0) - grid31_2(1, k, itype0)))
                     df13(num, :) = y2*dr13(:)
                  end if
               end do
               num_k13 = num

               ! dr23 = r(jat1) - r(jat2)
               call two_atom_distance(natom, xatom, AL, jat2, jat1, dr23, d)
               num = 0
               do k = 1, n3b2(itype0)
                  if (d .ge. grid32_2(1, k, itype0) .and. d .lt. grid32_2(2, k, itype0)) then
                     num = num + 1
                     x = (d - grid32_2(1, k, itype0))/(grid32_2(2, k, itype0) - grid32_2(1, k, itype0))
                     y = (x - 0.5d0)*pi2
                     f1 = 0.5d0*(cos(y) + 1)
                     f32(num) = f1
                     ind_f32(num) = k
                     y2 = -pi*sin(y)/(d*(grid32_2(2, k, itype0) - grid32_2(1, k, itype0)))
                     df32(num, 1, :) = y2*dr23(:)
                     df32(num, 2, :) = -df32(num, 1, :)
                  end if
               end do
               num_k23 = num

               !cccccccccccccccccccccccc
               !   Each R has two k features, so for the three R, we have the following

               jtype12 = jtype1 + ((jtype2 - 1)*jtype2)/2

               do i1 = 1, num_k12
                  do i2 = 1, num_k13
                     do j12 = 1, num_k23
                        k1 = ind_f12(i1)
                        k2 = ind_f13(i2)
                        k12 = ind_f32(j12)

                        ii_f = 0
                        if (jtype1 .ne. jtype2) then
                           ii_f = k1 + (k2 - 1)*n3b1(itype0) + (k12 - 1)*n3b1(itype0)**2
                        end if
                        if (jtype1 .eq. jtype2 .and. k1 .le. k2) then
                           ii_f = k1 + ((k2 - 1)*k2)/2 + (k12 - 1)*(n3b1(itype0)*(n3b1(itype0) + 1))/2
                        end if

                        if (ii_f .ne. 0) then
                           feat3(ii_f, jtype12, iat) = feat3(ii_f, jtype12, iat) + f12(i1)*f13(i2)*f32(j12)

                           dfeat3(ii_f, jtype12, iat, jj1, :) = dfeat3(ii_f, jtype12, iat, jj1, :) + &
                              df12(i1, :)*f13(i2)*f32(j12) + f12(i1)*f13(i2)*df32(j12, 1, :)

                           dfeat3(ii_f, jtype12, iat, jj2, :) = dfeat3(ii_f, jtype12, iat, jj2, :) + &
                              f12(i1)*df13(i2, :)*f32(j12) + f12(i1)*f13(i2)*df32(j12, 2, :)

                           dfeat3(ii_f, jtype12, iat, 1, :) = dfeat3(ii_f, jtype12, iat, 1, :) - &
                              df12(i1, :)*f13(i2)*f32(j12) - f12(i1)*df13(i2, :)*f32(j12)
                           !cccc (ii_f,jtype12) is the feature index
                        end if

                     end do
                  end do
               end do
            end do ! bid

         else
            ! iat and jat1 in the same mol (center mol)

            do jat2 = 1, natom

               if ((jat2 .eq. iat) .or. (jat2 .eq. jat1)) cycle

               call two_atom_distance(natom, xatom, AL, iat, jat2, dr13, d)
               ! the second neigh should within neigh_cutoff
               if (d .gt. neigh_cutoff) cycle

               call Add_Neigh(m_neigh, natom, neigh_numR, neigh_listR, iat, jat2, jj2)
               call Add_Neigh(m_neigh, natom, neigh_numR, neigh_listR, jat2, iat, ii2)

               jtype2 = itype_atom(jat2)

               num = 0
               do k = 1, n3b1(itype0)
                  if (d .ge. grid31_2(1, k, itype0) .and. d .lt. grid31_2(2, k, itype0)) then
                     num = num + 1
                     x = (d - grid31_2(1, k, itype0))/(grid31_2(2, k, itype0) - grid31_2(1, k, itype0))
                     y = (x - 0.5d0)*pi2
                     f1 = 0.5d0*(cos(y) + 1)
                     f13(num) = f1
                     ind_f13(num) = k
                     y2 = -pi*sin(y)/(d*(grid31_2(2, k, itype0) - grid31_2(1, k, itype0)))
                     df13(num, :) = y2*dr13(:)
                  end if
               end do
               num_k13 = num

               ! dr23 = r(jat1) - r(jat2)
               call two_atom_distance(natom, xatom, AL, jat2, jat1, dr23, d)
               num = 0
               do k = 1, n3b2(itype0)
                  if (d .ge. grid32_2(1, k, itype0) .and. d .lt. grid32_2(2, k, itype0)) then
                     num = num + 1
                     x = (d - grid32_2(1, k, itype0))/(grid32_2(2, k, itype0) - grid32_2(1, k, itype0))
                     y = (x - 0.5d0)*pi2
                     f1 = 0.5d0*(cos(y) + 1)
                     f32(num) = f1
                     ind_f32(num) = k
                     y2 = -pi*sin(y)/(d*(grid32_2(2, k, itype0) - grid32_2(1, k, itype0)))
                     df32(num, 1, :) = y2*dr23
                     df32(num, 2, :) = -df32(num, 1, :)
                  end if
               end do
               num_k23 = num

               !cccccccccccccccccccccccc
               !   Each R has two k features, so for the three R, we have the following

               jtype12 = jtype1 + ((jtype2 - 1)*jtype2)/2

               do i1 = 1, num_k12
                  do i2 = 1, num_k13
                     do j12 = 1, num_k23
                        k1 = ind_f12(i1)
                        k2 = ind_f13(i2)
                        k12 = ind_f32(j12)

                        ii_f = 0
                        if (jtype1 .ne. jtype2) then
                           ii_f = k1 + (k2 - 1)*n3b1(itype0) + (k12 - 1)*n3b1(itype0)**2
                        end if
                        if (jtype1 .eq. jtype2 .and. k1 .le. k2) then
                           ii_f = k1 + ((k2 - 1)*k2)/2 + (k12 - 1)*(n3b1(itype0)*(n3b1(itype0) + 1))/2
                        end if

                        if (ii_f .ne. 0) then
                           feat3(ii_f, jtype12, iat) = feat3(ii_f, jtype12, iat) + f12(i1)*f13(i2)*f32(j12)

                           dfeat3(ii_f, jtype12, iat, jj1, :) = dfeat3(ii_f, jtype12, iat, jj1, :) + &
                              df12(i1, :)*f13(i2)*f32(j12) + f12(i1)*f13(i2)*df32(j12, 1, :)

                           dfeat3(ii_f, jtype12, iat, jj2, :) = dfeat3(ii_f, jtype12, iat, jj2, :) + &
                              f12(i1)*df13(i2, :)*f32(j12) + f12(i1)*f13(i2)*df32(j12, 2, :)

                           dfeat3(ii_f, jtype12, iat, 1, :) = dfeat3(ii_f, jtype12, iat, 1, :) - &
                              df12(i1, :)*f13(i2)*f32(j12) - f12(i1)*df13(i2, :)*f32(j12)
                           !cccc (ii_f,jtype12) is the feature index
                        end if

                     end do
                  end do
               end do
            end do ! jat2

         endif

2000  end do ! jat1

3000 end do ! iat

   !cccccccccccccccccccccccccccccccccccccccccccccccc
   !cccccccccccccccccccccccccccccccccccccccccccccccc
   !   Now, we collect everything together, collapse the index (k,jtype)
   !   and feat2,feat3, into a single feature.

   ! for 2b feature
   do iat = 1, natom
      itype0 = itype_atom(iat)
      nneigh = neigh_numR(iat)

      num = 0
      do jtype = 1, ntype
         do k = 1, n2b(itype0)
            num = num + 1
            feat2_all(num, iat) = feat2(k, jtype, iat)
            dfeat2_all(num, iat, 1:nneigh, :) = dfeat2(k, jtype, iat, 1:nneigh, :)
         end do
      end do

      nfeat2_atom(iat) = num

      if (num .gt. nfeat0m2) then
         write (6, *) "num.gt.nfeat0m2,stop", nfeat0m2, num
         stop
      end if
   enddo


   ! for 3b feature
   do iat = 1, natom
      itype0 = itype_atom(iat)
      nneigh = neigh_numR(iat)

      num = 0
      do jtype2 = 1, ntype
         do jtype1 = 1, jtype2

            jtype12 = jtype1 + ((jtype2 - 1)*jtype2)/2

            do k1 = 1, n3b1(itype0)
               do k2 = 1, n3b1(itype0)
                  do k12 = 1, n3b2(itype0)

                     ii_f = 0
                     if (jtype1 .ne. jtype2) then
                        ii_f = k1 + (k2 - 1)*n3b1(itype0) + (k12 - 1)*n3b1(itype0)**2
                     end if
                     if (jtype1 .eq. jtype2 .and. k1 .le. k2) then
                        ii_f = k1 + ((k2 - 1)*k2)/2 + (k12 - 1)*(n3b1(itype0)*(n3b1(itype0) + 1))/2
                     end if

                     if (ii_f .gt. 0) then
                        num = num + 1
                        feat3_all(num, iat) = feat3(ii_f, jtype12, iat)
                        dfeat3_all(num, iat, 1:nneigh, :) = dfeat3(ii_f, jtype12, iat, 1:nneigh, :)
                     end if

                  end do
               end do
            end do

         end do
      end do

      nfeat3_atom(iat) = num

      if (num .gt. nfeat0m3) then
         write (6, *) "num.gt.nfeat0m3,stop", num, nfeat0m3
         stop
      end if
   enddo

   !ccccccccccccccccccccccccccccccccccc
   !  Now, we have to redefine the dfeat_all in another way.
   !  dfeat_all(nfeat,iat,jneigh,3) means:
   !  d_jth_feat_of_iat/d_R(jth_neigh_of_iat)
   !  dfeat_allR(nfeat,iat,jneigh,3) means:
   !  d_jth_feat_of_jth_neigh/d_R(iat)
   !cccccccccccccccccccccccccccccccccccccc

   dfeat2_allR = 0.d0
   dfeat3_allR = 0.d0

   do iat = 1, natom
      do j = 1, neigh_numR(iat)
         ! include the one which is itself, j=1
         jat = neigh_listR(j, iat)
         call Get_Neigh_Ind(m_neigh, natom, neigh_numR, neigh_listR, jat, iat, jj)
         if (jj .eq. - 1) then
            write(6,*) "@@@ transfer effor ", iat, j, jat, jj
         endif
         do ii_f = 1, nfeat2_atom(iat)
            dfeat2_allR(ii_f, jat, jj, :) = dfeat2_all(ii_f, iat, j, :)
            !ccc Note, dfeat2_allR(i,jat,jj,3), it can have more i then nfeat2_atom(jat),
            ! since it is the nfeat of j2_neighbore
         end do
         do ii_f = 1, nfeat3_atom(iat)
            dfeat3_allR(ii_f, jat, jj, :) = dfeat3_all(ii_f, iat, j, :)
            !ccc Note, dfeat3_allR(i,jat,jj,3), it can have more i then nfeat3_atom(jat),
            ! since it is the nfeat of j2_neighbore
         end do
      end do
   end do

end subroutine find_feature_2b_3b_type3


