subroutine calc_coul_xc_kin_V(rho, rho_z, vxc, vcoul, n1, n2, n3, &
   AL, E_coul, E_xc, E_kin1, E_kin2, fact_kin2, vxc2)

   implicit double precision(a - h, o - z)

   parameter(matom=2000)

   real*8 ALI(3, 3), AL(3, 3), vol, vins
   complex*16 cci, cc, c_zero, c_one

   real*8 rho(n1, n2, n3), rho_z(n1, n2, n3)
   real*8 vcoul(n1, n2, n3), vxc(n1, n2, n3)
   real*8 vxc2(n1, n2, n3)

   real*8, allocatable, dimension(:, :, :) :: workr, worki
   real*8, allocatable, dimension(:, :, :) :: psi
   complex*16, allocatable, dimension(:, :, :) :: rho_q, psi_q

   real*8, allocatable, dimension(:) :: wrk

   vol = al(3, 1)*(al(1, 2)*al(2, 3) - al(1, 3)*al(2, 2)) &
       + al(3, 2)*(al(1, 3)*al(2, 1) - al(1, 1)*al(2, 3)) &
       + al(3, 3)*(al(1, 1)*al(2, 2) - al(1, 2)*al(2, 1))

   vol = abs(vol)
   vins = 1/vol

   call get_ALI(AL, ALI)  ! strange compiler error

   sum1 = 0.d0
   sum2 = 0.d0
   do k = 1, n3
       do j = 1, n2
           do i = 1, n1
               sum1 = sum1 + rho(i, j, k)
               sum2 = sum2 + rho_z(i, j, k)
           end do
       end do
   end do
   sum1 = sum1*vol/(n1*n2*n3)
   sum2 = sum2*vol/(n1*n2*n3)
   !      write(6,*) "input rho,rho_z charge",sum1,sum2

   allocate (workr(n1, n2, n3))
   allocate (worki(n1, n2, n3))
   allocate (rho_q(n1, n2, n3))
   allocate (psi(n1, n2, n3))
   allocate (psi_q(n1, n2, n3))

   lwrk = 6*(n1 + n2 + n3) + 15
   allocate (wrk(lwrk))

   workr = rho_z
   worki = 0.d0

   call cfft(n1, n2, n3, workr, worki, wrk, lwrk, -1)

   rho_q = dcmplx(workr, worki)

   psi = dsqrt(abs(rho)) ! or set to zero for negative rho

   workr = psi
   worki = 0.d0

   call cfft(n1, n2, n3, workr, worki, wrk, lwrk, -1)

   psi_q = dcmplx(workr, worki)

   cci = dcmplx(0.d0, 1.d0)

   vins = 1.d0/vol
   pi = 4*datan(1.d0)
   gkx = dcmplx(0.d0, 0.d0)
   gky = dcmplx(0.d0, 0.d0)
   gkz = dcmplx(0.d0, 0.d0)

   !ccccccccccccc perhaps we should precalculate gkx,gky,gkz?
   q2m = 1000.d0

   do i = 1, n1
       i1 = i - 1
       if (i .gt. n1/2) i1 = i - n1 - 1
       do j = 1, n2
           j1 = j - 1
           if (j .gt. n2/2) j1 = j - n2 - 1
           do k = 1, n3
               k1 = k - 1
               if (k .gt. n3/2) k1 = k - n3 - 1

               gkx = 2*pi*(ALI(1, 1)*i1 + ALI(1, 2)*j1 + ALI(1, 3)*k1)
               gky = 2*pi*(ALI(2, 1)*i1 + ALI(2, 2)*j1 + ALI(2, 3)*k1)
               gkz = 2*pi*(ALI(3, 1)*i1 + ALI(3, 2)*j1 + ALI(3, 3)*k1)
               q2 = gkx**2 + gky**2 + gkz**2

               if (q2 .gt. 1.D-20) then
                   rho_q(i, j, k) = 4*pi/q2*rho_q(i, j, k)
               else
                   rho_q(i, j, k) = 0.d0
               end if

               psi_q(i, j, k) = 0.5*q2*psi_q(i, j, k)

           end do
       end do
   end do

   workr = real(rho_q)
   worki = aimag(rho_q)

   call cfft(n1, n2, n3, workr, worki, wrk, lwrk, 1)

   vcoul = workr

   workr = real(psi_q)
   worki = aimag(psi_q)

   call cfft(n1, n2, n3, workr, worki, wrk, lwrk, 1)

   psi = workr

   f53 = 5.d0/3.d0
   f23 = 2.d0/3.d0
   pi = 4*datan(1.d0)
   fact = 3/10.d0*(3*pi**2)**(2.d0/3)
   E_coul = 0.d0
   E_xc = 0.d0
   E_kin = 0.d0
   E_kin2 = 0.d0
   vxc = 0.d0
   vxc2 = 0.d0

   do k = 1, n3
       do j = 1, n2
           do i = 1, n1
               E_coul = E_coul + vcoul(i, j, k)*rho_z(i, j, k)
               vxc(i, j, k) = UxcCA(rho(i, j, k), uxc)
               E_xc = E_xc + uxc*rho(i, j, k)
               E_kin = E_kin + fact*abs(rho(i, j, k))**f53
               E_kin2 = E_kin2 + dsqrt(abs(rho(i, j, k)))*psi(i, j, k)
               vxc(i, j, k) = vxc(i, j, k) + fact*f53*abs(rho(i, j, k))**f23

               rhom = abs(rho(i, j, k))
               if (rhom .lt. 0.0001) rhom = 0.0001
               vxc2(i, j, k) = vxc(i, j, k)   ! vxc2 is used for polarization
               vxc(i, j, k) = vxc(i, j, k) + fact_kin2*psi(i, j, k)/dsqrt(rhom)
           end do
       end do
   end do
   E_coul = 0.5*E_coul*vol/(n1*n2*n3)
   E_xc = E_xc*vol/(n1*n2*n3)
   E_kin = E_kin*vol/(n1*n2*n3)
   E_kin2 = E_kin2*vol/(n1*n2*n3)

   E_kin1 = E_kin
   E_kin2 = E_kin2

   deallocate (workr)
   deallocate (worki)
   deallocate (rho_q)
   deallocate (wrk)
   deallocate (psi)
   deallocate (psi_q)

2000 continue

    return
end subroutine calc_coul_xc_kin_V

