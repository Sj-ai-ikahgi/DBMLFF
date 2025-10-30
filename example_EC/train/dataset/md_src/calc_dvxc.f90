      subroutine calc_dvxc(rho,dvxc,n1,n2,n3)

      implicit double precision (a-h,o-z)


      real*8 rho(n1,n2,n3)
      real*8 dvxc(n1,n2,n3)


      f53=5.d0/3.d0
      f23=2.d0/3.d0
      pi=4*datan(1.d0)
      fact=3/10.d0*(3*pi**2)**(2.d0/3)
      dvxc=0.d0

      do k=1,n3
      do j=1,n2
      do i=1,n1
      vxc_tmp1=UxcCA(rho(i,j,k),uxc)
      vxc_tmp1=vxc_tmp1+fact*f53*abs(rho(i,j,k))**f23
      vxc_tmp2=UxcCA(rho(i,j,k)+1.E-4,uxc_tmp)
      vxc_tmp2=vxc_tmp2+fact*f53*abs(rho(i,j,k)+1.E-4)**f23
      dvxc(i,j,k)=(vxc_tmp2-vxc_tmp1)/1.E-4 
      enddo
      enddo
      enddo

      return
      end subroutine calc_dvxc


