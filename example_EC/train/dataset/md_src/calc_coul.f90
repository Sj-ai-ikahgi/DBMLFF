      subroutine calc_coul(rho,vcoul,n1,n2,n3,AL)

      implicit double precision (a-h,o-z)

!      parameter (n1=328,n2=328,n3=328,natom=216)
!      integer, parameter :: n1=328, n2=328, n3=328, natom=216

      real*8 ALI(3,3),AL(3,3),vol,vins
      complex*16 cci,cc,c_zero,c_one

      real*8 rho(n1,n2,n3)
      real*8 vcoul(n1,n2,n3)

      real*8, allocatable,dimension(:,:,:) :: workr,worki
      real*8, allocatable,dimension(:,:,:) :: psi
      complex*16, allocatable,dimension(:,:,:) :: rho_q,psi_q

      real*8,allocatable,dimension(:) :: wrk


      vol=al(3,1)*(al(1,2)*al(2,3)-al(1,3)*al(2,2))   &
          +al(3,2)*(al(1,3)*al(2,1)-al(1,1)*al(2,3))  &
         +al(3,3)*(al(1,1)*al(2,2)-al(1,2)*al(2,1))

      vol=abs(vol)
      vins=1/vol

      call get_ALI(AL,ALI)  ! strange compiler error



      allocate(workr(n1,n2,n3))
      allocate(worki(n1,n2,n3))
      allocate(rho_q(n1,n2,n3))

      lwrk=6*(n1+n2+n3)+15
      allocate(wrk(lwrk))


      workr=rho
      worki=0.d0

      call cfft(n1,n2,n3,workr,worki,wrk,lwrk,-1)

      rho_q=cmplx(workr,worki)


      cci=dcmplx(0.d0,1.d0)

      vins=1.d0/vol
      pi=4*datan(1.d0)   
      gkx=dcmplx(0.d0,0.d0)
      gky=dcmplx(0.d0,0.d0)
      gkz=dcmplx(0.d0,0.d0)


!ccccccccccccc perhaps we should precalculate gkx,gky,gkz?
      q2m=1000.d0

      do i=1,n1
       i1=i-1
       if(i .gt. n1/2) i1=i-n1-1
      do j=1,n2
       j1=j-1
       if(j .gt. n2/2) j1=j-n2-1
      do k=1,n3
       k1=k-1
       if(k .gt. n3/2) k1=k-n3-1

         gkx=2*pi*(ALI(1,1)*i1+ALI(1,2)*j1+ALI(1,3)*k1)     
         gky=2*pi*(ALI(2,1)*i1+ALI(2,2)*j1+ALI(2,3)*k1)     
         gkz=2*pi*(ALI(3,1)*i1+ALI(3,2)*j1+ALI(3,3)*k1)   
         q2=gkx**2+gky**2+gkz**2

         if(q2.gt.1.D-20) then
         rho_q(i,j,k)=4*pi/q2*rho_q(i,j,k)
         else
         rho_q(i,j,k)=0.d0
         endif

      enddo
      enddo
      enddo 

      workr=real(rho_q)
      worki=aimag(rho_q)
 
      call cfft(n1,n2,n3,workr,worki,wrk,lwrk,1)

      vcoul=workr


      deallocate(workr)
      deallocate(worki)
      deallocate(rho_q)
      deallocate(wrk)


2000  continue

      return
      end subroutine calc_coul


