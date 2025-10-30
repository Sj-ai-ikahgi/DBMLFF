program gen_points
   implicit double precision (a-h,o-z)
   real*8 AL(3,3),AL2(3,3)
   real*8 ALI(3,3)
   real*8 xatom(3,200),xyz(3,200)
   real*8 xp(3),xyz0(3)
   real*8 iatom(200)
   real*8 px_at(200),py_at(200),pz_at(200)
   real*8,allocatable,dimension (:,:,:) :: vr,rho
   real*8,allocatable,dimension (:) :: vr_tmp

   open(10,file="xatom0.config")
   rewind(10)
   read(10,*) natom
   read(10,*)
   read(10,*) AL(1,1),AL(2,1),AL(3,1)
   read(10,*) AL(1,2),AL(2,2),AL(3,2)
   read(10,*) AL(1,3),AL(2,3),AL(3,3)
   read(10,*)
   do i=1,natom
      read(10,*) iatom(i),xatom(1,i),xatom(2,i),xatom(3,i)
   enddo
   close(10)

   call get_ALI(AL,ALI)

   !cccccccccccccccccccccccccccccccccccccccccccccccccc

   do i=1,natom
      x1=xatom(1,i)
      x2=xatom(2,i)
      x3=xatom(3,i)
      xyz(1,i)=AL(1,1)*x1+AL(1,2)*x2+AL(1,3)*x3
      xyz(2,i)=AL(2,1)*x1+AL(2,2)*x2+AL(2,3)*x3
      xyz(3,i)=AL(3,1)*x1+AL(3,2)*x2+AL(3,3)*x3
   enddo

   !cccccccccccccccccccccccccccccccccccccccccccccccccc
   xyz0=0.d0
   do i=1,natom
      xyz0(1)=xyz0(1)+xyz(1,i)
      xyz0(2)=xyz0(2)+xyz(2,i)
      xyz0(3)=xyz0(3)+xyz(3,i)
   enddo
   xyz0=xyz0/natom
   !cccccccccccccccccccccccccccccccccccccccccccccccccc
   pi=4*datan(1.d0)
   R0=20

   ntheta=8
   ! ntheta=10
   ! ntheta=6
   ! ntheta=4

   write(6,*) "input dcut"
   read(5,*) dcut

   open(12,file="point.temp")
   rewind(12)
   num=0
   do 100 i=0,ntheta
      do 100 j=0,i*4
         do 100 isign=-1,1,2
            !       if(i.eq.10.and.isign.eq.1) goto 100
            theta=i*pi/2/ntheta
            if(i.eq.0) then
               phi=0
            else
               phi=j*2*pi/(i*4)
            endif
            x0=R0*sin(theta)*cos(phi)
            y0=R0*sin(theta)*sin(phi)
            z0=isign*R0*cos(theta)

            do ii=200,1,-1
               x=x0*ii/200.d0+xyz0(1)
               y=y0*ii/200.d0+xyz0(2)
               z=z0*ii/200.d0+xyz0(3)

               dd=(x-xyz0(1))**2+(y-xyz0(2))**2+(z-xyz0(3))**2
               !       if(dd.lt.2.0**2) goto 200
               if(dd.lt.dcut**2) goto 200
               !       if(dd.lt.8.0**2) goto 200

               do iat=1,natom
                  dd=(x-xyz(1,iat))**2+(y-xyz(2,iat))**2+(z-xyz(3,iat))**2
                  if(dd.lt.1.8**2) goto 200
               enddo
            enddo  ! ii

200         continue
            x1=ALI(1,1)*x+ALI(2,1)*y+ALI(3,1)*z
            x2=ALI(1,2)*x+ALI(2,2)*y+ALI(3,2)*z
            x3=ALI(1,3)*x+ALI(2,3)*y+ALI(3,3)*z

            num=num+1
            write(12,*) 1,num   ! strength
            write(12,"(2(E13.5,1x),2x,3(E20.13,1x))") &
               0.50d0,1.4,x1,x2,x3
            num=num+1
            write(12,*) 1,num   ! strength
            write(12,"(2(E13.5,1x),2x,3(E20.13,1x))") &
               0.50d0,3.4,x1,x2,x3
              !       write(12,"(2(E13.5,1x),2x,3(E20.13,1x))") &
              !                    1.00d0,4.4,x1,x2,x3
              !       num=num+1
              !       write(12,*) 1,num   ! strength
              !       write(12,"(2(E13.5,1x),2x,3(E20.13,1x))") &
              !                    0.50d0,3.4,x1,x2,x3
100 continue
   close(12)
   stop
end
