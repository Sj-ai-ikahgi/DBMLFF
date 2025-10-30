module mod_calc_pol_dir
    
    ! PLEASE pay special attention to the meaning of dxpp1_diat2, dxpp2_diat4, dxpp3_diat2, dxpp3_diat4
    ! dxpp3_diat2(1,2) is dxpp(1,3)/dx22
    ! dxpp3_diat4(2,3) is dxpp(2,3)/dx43
    ! Misusing the first two indexes cause trouble and waste time
    
    use mod_data, only: natom_mm,natom_m
    use mod_param_densityFF, only: numb,ipol,indb
    
    implicit none
    
    real*8,allocatable,dimension(:,:,:,:) :: dxpp
    real*8,allocatable,dimension(:,:,:,:) :: dxpp1_diat2
    real*8,allocatable,dimension(:,:,:,:) :: dxpp2_diat4
    real*8,allocatable,dimension(:,:,:,:) :: dxpp3_diat2
    real*8,allocatable,dimension(:,:,:,:) :: dxpp3_diat4
    
contains
    subroutine calc_pol_dir(itype_mol,xatom,AL)
        implicit none
        
        integer itype_mol
        real*8 xatom(3,natom_m(itype_mol))
        real*8 AL(3,3)
        
        integer max_natom_numb
        
        integer i,j
        integer iatom,ibond
        integer iat1,iat2,iat3,iat4
        
        real*8 dx1_frac(3),dx2_frac(3)
        real*8 dx1(3),dx2(3),dx3(3)
        real*8 d1,d2,d3
        
        max_natom_numb=max(natom_m(itype_mol),numb(itype_mol))
        
        if(allocated(dxpp)) then
            deallocate(dxpp)
            deallocate(dxpp1_diat2)
            deallocate(dxpp2_diat4)
            deallocate(dxpp3_diat2)
            deallocate(dxpp3_diat4)
        endif
        allocate(dxpp(3,3,max_natom_numb,2))
        allocate(dxpp1_diat2(3,3,max_natom_numb,2))
        allocate(dxpp2_diat4(3,3,max_natom_numb,2))
        allocate(dxpp3_diat2(3,3,max_natom_numb,2))
        allocate(dxpp3_diat4(3,3,max_natom_numb,2))
        dxpp=0.d0
        dxpp1_diat2=0.d0
        dxpp2_diat4=0.d0
        dxpp3_diat2=0.d0
        dxpp3_diat4=0.d0
        
        do iatom=1,natom_m(itype_mol)
            iat1=ipol(1,iatom,itype_mol)
            iat2=ipol(2,iatom,itype_mol)
            iat3=ipol(3,iatom,itype_mol)
            iat4=ipol(4,iatom,itype_mol)
            
            do i=1,3
                dx1_frac(i)=xatom(i,iat2)-xatom(i,iat1)
                dx2_frac(i)=xatom(i,iat4)-xatom(i,iat3)
                do while(abs(dx1_frac(i)+1).lt.abs(dx1_frac(i)))
                    dx1_frac(i)=dx1_frac(i)+1
                enddo
                do while(abs(dx1_frac(i)-1).lt.abs(dx1_frac(i)))
                    dx1_frac(i)=dx1_frac(i)-1
                enddo
                do while(abs(dx2_frac(i)+1).lt.abs(dx2_frac(i)))
                    dx2_frac(i)=dx2_frac(i)+1
                enddo
                do while(abs(dx2_frac(i)-1).lt.abs(dx2_frac(i)))
                    dx2_frac(i)=dx2_frac(i)-1
                enddo
            enddo
            
            dx1(:)=AL(:,1)*dx1_frac(1)+AL(:,2)*dx1_frac(2)+AL(:,3)*dx1_frac(3)
            dx2(:)=AL(:,1)*dx2_frac(1)+AL(:,2)*dx2_frac(2)+AL(:,3)*dx2_frac(3)
            
            dx3(1)=dx1(2)*dx2(3)-dx1(3)*dx2(2)
            dx3(2)=dx1(3)*dx2(1)-dx1(1)*dx2(3)   
            dx3(3)=dx1(1)*dx2(2)-dx1(2)*dx2(1)   
            
            d1=dsqrt(dx1(1)**2+dx1(2)**2+dx1(3)**2) 
            do i=1,3
                do j=1,3
                    dxpp1_diat2(i,j,iatom,1)=-dx1(i)*dx1(j)/d1**3
                enddo
                dxpp1_diat2(i,i,iatom,1)=dxpp1_diat2(i,i,iatom,1)+1/d1
            enddo
            
            d2=dsqrt(dx2(1)**2+dx2(2)**2+dx2(3)**2) 
            do i=1,3
                do j=1,3
                    dxpp2_diat4(i,j,iatom,1)=-dx2(i)*dx2(j)/d2**3
                enddo
                dxpp2_diat4(i,i,iatom,1)=dxpp2_diat4(i,i,iatom,1)+1/d2
            enddo
            
            d3=dsqrt(dx3(1)**2+dx3(2)**2+dx3(3)**2) 
            do i=1,3
                dxpp3_diat2(i,1,iatom,1)=-dx3(i)/d3**3*(-dx3(2)*dx2(3)+dx3(3)*dx2(2))  
                dxpp3_diat2(i,2,iatom,1)=-dx3(i)/d3**3*(dx3(1)*dx2(3)-dx3(3)*dx2(1))
                dxpp3_diat2(i,3,iatom,1)=-dx3(i)/d3**3*(-dx3(1)*dx2(2)+dx3(2)*dx2(1))
            
                dxpp3_diat4(i,1,iatom,1)=-dx3(i)/d3**3*(dx3(2)*dx1(3)-dx3(3)*dx1(2))
                dxpp3_diat4(i,2,iatom,1)=-dx3(i)/d3**3*(-dx3(1)*dx1(3)+dx3(3)*dx1(1))
                dxpp3_diat4(i,3,iatom,1)=-dx3(i)/d3**3*(dx3(1)*dx1(2)-dx3(2)*dx1(1))
            enddo
            
            dxpp3_diat2(1,2,iatom,1)=dxpp3_diat2(1,2,iatom,1)+dx2(3)/d3
            dxpp3_diat2(1,3,iatom,1)=dxpp3_diat2(1,3,iatom,1)-dx2(2)/d3
            dxpp3_diat2(2,1,iatom,1)=dxpp3_diat2(2,1,iatom,1)-dx2(3)/d3
            dxpp3_diat2(2,3,iatom,1)=dxpp3_diat2(2,3,iatom,1)+dx2(1)/d3
            dxpp3_diat2(3,1,iatom,1)=dxpp3_diat2(3,1,iatom,1)+dx2(2)/d3
            dxpp3_diat2(3,2,iatom,1)=dxpp3_diat2(3,2,iatom,1)-dx2(1)/d3
            
            dxpp3_diat4(1,2,iatom,1)=dxpp3_diat4(1,2,iatom,1)-dx1(3)/d3
            dxpp3_diat4(1,3,iatom,1)=dxpp3_diat4(1,3,iatom,1)+dx1(2)/d3
            dxpp3_diat4(2,1,iatom,1)=dxpp3_diat4(2,1,iatom,1)+dx1(3)/d3
            dxpp3_diat4(2,3,iatom,1)=dxpp3_diat4(2,3,iatom,1)-dx1(1)/d3
            dxpp3_diat4(3,1,iatom,1)=dxpp3_diat4(3,1,iatom,1)-dx1(2)/d3
            dxpp3_diat4(3,2,iatom,1)=dxpp3_diat4(3,2,iatom,1)+dx1(1)/d3
            
            dx1=dx1/d1
            dx2=dx2/d2
            dx3=dx3/d3
            
            dxpp(:,1,iatom,1)=dx1
            dxpp(:,2,iatom,1)=dx2
            dxpp(:,3,iatom,1)=dx3
        enddo
        
        do ibond=1,numb(itype_mol)
            iat1=indb(1,ibond,itype_mol)
            iat2=indb(2,ibond,itype_mol)
            
            do i=1,3
                dx1_frac(i)=xatom(i,iat2)-xatom(i,iat1)
                do while(abs(dx1_frac(i)+1).lt.abs(dx1_frac(i)))
                    dx1_frac(i)=dx1_frac(i)+1
                enddo
                do while(abs(dx1_frac(i)-1).lt.abs(dx1_frac(i)))
                    dx1_frac(i)=dx1_frac(i)-1
                enddo
            enddo
            
            dx1(:)=AL(:,1)*dx1_frac(1)+AL(:,2)*dx1_frac(2)+AL(:,3)*dx1_frac(3) 
            
            d1=dsqrt(dx1(1)**2+dx1(2)**2+dx1(3)**2) 
            do i=1,3
                do j=1,3
                    dxpp1_diat2(i,j,ibond,2)=-dx1(i)*dx1(j)/d1**3
                enddo
                dxpp1_diat2(i,i,ibond,2)=dxpp1_diat2(i,i,ibond,2)+1/d1
            enddo
            
            dx1=dx1/d1
            
            dxpp(:,1,ibond,2)=dx1
        enddo
    end subroutine calc_pol_dir
    
end module mod_calc_pol_dir