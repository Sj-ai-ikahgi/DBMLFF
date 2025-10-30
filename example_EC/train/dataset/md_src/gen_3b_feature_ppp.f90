module calc_ftype2_ppp

    ! ---------- USAGE!!! ----------
    ! call load_model_type2_ppp_from_mem (mod_ppp.f90)
    !     -> call set_image_info_type2_ppp
    ! call calc_pol_dir (mod_calc_pol_dir.f90)
    ! call load_neighbore (find_neighbore_ppp.f90) or load_neighbore_from_mem (mod_ppp.f90)
    !     -> call calc_dR_neigh(find_neighbore_ppp.f90)

    use mod_mpi
    use mod_md, only: A_AU_1
    use mod_param_densityFF, only: indb,ipol
    use mod_calc_pol_dir, only: dxpp,dxpp1_diat2,dxpp2_diat4,dxpp3_diat2,dxpp3_diat4
    use find_neighbore_ppp, only: num_neigh,list_neigh,dR_neigh,num_neigh_alltype, &
        list_neigh_alltype,ind_all_neigh
    
    implicit none
    
    ! ---------- load module ----------
    ! loaded from mem in load_model_type2_ppp_from_mem (mod_ppp.f90)
    ! one should call load_model_type2_ppp_from_mem (mod_ppp.f90) and then use this module
    real*8 Rc_M
    integer m_neigh
    integer ntype
    integer iat_type(100)
    real*8 Rc_type(100)
    real*8 Rc2_type(100)
    real*8 n3b1_type(100)
    real*8 n3b2_type(100)
    integer iflag_ftype
    integer n3b1m
    integer n3b2m
    integer nfeat0m
    real*8,allocatable,dimension(:,:,:) :: grid31_2
    real*8,allocatable,dimension(:,:,:) :: grid32_2
    
    ! loaded from set_image_info_type2_ppp
    integer natom
    integer,allocatable,dimension(:) :: iatom
    
    ! ---------- public ----------
    real*8,allocatable,dimension(:,:) :: feat_M2
    real*8,allocatable,dimension(:,:,:,:) :: dfeat_M2
    integer nfeat0M2
    
    real*8,allocatable,dimension(:,:,:,:) :: dfeat_iat0                 ! dfeat/dR_self=d(2b01*2b02*2b12)/dR0=(d2b01/dR0)*2b02*2b12+2b01*(d2b02/dR0)*2b12
    real*8,allocatable,dimension(:,:,:,:,:,:) :: ddfeat                 ! d(dfeat/dR_self)/dR
    
    integer,allocatable,dimension(:,:) :: num_nonzero_feat
    integer,allocatable,dimension(:,:,:) :: ind_nonzero_feat
    
    ! ---------- private ----------
    integer,allocatable,dimension(:) :: itype_atom

    private :: itype_atom
    
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
contains
    subroutine set_image_info_type2_ppp(iatom_in,natom_in)
        integer,intent(in) :: iatom_in(natom_in)
        integer,intent(in) :: natom_in
        
        integer i,j
        
        ! ----------
        
        if(allocated(iatom)) then
            deallocate(iatom)
            deallocate(itype_atom)
        endif
        allocate(iatom(natom_in))
        allocate(itype_atom(natom_in))
        
        ! ----------
        
        natom=natom_in
        iatom(1:natom)=iatom_in(1:natom)
        
        ! ----------
        
        do i=1,natom
            do j=1,ntype
                if(iatom(i).eq.iat_type(j)) then
                    itype_atom(i)=j
                endif
            enddo
            if (itype_atom(i).eq.0) then
                write(6,*) "this atom type didn't found",itype_atom(i)
                stop
            endif
        enddo
    end subroutine set_image_info_type2_ppp
    
    
    subroutine calc_orig_dfeat_and_ddfeat_3b_ppp()
        ! --------------------------------------------------
        ! added_by_fx
        ! --------------------------------------------------
        ! O,M   dfeat_iat0              ! dfeat/dR_self=d(2b01*2b02*2b12)/dR0=(d2b01/dR0)*2b02*2b12+2b01*(d2b02/dR0)*2b12
        ! O,M   ddfeat
        ! O,M   num_nonzero_feat
        ! O,M   ind_nonzero_feat
        ! --------------------------------------------------
        implicit none
        
        integer i,j,k
        integer iat,itype
        integer itype0
        integer num
        integer itype2,itype1,itype12,j1,j2,jj1,jj2
        integer i1,i2,j12,k1,k2,k12,ii_f
        
        real*8 pi,pi2
        
        real*8 d,dd
        real*8 drange,x,y,f1,y2
        
        integer num_2b0x(m_neigh,ntype,natom)
        integer ind_2b0x(n3b1m,m_neigh,ntype,natom)
        real*8 feat_2b0x(n3b1m,m_neigh,ntype)
        real*8 y2_2b0x(n3b1m,m_neigh,ntype)
        real*8 dfeat_2b0x(n3b1m,m_neigh,ntype,3)
        real*8 temp_2b0x(n3b1m,m_neigh,ntype)
        
        integer num_2b12(m_neigh,m_neigh,ntype*(ntype+1)/2,natom)
        integer ind_2b12(n3b2m,m_neigh,m_neigh,ntype*(ntype+1)/2,natom)
        real*8 feat_2b12(n3b2m)
        real*8 dfeat_2b12(n3b2m,3)
        
        integer counter
        logical is_nonzero_feat(n3b1m*n3b1m*n3b2m)
        logical is_nonzero_feat_all(n3b1m*n3b1m*n3b2m)
        
        ! real*8 dfeat_iat0              ! dfeat/dR_self=d(2b01*2b02*2b12)/dR0=(d2b01/dR0)*2b02*2b12+2b01*(d2b02/dR0)*2b12
        real*8 ddfeat_iat1_part1(3,3)    ! (d(d2b01/dR0)/dR1)*2b02*2b12+(d2b01/dR1)*(d2b02/dR0)*2b12
        real*8 ddfeat_iat1_part2(3,3)    ! (d2b01/dR0)*2b02*(d2b12/dR1)+2b01*(d2b02/dR0)*(d2b12/dR1)
        real*8 ddfeat_iat2_part1(3,3)    ! (d2b01/dR0)*(d2b02/dR2)*2b12+2b01*(d(d2b02/dR0)/dR2)*2b12
        
        ! ----------
        
        if(allocated(dfeat_iat0)) then
            deallocate(dfeat_iat0)
            deallocate(ddfeat)
            deallocate(num_nonzero_feat)
            deallocate(ind_nonzero_feat)
        endif
        allocate(dfeat_iat0(3,n3b1m*n3b1m*n3b2m,ntype*(ntype+1)/2,natom))
        allocate(ddfeat(3,m_neigh,3,n3b1m*n3b1m*n3b2m,ntype*(ntype+1)/2,natom))
        allocate(num_nonzero_feat(ntype*(ntype+1)/2,natom))
        allocate(ind_nonzero_feat(n3b1m*n3b1m*n3b2m,ntype*(ntype+1)/2,natom))
        
        ! ----------
        
        pi=4*datan(1.d0)
        pi2=2*pi
        
        dfeat_iat0=0.d0
        ddfeat=0.d0
        do 3000 iat=1,natom
          itype0=itype_atom(iat)
          
          ! ---------- precalculate 2b feature of iat (to be used in the cent-neigh1 and cent-neigh2 part of 3b feature)
          
          do 1000 itype=1,ntype
            do 1000 j=1,num_neigh(itype,iat)
              dd=dR_neigh(1,j,itype,iat)**2+dR_neigh(2,j,itype,iat)**2+dR_neigh(3,j,itype,iat)**2
              d=dsqrt(dd)
              
              num=0
              do k=1,n3b1_type(itype0)
                
                if(d.ge.grid31_2(1,k,itype0).and.d.lt.grid31_2(2,k,itype0)) then
                  num=num+1
                  
                  ind_2b0x(num,j,itype,iat)=k
                  
                  drange=grid31_2(2,k,itype0)-grid31_2(1,k,itype0)
                  x=(d-grid31_2(1,k,itype0))/drange
                  y=(x-0.5d0)*pi2
                  f1=0.5d0*(cos(y)+1)
                  feat_2b0x(num,j,itype)=f1
                  
                  y2=-pi*sin(y)/(d*drange)
                  y2_2b0x(num,j,itype)=y2
                  
                  temp_2b0x(num,j,itype)=pi*(drange*sin(y)-pi2*d*cos(y))/(dd*d*drange**2)
                  
                  dfeat_2b0x(num,j,itype,:)=y2*dR_neigh(:,j,itype,iat)    ! dfeat_2b0x/dRx
                endif
              enddo
              num_2b0x(j,itype,iat)=num
          1000 continue
          
          ! ----------
          
          do 2000 itype2=1,ntype
            do 2000 itype1=1,itype2
              
              itype12=itype1+((itype2-1)*itype2)/2
              
              is_nonzero_feat_all=.false.
              do j2=1,num_neigh(itype2,iat)
                do j1=1,num_neigh(itype1,iat)
                  
                  if(itype1.eq.itype2.and.j1.eq.j2) cycle
                  
                  ! ---------- a case of combination of neigh1 and neigh2
                  jj1=ind_all_neigh(j1,itype1,iat)
                  jj2=ind_all_neigh(j2,itype2,iat)
                  
                  ! ---------- calc 2b feature of neigh1-neigh2
                  dd=(dR_neigh(1,j1,itype1,iat)-dR_neigh(1,j2,itype2,iat))**2+ &
                    (dR_neigh(2,j1,itype1,iat)-dR_neigh(2,j2,itype2,iat))**2+ &
                    (dR_neigh(3,j1,itype1,iat)-dR_neigh(3,j2,itype2,iat))**2
                  d=dsqrt(dd)
                  
                  if(d.gt.Rc2_type(itype0).or.d.lt.1.D-4) cycle
                  
                  num=0
                  do k=1,n3b2_type(itype0)
                    if(d.ge.grid32_2(1,k,itype0).and.d.lt.grid32_2(2,k,itype0)) then
                      num=num+1
                      
                      ind_2b12(num,j1,j2,itype12,iat)=k
                      
                      drange=grid32_2(2,k,itype0)-grid32_2(1,k,itype0)
                      x=(d-grid32_2(1,k,itype0))/drange
                      y=(x-0.5d0)*pi2
                      f1=0.5d0*(cos(y)+1)
                      feat_2b12(num)=f1
                      
                      y2=-pi*sin(y)/(d*drange)
                      
                      dfeat_2b12(num,:)=y2*(dR_neigh(:,j1,itype1,iat)-dR_neigh(:,j2,itype2,iat))    ! dfeat_2b12/dR1
                    endif
                  enddo
                  num_2b12(j2,j1,itype12,iat)=num
                  
                  ! ---------- combine 2b features of cent-neigh1, cent-neigh2, neigh1-neigh2, and get 3b feature of cent-neigh1-neigh2
                  
                  is_nonzero_feat=.false.
                  do i1=1,num_2b0x(j1,itype1,iat)
                    do i2=1,num_2b0x(j2,itype2,iat)
                      do j12=1,num_2b12(j2,j1,itype12,iat)
                        k1=ind_2b0x(i1,j1,itype1,iat)
                        k2=ind_2b0x(i2,j2,itype2,iat)
                        k12=ind_2b12(j12,j1,j2,itype12,iat)
                        
                        ii_f=0
                        if(itype1.ne.itype2) then
                          ii_f=k1+(k2-1)*n3b1_type(itype0)+(k12-1)*n3b1_type(itype0)**2
                        endif
                        if(itype1.eq.itype2.and.k1.le.k2) then
                          ii_f=k1+((k2-1)*k2)/2+(k12-1)*(n3b1_type(itype0)*(n3b1_type(itype0)+1))/2
                        endif
                        
                        if(ii_f.ne.0) then
                          is_nonzero_feat(ii_f)=.true.
                          is_nonzero_feat_all(ii_f)=.true.
                          
                          ! (d2b01/dR0)*2b02*2b12+2b01*(d2b02/dR0)*2b12
                          dfeat_iat0(:,ii_f,itype12,iat)=dfeat_iat0(:,ii_f,itype12,iat)- &
                            dfeat_2b0x(i1,j1,itype1,:)*feat_2b0x(i2,j2,itype2)*feat_2b12(j12)- &
                            feat_2b0x(i1,j1,itype1)*dfeat_2b0x(i2,j2,itype2,:)*feat_2b12(j12)
                          
                          ! (d(d2b01/dR0)/dR1)*2b02*2b12+(d2b01/dR1)*(d2b02/dR0)*2b12
                          ddfeat_iat1_part1(1,1)=-(temp_2b0x(i1,j1,itype1)*dR_neigh(1,j1,itype1,iat)*dR_neigh(1,j1,itype1,iat)+y2_2b0x(i1,j1,itype1))*feat_2b0x(i2,j2,itype2)*feat_2b12(j12)- &
                            dfeat_2b0x(i1,j1,itype1,1)*dfeat_2b0x(i2,j2,itype2,1)*feat_2b12(j12)
                          ddfeat_iat1_part1(2,1)=-(temp_2b0x(i1,j1,itype1)*dR_neigh(2,j1,itype1,iat)*dR_neigh(1,j1,itype1,iat))*feat_2b0x(i2,j2,itype2)*feat_2b12(j12)- &
                            dfeat_2b0x(i1,j1,itype1,2)*dfeat_2b0x(i2,j2,itype2,1)*feat_2b12(j12)
                          ddfeat_iat1_part1(3,1)=-(temp_2b0x(i1,j1,itype1)*dR_neigh(3,j1,itype1,iat)*dR_neigh(1,j1,itype1,iat))*feat_2b0x(i2,j2,itype2)*feat_2b12(j12)- &
                            dfeat_2b0x(i1,j1,itype1,3)*dfeat_2b0x(i2,j2,itype2,1)*feat_2b12(j12)
                          ddfeat_iat1_part1(1,2)=-(temp_2b0x(i1,j1,itype1)*dR_neigh(1,j1,itype1,iat)*dR_neigh(2,j1,itype1,iat))*feat_2b0x(i2,j2,itype2)*feat_2b12(j12)- &
                            dfeat_2b0x(i1,j1,itype1,1)*dfeat_2b0x(i2,j2,itype2,2)*feat_2b12(j12)
                          ddfeat_iat1_part1(2,2)=-(temp_2b0x(i1,j1,itype1)*dR_neigh(2,j1,itype1,iat)*dR_neigh(2,j1,itype1,iat)+y2_2b0x(i1,j1,itype1))*feat_2b0x(i2,j2,itype2)*feat_2b12(j12)- &
                            dfeat_2b0x(i1,j1,itype1,2)*dfeat_2b0x(i2,j2,itype2,2)*feat_2b12(j12)
                          ddfeat_iat1_part1(3,2)=-(temp_2b0x(i1,j1,itype1)*dR_neigh(3,j1,itype1,iat)*dR_neigh(2,j1,itype1,iat))*feat_2b0x(i2,j2,itype2)*feat_2b12(j12)- &
                            dfeat_2b0x(i1,j1,itype1,3)*dfeat_2b0x(i2,j2,itype2,2)*feat_2b12(j12)
                          ddfeat_iat1_part1(1,3)=-(temp_2b0x(i1,j1,itype1)*dR_neigh(1,j1,itype1,iat)*dR_neigh(3,j1,itype1,iat))*feat_2b0x(i2,j2,itype2)*feat_2b12(j12)- &
                            dfeat_2b0x(i1,j1,itype1,1)*dfeat_2b0x(i2,j2,itype2,3)*feat_2b12(j12)
                          ddfeat_iat1_part1(2,3)=-(temp_2b0x(i1,j1,itype1)*dR_neigh(2,j1,itype1,iat)*dR_neigh(3,j1,itype1,iat))*feat_2b0x(i2,j2,itype2)*feat_2b12(j12)- &
                            dfeat_2b0x(i1,j1,itype1,2)*dfeat_2b0x(i2,j2,itype2,3)*feat_2b12(j12)
                          ddfeat_iat1_part1(3,3)=-(temp_2b0x(i1,j1,itype1)*dR_neigh(3,j1,itype1,iat)*dR_neigh(3,j1,itype1,iat)+y2_2b0x(i1,j1,itype1))*feat_2b0x(i2,j2,itype2)*feat_2b12(j12)- &
                            dfeat_2b0x(i1,j1,itype1,3)*dfeat_2b0x(i2,j2,itype2,3)*feat_2b12(j12)
                          
                          ! (d2b01/dR0)*2b02*(d2b12/dR1)+2b01*(d2b02/dR0)*(d2b12/dR1)
                          ddfeat_iat1_part2(1,1)=-dfeat_2b0x(i1,j1,itype1,1)*feat_2b0x(i2,j2,itype2)*dfeat_2b12(j12,1)-feat_2b0x(i1,j1,itype1)*dfeat_2b0x(i2,j2,itype2,1)*dfeat_2b12(j12,1)
                          ddfeat_iat1_part2(2,1)=-dfeat_2b0x(i1,j1,itype1,1)*feat_2b0x(i2,j2,itype2)*dfeat_2b12(j12,2)-feat_2b0x(i1,j1,itype1)*dfeat_2b0x(i2,j2,itype2,1)*dfeat_2b12(j12,2)
                          ddfeat_iat1_part2(3,1)=-dfeat_2b0x(i1,j1,itype1,1)*feat_2b0x(i2,j2,itype2)*dfeat_2b12(j12,3)-feat_2b0x(i1,j1,itype1)*dfeat_2b0x(i2,j2,itype2,1)*dfeat_2b12(j12,3)
                          ddfeat_iat1_part2(1,2)=-dfeat_2b0x(i1,j1,itype1,2)*feat_2b0x(i2,j2,itype2)*dfeat_2b12(j12,1)-feat_2b0x(i1,j1,itype1)*dfeat_2b0x(i2,j2,itype2,2)*dfeat_2b12(j12,1)
                          ddfeat_iat1_part2(2,2)=-dfeat_2b0x(i1,j1,itype1,2)*feat_2b0x(i2,j2,itype2)*dfeat_2b12(j12,2)-feat_2b0x(i1,j1,itype1)*dfeat_2b0x(i2,j2,itype2,2)*dfeat_2b12(j12,2)
                          ddfeat_iat1_part2(3,2)=-dfeat_2b0x(i1,j1,itype1,2)*feat_2b0x(i2,j2,itype2)*dfeat_2b12(j12,3)-feat_2b0x(i1,j1,itype1)*dfeat_2b0x(i2,j2,itype2,2)*dfeat_2b12(j12,3)
                          ddfeat_iat1_part2(1,3)=-dfeat_2b0x(i1,j1,itype1,3)*feat_2b0x(i2,j2,itype2)*dfeat_2b12(j12,1)-feat_2b0x(i1,j1,itype1)*dfeat_2b0x(i2,j2,itype2,3)*dfeat_2b12(j12,1)
                          ddfeat_iat1_part2(2,3)=-dfeat_2b0x(i1,j1,itype1,3)*feat_2b0x(i2,j2,itype2)*dfeat_2b12(j12,2)-feat_2b0x(i1,j1,itype1)*dfeat_2b0x(i2,j2,itype2,3)*dfeat_2b12(j12,2)
                          ddfeat_iat1_part2(3,3)=-dfeat_2b0x(i1,j1,itype1,3)*feat_2b0x(i2,j2,itype2)*dfeat_2b12(j12,3)-feat_2b0x(i1,j1,itype1)*dfeat_2b0x(i2,j2,itype2,3)*dfeat_2b12(j12,3)
                          
                          ! (d2b01/dR0)*(d2b02/dR2)*2b12+2b01*(d(d2b02/dR0)/dR2)*2b12
                          ddfeat_iat2_part1(1,1)=-dfeat_2b0x(i1,j1,itype1,1)*dfeat_2b0x(i2,j2,itype2,1)*feat_2b12(j12)- &
                            feat_2b0x(i1,j1,itype1)*(temp_2b0x(i2,j2,itype2)*dR_neigh(1,j2,itype2,iat)*dR_neigh(1,j2,itype2,iat)+y2_2b0x(i2,j2,itype2))*feat_2b12(j12)
                          ddfeat_iat2_part1(2,1)=-dfeat_2b0x(i1,j1,itype1,1)*dfeat_2b0x(i2,j2,itype2,2)*feat_2b12(j12)- &
                            feat_2b0x(i1,j1,itype1)*(temp_2b0x(i2,j2,itype2)*dR_neigh(2,j2,itype2,iat)*dR_neigh(1,j2,itype2,iat))*feat_2b12(j12)
                          ddfeat_iat2_part1(3,1)=-dfeat_2b0x(i1,j1,itype1,1)*dfeat_2b0x(i2,j2,itype2,3)*feat_2b12(j12)- &
                            feat_2b0x(i1,j1,itype1)*(temp_2b0x(i2,j2,itype2)*dR_neigh(3,j2,itype2,iat)*dR_neigh(1,j2,itype2,iat))*feat_2b12(j12)
                          ddfeat_iat2_part1(1,2)=-dfeat_2b0x(i1,j1,itype1,2)*dfeat_2b0x(i2,j2,itype2,1)*feat_2b12(j12)- &
                            feat_2b0x(i1,j1,itype1)*(temp_2b0x(i2,j2,itype2)*dR_neigh(1,j2,itype2,iat)*dR_neigh(2,j2,itype2,iat))*feat_2b12(j12)
                          ddfeat_iat2_part1(2,2)=-dfeat_2b0x(i1,j1,itype1,2)*dfeat_2b0x(i2,j2,itype2,2)*feat_2b12(j12)- &
                            feat_2b0x(i1,j1,itype1)*(temp_2b0x(i2,j2,itype2)*dR_neigh(2,j2,itype2,iat)*dR_neigh(2,j2,itype2,iat)+y2_2b0x(i2,j2,itype2))*feat_2b12(j12)
                          ddfeat_iat2_part1(3,2)=-dfeat_2b0x(i1,j1,itype1,2)*dfeat_2b0x(i2,j2,itype2,3)*feat_2b12(j12)- &
                            feat_2b0x(i1,j1,itype1)*(temp_2b0x(i2,j2,itype2)*dR_neigh(3,j2,itype2,iat)*dR_neigh(2,j2,itype2,iat))*feat_2b12(j12)
                          ddfeat_iat2_part1(1,3)=-dfeat_2b0x(i1,j1,itype1,3)*dfeat_2b0x(i2,j2,itype2,1)*feat_2b12(j12)- &
                            feat_2b0x(i1,j1,itype1)*(temp_2b0x(i2,j2,itype2)*dR_neigh(1,j2,itype2,iat)*dR_neigh(3,j2,itype2,iat))*feat_2b12(j12)
                          ddfeat_iat2_part1(2,3)=-dfeat_2b0x(i1,j1,itype1,3)*dfeat_2b0x(i2,j2,itype2,2)*feat_2b12(j12)- &
                            feat_2b0x(i1,j1,itype1)*(temp_2b0x(i2,j2,itype2)*dR_neigh(2,j2,itype2,iat)*dR_neigh(3,j2,itype2,iat))*feat_2b12(j12)
                          ddfeat_iat2_part1(3,3)=-dfeat_2b0x(i1,j1,itype1,3)*dfeat_2b0x(i2,j2,itype2,3)*feat_2b12(j12)- &
                            feat_2b0x(i1,j1,itype1)*(temp_2b0x(i2,j2,itype2)*dR_neigh(3,j2,itype2,iat)*dR_neigh(3,j2,itype2,iat)+y2_2b0x(i2,j2,itype2))*feat_2b12(j12)
                          
                          ddfeat(:,jj1,:,ii_f,itype12,iat)=ddfeat(:,jj1,:,ii_f,itype12,iat)+ddfeat_iat1_part1(:,:)+ddfeat_iat1_part2(:,:)
                          ddfeat(:,jj2,:,ii_f,itype12,iat)=ddfeat(:,jj2,:,ii_f,itype12,iat)+ddfeat_iat2_part1(:,:)-ddfeat_iat1_part2(:,:)
                          ddfeat(:,1,:,ii_f,itype12,iat)=ddfeat(:,1,:,ii_f,itype12,iat)-ddfeat_iat1_part1(:,:)-ddfeat_iat2_part1(:,:)
                          
                        endif
                      enddo
                    enddo
                  enddo
                enddo
              enddo
              
              counter=0
              do i=1,n3b1m*n3b1m*n3b2m
                if(is_nonzero_feat_all(i)) then
                  counter=counter+1
                  ind_nonzero_feat(counter,itype12,iat)=i
                endif
              enddo
              num_nonzero_feat(itype12,iat)=counter
              
          2000 continue
        3000 continue
    end subroutine calc_orig_dfeat_and_ddfeat_3b_ppp
    
    
    subroutine calc_new_feat_and_dfeat_3b_ppp(itype_mol,proj_type,iatom_or_ibond,ippp)
        implicit none
        
        integer,intent(in) :: itype_mol
        integer,intent(in) :: proj_type  ! 1, atom; 2, bond
        integer,intent(in) :: iatom_or_ibond
        integer,intent(in) :: ippp
        
        integer i,j,k
        integer iat,itype1,itype2,itype12,ifeat,ii_f,ineigh1,ineigh2
        integer iatom_nonzero_feat,idir
        integer j1,j2
        integer iat1,iat2,iat3,iat4
        integer jj1,jj2,jj3,jj4
        integer num0,num
        integer itype0
        integer k1,k2,k12
        
        real*8 real_bin(3)
        real*8 real_bin2(3)
        real*8 real_bin3(3)
        
        integer map1(n3b1m*n3b1m*n3b2m)
        integer map2(n3b1m*n3b1m*n3b2m)
        
        ! ----------
        
        if(allocated(feat_M2)) then
            deallocate(feat_M2)
            deallocate(dfeat_M2)
        endif
        allocate(feat_M2(nfeat0m,natom))
        allocate(dfeat_M2(nfeat0m,natom,m_neigh,3))
        
        ! ----------
        
        iat1=0
        iat2=0
        iat3=0
        iat4=0
        if(proj_type.eq.1) then
            iat1=ipol(1,iatom_or_ibond,itype_mol)
            iat2=ipol(2,iatom_or_ibond,itype_mol)
            iat3=ipol(3,iatom_or_ibond,itype_mol)
            iat4=ipol(4,iatom_or_ibond,itype_mol)
            
            iatom_nonzero_feat=iatom_or_ibond
            idir=(ippp-1)/2+1
        elseif(proj_type.eq.2) then
            iat1=indb(1,iatom_or_ibond,itype_mol)
            iat2=indb(2,iatom_or_ibond,itype_mol)
            
            if((ippp.eq.1).or.(ippp.eq.2)) iatom_nonzero_feat=iat1
            if((ippp.eq.3).or.(ippp.eq.4)) iatom_nonzero_feat=iat2
            idir=1
        endif
        
        jj1=0
        jj2=0
        jj3=0
        jj4=0
        do i=1,num_neigh_alltype(iatom_nonzero_feat)
            if(list_neigh_alltype(i,iatom_nonzero_feat).eq.iat1) jj1=i
            if(list_neigh_alltype(i,iatom_nonzero_feat).eq.iat2) jj2=i
            if(list_neigh_alltype(i,iatom_nonzero_feat).eq.iat3) jj3=i
            if(list_neigh_alltype(i,iatom_nonzero_feat).eq.iat4) jj4=i
        enddo
        
        if((idir.eq.1).and.((jj1.eq.0).or.(jj2.eq.0))) then
            write(6,*) "(idir.eq.1).and.((jj1.eq.0).or.(jj2.eq.0)),stop"
            stop
        endif
        if((idir.eq.2).and.((jj3.eq.0).or.(jj4.eq.0))) then
            write(6,*) "(idir.eq.2).and.((jj3.eq.0).or.(jj4.eq.0)),stop"
            stop
        endif
        if((idir.eq.3).and.((jj1.eq.0).or.(jj2.eq.0).or.(jj3.eq.0).or.(jj4.eq.0))) then
            write(6,*) "(idir.eq.3).and.((jj1.eq.0).or.(jj2.eq.0).or.(jj3.eq.0).or.(jj4.eq.0)),stop"
            stop
        endif
        
        ! ----------
        
        itype0=itype_atom(iatom_nonzero_feat)
        
        map1=0
        map2=0
        
        num=0
        do k1=1,n3b1_type(itype0)
            do k2=1,n3b1_type(itype0)
                do k12=1,n3b2_type(itype0)
                    ii_f=k1+(k2-1)*n3b1_type(itype0)+(k12-1)*n3b1_type(itype0)**2
                    
                    num=num+1
                    
                    map1(ii_f)=num
                enddo
            enddo
        enddo
        
        num=0
        do k1=1,n3b1_type(itype0)
            do k2=1,n3b1_type(itype0)
                do k12=1,n3b2_type(itype0)
                    ii_f=0
                    if(k1.le.k2) ii_f=k1+((k2-1)*k2)/2+(k12-1)*(n3b1_type(itype0)*(n3b1_type(itype0)+1))/2
                    
                    if(ii_f.gt.0) then
                        num=num+1
                        map2(ii_f)=num
                    endif
                enddo
            enddo
        enddo
        
        ! ----------
                
        feat_M2(:,iatom_nonzero_feat)=0.d0
        dfeat_M2(:,iatom_nonzero_feat,:,:)=0.d0
        
        num0=0
        do itype2=1,ntype
          do itype1=1,itype2
            
            itype12=itype1+((itype2-1)*itype2)/2
            
            do ifeat=1,num_nonzero_feat(itype12,iatom_nonzero_feat)
              
              ii_f=ind_nonzero_feat(ifeat,itype12,iatom_nonzero_feat)
              if(itype1.ne.itype2) num=map1(ii_f)
              if(itype1.eq.itype2) num=map2(ii_f)
              
              ! feat part
              feat_M2(num0+num,iatom_nonzero_feat)=dfeat_iat0(1,ii_f,itype12,iatom_nonzero_feat)*dxpp(1,idir,iatom_or_ibond,proj_type)+ &
                dfeat_iat0(2,ii_f,itype12,iatom_nonzero_feat)*dxpp(2,idir,iatom_or_ibond,proj_type)+ &
                dfeat_iat0(3,ii_f,itype12,iatom_nonzero_feat)*dxpp(3,idir,iatom_or_ibond,proj_type)
            
              ! dfeat part
              do i=1,m_neigh
                dfeat_M2(num0+num,iatom_nonzero_feat,i,:)=dfeat_M2(num0+num,iatom_nonzero_feat,i,:)+ &
                  ddfeat(:,i,1,ii_f,itype12,iatom_nonzero_feat)*dxpp(1,idir,iatom_or_ibond,proj_type)+ &
                  ddfeat(:,i,2,ii_f,itype12,iatom_nonzero_feat)*dxpp(2,idir,iatom_or_ibond,proj_type)+ &
                  ddfeat(:,i,3,ii_f,itype12,iatom_nonzero_feat)*dxpp(3,idir,iatom_or_ibond,proj_type)
              enddo
              
              if(idir.eq.1) then
                  real_bin(:)=dfeat_iat0(1,ii_f,itype12,iatom_nonzero_feat)*dxpp1_diat2(1,:,iatom_or_ibond,proj_type)/A_AU_1+ &
                      dfeat_iat0(2,ii_f,itype12,iatom_nonzero_feat)*dxpp1_diat2(2,:,iatom_or_ibond,proj_type)/A_AU_1+ &
                      dfeat_iat0(3,ii_f,itype12,iatom_nonzero_feat)*dxpp1_diat2(3,:,iatom_or_ibond,proj_type)/A_AU_1
                  
                  dfeat_M2(num0+num,iatom_nonzero_feat,jj1,:)=dfeat_M2(num0+num,iatom_nonzero_feat,jj1,:)-real_bin(:)
                  dfeat_M2(num0+num,iatom_nonzero_feat,jj2,:)=dfeat_M2(num0+num,iatom_nonzero_feat,jj2,:)+real_bin(:)
              endif
              if(idir.eq.2) then
                  real_bin(:)=dfeat_iat0(1,ii_f,itype12,iatom_nonzero_feat)*dxpp2_diat4(1,:,iatom_or_ibond,proj_type)/A_AU_1+ &
                      dfeat_iat0(2,ii_f,itype12,iatom_nonzero_feat)*dxpp2_diat4(2,:,iatom_or_ibond,proj_type)/A_AU_1+ &
                      dfeat_iat0(3,ii_f,itype12,iatom_nonzero_feat)*dxpp2_diat4(3,:,iatom_or_ibond,proj_type)/A_AU_1
                  
                  dfeat_M2(num0+num,iatom_nonzero_feat,jj3,:)=dfeat_M2(num0+num,iatom_nonzero_feat,jj3,:)-real_bin(:)
                  dfeat_M2(num0+num,iatom_nonzero_feat,jj4,:)=dfeat_M2(num0+num,iatom_nonzero_feat,jj4,:)+real_bin(:)
              endif
              if(idir.eq.3) then
                  ! please pay special attention to the meaning of dxpp1_diat2, dxpp2_diat4, dxpp3_diat2, dxpp3_diat4
                  ! for example
                  ! dxpp3_diat2(1,2) is dxpp(1,3)/dx22
                  ! dxpp3_diat4(2,3) is dxpp(2,3)/dx43
                  real_bin(:)=dfeat_iat0(1,ii_f,itype12,iatom_nonzero_feat)*dxpp3_diat2(1,:,iatom_or_ibond,proj_type)/A_AU_1+ &
                      dfeat_iat0(2,ii_f,itype12,iatom_nonzero_feat)*dxpp3_diat2(2,:,iatom_or_ibond,proj_type)/A_AU_1+ &
                      dfeat_iat0(3,ii_f,itype12,iatom_nonzero_feat)*dxpp3_diat2(3,:,iatom_or_ibond,proj_type)/A_AU_1
                  
                  dfeat_M2(num0+num,iatom_nonzero_feat,jj1,:)=dfeat_M2(num0+num,iatom_nonzero_feat,jj1,:)-real_bin(:)
                  dfeat_M2(num0+num,iatom_nonzero_feat,jj2,:)=dfeat_M2(num0+num,iatom_nonzero_feat,jj2,:)+real_bin(:)
                  
                  real_bin(:)=dfeat_iat0(1,ii_f,itype12,iatom_nonzero_feat)*dxpp3_diat4(1,:,iatom_or_ibond,proj_type)/A_AU_1+ &
                      dfeat_iat0(2,ii_f,itype12,iatom_nonzero_feat)*dxpp3_diat4(2,:,iatom_or_ibond,proj_type)/A_AU_1+ &
                      dfeat_iat0(3,ii_f,itype12,iatom_nonzero_feat)*dxpp3_diat4(3,:,iatom_or_ibond,proj_type)/A_AU_1
                  
                  dfeat_M2(num0+num,iatom_nonzero_feat,jj3,:)=dfeat_M2(num0+num,iatom_nonzero_feat,jj3,:)-real_bin(:)
                  dfeat_M2(num0+num,iatom_nonzero_feat,jj4,:)=dfeat_M2(num0+num,iatom_nonzero_feat,jj4,:)+real_bin(:)
              endif
            
            enddo

            if(itype1.eq.itype2) then
              num0=num0+(1+n3b1_type(itype0))*n3b1_type(itype0)/2*n3b2_type(itype0)
            else
              num0=num0+n3b1_type(itype0)*n3b1_type(itype0)*n3b2_type(itype0)
            endif
            
          enddo
        enddo
        
        nfeat0M2=nfeat0m
    end subroutine calc_new_feat_and_dfeat_3b_ppp
    
end module calc_ftype2_ppp
