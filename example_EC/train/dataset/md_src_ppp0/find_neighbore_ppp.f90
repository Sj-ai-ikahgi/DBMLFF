module find_neighbore_ppp

    implicit none

    integer natom
    integer ntype
    integer m_neigh
        
    integer,allocatable,dimension(:,:) :: num_neigh
    integer,allocatable,dimension(:) :: num_neigh_alltype
    integer,allocatable,dimension(:,:,:) :: iat_neigh
    integer,allocatable,dimension(:,:) :: iat_neigh_alltype
    integer,allocatable,dimension(:,:,:) :: list_neigh
    integer,allocatable,dimension(:,:) :: list_neigh_alltype
    real*8,allocatable,dimension(:,:,:,:) :: dR_neigh
    real*8,allocatable,dimension(:,:,:) :: dR_neigh_alltype
    integer,allocatable,dimension(:,:,:) :: ind_all_neigh
    
    integer,allocatable,dimension(:,:) :: num_neigh_M
    integer,allocatable,dimension(:) :: num_neigh_M_alltype
    integer,allocatable,dimension(:,:,:) :: iat_neigh_M
    integer,allocatable,dimension(:,:) :: iat_neigh_M_alltype
    integer,allocatable,dimension(:,:,:) :: list_neigh_M
    integer,allocatable,dimension(:,:) :: list_neigh_M_alltype
    integer,allocatable,dimension(:,:,:) :: map2neigh_M
    integer,allocatable,dimension(:,:) :: map2neigh_M_alltype
    integer,allocatable,dimension(:,:,:) :: ind_all_neigh_M

contains
    subroutine load_neighbore(path_IN_NEIGHBORE,natom_in,ntype_in,m_neigh_in)
    ! ------------------------------------------------------------
    ! added_by_fx
    ! ------------------------------------------------------------
    ! Users can specify the neighbores of each atom in the file
    ! IN.NEIGHBORE. This subroutine is designed to read the file
    ! IN.NEIGHBORE and save data in mem.
    ! ------------------------------------------------------------
    ! I     path_IN_NEIGHBORE
    ! I     natom_in
    ! I     ntype_in
    ! I     m_neigh_in
    !
    ! O,M   natom
    ! O,M   ntype
    ! O,M   m_neigh
    ! O,M   num_neigh
    ! O,M   num_neigh_alltype
    ! O,M   iat_neigh
    ! O,M   iat_neigh_alltype
    ! O,M   list_neigh
    ! O,M   list_neigh_alltype
    ! O,M   ind_all_neigh
    ! ------------------------------------------------------------
        implicit none
        
        character(len=1000),intent(in) :: path_IN_NEIGHBORE
        integer,intent(in) :: natom_in
        integer,intent(in) :: ntype_in
        integer,intent(in) :: m_neigh_in
        
        integer int_bin
        integer int_bin2
        
        integer i,j,k
        integer iat,itype
        integer num
        
        ! ----------
        
        if(allocated(num_neigh)) then
            deallocate(num_neigh)
            deallocate(num_neigh_alltype)
            deallocate(iat_neigh)
            deallocate(iat_neigh_alltype)
            deallocate(list_neigh)
            deallocate(list_neigh_alltype)
            deallocate(ind_all_neigh)
        endif
        allocate(num_neigh(ntype,natom))
        allocate(num_neigh_alltype(natom))
        allocate(iat_neigh(m_neigh,ntype,natom))
        allocate(iat_neigh_alltype(m_neigh,natom))
        allocate(list_neigh(m_neigh,ntype,natom))
        allocate(list_neigh_alltype(m_neigh,natom))
        allocate(ind_all_neigh(m_neigh,ntype,natom))
        
        ! ----------
        
        natom=natom_in
        ntype=ntype_in
        m_neigh=m_neigh_in
        
        ! ----------
        
        open(400,file=trim(path_IN_NEIGHBORE))
        rewind(400)
        
        read(400,*) int_bin,int_bin2
        if(int_bin.ne.natom) then
            write(6,*) "int_bin.ne.natom,stop",int_bin,natom
            stop
        endif
        if(int_bin2.ne.ntype) then
            write(6,*) "int_bin2.ne.ntype,stop",int_bin2,ntype
            stop
        endif
        
        num_neigh=0
        iat_neigh=0
        list_neigh=0
        do i=1,natom
            read(400,*) int_bin,(num_neigh(j,i),j=1,ntype)
        enddo
        do i=1,natom
            do j=1,ntype
                do k=1,num_neigh(j,i)
                    read(400,*) int_bin,int_bin2,iat_neigh(k,j,i),list_neigh(k,j,i)
                enddo
            enddo
        enddo
        
        close(400)
        
        ! ----------
        
        num_neigh_alltype=0
        iat_neigh_alltype=0
        list_neigh_alltype=0
        ind_all_neigh=0
        
        do iat=1,natom
            num=1
            list_neigh_alltype(1,iat)=iat    ! the first neighbore is itself
            
            do itype=1,ntype
                do j=1,num_neigh(itype,iat)
                    num=num+1
                    if(num.gt.m_neigh) then
                        write(6,*) "total num_neigh.gt.m_neigh,stop",m_neigh
                        stop
                    endif
                    ind_all_neigh(j,itype,iat)=num
                    iat_neigh_alltype(num,iat)=iat_neigh(j,itype,iat)    ! iat_neigh_alltype(1,iat) is meaningless, do not know the type of the center
                    list_neigh_alltype(num,iat)=list_neigh(j,itype,iat)
                enddo
            enddo
            num_neigh_alltype(iat)=num
        enddo
    end subroutine load_neighbore
    
    
    subroutine find_neighbore_within_Rc_M(Rc_M,iatom,xatom,AL,iat_type,is_strict)
    ! ------------------------------------------------------------
    ! added_by_fx
    ! ------------------------------------------------------------
    ! After loading the neighbore list, find the neighbores 
    ! within Rc_M and calc the mapping from the loaded neighbores 
    ! to the neighbores within Rc_M.
    ! ------------------------------------------------------------
    ! I     Rc_M
    ! I     iatom: The type of each atom in the current configuration
    ! I     xatom
    ! I     AL
    ! I     iat_type: The absolute order of types
    ! I     is_strict
    !
    ! I,M   natom
    ! I,M   ntype
    ! I,M   m_neigh
    ! I,M   num_neigh
    ! I,M   list_neigh
    ! O,M   num_neigh_M
    ! O,M   num_neigh_M_alltype
    ! O,M   iat_neigh_M
    ! O,M   iat_neigh_M_alltype
    ! O,M   list_neigh_M
    ! O,M   list_neigh_M_alltype
    ! O,M   map2neigh_M
    ! O,M   map2neigh_M_alltype
    ! O,M   ind_all_neigh_M
    ! ------------------------------------------------------------
        implicit none
        
        real*8,intent(in) :: Rc_M
        integer,intent(in) :: iatom(natom)
        real*8,intent(in) :: xatom(3,natom)
        real*8,intent(in) :: AL(3,3)
        integer,intent(in) :: iat_type(natom)
        logical,intent(in) :: is_strict
        
        integer itype_atom(natom) ! Absolute index of the type of each atom in the current configuration 
        integer num_type_M(ntype)
        real*8 dx_frac(3)
        real*8 dx(3)
        real*8 dd
        
        integer i,j,k
        integer ixyz,itype,iat
        integer num
        
        ! ----------
        
        if(allocated(num_neigh_M)) then
            deallocate(num_neigh_M)
            deallocate(num_neigh_M_alltype)
            deallocate(iat_neigh_M)
            deallocate(iat_neigh_M_alltype)
            deallocate(list_neigh_M)
            deallocate(list_neigh_M_alltype)
            deallocate(map2neigh_M)
            deallocate(map2neigh_M_alltype)
            deallocate(ind_all_neigh_M)
        endif
        allocate(num_neigh_M(ntype,natom))
        allocate(num_neigh_M_alltype(natom))
        allocate(iat_neigh_M(m_neigh,ntype,natom))
        allocate(iat_neigh_M_alltype(m_neigh,natom))
        allocate(list_neigh_M(m_neigh,ntype,natom))
        allocate(list_neigh_M_alltype(m_neigh,natom))
        allocate(map2neigh_M(m_neigh,ntype,natom))
        allocate(map2neigh_M_alltype(m_neigh,natom))
        allocate(ind_all_neigh_M(m_neigh,ntype,natom))
        
        ! ----------
        
        itype_atom=0
        do i=1,natom
            do j=1,ntype
                if(iatom(i).eq.iat_type(j)) then
                    itype_atom(i)=j
                endif
            enddo
            if(itype_atom(i).eq.0) then
                write(6,*) "this atom type didn't found,stop", itype_atom(i)
                stop
            endif
        enddo
        
        ! ----------
        
        num_type_M=0
        do i=1,natom
            do j=1,natom
                dx_frac(:)=xatom(:,j)-xatom(:,i)
                do ixyz=1,3
                    do while(abs(dx_frac(ixyz)+1).lt.abs(dx_frac(ixyz)))
                        dx_frac(ixyz)=dx_frac(ixyz)+1
                    enddo
                    do while(abs(dx_frac(ixyz)-1).lt.abs(dx_frac(ixyz)))
                        dx_frac(ixyz)=dx_frac(ixyz)-1
                    enddo
                enddo
                dx(:)=AL(:,1)*dx_frac(1)+AL(:,2)*dx_frac(2)+AL(:,3)*dx_frac(3)
                dd=dx(1)**2+dx(2)**2+dx(3)**2
                
                if((dd.lt.Rc_M**2).and.(dd.gt.1.D-8)) then
                    itype=itype_atom(j)
                    num_type_M(itype)=num_type_M(itype)+1
                    
                    if(num_type_M(itype).gt.m_neigh) then
                        write(6,*) "num.gt.m_neigh, stop",m_neigh
                        stop
                    endif
                    list_neigh_M(num_type_M(itype),itype,i)=j
                    iat_neigh_M(num_type_M(itype),itype,i)=iatom(j)
                    
                    do k=1,num_neigh(itype,i)
                        if(list_neigh(k,itype,i).eq.j) map2neigh_M(k,itype,i)=num_type_M(itype)
                    enddo
                endif
            enddo
            
            do itype=1,ntype
                do j=1,num_neigh(itype,i)
                    if(map2neigh_M(j,itype,i).eq.0) then
                        if(is_strict) then
                            write(6,*) 'map2neigh_M(j,itype,i).eq.0! There is a specified neighbore with a distance greater than Rc_M from the central atom! stop'
                            stop
                        else
                            write(6,*) 'map2neigh_M(j,itype,i).eq.0! There is a specified neighbore with a distance greater than Rc_M from the central atom!'
                        endif
                    endif
                enddo
            enddo
            
            num_neigh_M(:,i)=num_type_M(:)
        enddo
        
        ! ----------
        
        num_neigh_M_alltype=0
        iat_neigh_M_alltype=0
        list_neigh_M_alltype=0
        map2neigh_M_alltype=0
        ind_all_neigh_M=0
        
        do iat=1,natom
            num=1
            list_neigh_M_alltype(1,iat)=iat   ! the first neighbore is itself
            
            do itype=1,ntype
                do j=1,num_neigh_M(itype,iat)
                    num=num+1
                    if(num.gt.m_neigh) then
                        write(6,*) "total num_neigh.gt.m_neigh,stop",m_neigh
                        stop
                    endif
                    ind_all_neigh_M(j,itype,iat)=num
                    iat_neigh_M_alltype(num,iat)=iat_neigh_M(j,itype,iat)    ! iat_neigh_M_alltype(1,iat) is meaningless, do not know the type of the center
                    list_neigh_M_alltype(num,iat)=list_neigh_M(j,itype,iat)
                enddo
            enddo
            
            do i=1,num_neigh_alltype(iat)
                do j=1,num_neigh_M_alltype(iat)
                    if(list_neigh_alltype(i,iat).eq.list_neigh_M_alltype(j,iat)) map2neigh_M_alltype(i,iat)=j
                enddo
            enddo
            
            num_neigh_M_alltype(iat)=num
        enddo       
    end subroutine find_neighbore_within_Rc_M
    

    subroutine calc_dR_neigh(xatom,AL)
    ! ------------------------------------------------------------
    ! added_by_fx
    ! ------------------------------------------------------------
    ! After loading the neighbore list, calc the distance between 
    ! each atom and its neighbores.
    ! ------------------------------------------------------------
    ! I     xatom
    ! I     AL
    !
    ! I,M   natom
    ! I,M   ntype
    ! I,M   m_neigh
    ! I,M   num_neigh
    ! I,M   list_neigh
    ! O,M   dR_neigh
    ! O,M   dR_neigh_alltype
    ! ------------------------------------------------------------
        implicit none
        
        real*8,intent(in) :: xatom(3,natom)
        real*8,intent(in) :: AL(3,3)
        
        real*8 dx_frac(3)
        real*8 dx(3)
        
        integer i,j,k
        integer ixyz
        integer num
        
        if(allocated(dR_neigh)) then
            deallocate(dR_neigh)
            deallocate(dR_neigh_alltype)
        endif
        allocate(dR_neigh(3,m_neigh,ntype,natom))
        allocate(dR_neigh_alltype(3,m_neigh,natom))
        
        do i=1,natom
            num=1
            dR_neigh_alltype(:,1,i)=0.d0    ! the first neighbore is itself
            do j=1,ntype
                do k=1,num_neigh(j,i)
                    num=num+1
                    
                    dx_frac(:)=xatom(:,list_neigh(k,j,i))-xatom(:,i)
                    do ixyz=1,3
                        do while(abs(dx_frac(ixyz)+1).lt.abs(dx_frac(ixyz)))
                            dx_frac(ixyz)=dx_frac(ixyz)+1
                        enddo
                        do while(abs(dx_frac(ixyz)-1).lt.abs(dx_frac(ixyz)))
                            dx_frac(ixyz)=dx_frac(ixyz)-1
                        enddo
                    enddo
                    
                    dx(:)=AL(:,1)*dx_frac(1)+AL(:,2)*dx_frac(2)+AL(:,3)*dx_frac(3)
        
                    dR_neigh(:,k,j,i)=dx(:)
                    dR_neigh_alltype(:,num,i)=dx(:)
                enddo
            enddo
        enddo
    end subroutine calc_dR_neigh
    
end module find_neighbore_ppp
