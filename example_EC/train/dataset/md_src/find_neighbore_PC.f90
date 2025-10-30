subroutine find_neighbore_PC(iatom,natom,xatom,AL,Rc_type,num_neigh,list_neigh, &
    dR_neigh,iat_neigh,ntype,iat_type,m_neigh,Rc_neigh,map2neigh_M,list_neigh_M, &
    num_neigh_M)
    ! num_neigh_M,iat_neigh_M)
    
    implicit none
    integer natom,ntype
    real*8 Rc_type(100),Rc_neigh
    real*8 xatom(3,natom),AL(3,3)
    integer iatom(natom)
    real*8 dR_neigh(3,m_neigh,ntype,natom)
    integer iat_neigh(m_neigh,ntype,natom),list_neigh(m_neigh,ntype,natom)
    integer iat_neigh_M(m_neigh,ntype,natom)
    integer list_neigh_M(m_neigh,ntype,natom),map2neigh_M(m_neigh,ntype,natom)
    integer num_neigh(ntype,natom),num_neigh_M(ntype,natom)
    integer num_type(ntype),num_type_M(ntype)
    integer nperiod(3)
    integer iflag,i,j,k,num
    integer i1,i2,i3,itype
    integer iat_type(100)
    integer itype_atom(natom)
    real*8 d,Rc2,dx1,dx2,dx3,dx,dy,dz,dd
    integer m_neigh
    real*8 Rc_neigh2

    Rc_neigh2 = Rc_neigh*Rc_neigh

    map2neigh_M = 0
    list_neigh_M = 0
    num_neigh_M = 0
    list_neigh = 0
    dR_neigh = 0

        
    itype_atom=0
    
    do i=1,natom
        do j=1,ntype
            if(iatom(i).eq.iat_type(j)) then
                itype_atom(i)=j
            endif
        enddo
        if(itype_atom(i).eq.0) then
            write(6,*) "this atom type didn't found", iatom(i),iat_type(1:ntype)
            stop
        endif
    enddo


    do i = 1, natom
        
        num_type=0
        num_type_M=0
        
        do  j=1,natom

            if (i .eq. j) cycle

            dx1=xatom(1,j)-xatom(1,i)
            ! if(abs(dx1+1).lt.abs(dx1)) dx1=dx1+1
            ! if(abs(dx1-1).lt.abs(dx1)) dx1=dx1-1
            
            dx2=xatom(2,j)-xatom(2,i)
            ! if(abs(dx2+1).lt.abs(dx2)) dx2=dx2+1
            ! if(abs(dx2-1).lt.abs(dx2)) dx2=dx2-1
            
            dx3=xatom(3,j)-xatom(3,i)
            ! if(abs(dx3+1).lt.abs(dx3)) dx3=dx3+1
            ! if(abs(dx3-1).lt.abs(dx3)) dx3=dx3-1

            dx=AL(1,1)*dx1+AL(1,2)*dx2+AL(1,3)*dx3
            dy=AL(2,1)*dx1+AL(2,2)*dx2+AL(2,3)*dx3
            dz=AL(3,1)*dx1+AL(3,2)*dx2+AL(3,3)*dx3

            dd=dx**2+dy**2+dz**2
            
            ! if((dd < Rc_neigh2) .and. (dd.gt.1.D-8)) then
                itype=itype_atom(j)
                num_type_M(itype)=num_type_M(itype)+1

                if(num_type_M(itype).gt.m_neigh) then
                    write(6,*) "Error! maxNeighborNum too small", m_neigh
                    stop
                endif
                
                list_neigh_M(num_type_M(itype),itype,i)=j
                iat_neigh_M(num_type_M(itype),itype,i)=iatom(j)
                
                num_type(itype)=num_type(itype)+1
                
                list_neigh(num_type(itype),itype,i)=j
                iat_neigh(num_type(itype),itype,i)=iatom(j)
                map2neigh_M(num_type(itype),itype,i)=num_type_M(itype)

                dR_neigh(1,num_type(itype),itype,i)=dx
                dR_neigh(2,num_type(itype),itype,i)=dy
                dR_neigh(3,num_type(itype),itype,i)=dz
            ! endif
        
        enddo
        num_neigh(:,i)=num_type(:)
        num_neigh_M(:,i)=num_type_M(:)
        ! write(6,*) "@@@ num neigh ", i, SUM(num_neigh(:,i)), SUM(num_neigh_M(:,i))

    enddo 
    

    ! write(6,*) "check md list neigh"
    ! do i = 1, natom
    !     do itype = 1, ntype
    !         do j = 1, num_neigh(itype,i)
    !             write(6,*) i, itype, j, list_neigh(j,itype,i) 
    !         enddo
    !     enddo
    ! enddo




    return
end subroutine find_neighbore_PC
