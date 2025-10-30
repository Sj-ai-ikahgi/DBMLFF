subroutine find_neighbore(iatom,natom,xatom,AL,Rc_type,num_neigh,list_neigh, &
    dR_neigh,iat_neigh,ntype,iat_type,m_neigh,Rc_M,map2neigh_M,list_neigh_M, &
    num_neigh_M,iat_neigh_M)
     
    implicit none
    integer natom,ntype
    real*8 Rc_type(100),Rc_M
    real*8 xatom(3,natom),AL(3,3)
    integer iatom(natom)
    real*8 dR_neigh(3,m_neigh,ntype,natom)
    integer iat_neigh(m_neigh,ntype,natom),list_neigh(m_neigh,ntype,natom)
    integer iat_neigh_M(m_neigh,ntype,natom)
    integer list_neigh_M(m_neigh,ntype,natom),map2neigh_M(m_neigh,ntype,natom)
    integer num_neigh(ntype,natom),num_neigh_M(ntype,natom)
    integer num_type(ntype),num_type_M(ntype)
    integer nperiod(3)
    integer iflag,i,j,k,l,num
    integer i1,i2,i3,itype
    integer iat_type(100)
    integer itype_atom(natom)
    real*8 d,Rc2,dx1,dx2,dx3,dx,dy,dz,dd
    integer m_neigh
    integer num_test
    integer natom_read,ntype_read
    integer i_read,j_read

    iflag=0
    do i=1,3
        d=dsqrt(AL(1,i)**2+AL(2,i)**2+AL(3,i)**2)
        nperiod(i)=int(Rc_M/d)+1
        if(d.lt.2*Rc_M) iflag=1
    enddo

    ! ----------------------------------------------------------------------
    ! Modified_by_fx_20230705
    ! ----------------------------------------------------------------------
    ! Here, we add two options for iflag. 
    ! The first one is iflag.eq.2, which indicates num_neigh_M, 
    ! iat_neigh_M, list_neigh_M, num_neigh, iat_neigh, list_neigh, 
    ! map2neigh_M will be stored in a file named OUT.NEIGH, and then the 
    ! program exit.
    ! The second one is iflag.eq.3, which indicates all the above variables 
    ! will be read from OUT.NEIGH (dR_neigh will be calculated according to 
    ! the above variables) and return to the code calling this subroutine.
    ! TODO: We can even add the third one iflag.eq.4, which indicates only 
    ! inter-molecular neighbores will be returned. This option can be used 
    ! for pair correction model fitting. (For Li+ - molecule correction, 
    ! one can simply set Rc for types other than Li to 0. But for molecule
    ! - molecule correction, this function is necessary.)
    ! ----------------------------------------------------------------------
    ! !ccccccccccccccccccccccccccccccccccc
    ! !  for molecular system, always set iflag=0
    ! iflag=0
    ! ----------------------------------------------------------------------
    open(400,file="input/find_neighbore.in")
    rewind(400)
    read(400,*) iflag
    close(400)
    ! ----------------------------------------------------------------------

    do i=1,natom
        xatom(1,i)=mod(xatom(1,i)+10.d0,1.d0)
        xatom(2,i)=mod(xatom(2,i)+10.d0,1.d0)
        xatom(3,i)=mod(xatom(3,i)+10.d0,1.d0)
    enddo

    itype_atom=0
    do i=1,natom
        do j=1,ntype
            if(iatom(i).eq.iat_type(j)) then
                itype_atom(i)=j
            endif
        enddo
        if(itype_atom(i).eq.0) then
            write(6,*) "this atom type didn't found(neigh)", &
                iatom(i),(iat_type(j),j=1,ntype) 
            stop
        endif
    enddo


    Rc2=Rc_M**2

    do 2000 i=1,natom

        if(iflag.eq.1) then
            num_type=0
            num_type_M=0
            num_test=0
            do  j=1,natom
                do i1=-nperiod(1),nperiod(1)
                    do i2=-nperiod(2),nperiod(2)
                        do i3=-nperiod(3),nperiod(3) 
                            dx1=xatom(1,j)-xatom(1,i)+i1
                            dx2=xatom(2,j)-xatom(2,i)+i2
                            dx3=xatom(3,j)-xatom(3,i)+i3
                            dx=AL(1,1)*dx1+AL(1,2)*dx2+AL(1,3)*dx3
                            dy=AL(2,1)*dx1+AL(2,2)*dx2+AL(2,3)*dx3
                            dz=AL(3,1)*dx1+AL(3,2)*dx2+AL(3,3)*dx3
                            dd=dx**2+dy**2+dz**2
                            if(dd.lt.Rc2.and.dd.gt.1.D-8) then
                                num_test=num_test+1
                                itype=itype_atom(j)
                                num_type_M(itype)=num_type_M(itype)+1
                                if(num_type_M(itype).gt.m_neigh) then
                                    write(6,*) "num.gt.m_neigh, stop", m_neigh
                                    stop
                                endif
                                list_neigh_M(num_type_M(itype),itype,i)=j
                                iat_neigh_M(num_type_M(itype),itype,i)=iatom(j)

                                if(dd.lt.Rc_type(itype_atom(i))**2.and.dd.gt.1.D-8) then
                                    num_type(itype)=num_type(itype)+1
                                    list_neigh(num_type(itype),itype,i)=j
                                    iat_neigh(num_type(itype),itype,i)=iatom(j)
                                    map2neigh_M(num_type(itype),itype,i)=num_type_M(itype)
                                    dR_neigh(1,num_type(itype),itype,i)=dx
                                    dR_neigh(2,num_type(itype),itype,i)=dy
                                    dR_neigh(3,num_type(itype),itype,i)=dz
                                endif
                            endif
                        enddo
                    enddo
                enddo
            enddo
            num_neigh(:,i)=num_type(:)
            num_neigh_M(:,i)=num_type_M(:)
        endif

        ! ----------------------------------------------------------------------
        ! Modified_by_fx_20230705
        ! ----------------------------------------------------------------------
        ! if(iflag.eq.0) then
        ! ----------------------------------------------------------------------
        if((iflag.eq.0).or.(iflag.eq.2)) then
        ! ----------------------------------------------------------------------
            num_type=0
            num_type_M=0
            do  j=1,natom
                dx1=xatom(1,j)-xatom(1,i)
                if(abs(dx1+1).lt.abs(dx1)) dx1=dx1+1
                if(abs(dx1-1).lt.abs(dx1)) dx1=dx1-1
                dx2=xatom(2,j)-xatom(2,i)
                if(abs(dx2+1).lt.abs(dx2)) dx2=dx2+1
                if(abs(dx2-1).lt.abs(dx2)) dx2=dx2-1
                dx3=xatom(3,j)-xatom(3,i)
                if(abs(dx3+1).lt.abs(dx3)) dx3=dx3+1
                if(abs(dx3-1).lt.abs(dx3)) dx3=dx3-1

                dx=AL(1,1)*dx1+AL(1,2)*dx2+AL(1,3)*dx3
                dy=AL(2,1)*dx1+AL(2,2)*dx2+AL(2,3)*dx3
                dz=AL(3,1)*dx1+AL(3,2)*dx2+AL(3,3)*dx3
                dd=dx**2+dy**2+dz**2
                if(dd.lt.Rc2.and.dd.gt.1.D-8) then
                    itype=itype_atom(j)
                    num_type_M(itype)=num_type_M(itype)+1

                    if(num_type_M(itype).gt.m_neigh) then
                        write(6,*) "num.gt.m_neigh, stop",m_neigh
                        stop
                    endif
                    list_neigh_M(num_type_M(itype),itype,i)=j
                    iat_neigh_M(num_type_M(itype),itype,i)=iatom(j)

                    if(dd.lt.Rc_type(itype_atom(i))**2.and.dd.gt.1.D-8) then
                        num_type(itype)=num_type(itype)+1

                        list_neigh(num_type(itype),itype,i)=j
                        iat_neigh(num_type(itype),itype,i)=iatom(j)
                        map2neigh_M(num_type(itype),itype,i)=num_type_M(itype)
                        dR_neigh(1,num_type(itype),itype,i)=dx
                        dR_neigh(2,num_type(itype),itype,i)=dy
                        dR_neigh(3,num_type(itype),itype,i)=dz
                    endif
                endif
            enddo
            num_neigh(:,i)=num_type(:)
            num_neigh_M(:,i)=num_type_M(:)
        endif

2000 continue
    
    ! ----------------------------------------------------------------------
    ! Modified_by_fx_20230705
    ! ----------------------------------------------------------------------
    if(iflag.eq.2) then
        open(400,file="input/IN.NEIGHBORE")
        rewind(400)
        write(400,*) natom,ntype
        do i=1,natom
            write(400,"(i5,100(i5))") i,(num_neigh(j,i),j=1,ntype)
        enddo
        do i=1,natom
            do j=1,ntype
                do k=1,num_neigh(j,i)
                    ! write(400,"(2(i5),3(i5),3(es25.16))") i,j,iat_neigh(k,j,i),list_neigh(k,j,i),map2neigh_M(k,j,i),(dR_neigh(l,k,j,i),l=1,3)
                    write(400,"(2(i5),3(i5))") i,j,iat_neigh(k,j,i),list_neigh(k,j,i),map2neigh_M(k,j,i)
                enddo
            enddo
        enddo
        do i=1,natom
            write(400,"(i5,100(i5))") i,(num_neigh_M(j,i),j=1,ntype)
        enddo
        do i=1,natom
            do j=1,ntype
                do k=1,num_neigh_M(j,i)
                    write(400,"(2(i5),2(i5))") i,j,iat_neigh_M(k,j,i),list_neigh_M(k,j,i)
                enddo
            enddo
        enddo
        close(400)
        stop
    endif
    
    if(iflag.eq.3) then
        open(400,file="input/IN.NEIGHBORE")
        rewind(400)
        read(400,*) natom_read,ntype_read
        if(natom_read.ne.natom) then
            write(6,*) "natom_read.ne.natom,stop",natom_read,natom
            stop
        endif
        if(ntype_read.ne.ntype) then
            write(6,*) "ntype_read.ne.ntype,stop",ntype_read,ntype
            stop
        endif
        
        num_neigh=0
        iat_neigh=0
        list_neigh=0
        map2neigh_M=0
        dR_neigh=0.d0
        do i=1,natom
            read(400,*) i_read,(num_neigh(j,i),j=1,ntype)
        enddo
        do i=1,natom
            do j=1,ntype
                do k=1,num_neigh(j,i)
                    read(400,*) i_read,j_read,iat_neigh(k,j,i),list_neigh(k,j,i),map2neigh_M(k,j,i)
                    
                    dx1=xatom(1,list_neigh(k,j,i))-xatom(1,i)
                    if(abs(dx1+1).lt.abs(dx1)) dx1=dx1+1
                    if(abs(dx1-1).lt.abs(dx1)) dx1=dx1-1
                    dx2=xatom(2,list_neigh(k,j,i))-xatom(2,i)
                    if(abs(dx2+1).lt.abs(dx2)) dx2=dx2+1
                    if(abs(dx2-1).lt.abs(dx2)) dx2=dx2-1
                    dx3=xatom(3,list_neigh(k,j,i))-xatom(3,i)
                    if(abs(dx3+1).lt.abs(dx3)) dx3=dx3+1
                    if(abs(dx3-1).lt.abs(dx3)) dx3=dx3-1
                    
                    dx=AL(1,1)*dx1+AL(1,2)*dx2+AL(1,3)*dx3
                    dy=AL(2,1)*dx1+AL(2,2)*dx2+AL(2,3)*dx3
                    dz=AL(3,1)*dx1+AL(3,2)*dx2+AL(3,3)*dx3

                    dR_neigh(1,k,j,i)=dx
                    dR_neigh(2,k,j,i)=dy
                    dR_neigh(3,k,j,i)=dz
                    
                    ! open(401,file="input/OUT.NEIGHBORE.CHECK",position="append")
                    ! write(401,"(3(es25.16))") dR_neigh(:,k,j,i)
                    ! close(401)
                enddo
            enddo
        enddo
        
        num_neigh_M=0
        iat_neigh_M=0
        list_neigh_M=0
        do i=1,natom
            read(400,*) i_read,(num_neigh_M(j,i),j=1,ntype)
        enddo
        do i=1,natom
            do j=1,ntype
                do k=1,num_neigh_M(j,i)
                    read(400,*) i_read,j_read,iat_neigh_M(k,j,i),list_neigh_M(k,j,i)
                enddo
            enddo
        enddo
        close(400)
    endif
    ! ----------------------------------------------------------------------
    
    return
end subroutine find_neighbore
      
