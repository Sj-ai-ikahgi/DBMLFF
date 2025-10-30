program funcr_fit_Ecut
    implicit double precision (a-h,o-z)
    parameter (mq=200)
    real*8, allocatable, dimension (:,:) :: funcr2,funcr2_new
    real*8, allocatable, dimension (:) :: Q_type,Q_type_new
    real*8, allocatable, dimension (:) :: rr
    real*8, allocatable, dimension (:) :: tmm

    real*8 q(mq),v(mq)
    real*8 aI(0:100),aM(0:100,mq)
    real*8 aM0(0:100),weight(0:100)

    real*8 AA(mq,mq),BB(mq)
    real*8 AA2(mq,mq),BB2(mq)
    real*8 SS(mq),work(30*mq+60)

    real*8 rc_list(50)


    character*50 fname

    pi=4*datan(1.d0)

    ! rc1 = 6.d0
    ! rc2 = 9.d0

    rc1 = 2.d0/0.529177208
    rc2 = 4.d0/0.529177208

    rc_list(:) = 0.3
    

    rc_list(1:5) = 0.3
    rc_list(6:7) = 0.2
    rc_list(8) = 0.3
    rc_list(9) = 0.5
    rc_list(10) = 0.3
    rc_list(11) = 0.5
    rc_list(12:15) = 0.2


    write(6,*) "input the funcr_atom.fit.bin file name"
    ! read(5,*) fname
    fname = "funcr_atom.fit.bin_origin"

    open(13,file=fname,form="unformatted",action="read")
    rewind(13)
    read(13) nr,Rm2,ntype_cent

    allocate(Q_type(ntype_cent))
    allocate(Q_type_new(ntype_cent))
    allocate(funcr2(nr,ntype_cent))
    allocate(funcr2_new(nr,ntype_cent))
    allocate(rr(nr))
    allocate(tmm(nr))

    read(13) Q_type
    read(13) funcr2
    close(13)

    write(6,*) "input rc(a.u) and Ecut2(Ryd),mc(6),weight0(10),mm(3)"
    ! read(5,*) rc0,Ecut2,mc,weight0,mm
    ! read(5,*) Ecut2,mc,weight0,mm
    ! Ecut2 = 200
    ! mc = 6
    ! weight = 10
    ! mm = 3

    Ecut2 = 150
    mc = 6
    weight0 = 10
    mm = 3

    ! rc0: the cut of point
    ! Ecut2(Ryd): the Ecut2 to be used in the code
    ! mc: the max moment to fit (don't be too high)
    ! weight0: weight for the moment
    ! mm: the r**mm weight for r>rc
    !  This code refit the functionrc using reciprocal q within a Ecut2.
    !  So, it can be used later in the grid calculations.
    ! For r<rc, it fits the moments upto \rho(r)*r**(2+mc)*dr (so mc=0 is the
    ! charge).
    ! For r>rc, it requirs all the values are the same.

    dr=Rm2/nr
    do i=1,nr
        rr(i)=(i-1)*Rm2/nr
    enddo
    qmn=dsqrt(Ecut2)
    do iq=1,mq
        q(iq)=iq*qmn/mq    ! no zero q, contribution is zero
    enddo

    weight=weight0
    weight(0)=weight0*10

    write(6,*) "ntype_cent=",ntype_cent

    do 100 jj=1,ntype_cent
    ! do 100 9=1,9

        rc0 = rc_list(jj)

        write(6,*) "jj=",jj, rc0


        irc=rc0*nr/Rm2+1
        rc=(irc-1)*Rm2/nr

        do iq=1,mq
            x=rc*q(iq)
            aI(1)=-x*dcos(x)+dsin(x)
            aI(2)=(-x**2+2.d0)*dcos(x)+2*x*dsin(x)-2.d0
            do im=3,mc+1
                aI(im)=x**(im-1)*(im*dsin(x)-x*dcos(x))-im*(im-1)*aI(im-2)
            enddo

            do im=0,mc
                aM(im,iq)=aI(im+1)/q(iq)**(im+1)
            enddo
        enddo

        do im=0,mc
            sumr=0.d0
            do ir=2,irc
                sumr=sumr+0.5*(funcr2(ir,jj)+funcr2(ir-1,jj))*((rr(ir-1)+rr(ir))/2)**(2+im)
            enddo
            aM0(im)=sumr*dr
        enddo

        !    (sum_iq Vq(iq)*aM(im,iq)-aM0(im))**2

        AA=0.d0
        BB=0.d0
        do iq1=1,mq
            do iq2=1,mq
                do im=0,mc
                    AA(iq1,iq2)=AA(iq1,iq2)+aM(im,iq1)*aM(im,iq2)*weight(im)
                enddo
            enddo
        enddo

        do iq=1,mq
            do im=0,mc
                BB(iq)=BB(iq)+aM(im,iq)*aM0(im)*weight(im)
            enddo
        enddo


        AA2=0.d0

        do iq1=1,mq
            q1=q(iq1)
            do iq2=1,iq1
                q2=q(iq2)
                sumr=0.d0
                do ir=irc,nr
                    sumr=sumr+dsin(q1*rr(ir))*dsin(q2*rr(ir))*rr(ir)**(mm-2)
                enddo
                AA2(iq1,iq2)=sumr*dr*q1*q2
                AA2(iq2,iq1)=sumr*dr*q1*q2
            enddo
        enddo

        BB2=0.d0
        do iq=1,mq
            q1=q(iq)
            sumr=0.d0
            do ir=irc,nr
                sumr=sumr+dsin(q1*rr(ir))*funcr2(ir,jj)*rr(ir)**(mm-1)
            enddo
            BB2(iq)=sumr*dr*q1
        enddo

        AA=AA+AA2
        BB=BB+BB2
        !cccc solve line

        lwork=30*mq+60
        call dgelss(mq,mq,1,AA,mq,BB,mq,SS,-0.1,irank,work,lwork,info)


        do ir=2,nr
            sumr=0.d0
            r=rr(ir)
            do iq=1,mq
                sumr=sumr+BB(iq)*q(iq)*dsin(q(iq)*r)
            enddo
            funcr2_new(ir,jj)=sumr/r
        enddo
        funcr2_new(1,jj)=funcr2_new(2,jj)


        ! change new to old from rc1 --- begin
        do ir=2,nr
            r = rr(ir)
            if (r < rc1) then
                tmm(ir) = 1.0
            elseif (r < rc2) then
                tmm(ir) = 0.5*(cos((r-rc1)/(rc2-rc1)*PI) + 1.0d0)
            else
                tmm(ir) = 0.0
            endif
        enddo

        do ir=2, nr
            funcr2_new(ir,jj) = funcr2_new(ir,jj)*tmm(ir) + funcr2(ir,jj)*(1.0-tmm(ir))
        enddo
        ! change new to old from rc1 --- end


        do im=0,mc
            sumr=0.d0
            do ir=2,irc
                sumr=sumr+0.5*(funcr2_new(ir,jj)+funcr2_new(ir-1,jj))*((rr(ir-1)+rr(ir))/2)**(2+im)
            enddo
            aM00=sumr*dr
            write(6,*) "aM0",im,aM0(im),aM00
        enddo

        sumr=0.d0
        do ir=2,nr
            sumr=sumr+0.5*(funcr2_new(ir,jj)+funcr2_new(ir-1,jj))*((rr(ir-1)+rr(ir))/2)**2
        enddo

        Q_type_new(jj)=sumr*4*pi*dr

        write(6,*) "Q_type,Q_type_new",Q_type(jj),Q_type_new(jj)
        write(6,*) "rescale", Q_type(jj)/Q_type_new(jj)

         funcr2_new(:,jj)=funcr2_new(:,jj)*Q_type(jj)/Q_type_new(jj)

100 continue


    open(13,file="funcr_atom.fit.bin_new",form="unformatted")
    rewind(13)
    write(13) nr,Rm2,ntype_cent
    write(13) Q_type
    write(13) funcr2_new
    close(13)

    open(11,file="plot.funcr")
    rewind(11)

    do i=1,nr
        r=(i-1)*Rm2/nr
        write(11,"(50(E15.7,1x))") r, (funcr2(i,j),j=1,ntype_cent)
    enddo

    close(11)

    open(11,file="plot.funcr_new")
    rewind(11)

    do i=1,nr
        r=(i-1)*Rm2/nr
        write(11,"(50(E15.7,1x))") r, (funcr2_new(i,j),j=1,ntype_cent)
    enddo

    close(11)

    stop
end


