       program choose_two
       implicit double precision (a-h,o-z)
       real*8 xp1(3,1000),strength1(1000),wcore1(1000)
       real*8 xp2(3,1000),strength2(1000),wcore2(1000)
       character*50 filename1,filename2

       write(6,*) "input filename1,filename2"
       read(5,*) filename1,filename2
       write(6,*) "iran"
       read(5,*) iran, nump1, nump2
       write(6,*) iran, nump1, nump2

       open(10,file=filename1)
       rewind(10)
       ! nump1=364
       do ii=1,nump1
       read(10,*) jjp
       read(10,*) strength1(ii),wcore1(ii),xp1(1,ii),xp1(2,ii),xp1(3,ii)
       enddo
       close(10)

       open(10,file=filename2)
       rewind(10)
       ! nump2=364
       do ii=1,nump2
       read(10,*) jjp
       read(10,*) strength2(ii),wcore2(ii),xp2(1,ii),xp2(2,ii),xp2(3,ii)
       enddo
       close(10)

!       iran=-2989
       do ii=1,10
       x=ran1(iran)
       enddo


       open(11,file="point.all.2pt")
       rewind(11)

       numpp=0
       do ii=1,1000

       x1=ran1(iran)
       do jj=1,10
       x2=ran1(iran)
       enddo

       ii1=x1*(nump1-1)+1
       ii2=x2*(nump2-1)+1

       dd=(xp1(1,ii1)-xp2(1,ii2))**2+(xp1(2,ii1)-xp2(2,ii2))**2+(xp1(3,ii1)-xp2(3,ii2))**2
       dd=dsqrt(dd)*10.d0
       if(dd.lt.1.5) goto 2000 

       numpp=numpp+1
!       write(11,"(i4,2x,2(3(f12.8,1x),2x))") ii, xp(1,ii1),xp(2,ii1),xp(3,ii1),&
!            xp(1,ii2),xp(2,ii2),xp(3,ii2)
       write(11,*) 2,numpp
       write(11,"(2(E13.5,1x),2x,3(f12.8,1x))") strength1(ii1),wcore1(ii1),xp1(1,ii1),xp1(2,ii1),xp1(3,ii1)
       write(11,"(2(E13.5,1x),2x,3(f12.8,1x))") strength2(ii2),wcore2(ii2),xp2(1,ii2),xp2(2,ii2),xp2(3,ii2)
2000   continue
       enddo
       close(11)

       stop
       end






