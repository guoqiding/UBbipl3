C

      SUBROUTINE abagplot(n,alpha1,x,y,whisk,tukm,interpx, interpy,num,
     +    datatyp,indoutl,datatyp2,pxpy,boxpl,nointer)
C
C  This subroutine computes a alphabagplot for any bivariate data set.
C
C  There is a choice between computing all the whiskers (take whisk equal
C  to 1), or few whiskers (only one whisker on each edge of the bag - take
C  whisk equal to 2), or no whiskers (take whisk equal to 3), or star-shaped 
C  whiskers (take whisk equal to 4). 
C  If the data set is not in general position, dithering is used.
C  Some adaptations are made for small n, and for large n (take subsets).
C  If all the points are (nearly) collinear, the bagplot reduces to the
C  usual boxplot for univariate data. 
C  In this case, boxplot is set equal to 1 or 2.
C
C  Version: 15 December 1999
C
  
      implicit integer(i-n), double precision(a-h,o-z)
      PARAMETER(MAXN=2500, MAXM=MAXN*(MAXN-1)/2, MAXNUM=500000)
      INTEGER NCIRQ(MAXN),MCIRQ(MAXN),NRANK(MAXN),F(MAXN)
      integer JLV(MAXM),JRV(MAXM),kstar
      INTEGER IND1(MAXM),IND2(MAXM),ind(maxn),jnd(maxn)
      INTEGER KAND1(MAXNUM),KAND2(MAXNUM),KORNR(MAXNUM,4)
      INTEGER LEFT,KOUNT,kount1,NUM,num1,num2,num3,hdep
      INTEGER I,J,L,M,N,empty,ib,ie,le,tm,nc,jk,tel,tel2
      integer numdep(maxn),ja,jb,index(maxn),typ(maxn)
      integer indhulp,start,ii,istarx,nointer,indoutl(n)
      integer notingp,seed,starttel,a(maxn),ntot,nsub
      integer boxpl,whisk,dattm(maxn),numdattm,angi
      integer numk,numk1
	integer alpha1
      double precision X(N),Y(N),WX(MAXN),WY(MAXN),dpf(maxn)
	double precision interpx(n*2), interpy(n*2)
      double precision ANGLE(MAXM),D(MAXM),alpha(maxnum)
      double precision beta(maxn),wx1(maxn),wy1(maxn)
      double precision angz(maxn),angy1(maxn),angy2(maxn)
      double precision px(maxn,4),py(maxn,4),gamma(maxn)
      double precision lambda(maxn),lambdanc
      double precision wxjb1,wyjb1,wxja1,wyja1,hulp
      double precision PI,PI2,EPS,xcord1,ycord1,xsum,ysum
      double precision xcordp,ycordp,fac,xdev,ydev
      double precision XCORD,YCORD,E1,E2,F1,F2,G1,G2,ANG1
      double precision sum,tukmed(2),tukm(2),dist,c
      double precision star(maxn*3,2),ang,xmean,ymean
      double precision zori(maxn,2),datatyp(n,3)
      double precision interpol(MAXN*2,2),pxpy(n,3)
      double precision datatyp2(n,2),xrand(maxn),yrand(2)


      data seed/256/

	
      PI=DACOS(DBLE(-1.0))
      PI2=PI/2.0
      EPS=0.0000001
      fac=1000000000.0
      do 11 i=1,maxn
         ncirq(i)=0
         mcirq(i)=0
         nrank(i)=0
         f(i)=0
         ind(i)=0
         jnd(i)=0
         numdep(i)=0
         index(i)=0
         typ(i)=0
         a(i)=0
         dattm(i)=0
11    continue
      do 12 i=1,maxm
         jlv(i)=0
         jrv(i)=0
         ind1(i)=0
         ind2(i)=0
12    continue
      do 13 i=1,maxnum
         kand1(i)=0
         kand2(i)=0
         kornr(i,1)=0
         kornr(i,2)=0
         kornr(i,3)=0
         kornr(i,4)=0
13    continue
      left=0
      kount=0
      kount1=0
      num=0
      num1=0
      num2=0
      num3=0
      hdep=0
      j=0
      l=0
      empty=0
      ib=0
      ie=0
      le=0
      tm=0
      nc=0
      jk=0
      tel=0
      tel2=0
      ja=0
      jb=0
      indhulp=0
      start=0
      ii=0
      istarx=0
      notingp=0
      starttel=0
      ntot=0
      nsub=0
      boxpl=0
      numdattm=0
      angi=0
      nointer=0
      do 14 i=1,maxn
         wx(i)=0.0
         wy(i)=0.0
         dpf(i)=0.0
         beta(i)=0.0
         angz(i)=0.0
         angy1(i)=0.0
         angy2(i)=0.0
         px(i,1)=0.0
         px(i,2)=0.0
         px(i,3)=0.0
         px(i,4)=0.0
         py(i,1)=0.0
         py(i,2)=0.0
         py(i,3)=0.0
         py(i,4)=0.0
         wx1(i)=0.0
         wy1(i)=0.0
         lambda(i)=0.0
         gamma(i)=0.0
         zori(i,1)=0.0
         zori(i,2)=0.0
         xrand(i)=0.0
14    continue
      do 16 i=1,maxm
         angle(i)=0.0
         d(i)=0.0
16    continue
      do 17 i=1,maxnum
         alpha(i)=0.0
17    continue
      do 18 i=1,n*3
         star(i,1)=0.0
         star(i,2)=0.0
18    continue
      kstar=0
      lambdanc=0.0
      wxjb1=0.0
      wyjb1=0.0
      wxja1=0.0
      wyja1=0.0
      hulp=0.0
      xcord=0.0
      ycord=0.0
      xcord1=0.0
      ycord1=0.0
      xcordp=0.0
      ycordp=0.0
      xsum=0.0
      ysum=0.0
      e1=0.0
      e2=0.0
      f1=0.0
      f2=0.0
      g1=0.0
      g2=0.0
      ang1=0.0
      ang=0.0
      sum=0.0
      tukmed(1)=0.0
      tukmed(2)=0.0
      tukm(1)=0.0
      tukm(2)=0.0
      dist=0.0
      c=0.0
      yrand(1)=0.0
      yrand(2)=0.0
      if (n.lt.20) whisk=1
C
C     Standardization of the data set.
C
      xmean=0.0
      ymean=0.0
      xdev=0.0
      ydev=0.0
      DO 2 I=1,N
         xmean=xmean+x(i)
         ymean=ymean+y(i)
2     CONTINUE
      xmean=xmean/n
      ymean=ymean/n
      DO 3 I=1,n
         xdev=xdev+((x(i)-xmean)*(x(i)-xmean))
         ydev=ydev+((y(i)-ymean)*(y(i)-ymean))
3     continue
      xdev=xdev/(n-1)
      ydev=ydev/(n-1)
      xdev=dsqrt(xdev)
      ydev=dsqrt(ydev)
      DO 5 I=1,N
         if (xdev.gt.eps) x(i)=(x(i)-xmean)/xdev
         if (ydev.gt.eps) y(i)=(y(i)-ymean)/ydev
         zori(i,1)=x(i)
         zori(i,2)=y(i)
	 WX(I)=X(I)
	 WY(I)=Y(I)
         wx1(i)=x(i)
         wy1(i)=y(i)
	 NCIRQ(I)=I
	 MCIRQ(I)=I
 5    CONTINUE
C
C     If n is large, take a subset.
C
      nsub=2500
      if (n.gt.nsub) then
         ntot=n
         n=nsub
         call bp_rdraw(a,ntot,seed,n)
         do 41 i=1,n
            x(i)=wx1(a(i))
            y(i)=wy1(a(i))
            wx(i)=x(i)
            wy(i)=y(i)
41       continue
      endif
C      
C     Test whether more than half of the points lie on a vertical line.
C
      boxpl=0
      numdep(1)=1
      tel=1
      i=2
      j=1
42    if (dabs(x(i)-x(j)).lt.eps) then
         numdep(j)=numdep(j)+1
         if (numdep(j).gt.idnint(dble(n/2))) then
            boxpl=1
            if (xdev.gt.eps) then
            interpol(1,1)=x(j)*xdev+xmean
            interpol(2,1)=x(i)*xdev+xmean
            else
            interpol(1,1)=x(j)
            interpol(2,1)=x(i)
            endif
            if (ydev.gt.eps) then
            interpol(1,2)=y(j)*ydev+ymean
            interpol(2,2)=y(i)*ydev+ymean
            else
            interpol(1,2)=y(j)
            interpol(2,2)=y(i)
            endif
            goto 610
         endif
         i=i+1
         j=1
         goto 42
      else
         j=j+1
         if (j.eq.i) then
            numdep(i)=1
            i=i+1
            j=1
            tel=tel+1
            if (tel.gt.idnint(dble(n/2))+1) goto 46
         endif
         goto 42
      endif 
46    continue
C
C     To test whether the data points are in general position,
C     we first check whether any two points coincide.
C     When we sort the points with respect to the x-coordinate, we
C     permute NCIRQ and MCIRQ in the same way, to initialize them.
C
      notingp=0
1      CALL BP_SORT(WX,NCIRQ,MCIRQ,WY,N,JLV,JRV)
      I=0
10    I=I+1
      IF (I+1.GT.N) GOTO 15
      J=I+1
      IF (WX(I).NE.WX(J)) GOTO 10
      IF (WY(I).EQ.WY(J)) THEN
C         The data are not in general position:
C         The points NCIRQ(I) and NCIRQ(J) coincide.
         GOTO 15
      ELSE
         if (wy(i).lt.wy(j)) then
            iv=mcirq(j)
            mcirq(j)=mcirq(i)
            mcirq(i)=iv
         endif
         IF (J+1.LE.N) THEN
            L=J+1
            IF (WX(I).EQ.WX(L)) THEN
C        The data are not in general position:
C        The points NCIRQ(I), NCIRQ(J), and NCIRQ(L) lie on a vertical line.
               if ((i.eq.1).and.(j.eq.2).and.(l.eq.3)) goto 15
               GOTO 210
            ELSE
               GOTO 10
            ENDIF
         ENDIF
      ENDIF
C    
C     Compute all the angles formed by pairs of data points.
C
15    M=((N*(N-1))/2)
      L=1
      DO 20 I=1,N
	 DO 25 J=I+1,N
	    IF (X(I).EQ.X(J)) THEN
	       ANGLE(L)=PI2
            ELSE
	       ANGLE(L)=DATAN((Y(I)-Y(J))/(X(I)-X(J)))
	       IF (ANGLE(L).LE.0.0) ANGLE(L)=ANGLE(L)+PI
            ENDIF
	    IND1(L)=I
	    IND2(L)=J
	    L=L+1
25    CONTINUE
20    CONTINUE
C 
C     Sort all the angles and permute IND1 and IND2 in the same way.
C     To avoid using several SORT-routines, we will always permute
C     two integer arrays and one double precision array.
C
      CALL BP_SORT(ANGLE,IND1,IND2,D,M,JLV,JRV)
C     
C     If all points are collinear, use a univariate boxplot.
C
      boxpl=0
      angi=1
      if (angle(m)-angle(1).le.eps) then
         boxpl=1
         goto 27
      endif
      if ((angle(1).le.eps).and.(angle(m).ge.pi*2-eps)) then
         do 26 i=2,m-1
            if ((angle(i).le.eps).or.(angle(i).ge.pi*2-eps)) then
               boxpl=1
            else
               boxpl=0
               goto 27
            endif
26       continue
      endif
C
C     Count whether 50% of the data points are collinear
C
      boxpl=0
      ang=angle(1)
      angi=1
      tel=1
      do 126 i=2,m
         if (dabs(angle(i)-ang).le.eps) then
            tel=tel+1
            if (tel.gt.(n*(n-2)/8)) then
               boxpl=1
               goto 27
            endif
         else
            ang=angle(i)
            angi=i 
            tel=1
         endif
126   continue

27    if (boxpl.eq.1) then
C     All points collinear
         if (xdev.gt.eps) then
            interpol(1,1)=x(ind1(angi))*xdev+xmean
            interpol(2,1)=x(ind2(angi))*xdev+xmean
         else
            interpol(1,1)=x(ind1(angi))
            interpol(2,1)=x(ind2(angi))
         endif
         if (ydev.gt.eps) then
            interpol(1,2)=y(ind1(angi))*ydev+ymean
            interpol(2,2)=y(ind2(angi))*ydev+ymean
         else
            interpol(1,2)=y(ind1(angi))
            interpol(2,2)=y(ind2(angi))
         endif
         goto 610
      endif
C 
C  Test whether any three points are collinear.
C
      LEFT=1
30    ANG1=ANGLE(LEFT)
      DO 35 J=LEFT+1,M
	 IF (ANGLE(J).GT.ANG1) THEN
	    LEFT=J
	    GOTO 30
         ELSE
            DO 36 I=LEFT,J-1
               IF ((IND1(I).EQ.IND1(J)).or. 
     +	          (IND1(I).EQ.IND2(J))) THEN 
C	       The data are not in general position:
C	       The points IND1(J), IND2(J), and IND2(I) are collinear.
	       GOTO 210
	       ENDIF
	       IF ((IND2(I).EQ.IND1(J)).or.
     +            (IND2(I).EQ.IND2(J))) THEN
C	       The data are not in general position:
C	       The points IND1(J), IND2(J), and IND1(I) are collinear.
	       GOTO 210
	       ENDIF
36          CONTINUE
         ENDIF
35    CONTINUE
      goto 37

C
C     If the data are not in general position, use dithering.
C
210   notingp=1
      fac=fac/10.0
      do 211 i=1,n
         ncirq(i)=i
         mcirq(i)=i
         call nbp_norrandp(2,seed,tukmed)
         wx(i)=wx1(i)+tukmed(1)/fac
         x(i)=wx(i)
         wy(i)=wy1(i)+tukmed(2)/fac
         y(i)=wy(i)
211   continue
      goto 1

C
C      The data are in general position.
C      
37    if (notingp.eq.1) then
      do 212 i=1,n
          numdep(i)=0
212   continue
      endif
      kstar=0
      do 138 i=1,n
         numdep(i)=0
138   continue
      do 38 i=1,n
         call bp_depth(x(i),y(i),n,x,y,beta,f,dpf,jlv,jrv,hdep)
         if (hdep.gt.kstar) kstar=hdep
         numdep(hdep)=numdep(hdep)+1
38    continue
C
C     Calculation of the Tukey median.
C     
      tm=0
      xsum=0
      ysum=0
      tukmed(1)=0
      tukmed(2)=0
      if (n.le.3) then
         do 40 i=1,n
         xsum=xsum+x(i)
         ysum=ysum+y(i)
      if (xdev.gt.eps) datatyp(i,1)=x(i)*xdev+xmean
      if (ydev.gt.eps) datatyp(i,2)=y(i)*ydev+ymean
40       continue
         tukm(1)=xsum/n
      if (xdev.gt.eps) tukm(1)=tukm(1)*xdev+xmean
         tukm(2)=ysum/n
      if (ydev.gt.eps) tukm(2)=tukm(2)*ydev+ymean
         goto 800 
      endif
      empty=0
      ib=kstar
      ie=int(dble(n)/2)
180     le=ie-ib
        if (le.lt.0) then
           nointer=1
           le=0
        endif
        if (le.eq.0) goto 185
      CALL BP_ISODEPTH(N,M,X,Y,MAXN,MAXM,MAXNUM,NRANK,D,F,BETA,
     +  KAND1,KAND2,ALPHA,IND1,IND2,NCIRQ,MCIRQ,ANGLE,KORNR,L,
     +  JRV,JLV,DPF,NUM,ib+nbp_nceil(le,2),EMPTY)
        if (empty.eq.1) ie=ib+nbp_nceil(le,2)
        if (empty.eq.0) ib=ib+nbp_nceil(le,2)
        if (le.eq.1) goto 185
        goto 180
185   CALL BP_ISODEPTH(N,M,X,Y,MAXN,MAXM,MAXNUM,NRANK,D,F,BETA,
     +  KAND1,KAND2,ALPHA,IND1,IND2,NCIRQ,MCIRQ,ANGLE,KORNR,L,
     +  JRV,JLV,DPF,NUM,ib,EMPTY)
C
C     Scan KORNR and compute coordinates of the vertices. 
C
186   KOUNT=0
      I=1
      E1=Y(KORNR(I,2))-Y(KORNR(I,1))
      F1=X(KORNR(I,1))-X(KORNR(I,2))
      G1=X(KORNR(I,1))*(Y(KORNR(I,2))-Y(KORNR(I,1)))
     +   -Y(KORNR(I,1))*(X(KORNR(I,2))-X(KORNR(I,1)))
      E2=Y(KORNR(I,4))-Y(KORNR(I,3))
      F2=X(KORNR(I,3))-X(KORNR(I,4))
      G2=X(KORNR(I,3))*(Y(KORNR(I,4))-Y(KORNR(I,3)))
     +   -Y(KORNR(I,3))*(X(KORNR(I,4))-X(KORNR(I,3)))
      XCORD=(-F2*G1+F1*G2)/(E2*F1-E1*F2)
      YCORD=(-E2*G1+E1*G2)/(E1*F2-E2*F1)
         wx(kount+1)=xcord
         wy(kount+1)=ycord
      if (tm.eq.0) then
         xsum=xcord
         ysum=ycord
      endif
      xcord1=xcord
      ycord1=ycord
      xcordp=xcord
      ycordp=ycord
      KOUNT=KOUNT+1
      I=I+1
      if (num.eq.1) goto 195 
190   IF ((KORNR(I,1).EQ.KORNR(I-1,1)).AND.(KORNR(I,2).EQ.KORNR(I-1,2))
     +.AND.(KORNR(I,3).EQ.KORNR(I-1,3).AND.KORNR(I,4).EQ.KORNR(I-1,4)))
     +THEN
	I=I+1
      ELSE
        IF ((KORNR(I,1).EQ.KORNR(1,1)).AND.(KORNR(I,2).EQ.KORNR(1,2))
     +    .AND.(KORNR(I,3).EQ.KORNR(1,3).AND.KORNR(I,4).EQ.KORNR(1,4))) 
     +	THEN
	  GOTO 195
        ELSE
          E1=Y(KORNR(I,2))-Y(KORNR(I,1))
          F1=X(KORNR(I,1))-X(KORNR(I,2))
          G1=X(KORNR(I,1))*(Y(KORNR(I,2))-Y(KORNR(I,1)))
     +       -Y(KORNR(I,1))*(X(KORNR(I,2))-X(KORNR(I,1)))
          E2=Y(KORNR(I,4))-Y(KORNR(I,3))
          F2=X(KORNR(I,3))-X(KORNR(I,4))
          G2=X(KORNR(I,3))*(Y(KORNR(I,4))-Y(KORNR(I,3)))
     +       -Y(KORNR(I,3))*(X(KORNR(I,4))-X(KORNR(I,3)))
          XCORD=(-F2*G1+F1*G2)/(E2*F1-E1*F2)
          YCORD=(-E2*G1+E1*G2)/(E1*F2-E2*F1)
          if (((dabs(xcord-xcordp).lt.eps).and.
     +        (dabs(ycord-ycordp).lt.eps)).or. 
     +       ((dabs(xcord-xcord1).lt.eps).and.
     +        (dabs(ycord-ycord1).lt.eps))) then
             i=i+1
          else
             xcordp=xcord
             ycordp=ycord
             wx(kount+1)=xcord
             wy(kount+1)=ycord
      if (tm.eq.0) then
             xsum=xsum+xcord
             ysum=ysum+ycord
      endif
	  KOUNT=KOUNT+1
	  I=I+1
          endif
        ENDIF
      ENDIF
      IF (I.NE.(NUM+1)) GOTO 190
195   if (tm.eq.2) goto 500
      if (tm.eq.1) goto 300
c    
c     Calculation of the center of gravity.
c     
      if (tm.eq.0) then
         if (kount.gt.1) then
         do 205 i=1,kount
         wx(i)=wx(i)-(xsum/kount)
         wy(i)=wy(i)-(ysum/kount)
205      continue
         sum=0
         tukmed(1)=0
         tukmed(2)=0
         do 206 i=1,kount-1
         sum=sum+dabs(wx(i)*wy(i+1)-wx(i+1)*wy(i))
         tukmed(1)=tukmed(1)+
     +  ((wx(i)+wx(i+1))*dabs(wx(i)*wy(i+1)-wx(i+1)*wy(i)))   
         tukmed(2)=tukmed(2)+
     +  ((wy(i)+wy(i+1))*dabs(wx(i)*wy(i+1)-wx(i+1)*wy(i)))   
206      continue
         sum=sum+dabs(wx(kount)*wy(1)-wx(1)*wy(kount))
         tukmed(1)=tukmed(1)+
     +  ((wx(kount)+wx(1))*dabs(wx(kount)*wy(1)-wx(1)*wy(kount)))   
         tukmed(2)=tukmed(2)+
     +  ((wy(kount)+wy(1))*dabs(wx(kount)*wy(1)-wx(1)*wy(kount)))   
         tukmed(1)=(tukmed(1)/(3*sum))+(xsum/kount)
         tukmed(2)=(tukmed(2)/(3*sum))+(ysum/kount)
         else
         tukmed(1)=xsum
         tukmed(2)=ysum
         endif
207      if (xdev.gt.eps) then
            xcord=tukmed(1)*xdev+xmean
         else
            xcord=tukmed(1)
         endif
         if (ydev.gt.eps) then
            ycord=tukmed(2)*ydev+ymean
         else
            ycord=tukmed(2)
         endif
            tukm(1)=xcord
            tukm(2)=ycord
      endif
      tm=1
      if (nointer.eq.1) goto 800
C
C     Calculation of correct value of k. 
C
      nc=int(n*alpha1/100)
      j=kstar+1
200   j=j-1
      if (numdep(kstar).le.nc) then
         numdep(kstar)=numdep(kstar)+numdep(j-1) 
         goto 200
      endif
      k=j+1
      numk1=numdep(kstar)
      numk=numk1-numdep(k-1)
      lambdanc=dble(nc-numk)/dble(numk1-numk)
C
C     Calculation of the vertices of Dk.
C
      CALL BP_ISODEPTH(N,M,X,Y,MAXN,MAXM,MAXNUM,NRANK,D,F,BETA,
     +  KAND1,KAND2,ALPHA,IND1,IND2,NCIRQ,MCIRQ,ANGLE,KORNR,L,
     +  JRV,JLV,DPF,NUM,k,EMPTY)
      if (empty.eq.0) goto 186
300   continue
      kount1=kount
      do 201 i=1,kount1
      wx1(i)=wx(i)
      wy1(i)=wy(i)
201   continue
C
C     Calculation of the vertices of Dk-1.
C 
      tm=2
      CALL BP_ISODEPTH(N,M,X,Y,MAXN,MAXM,MAXNUM,NRANK,D,F,BETA,
     +  KAND1,KAND2,ALPHA,IND1,IND2,NCIRQ,MCIRQ,ANGLE,KORNR,L,
     +  JRV,JLV,DPF,NUM,k-1,EMPTY)
      if (empty.eq.0) goto 186
500   continue
C
C     If Dk-1 is a line segment, use a univariate boxplot.
C 
      if (dabs(wx(2)-wx(1)).gt.eps) then
         hulp=(wy(2)-wy(1))/(wx(2)-wx(1))
         do 506 i=3,kount
         if (dabs(wy(i)-wy(1)-hulp*(wx(i)-wx(1))).gt.eps*10) goto 505
506      continue
         boxpl=2
         goto 509
      else
         do 504 i=3,kount
         if (dabs(wx(i)-wx(1)).gt.eps*10) goto 505
504      continue
         boxpl=2
         goto 509  
      endif
509   continue
      do 507 i=1,n
      if (xdev.gt.eps) datatyp(i,1)=x(i)*xdev+xmean
      if (ydev.gt.eps) datatyp(i,2)=y(i)*ydev+ymean
      datatyp(i,3)=4.0
507   continue
      do 508 i=1,kount
      if (xdev.gt.eps) interpol(i,1)=wx(i)*xdev+xmean
      if (ydev.gt.eps) interpol(i,2)=wy(i)*ydev+ymean
508   continue
      goto 610

505      if (n.ge.10) then
      do 202 i=1,kount1
      wx1(i)=wx1(i)-tukmed(1)
      wy1(i)=wy1(i)-tukmed(2)
202   continue
      do 203 i=1,kount
      wx(i)=wx(i)-tukmed(1)
      wy(i)=wy(i)-tukmed(2)
203   continue
      do 204 i=1,n
      x(i)=x(i)-tukmed(1)
      y(i)=y(i)-tukmed(2)
204   continue
      do 1204 i=1,ntot
      zori(i,1)=zori(i,1)-tukmed(1)
      zori(i,2)=zori(i,2)-tukmed(2)
1204   continue
C
C     Compute angles of the data points.
C
      numdattm=0
      do 400 i=1,n
      index(i)=i
      if ((dabs(x(i)).lt.eps).and.(dabs(y(i)).lt.eps)) then
         numdattm=numdattm+1
         dattm(numdattm)=i
         angz(i)=1000.0
      else
      dist=dsqrt(x(i)*x(i)+y(i)*y(i))
      xcord=x(i)/dist
      ycord=y(i)/dist 
      if (dabs(xcord).gt.dabs(ycord)) then
         if (xcord.ge.0.0) then
            angz(i)=dasin(ycord)
            if (angz(i).lt.0.0) angz(i)=angz(i)+pi*2
         else
            angz(i)=pi-dasin(ycord)
         endif
      else
         if (ycord.ge.0.0) then
            angz(i)=dacos(xcord)
         else
            angz(i)=pi*2-dacos(xcord)
         endif
      endif
         if (angz(i).ge.(pi*2-eps)) angz(i)=0.0
      endif
400   continue
C     numdattm is the number of data points equal to the Tukey median.
      call bp_sort(angz,index,index,dpf,n,jlv,jrv)
      do 401 i=1,n
         indoutl(i)=index(i)
401   continue
      n=n-numdattm
C
C     Compute angles of the vertices of the depth region Dk.
C
      if (kount1.eq.1) goto 412
      do 410 i=1,kount1
      ind(i)=i
      dist=dsqrt(wx1(i)*wx1(i)+wy1(i)*wy1(i))
      if (dist.gt.eps) then
      xcord=wx1(i)/dist
      ycord=wy1(i)/dist 
      if (dabs(xcord).gt.dabs(ycord)) then
         if (xcord.ge.0.0) then
            angy1(i)=dasin(ycord)
            if (angy1(i).lt.0.0) angy1(i)=angy1(i)+pi*2
         else
            angy1(i)=pi-dasin(ycord)
         endif
      else
         if (ycord.ge.0.0) then
            angy1(i)=dacos(xcord)
         else
            angy1(i)=pi*2-dacos(xcord)
         endif
      endif
         if (angy1(i).ge.(pi*2-eps)) angy1(i)=0.0
      else
      nointer=1
      endif
410   continue
      call bp_sort(angy1,ind,ind,dpf,kount1,jlv,jrv)
C
C     Compute angles of the vertices of the depth region Dk-1.
C
412      do 510 i=1,kount
      jnd(i)=i
      dist=dsqrt(wx(i)*wx(i)+wy(i)*wy(i))
      if (dist.gt.eps) then
      xcord=wx(i)/dist
      ycord=wy(i)/dist 
      if (dabs(xcord).gt.dabs(ycord)) then
         if (xcord.ge.0.0) then
            angy2(i)=dasin(ycord)
            if (angy2(i).lt.0.0) angy2(i)=angy2(i)+pi*2
         else
            angy2(i)=pi-dasin(ycord)
         endif
      else
         if (ycord.ge.0.0) then
            angy2(i)=dacos(xcord)
         else
            angy2(i)=pi*2-dacos(xcord)
         endif
      endif
         if (angy2(i).ge.(pi*2-eps)) angy2(i)=0.0
      else
      nointer=1
      endif
510   continue
      call bp_sort(angy2,jnd,jnd,dpf,kount,jlv,jrv)
C
C     Calculation of arrays px and py for Dk.
C
      if (kount1.eq.1) goto 431
      jk=0
      wx1(kount1+1)=wx1(1)
      wy1(kount1+1)=wy1(1)
      angy1(kount1+1)=angy1(1)
      ind(kount1+1)=ind(1)
      if (angz(1).lt.angy1(1)) j=kount1
      if (angz(1).ge.angy1(1)-eps) j=1
      do 429 i=1,n  
420   if ((angz(i).ge.angy1(j+1)-eps).and.(jk.eq.0)) then
         j=j+1
         if (j.eq.kount1) jk=1
         if (j.eq.kount1+1) j=1
         goto 420
      endif
      if ((dabs(wx1(ind(j+1))-wx1(ind(j))).gt.eps).and.
     +    (dabs(x(index(i))).gt.eps)) then
      dist=(y(index(i))/x(index(i)))-
     + ((wy1(ind(j+1))-wy1(ind(j)))/(wx1(ind(j+1))-wx1(ind(j))))
      px(i,1)=(wy1(ind(j))*wx1(ind(j+1))-
     + wx1(ind(j))*wy1(ind(j+1)))/(wx1(ind(j+1))-wx1(ind(j)))
      px(i,1)=px(i,1)/dist
      py(i,1)=px(i,1)*y(index(i))/x(index(i))
      else
         if (dabs(wx1(ind(j+1))-wx1(ind(j))).le.eps) then
      px(i,1)=wx1(ind(j))
      py(i,1)=px(i,1)*y(index(i))/x(index(i))
         else
      px(i,1)=0.0
      py(i,1)=((wy1(ind(j))-wy1(ind(j+1)))*wx1(ind(j))/
     +        (wx1(ind(j+1))-wx1(ind(j))))+wy1(ind(j))
         endif
      endif
429   continue
C
C     Calculation of arrays px and py for Dk-1.
C
431      jk=0
      wx(kount+1)=wx(1)
      wy(kount+1)=wy(1)
      angy2(kount+1)=angy2(1)
      jnd(kount+1)=jnd(1)
      if (angz(1).lt.angy2(1)) j=kount
      if (angz(1).ge.angy2(1)-eps) j=1
      do 529 i=1,n  
520   if ((angz(i).ge.angy2(j+1)-eps).and.(jk.eq.0)) then
         j=j+1
         if (j.eq.kount) jk=1
         if (j.eq.kount+1) j=1
         goto 520
      endif
      if ((dabs(wx(jnd(j+1))-wx(jnd(j))).gt.eps).and.
     +    (dabs(x(index(i))).gt.eps)) then
      dist=(y(index(i))/x(index(i)))-
     + ((wy(jnd(j+1))-wy(jnd(j)))/(wx(jnd(j+1))-wx(jnd(j))))
      px(i,2)=(wy(jnd(j))*wx(jnd(j+1))-
     + wx(jnd(j))*wy(jnd(j+1)))/(wx(jnd(j+1))-wx(jnd(j)))
      px(i,2)=px(i,2)/dist
      py(i,2)=px(i,2)*y(index(i))/x(index(i))
      else
         if (dabs(wx(jnd(j+1))-wx(jnd(j))).le.eps) then
      px(i,2)=wx(jnd(j))
      py(i,2)=px(i,2)*y(index(i))/x(index(i))
         else
      px(i,2)=0.0
      py(i,2)=((wy(jnd(j))-wy(jnd(j+1)))*wx(jnd(j))/
     +  (wx(jnd(j+1))-wx(jnd(j))))+wy(jnd(j))
         endif
      endif
529   continue
      do 540 i=1,n
      if (kount1.eq.1) then
         px(i,1)=wx1(1)
         py(i,1)=wy1(1)
      endif
540   continue
C
C     Mergesort of angy1 and angy2 to obtain gamma
C     and calculation of arrays px and py.
C
      if (kount1.eq.1) then
         do 599 i=1,kount
            gamma(i)=angy2(i)
            px(i,3)=wx1(i)
            py(i,3)=wy1(i)
            px(i,4)=wx(i)
            py(i,4)=wy(i)
599      continue
         goto 602
      endif
      ja=1
      jb=1
      i=1
600   if (angy1(ja).le.angy2(jb)) then
         if (ja.le.kount1) then
            gamma(i)=angy1(ja)
            px(i,3)=wx1(ind(ja))
            py(i,3)=wy1(ind(ja))
            if (jb.eq.1) then
               wxjb1=wx(jnd(kount))
               wyjb1=wy(jnd(kount))
            else
               wxjb1=wx(jnd(jb-1))
               wyjb1=wy(jnd(jb-1))
            endif
       if ((dabs(wx(jnd(jb))-wxjb1).gt.eps).and.
     +     (dabs(wx1(ind(ja))).gt.eps)) then
       dist=(wy1(ind(ja))/wx1(ind(ja)))-
     + ((wy(jnd(jb))-wyjb1)/(wx(jnd(jb))-wxjb1))
            px(i,4)=(wyjb1*wx(jnd(jb))-
     + wxjb1*wy(jnd(jb)))/(wx(jnd(jb))-wxjb1)
            px(i,4)=px(i,4)/dist
            py(i,4)=px(i,4)*wy1(ind(ja))/wx1(ind(ja))
        else
           if (dabs(wx(jnd(jb))-wxjb1).le.eps) then
               px(i,4)=wxjb1
               py(i,4)=px(i,4)*wy1(ind(ja))/wx1(ind(ja))
           else
               px(i,4)=0.0
               py(i,4)=((wyjb1-wy(jnd(jb)))*wxjb1/
     +       (wx(jnd(jb))-wxjb1))+wyjb1
           endif
        endif
            ja=ja+1
            i=i+1
         else
            angy1(ja)=angy1(ja)+10
         endif
      else
         if (jb.le.kount) then
            gamma(i)=angy2(jb)
            px(i,4)=wx(jnd(jb))
            py(i,4)=wy(jnd(jb))
            if (ja.eq.1) then
               wxja1=wx1(ind(kount1))
               wyja1=wy1(ind(kount1))
            else
               wxja1=wx1(ind(ja-1))
               wyja1=wy1(ind(ja-1))
            endif
       if ((dabs(wx1(ind(ja))-wxja1).gt.eps).and.
     +     (dabs(wx(jnd(jb))).gt.eps)) then
       dist=(wy(jnd(jb))/wx(jnd(jb)))-
     + ((wy1(ind(ja))-wyja1)/(wx1(ind(ja))-wxja1))
            px(i,3)=(wyja1*wx1(ind(ja))-
     + wxja1*wy1(ind(ja)))/(wx1(ind(ja))-wxja1)
            px(i,3)=px(i,3)/dist
            py(i,3)=px(i,3)*wy(jnd(jb))/wx(jnd(jb))
       else
          if (dabs(wx1(ind(ja))-wxja1).le.eps) then
             px(i,3)=wxja1
             py(i,3)=px(i,3)*wy(jnd(jb))/wx(jnd(jb))
          else
             px(i,3)=0.0
             py(i,3)=((wyja1-wy1(ind(ja)))*wxja1/
     +    (wx1(ind(ja))-wxja1))+wyja1
          endif
      endif
            jb=jb+1
            i=i+1
         else
            angy2(jb)=angy2(jb)+10
         endif
      endif
      if (i.le.kount1+kount) goto 600
C
C     Interpolation of two arrays px and py.
C      
602      c=3.0
       if (nointer.eq.1) then
          kount1=0
       endif
       do 605 i=1,kount1+kount
         if (nointer.eq.0) then
       wx(i)=lambdanc*px(i,4)+(1-lambdanc)*px(i,3)
       wy(i)=lambdanc*py(i,4)+(1-lambdanc)*py(i,3)
         endif
       if (xdev.gt.eps) then
          xcord=(wx(i)+tukmed(1))*xdev+xmean
       else
          xcord=wx(i)+tukmed(1)
       endif
       if (ydev.gt.eps) then
          ycord=(wy(i)+tukmed(2))*ydev+ymean
       else
          ycord=wy(i)+tukmed(2)
       endif
   
       interpol(i,1)=xcord
       interpol(i,2)=ycord 
	  interpx(i) = xcord
	  interpy(i) = ycord
605    continue
       if (nointer.eq.1) then
C      No interpolation, too many points coincide.
           do 1605 i=1,n
           if (xdev.gt.eps) then
              xcord=(x(index(i))+tukmed(1))*xdev+xmean
           else
              xcord=x(index(i))+tukmed(1)
           endif
           if (ydev.gt.eps) then
              ycord=(y(index(i))+tukmed(2))*ydev+ymean
           else
              ycord=y(index(i))+tukmed(2)
           endif
              datatyp(i,1)=xcord
              datatyp(i,2)=ycord
              datatyp(i,3)=1.0
1605       continue
           do 1606 i=1,numdattm
           if (xdev.gt.eps) then
              xcord=(x(dattm(i))+tukmed(1))*xdev+xmean
           else
              xcord=x(dattm(i))+tukmed(1)
           endif
           if (ydev.gt.eps) then
              ycord=(y(dattm(i))+tukmed(2))*ydev+ymean
           else
              ycord=y(dattm(i))+tukmed(2)
           endif
               datatyp(n+i,1)=xcord
               datatyp(n+i,2)=ycord
               datatyp(n+i,3)=0.0
1606      continue
          goto 610
       endif
C
C      Repeat some calculations for the whole data set.
C
      if (ntot.gt.nsub) then
         n=ntot
         do 1399 i=1,n
         x(i)=zori(i,1)
         y(i)=zori(i,2)
1399      continue
      else
         goto 1531
      endif
C
C     Compute angles of the data points.
C
      do 1300 i=1,n
         x(i)=x(i)-tukmed(1)
         y(i)=y(i)-tukmed(2)
1300  continue
      numdattm=0
      do 1400 i=1,n
      index(i)=i
      if ((dabs(x(i)).lt.eps).and.(dabs(y(i)).lt.eps)) then
         angz(i)=1000.0
         numdattm=numdattm+1
         dattm(numdattm)=i
      else
      dist=dsqrt(x(i)*x(i)+y(i)*y(i))
      xcord=x(i)/dist
      ycord=y(i)/dist 
      if (dabs(xcord).gt.dabs(ycord)) then
         if (xcord.ge.0.0) then
            angz(i)=dasin(ycord)
            if (angz(i).lt.0.0) angz(i)=angz(i)+pi*2
         else
            angz(i)=pi-dasin(ycord)
         endif
      else
         if (ycord.ge.0.0) then
            angz(i)=dacos(xcord)
         else
            angz(i)=pi*2-dacos(xcord)
         endif
      endif
         if (angz(i).ge.(pi*2-eps)) angz(i)=0.0
      endif
1400   continue
      call bp_sort(angz,index,index,dpf,n,jlv,jrv)
      do 1401 i=1,n
         indoutl(i)=index(i)
1401   continue
      n=n-numdattm
C
C     Calculation of arrays px and py for B.
C
1431      jk=0
      wx(kount+kount1+1)=wx(1)
      wy(kount+kount1+1)=wy(1)
      gamma(kount+kount1+1)=gamma(1)
      if (angz(1).lt.gamma(1)) j=kount+kount1
      if (angz(1).ge.gamma(1)-eps) j=1
      do 1529 i=1,n  
1520   if ((angz(i).ge.gamma(j+1)-eps).and.(jk.eq.0)) then
         j=j+1
         if (j.eq.kount+kount1) jk=1
         if (j.eq.kount+kount1+1) j=1
         goto 1520
      endif
      if ((dabs(wx(j+1)-wx(j)).gt.eps).and.
     +    (dabs(x(index(i))).gt.eps)) then
      dist=(y(index(i))/x(index(i)))-
     + ((wy(j+1)-wy(j))/(wx(j+1)-wx(j)))
      px(i,1)=(wy(j)*wx(j+1)-
     + wx(j)*wy(j+1))/(wx(j+1)-wx(j))
      px(i,1)=px(i,1)/dist
      py(i,1)=px(i,1)*y(index(i))/x(index(i))
      else
         if (dabs(wx(j+1)-wx(j)).le.eps) then
      px(i,1)=wx(j)
      py(i,1)=px(i,1)*y(index(i))/x(index(i))
         else
      px(i,1)=0.0
      py(i,1)=((wy(j)-wy(j+1))*wx(j)/
     +  (wx(j+1)-wx(j)))+wy(j)
         endif
      endif
1529   continue
C
C      Decide on the type of each data point.
C
1531   c=3.0
       num=0
       num1=0
       num2=0
       num3=0
       do 1639 i=1,n
       if (ntot.le.nsub) then
          px(i,1)=lambdanc*px(i,2)+(1-lambdanc)*px(i,1)
          py(i,1)=lambdanc*py(i,2)+(1-lambdanc)*py(i,1)
       endif
       if (px(i,1)*px(i,1)+py(i,1)*py(i,1).gt.eps) then
          lambda(i)=dsqrt(x(index(i))*x(index(i))+
     + y(index(i))*y(index(i)))/
     + dsqrt(px(i,1)*px(i,1)+py(i,1)*py(i,1))
       else
          if ((dabs(x(index(i))).lt.eps).and.
     +        (dabs(y(index(i))).lt.eps)) then
             lambda(i)=0.0
          else
             lambda(i)=c+1.0
          endif
       endif
       if (lambda(i).lt.eps) then
         num=num+1
         typ(i)=0
         goto 639
       endif
       if (lambda(i).le.1+eps) then
         num1=num1+1
         typ(i)=1
         goto 639
       endif
       if (lambda(i).le.c+eps) then
         num2=num2+1
         typ(i)=2
         goto 639
       endif
       if (lambda(i).gt.c) then
         num3=num3+1
         typ(i)=3
c        indoutl(i)=index(i)
       endif
639    if (xdev.gt.eps) then
          xcord=(x(index(i))+tukmed(1))*xdev+xmean
       else
          xcord=x(index(i))+tukmed(1)
       endif
       if (ydev.gt.eps) then
          ycord=(y(index(i))+tukmed(2))*ydev+ymean
       else
          ycord=y(index(i))+tukmed(2)
       endif
           datatyp(i,1)=xcord
           datatyp(i,2)=ycord
           datatyp(i,3)=dble(typ(i))
1639       continue 
       if (numdattm.ne.0) then
             num=num+numdattm
       do 1644 i=1,numdattm
             typ(n+i)=0
             lambda(n+i)=0.0
             px(n+i,1)=x(dattm(i))
             py(n+i,1)=y(dattm(i))
       if (xdev.gt.eps) then
          xcord=(x(dattm(i))+tukmed(1))*xdev+xmean
       else
          xcord=x(dattm(i))+tukmed(1)
       endif
       if (ydev.gt.eps) then
          ycord=(y(dattm(i))+tukmed(2))*ydev+ymean
       else
          ycord=y(dattm(i))+tukmed(2)
       endif
           datatyp(n+i,1)=xcord
           datatyp(n+i,2)=ycord
           datatyp(n+i,3)=0.0
1644   continue
       endif

      if (whisk.eq.1) goto 1640
      if (whisk.eq.2) goto 1641
1640  do 640 i=1,n
      if (xdev.gt.eps) then
         xcord=(px(i,1)+tukmed(1))*xdev+xmean
      else
         xcord=px(i,1)+tukmed(1)
      endif
      if (ydev.gt.eps) then
         ycord=(py(i,1)+tukmed(2))*ydev+ymean
      else
         ycord=py(i,1)+tukmed(2)
      endif
          pxpy(i,1)=xcord
          pxpy(i,2)=ycord
          pxpy(i,3)=0.0
640   continue
      if (whisk.eq.1) goto 610
      if (whisk.eq.3) goto 610
C
C     Retain only one whisker for each edge.
C    
1641  if (ntot.gt.1000) then
        n=500
        call bp_rdraw(a,ntot,seed,n)
C       do 2642 i=1,n
C       if (xdev.gt.eps) then
C            xcord=(x(index(a(i)))+tukmed(1))*xdev+xmean
C        else
C           xcord=x(index(a(i)))+tukmed(1)
C         endif
C         if (ydev.gt.eps) then
C            ycord=(y(index(a(i)))+tukmed(2))*ydev+ymean
C         else
C            ycord=y(index(a(i)))+tukmed(2)
C         endif
C          datatyp2(i,1)=xcord
C          datatyp2(i,2)=ycord
C2642     continue
C         do 1642 i=1,n
C            beta(i)=angz(a(i))
C            dattm(i)=typ(a(i))
C            dpf(i)=lambda(a(i))
C            px(i,2)=px(a(i),1)
C            py(i,2)=py(a(i),1)
C1642     continue
C         do 1643 i=1,n
C            angz(i)=beta(i)
C            typ(i)=dattm(i)
C            lambda(i)=dpf(i)
C            px(i,1)=px(i,2)
C            py(i,1)=py(i,2)
C1643     continue
C      endif
C      if (num2.ne.0) then
C      tel=1
C      tel2=1
C      gamma(kount1+kount+1)=gamma(1)+(pi*2)
C      i=1
C      ii=1
C649   if (typ(i).ne.2) then
C         i=i+1
C         goto 649
C      endif
C      if (angz(i).lt.gamma(1)) then
C         i=i+1
C         goto 649
C      endif
C      start=i
C648   if (angz(i).ge.gamma(tel+1)) then
C         tel=tel+1
C      goto 648
C      endif
C      hulp=lambda(i)
C      indhulp=i
C      do 650 i=1,start
C         angz(i)=angz(i)+(pi*2)
C650   continue
C      i=start+1
C651   if (typ(i).ne.2) then
C         i=i+1
C         if (i.eq.n+1) i=1
C         goto 651
C      endif
C      if (i.eq.start) then
C         if (tel2.gt.0) then
C
C      if (xdev.gt.eps) then
C         xcord=(px(indhulp,1)+tukmed(1))*xdev+xmean
C      else
C         xcord=px(indhulp,1)+tukmed(1)
C      endif
C      if (ydev.gt.eps) then
C         ycord=(py(indhulp,1)+tukmed(2))*ydev+ymean
C      else
C         ycord=py(indhulp,1)+tukmed(2)
C      endif
C          pxpy(ii,1)=xcord
C          pxpy(ii,2)=ycord
C          pxpy(ii,3)=indhulp
C          ii=ii+1
C         endif
C         goto 660
C      endif
C      if (angz(i).lt.gamma(tel+1)) then
C         tel2=tel2+1
C         if (lambda(i).gt.hulp) then
C            hulp=lambda(i)
C            indhulp=i
C         endif
C         i=i+1
C         if (i.eq.n+1) i=1
C         goto 651
C      else
C      if (tel2.gt.0) then
C
C      if (xdev.gt.eps) then
C         xcord=(px(indhulp,1)+tukmed(1))*xdev+xmean
C      else
C         xcord=px(indhulp,1)+tukmed(1)
C      endif
C      if (ydev.gt.eps) then
C         ycord=(py(indhulp,1)+tukmed(2))*ydev+ymean
C      else
C         ycord=py(indhulp,1)+tukmed(2)
C      endif
C          pxpy(ii,1)=xcord
C          pxpy(ii,2)=ycord
C          pxpy(ii,3)=indhulp
C          ii=ii+1
C      endif
C        tel=tel+1
C        if (tel.eq.kount1+kount+1) tel=1
C        if (angz(i).ge.gamma(tel+1)) then
C           tel=tel+1
C           if (tel.eq.kount1+kount+1) tel=1
C        endif
C        hulp=lambda(i)
C        indhulp=i
C        tel2=0
C        goto 651
C      endif
C660   continue
C      do 661 i=1,start
C         angz(i)=angz(i)-(pi*2)
C661   continue
C 
C     Calculation of star-shaped whiskers.
C
      if (ntot.gt.nsub) goto 800
      if (whisk.ne.4) goto 800

      angz(n+1)=angz(1)
      do 669 i=1,n
         dattm(i)=0
669   continue
      do 670 i=1,n
         if ((dabs(angz(i+1)-angz(i)).lt.eps).and.
     +       (typ(i+1).eq.2).and.(typ(i).eq.2)) then
            if (lambda(i+1).gt.lambda(i)) dattm(i)=1
            if (lambda(i+1).le.lambda(i)) dattm(i+1)=1
         endif
670   continue
      dattm(1)=dattm(n+1)
      do 675 i=1,start
         angz(i)=angz(i)+(pi*2)
675   continue
      tel=1
      j=1
      if (dattm(start).eq.0) then
         i=start
      else
676      i=start+1
         if (i.eq.n+1) i=1
         if (typ(i).ne.2) goto 676
         if (dattm(i).eq.1) goto 676
      endif
      start=i
      ii=i
680   if (angz(i)-(pi*2).ge.gamma(tel+1)) then
         tel=tel+1
      goto 680
      endif
      do 690 l=1,tel
         gamma(l)=gamma(l)+(pi*2)
690   continue
      starttel=tel
      star(j,1)=x(index(i))
      star(j,2)=y(index(i))
      istarx=1
      if (xdev.gt.eps) then
         xcord=(star(j,1)+tukmed(1))*xdev+xmean
      else
         xcord=star(j,1)+tukmed(1)
      endif
      if (ydev.gt.eps) then
         ycord=(star(j,2)+tukmed(2))*ydev+ymean
      else
         ycord=star(j,2)+tukmed(2)
      endif
         star(j,1)=xcord
         star(j,2)=ycord
700   ii=ii+1
      if (ii.eq.n+1) ii=1
      if (typ(ii).ne.2) goto 700
      if (dattm(ii).eq.1) goto 700
      if (angz(ii).lt.gamma(tel+1)) then
      j=j+1
      star(j,1)=x(index(ii))
      star(j,2)=y(index(ii))
      istarx=1
      if (xdev.gt.eps) then
         xcord=(star(j,1)+tukmed(1))*xdev+xmean
      else
         xcord=star(j,1)+tukmed(1)
      endif
      if (ydev.gt.eps) then
         ycord=(star(j,2)+tukmed(2))*ydev+ymean
      else
         ycord=star(j,2)+tukmed(2)
      endif
         star(j,1)=xcord
         star(j,2)=ycord
      i=ii
      if (i.eq.start) goto 800
      goto 700
      else
730   tel=tel+1
      if (tel.eq.kount+kount1+1) tel=1
      IF (star(j,1).EQ.X(index(ii))) THEN
      ANG=PI2
      ELSE
      ANG=DATAN((star(j,2)-Y(index(ii)))/
     +(star(j,1)-X(index(ii))))
      IF (ANG.LE.0.0) ANG=ANG+PI
      ENDIF
      if (ang.lt.pi2) then
         if (x(index(ii)).lt.star(j,1)) ang=ang+pi
      else
         if (x(index(ii)).gt.star(j,1)) ang=ang+pi
      endif
720   IF (DSIN(ang)*wx(tel)-DCOS(ang)*wy(tel)
     +    .le.dsin(ang)*star(j,1)-dcos(ang)*star(j,2)) THEN
         if ((ii.eq.start).and.(tel.eq.starttel)) 
     +         gamma(tel+1)=gamma(tel+1)+(2*pi)
         if ((angz(ii).ge.gamma(tel+1))) then
            tel=tel+1
      if (tel.eq.kount+kount1+1) tel=1
            goto 720
         endif
      j=j+1
      star(j,1)=x(index(ii))
      star(j,2)=y(index(ii))
      istarx=1
      if (xdev.gt.eps) then
         xcord=(star(j,1)+tukmed(1))*xdev+xmean
      else
         xcord=star(j,1)+tukmed(1)
      endif
      if (ydev.gt.eps) then
         ycord=(star(j,2)+tukmed(2))*ydev+ymean
      else
         ycord=star(j,2)+tukmed(2)
      endif
         star(j,1)=xcord
         star(j,2)=ycord
      i=ii
      if (i.eq.start) goto 800
      goto 700
      ELSE
      j=j+1
740   star(j,1)=wx(tel)
      star(j,2)=wy(tel)
      if ((istarx.eq.1).and.(gamma(tel+1).lt.angz(ii))) then
         IF (x(index(i)).EQ.wx(tel)) THEN
         ANG=PI2
         ELSE
         ANG=DATAN((y(index(i))-wy(tel))/(x(index(i))-wx(tel)))
         IF (ANG.LE.0.0) ANG=ANG+PI
         ENDIF
         if (ang.lt.pi2) then
            if (wx(tel).lt.x(index(i))) ang=ang+pi
         else
            if (wx(tel).gt.x(index(i))) ang=ang+pi
         endif
         IF (DSIN(ang)*wx(tel+1)-DCOS(ang)*wy(tel+1)
     +    .gt.dsin(ang)*wx(tel)-dcos(ang)*wy(tel)) THEN
         tel=tel+1
         goto 740
         ELSE
         ENDIF
      endif
      istarx=0
      if (xdev.gt.eps) then
         xcord=(star(j,1)+tukmed(1))*xdev+xmean
      else
         xcord=star(j,1)+tukmed(1)
      endif
      if (ydev.gt.eps) then
         ycord=(star(j,2)+tukmed(2))*ydev+ymean
      else
         ycord=star(j,2)+tukmed(2)
      endif
         star(j,1)=xcord
         star(j,2)=ycord
      goto 730
      ENDIF
      endif
      else
      nointer=2
      endif
      else
      do 801 i=1,n
      datatyp(i,1)=(x(i)+tukmed(1))*xdev+xmean
      datatyp(i,2)=(y(i)+tukmed(2))*ydev+ymean
      datatyp(i,3)=0.0
801   continue
      endif
800   continue
610   num=kount1+kount   

	
      END
      
      INTEGER FUNCTION NBP_NCEIL(M,J)
      IF (MOD(M,J).EQ.0) THEN
         NBP_NCEIL=INT(dble(M)/J)
      ELSE
         NBP_NCEIL=NINT(dble(M)/J+0.5)
      ENDIF
      RETURN
      END


      SUBROUTINE BP_ISODEPTH(N,M,X,Y,MAXN,MAXM,MAXNUM,NRANK,D,F,BETA,
     +  KAND1,KAND2,ALPHA,IND1,IND2,NCIRQ,MCIRQ,ANGLE,KORNR,L,
     +  JRV,JLV,DPF,NUM,K,EMPTY)
C
C     Computes the depth contour of depth k. This subroutine was described
C     in: Ruts, I. and Rousseeuw, P.J. (1996). Computing depth contours of
C     bivariate point clouds. CSDA 23, 153-168.
C
      INTEGER NCIRQ(N),MCIRQ(N),NRANK(N),F(N)
      integer JLV(M),JRV(M)
      INTEGER IND1(M),IND2(M)
      INTEGER KAND1(MAXNUM),KAND2(MAXNUM),KORNR(MAXNUM,4)
      INTEGER KON,KONTROL,NDATA,NDK,HALT,halt2,jj,JFULL,EMPTY
      INTEGER IV,IW1,IW2,NEXT,JFLAG,KOUNT,NUM,tel
      INTEGER HDEP1,HDEP2,HDEP3,HDEP4,HDEP5,I,J,K,L,M,N
      double precision X(N),Y(N),BETA(N)
      double precision ANGLE(M),D(M),ALPHA(MAXNUM),DPF(N)
      double precision PI,PI2,EPS
      double precision XCORD,YCORD,ANG1,xcord1,ycord1,m1,m2
      PI=DACOS(DBLE(-1.0))
      PI2=PI/2.0
      EPS=0.0000001
      empty=0
C
C   (Re)initialize NCIRQ and NRANK.
C
      DO 45 I=1,N
	 NCIRQ(I)=MCIRQ(I)
45    CONTINUE
      DO 50 I=1,N
	 IV=NCIRQ(I)
	 NRANK(IV)=I
 50   CONTINUE
C
C  Let the line rotate from zero to ANGLE(1).
C
      KOUNT=1
      HALT=0
      if (angle(1).gt.pi2) then
         l=1
	 CALL BP_ADJUST(IND1,IND2,L,NRANK,NCIRQ,KOUNT,ALPHA,ANGLE,
     +	                  K,N,M,MAXNUM,KAND1,KAND2,D,X,Y)
         halt=1
      endif
      L=2
 60   KONTROL=0
      IF ((PI.LE.(ANGLE(L)+PI2)).AND.((ANGLE(L)-PI2).LT.ANGLE(1))) THEN
	 CALL BP_ADJUST(IND1,IND2,L,NRANK,NCIRQ,KOUNT,ALPHA,ANGLE,
     +	                  K,N,M,MAXNUM,KAND1,KAND2,D,X,Y)
	 KONTROL=1
      ENDIF
      L=L+1
      IF (KONTROL.EQ.1) HALT=1
      IF ((L.EQ.M+1).AND.(KONTROL.EQ.1)) THEN
	 JFLAG=1
	 GOTO 79
      ENDIF
      IF (((HALT.EQ.1).AND.(KONTROL.EQ.0)).OR.(L.EQ.M+1)) THEN
	 GOTO 70
      ELSE
	 GOTO 60
      ENDIF
 70   if (l.gt.1) then
         JFLAG=L-1
      else
         jflag=m
      endif
      J=0
C
C  In case the first switch didn't occur between zero and ANGLE(1),
C  look for it between the following angles.
C
      IF ((L.EQ.M+1).AND.(KONTROL.EQ.0)) THEN
	 HALT=0
         halt2=0
 73      J=J+1
         if (j.eq.m+1) j=1
	 L=J+1
         if (l.eq.m+1) l=1
 75      KONTROL=0
	 IF ((ANGLE(L)+PI2).LT.PI) THEN
	    ANG1=ANGLE(L)+PI2
         ELSE
	    ANG1=ANGLE(L)-PI2
         ENDIF
         if (j.eq.m) then
            jj=1
            if (halt2.eq.0) angle(1)=angle(1)+pi
         else
            jj=j+1
         endif
	 IF ((ANGLE(J).LE.ANG1).AND.(ANG1.LT.ANGLE(jj))) THEN
            if (angle(1).gt.pi) angle(1)=angle(1)-pi
	    CALL BP_ADJUST(IND1,IND2,L,NRANK,NCIRQ,KOUNT,ALPHA,ANGLE,
     +                       K,N,M,MAXNUM,KAND1,KAND2,D,X,Y)
	    KONTROL=1
         ENDIF
         if (angle(1).gt.pi) angle(1)=angle(1)-pi
	 IF (L.NE.M) THEN
	    L=L+1
         ELSE
	    L=1
         ENDIF
	 IF (KONTROL.EQ.1) HALT=1
	 IF ((HALT.EQ.1).AND.(KONTROL.EQ.0)) THEN
            if (halt2.eq.1) goto 101
            if (l.gt.1) then
               jflag=l-1
            else
               jflag=m
            endif
	    GOTO 79
         ELSE
            IF (L.EQ.jj) THEN
               if (jj.eq.1) halt2=1
               GOTO 73
            ELSE
	       GOTO 75
            ENDIF
         ENDIF
      ENDIF
C
C  The first switch has occurred. Now start looking for the next ones,
C  between the following angles.
C
79    DO 80 I=J+1,M-1
	 L=JFLAG
 90      KONTROL=0
	 IF ((ANGLE(L)+PI2).LT.PI) THEN
	    ANG1=ANGLE(L)+PI2
         ELSE
	    ANG1=ANGLE(L)-PI2
         ENDIF
	 IF ((ANGLE(I).LE.ANG1).AND.(ANG1.LT.ANGLE(I+1))) THEN
	    CALL BP_ADJUST(IND1,IND2,L,NRANK,NCIRQ,KOUNT,ALPHA,
     +                  ANGLE,K,N,M,MAXNUM,KAND1,KAND2,D,X,Y)
	    KONTROL=1
         ENDIF
	 IF (KONTROL.EQ.0) THEN
	    JFLAG=L
         ELSE
	    IF (L.NE.M) THEN
	       L=L+1
            ELSE
	       L=1
            ENDIF
	    GOTO 90
         ENDIF
 80   CONTINUE
      L=JFLAG
C
C  Finally, look for necessary switches between the last angle and zero.
C
100   KONTROL=0
      IF ((ANGLE(L)+PI2).LT.PI) THEN
	 ANG1=ANGLE(L)+PI2
      ELSE
	 ANG1=ANGLE(L)-PI2
      ENDIF
      IF ((ANGLE(M).LE.ANG1).AND.(ANG1.LT.PI)) THEN
	 CALL BP_ADJUST(IND1,IND2,L,NRANK,NCIRQ,KOUNT,ALPHA,
     +               ANGLE,K,N,M,MAXNUM,KAND1,KAND2,D,X,Y)
	 KONTROL=1
      ENDIF
      IF (KONTROL.EQ.1) THEN
         IF (L.NE.M) THEN
	     L=L+1
         ELSE
	     L=1
         ENDIF
         GOTO 100
      ENDIF 
101      NUM=KOUNT-1
C  
C  Sort the NUM special k-dividers. 
C  Permute KAND1, KAND2 and D in the same way.
C
      CALL BP_SORT(ALPHA,KAND1,KAND2,D,NUM,JLV,JRV)
      IW1=1
      IW2=2
      JFULL=0
      NDK=0
      tel=0

120   NDATA=0
C
C  Compute the intersection point.
C
      IF (DABS(-DSIN(ALPHA(IW2))*DCOS(ALPHA(IW1))
     +         +DSIN(ALPHA(IW1))*DCOS(ALPHA(IW2))).LT.EPS) THEN
	 IW2=IW2+1
	 IF (IW2.EQ.NUM+1) IW2=1
	 GOTO 120
      ENDIF
      XCORD=(DCOS(ALPHA(IW2))*D(IW1)-DCOS(ALPHA(IW1))*D(IW2))
     + /(-DSIN(ALPHA(IW2))*DCOS(ALPHA(IW1))
     +                   +DSIN(ALPHA(IW1))*DCOS(ALPHA(IW2)))
      YCORD=(-DSIN(ALPHA(IW2))*D(IW1)+DSIN(ALPHA(IW1))*D(IW2))
     + /(-DSIN(ALPHA(IW1))*DCOS(ALPHA(IW2))
     +                   +DSIN(ALPHA(IW2))*DCOS(ALPHA(IW1)))
C 
C  Test whether the intersection point is a data point. 
C  If so, adjust IW1 and IW2.
C
      IF ((KAND1(IW1).EQ.KAND1(IW2)).OR.(KAND1(IW1).EQ.KAND2(IW2)))
     +     NDATA=KAND1(IW1)
      IF ((KAND2(IW1).EQ.KAND1(IW2)).OR.(KAND2(IW1).EQ.KAND2(IW2)))
     +     NDATA=KAND2(IW1)
      IF (NDATA.NE.0) THEN
         iv=0
 125     NEXT=IW2+1
         iv=iv+1
         IF (NEXT.EQ.(NUM+1)) NEXT=1
         if (next.ne.iw1) then
            IF ((NDATA.EQ.KAND1(NEXT)).OR.(NDATA.EQ.KAND2(NEXT))) THEN
               IW2=IW2+1
               IF (IW2.EQ.(NUM+1)) IW2=1
               GOTO 125
            ENDIF
         endif
         if (iv.eq.(num-1)) then
            num=1
            KORNR(1,1)=KAND1(IW1)
            KORNR(1,2)=KAND2(IW1)
            KORNR(1,3)=KAND1(IW2)
            KORNR(1,4)=KAND2(IW2)
            return
         endif
      ENDIF
      IF (IW2.EQ.NUM) THEN
         KON=1
      ELSE
         KON=IW2+1
      ENDIF
      if (kon.eq.iw1) kon=kon+1
      if (kon.eq.num+1) kon=1
C
C  Test whether the intersection point lies to the left of the special 
C  k-divider which corresponds to ALPHA(KON). If so, compute its depth.
C
      IF ((DSIN(ALPHA(KON))*XCORD-DCOS(ALPHA(KON))*YCORD
     +     -D(KON)).le.eps) THEN
         
         CALL BP_DEPTH(XCORD,YCORD,N,X,Y,BETA,F,DPF,JLV,JRV,HDEP1)
         
         IF (HDEP1.EQ.K) NDK=1
         IF (HDEP1.NE.K) THEN
         CALL BP_DEPTH(XCORD-EPS*10,YCORD-EPS*10,N,X,Y,BETA,F,DPF,
     +           JLV,JRV,HDEP2)
         CALL BP_DEPTH(XCORD+EPS*10,YCORD+EPS*10,N,X,Y,BETA,F,DPF,
     +        JLV,JRV,HDEP3)
         CALL BP_DEPTH(XCORD-EPS*10,YCORD+EPS*10,N,X,Y,BETA,F,DPF,
     +        JLV,JRV,HDEP4)
         CALL BP_DEPTH(XCORD+EPS*10,YCORD-EPS*10,N,X,Y,BETA,F,DPF,
     +        JLV,JRV,HDEP5)
         IF ((NDK.EQ.0).AND.
     +        ((HDEP1.ge.K).OR.(HDEP2.ge.K).OR.(HDEP3.ge.K)
     +        .OR.(HDEP4.ge.K).OR.(HDEP5.ge.K))) THEN 
            NDK=1
         ENDIF
         IF ((HDEP1.LT.K).AND.(HDEP2.LT.K)
     +        .AND.(HDEP3.LT.K).AND.(HDEP4.LT.K)
     +        .AND.(HDEP5.LT.K).AND.(NDK.EQ.1)) THEN
C     
C  The intersection point is not the correct one, 
C  try the next special k-divider.
C
            IW2=IW2+1
            IF (IW2.EQ.(NUM+1)) IW2=1
            GOTO 120
         ENDIF
      ENDIF
C
C  Store IW1 and IW2 in KORNR. If KORNR has already been filled, check whether 
C  we have encountered this intersection point before.
C
      IF ((IW2.GT.IW1).AND.(JFULL.EQ.0)) THEN
         DO 130 I=IW1,IW2-1
            KORNR(I,1)=KAND1(IW1)
            KORNR(I,2)=KAND2(IW1)
            KORNR(I,3)=KAND1(IW2)
            KORNR(I,4)=KAND2(IW2)
 130     CONTINUE
      ELSE
         IF (IW2.GT.IW1) THEN
            DO 140 I=IW1,IW2-1
               IF ((KORNR(I,1).EQ.KAND1(IW1)).AND.
     +              (KORNR(I,2).EQ.KAND2(IW1)).AND.
     +              (KORNR(I,3).EQ.KAND1(IW2)).AND.
     +              (KORNR(I,4).EQ.KAND2(IW2)))
     +              THEN
		  GOTO 170
               ELSE
                  tel=tel+1
                  if (tel.gt.num*num) then
                     ndk=1
                     goto 170
                  endif
                  m1=(y(kornr(i,2))-y(kornr(i,1)))/
     +                 (x(kornr(i,2))-x(kornr(i,1)))
                  m2=(y(kornr(i,4))-y(kornr(i,3)))/
     +                 (x(kornr(i,4))-x(kornr(i,3)))
                  if (m1.ne.m2) then
                     xcord1=(m1*x(kornr(i,1))-y(kornr(i,1))-
     +                    m2*x(kornr(i,3))-y(kornr(i,3)))/(m1-m2)
                     ycord1=(m2*(m1*x(kornr(i,1))-y(kornr(i,1)))-
     +                    m1*(m2*x(kornr(i,3))-y(kornr(i,3))))/(m1-m2)
                  endif
                  if ((dabs(xcord1-xcord).le.eps).and.
     +                 (dabs(ycord1-ycord).le.eps)) then
                     goto 170
                  endif
                  
		  KORNR(I,1)=KAND1(IW1)
		  KORNR(I,2)=KAND2(IW1)
		  KORNR(I,3)=KAND1(IW2)
		  KORNR(I,4)=KAND2(IW2)
               ENDIF
 140        CONTINUE
         ELSE
            JFULL=1
            DO 150 I=IW1,NUM
               KORNR(I,1)=KAND1(IW1)
               KORNR(I,2)=KAND2(IW1)
               KORNR(I,3)=KAND1(IW2)
               KORNR(I,4)=KAND2(IW2)
 150        CONTINUE
            DO 160 I=1,IW2-1
               IF ((KORNR(I,1).EQ.KAND1(IW1)).AND.
     +              (KORNR(I,2).EQ.KAND2(IW1)).AND.
     +              (KORNR(I,3).EQ.KAND1(IW2)).AND.
     +              (KORNR(I,4).EQ.KAND2(IW2)))
     +              THEN
		  GOTO 170
               ELSE
                  tel=tel+1
                  if (tel.gt.num*num) then
                     ndk=1
                     goto 170
                  endif
                  m1=(y(kornr(i,2))-y(kornr(i,1)))/
     +                 (x(kornr(i,2))-x(kornr(i,1)))
                  m2=(y(kornr(i,4))-y(kornr(i,3)))/
     +                 (x(kornr(i,4))-x(kornr(i,3)))
                  if (m1.ne.m2) then
                     xcord1=(m1*x(kornr(i,1))-y(kornr(i,1))-
     +                    m2*x(kornr(i,3))-y(kornr(i,3)))/(m1-m2)
                     ycord1=(m2*(m1*x(kornr(i,1))-y(kornr(i,1)))-
     +                    m1*(m2*x(kornr(i,3))-y(kornr(i,3))))/(m1-m2)
                  endif
                  if ((dabs(xcord1-xcord).le.eps).and.
     +                 (dabs(ycord1-ycord).le.eps)) then
                     goto 170
                  endif
                  
		  KORNR(I,1)=KAND1(IW1)
		  KORNR(I,2)=KAND2(IW1)
		  KORNR(I,3)=KAND1(IW2)
		  KORNR(I,4)=KAND2(IW2)
               ENDIF
 160        CONTINUE
         ENDIF
      ENDIF
      ELSE
C     
C  The intersection point is not the correct one, 
C  try the next special k-divider.
C
	 IW2=IW2+1
	 IF (IW2.EQ.(NUM+1)) IW2=1
	 GOTO 120
      ENDIF
C
C  Look for the next vertex of the convex figure.
C
      IW1=IW2
      IW2=IW2+1
      IF (IW2.EQ.(NUM+1)) IW2=1
      GOTO 120
170   if (ndk.eq.0) empty=1
      RETURN
      END

      SUBROUTINE BP_SORT(B,I1,I2,R,N,JLV,JRV)
C
C  Sorts a double precision  array B of length N and permutes two integer 
C  arrays I1 and I2 and one double precision array R in the same way.
C
      INTEGER N,I1(N),I2(N),H1,H2
      double precision B(N),XX,AMM
      double precision R(N),H3
      integer JLV(N),JRV(N)
      JSS=1
      JLV(1)=1
      JRV(1)=N
 10   JNDL=JLV(JSS)
      JR=JRV(JSS)
      JSS=JSS-1
 20   JNC=JNDL
      J=JR
      JTWE=(JNDL+JR)/2
      XX=B(JTWE)
 30   IF (B(JNC).GE.XX) GOTO 40
      JNC=JNC+1
      GOTO 30
 40   IF (XX.GE.B(J)) GOTO 50
      J=J-1
      GOTO 40
 50   IF (JNC.GT.J) GOTO 60
      AMM=B(JNC)
      H1=I1(JNC)
      H2=I2(JNC)
      H3=R(JNC) 
      B(JNC)=B(J)
      I1(JNC)=I1(J)
      I2(JNC)=I2(J)
      R(JNC)=R(J)
      B(J)=AMM
      I1(J)=H1
      I2(J)=H2
      R(J)=H3
      JNC=JNC+1
      J=J-1
 60   IF (JNC.LE.J) GOTO 30
      IF ((J-JNDL).LT.(JR-JNC)) GOTO 80
      IF (JNDL.GE.J) GOTO 70
      JSS=JSS+1
      JLV(JSS)=JNDL
      JRV(JSS)=J
 70   JNDL=JNC
      GOTO 100
 80   IF (JNC.GE.JR) GOTO 90
      JSS=JSS+1
      JLV(JSS)=JNC
      JRV(JSS)=JR
 90   JR=J
100   IF (JNDL.LT.JR) GOTO 20
      IF (JSS.NE.0) GOTO 10
      RETURN
      END



      SUBROUTINE BP_ADJUST(IND1,IND2,L,NRANK,NCIRQ,KOUNT,
     +           ALPHA,ANGLE,K,N,M,MAXNUM,KAND1,KAND2,D,X,Y)
C
C  Updates NCIRQ and NRANK, detects the special k-dividers and stores 
C  their angles and the constant terms of their equations.
C
      INTEGER NCIRQ(N),NRANK(N),IND1(M),IND2(M)
      INTEGER KAND1(MAXNUM),KAND2(MAXNUM)
      INTEGER KOUNT,K,L,N,IV,IV1,IV2,D1,D2
      double precision X(N),Y(N),ANGLE(M),D(M)
      double precision ALPHA(MAXNUM),DUM,PI,PI2
      PI=DACOS(DBLE(-1.0))
      PI2=PI/2.0
      D1=IND1(L)
      IV1=NRANK(D1)
      D2=IND2(L)
      IV2=NRANK(D2)
      IV=NCIRQ(IV1)
      NCIRQ(IV1)=NCIRQ(IV2)
      NCIRQ(IV2)=IV
      IV=IV1
      NRANK(D1)=IV2
      NRANK(D2)=IV
	 IF (((IV1.EQ.K).AND.(IV2.EQ.(K+1)))
     +      .OR.((IV2.EQ.K).AND.(IV1.EQ.(K+1)))
     +      .OR.((IV1.EQ.(N-K)).AND.(IV2.EQ.(N-K+1))) 
     +      .OR.((IV2.EQ.(N-K)).AND.(IV1.EQ.(N-K+1)))) THEN
	    IF (ANGLE(L).LT.PI2) THEN
	       DUM=ANGLE(L)+PI2
            ELSE
	       DUM=ANGLE(L)-PI2
            ENDIF
            IF (((IV1.EQ.K).AND.(IV2.EQ.(K+1)))
     +         .OR.((IV2.EQ.K).AND.(IV1.EQ.(K+1)))) THEN
	       IF (DUM.LE.PI2) THEN
		  ALPHA(KOUNT)=ANGLE(L)+PI
               ELSE
		  ALPHA(KOUNT)=ANGLE(L)
               ENDIF
            ENDIF
	    IF (((IV1.EQ.(N-K)).AND.(IV2.EQ.(N-K+1)))
     +        .OR.((IV2.EQ.(N-K)).AND.(IV1.EQ.(N-K+1)))) THEN
	       IF (DUM.LE.PI2) THEN
		  ALPHA(KOUNT)=ANGLE(L)
               ELSE
		  ALPHA(KOUNT)=ANGLE(L)+PI
               ENDIF
            ENDIF
	    KAND1(KOUNT)=IND1(L)
	    KAND2(KOUNT)=IND2(L)
	    D(KOUNT)=DSIN(ALPHA(KOUNT))*X(IND1(L))
     +                -DCOS(ALPHA(KOUNT))*Y(IND1(L))
            KOUNT=KOUNT+1
         ENDIF
      RETURN
      END



      SUBROUTINE BP_DEPTH(U,V,N,X,Y,BETA,F,DPF,JLV,JRV,HDEP)
C
C  Computes the halfspace depth of a point. This subroutine was described
C  in: Rousseeuw, P.J. and Ruts, I. (1996). Algorithm AS 307: Bivariate 
C  location depth. Applied Statistics (JRSS-C) 45, 516-526.
C
      double precision U,V,BETA(N),X(N),Y(N),DPF(N)
      double precision P,P2,EPSI,D,XU,YU,ANG,ALPHK,BETAK
      INTEGER F(N),GI,HDEP
      integer JLV(N),JRV(N)
      NUMH=0
      HDEP=0
      IF (N.LT.1) RETURN
      P=DACOS(DBLE(-1.0))
      P2=P*2.0
      EPSI=0.000001
      NZ=0
C
C  Construct the array BETA.
C

      DO 10 I=1,N
          D=DSQRT((X(I)-U)*(X(I)-U)+(Y(I)-V)*(Y(I)-V))
          IF (D.LE.EPSI) THEN
              NZ=NZ+1
          ELSE
              XU=(X(I)-U)/D
              YU=(Y(I)-V)/D
              IF (DABS(XU).GT.DABS(YU)) THEN
                  IF (X(I).GE.U) THEN
                      BETA(I-NZ)=DASIN(YU)
                      IF(BETA(I-NZ).LT.0.0) THEN
                          BETA(I-NZ)=P2+BETA(I-NZ)
                      ENDIF
                  ELSE
                      BETA(I-NZ)=P-DASIN(YU)
                  ENDIF
              ELSE
                  IF (Y(I).GE.V) THEN
                      BETA(I-NZ)=DACOS(XU)
                  ELSE
                      BETA(I-NZ)=P2-DACOS(XU)
                  ENDIF
              ENDIF
              IF (BETA(I-NZ).GE.(P2-EPSI)) BETA(I-NZ)=0.0
          ENDIF
  10  CONTINUE
      NN=N-NZ
      IF (NN.LE.1) GOTO 60
C
C  Sort the array BETA.
C
      DO 15 I=1,NN
      DPF(I)=DBLE(F(I))
15    CONTINUE
      CALL BP_SORT(BETA,F,F,DPF,NN,JLV,JRV)
C
C  Check whether Z=(U,V) lies outside the data cloud.
C
      ANG=BETA(1)-BETA(NN)+P2
      DO 20 I=2,NN
          ANG=DMAX1(ANG,(BETA(I)-BETA(I-1)))
  20  CONTINUE
      IF (ANG.GT.(P+EPSI)) GOTO 60
C
C  Make smallest BETA equal to zero,
C  and compute NU = number of BETA < PI.
C
      ANG=BETA(1)
      NU=0
      DO 30 I=1,NN
          BETA(I)=BETA(I)-ANG
          IF (BETA(I).LT.(P-EPSI)) NU=NU+1
  30  CONTINUE
      IF (NU.GE.NN) GOTO 60
C
C  Mergesort the BETA with their antipodal angles,
C  and at the same time update I, F(I), and NBAD.
C
      JA=1
      JB=1
      ALPHK=BETA(1)
      BETAK=BETA(NU+1)-P
      NN2=NN*2
      NBAD=0
      I=NU
      NF=NN
      DO 40 J=1,NN2
          IF ((ALPHK+EPSI).LT.BETAK) THEN
              NF=NF+1
              IF (JA.LT.NN) THEN
                  JA=JA+1
                  ALPHK=BETA(JA)
              ELSE
                  ALPHK=P2+1.0
              ENDIF
          ELSE
              I=I+1
              IF (I.EQ.(NN+1)) THEN
                  I=1
                  NF=NF-NN
              ENDIF
              F(I)=NF
              NBAD=NBAD+NBP_K((NF-I),2)
              IF (JB.LT.NN) THEN
                  JB=JB+1
                  IF ((JB+NU).LE.NN) THEN
                      BETAK=BETA(JB+NU)-P
                  ELSE
                      BETAK=BETA(JB+NU-NN)+P
                  ENDIF
              ELSE
                  BETAK=P2+1.0
              ENDIF
          ENDIF
  40  CONTINUE
C
C  Computation of NUMH for halfspace depth.
C
      GI=0
      JA=1
      ANG=BETA(1)
      NUMH=MIN0(F(1),(NN-F(1)))
      DO 50 I=2,NN
          IF(BETA(I).LE.(ANG+EPSI)) THEN
              JA=JA+1
          ELSE
              GI=GI+JA
              JA=1
              ANG=BETA(I)
          ENDIF
          KI=F(I)-GI
          NUMH=MIN0(NUMH,MIN0(KI,(NN-KI)))
   50 CONTINUE
C
C  Adjust for the number NZ of data points equal to Z=(U,V).
C
   60 NUMH=NUMH+NZ
      HDEP=NUMH
      RETURN
      END


      INTEGER FUNCTION NBP_K(M,J)
      IF (M.LT.J) THEN
          NBP_K=0
      ELSE
          IF (J.EQ.1) NBP_K=M
          IF (J.EQ.2) NBP_K=(M*(M-1))/2
          IF (J.EQ.3) NBP_K=(M*(M-1)*(M-2))/6
      ENDIF
      RETURN
      END

      subroutine bp_rdraw(a,ntot,seed,n)
C
C     Draws n elements out of a dataset of size ntot, such that
C     the selected case numbers are uniformly distributed from 1 to ntot.
C	  Where bp_uniran(1, ...) has been changed to bp_uniran(2, ...)
C     and	int(urand*... to int(urand(1)*...

	integer a(n)
	integer seed,nrand
        double precision urand(2)
	jndex=0
	do 20 m=1,n
          call bp_uniran(2,seed,urand)
	  nrand=int(urand(2)*(ntot-jndex))+1 
	  jndex=jndex+1
	  if(jndex.eq.1) then
	    a(jndex)=nrand
	  else
	    a(jndex)=nrand+jndex-1
	    do 5 i=1,jndex-1
	      if(a(i).gt.nrand+i-1) then
	        do 6 j=jndex,i+1,-1
	          a(j)=a(j-1)
 6              continue
	        a(i)=nrand+i-1
	        goto 20
	      endif
 5          continue
	  endif
 20     continue
	return
	end

      subroutine nbp_norrandp(n,iseed,xrand)
C
C     This subroutine generates a random sample of size n
C     from the normal (Gaussian) distribution with mean = 0 and 
C     standard deviation = 1.
C     The generated random sample will be placed in the vector xrand.
C
      double precision xrand(n),u1,u2,arg1,arg2,sqrt1,z1,z2
      double precision yrand(2)
C
C     Generate n uniform (0,1) random numbers;
C     then generate 2 additional uniform (0,1) random numbers.
C
      ipr=6
      data pi/3.14159265359/
      call bp_uniran(n,iseed,xrand)
      call bp_uniran(2,iseed,yrand)
C
C     Generate n normal random numbers using the Box-Muller method.
C
      DO 9200 I=1,N,2
      IP1=I+1
      U1=xrand(I)
      IF(I.EQ.N) GOTO 9210
      U2=xrand(IP1)
      GOTO 9220
9210  U2=yrand(2)
9220  ARG1=-2.0*dLOG(U1)
      ARG2=2.0*PI*U2
      SQRT1=dSQRT(ARG1)
      Z1=SQRT1*dCOS(ARG2)
      Z2=SQRT1*dSIN(ARG2)
      xrand(I)=Z1
      IF(I.EQ.N) GOTO 9200
      xrand(IP1)=Z2
9200  CONTINUE
      RETURN
      END

      subroutine bp_uniran(n,iseed,xrand)
C
C     This subroutine generates a random sample of size n from the 
C     uniform (rectangular) distribution on the unit interval (0,1).
C     The generated random sample will be placed in the vector xrand.
C
      double precision xrand(n)
      integer M(17)
      save i,j,m,m1,m2
      DATA M(1),M(2),M(3),M(4),M(5),M(6),M(7),M(8),M(9),M(10),M(11),
     1     M(12),M(13),M(14),M(15),M(16),M(17)
     1/    30788,23052,2053,19346,10646,19427,23975,
     1     19049,10949,19693,29746,26748,2796,23890,
     1     29168,31924,16499/ 
      DATA M1,M2,I,J / 32767,256,5,17 / 
      IPR=6
      IF (ISEED.LE.0) GOTO 9290
      MDIG=32
      M1=2**(MDIG-2)+(2**(MDIG-2)-1)
      M2=2**(MDIG/2)
      ISEED3=IABS(ISEED)
      IF(M1.LT.IABS(ISEED))ISEED3=M1
      IF(MOD(ISEED3,2).EQ.0)ISEED3=ISEED3-1
      K0=MOD(9069,M2)
      K1=9069/M2
      J0=MOD(ISEED3,M2)
      J1=ISEED3/M2
      DO 9200 I=1,17
      ISEED3=J0*K0
      J1=MOD(ISEED3/M2+J0*K1+J1*K0,M2/2)
      J0=MOD(ISEED3,M2)
      M(I)=J0+M2*J1 
9200  CONTINUE
      I=5 
      J=17
9290  CONTINUE
C
C     Generate the n random numbers.  
C
      DO 9300 L=1,N
      K=M(I)-M(J)
      IF(K.LT.0)K=K+M1
      M(J)=K
      I=I-1
      IF(I.EQ.0)I=17
      J=J-1
      IF(J.EQ.0)J=17
      AK=K
      AM1=M1
      xrand(L)=dble(AK/AM1)
9300  CONTINUE
      ISEED=(-1)
9000  CONTINUE
      RETURN
      END 

