      module nr_mod
      contains

      FUNCTION gammq(a,x)
      REAL(8) a,gammq,x
!     USES gcf,gser
      REAL(8) gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.) then 
         write(*,*) 'bad arguments in gammq'
         stop
      endif
      if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammq=1.-gamser
      else
        call gcf(gammcf,a,x,gln)
        gammq=gammcf
      endif
      return
      END

      SUBROUTINE gser(gamser,a,x,gln)
      INTEGER ITMAX
      REAL(8) a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.e-7)
!     USES gammln
      INTEGER n
      REAL(8) ap,del,sum,gammln
      gln=gammln(a)
      if(x.le.0.)then
        if(x.lt.0.) then 
           write(*,*) 'x < 0 in gser'
        endif
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
 11         continue
      write(*,*) 'a too large, ITMAX too small in gser'
      stop
 1         gamser=sum*exp(-x+a*log(x)-gln)
      return
      END


      SUBROUTINE gcf(gammcf,a,x,gln)
      INTEGER ITMAX
      REAL(8) a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
!     USES gammln
      INTEGER i
      REAL(8) an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1.-a
      c=1./FPMIN
      d=1./b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS)goto 1
 11         continue
            write(*,*) 'a too large, ITMAX too small in gcf'
            stop
 1         gammcf=exp(-x+a*log(x)-gln)*h
      return
      END



      SUBROUTINE FRPRMN(p,n,ftol,iter,fret)
      use nr_common_data
      INTEGER iter,n,NMAX,ITMAX
      REAL(8) fret,ftol,p(n),EPS,func
      EXTERNAL func
      PARAMETER (NMAX=60000,ITMAX=16000000,EPS=1.d-10)
!     USES dfunc,func,linmin
      INTEGER its,j,nsteps
      REAL(8) dgg,fp,gam,gg,g(NMAX),h(NMAX),xi(NMAX)


      nsteps = 0

      fp=func(p)
      call dfunc(p,xi)
      do 11 j=1,n
        g(j)=-xi(j)
        h(j)=g(j)
        xi(j)=h(j)
 11         continue
      do 14 its=1,ITMAX
        iter=its
        call linmin(p,xi,n,fret)
        if(2.*abs(fret-fp).le.ftol*(abs(fret)+abs(fp)+EPS))return

        nsteps = nsteps + 1
        fp=func(p)
        call dfunc(p,xi)

        if(mod(nsteps,20) == 0) then
           write(*,*) nsteps,fp
!           write(90,*) nsteps,fp
!           call flush(90)
!           call print_coordinates(89)
!           call print_rigid_coordinates(88)
        endif
        if(mod(nsteps,nprintdata).eq.0) then 
           call print_data
        endif

        gg=0.
        dgg=0.
        do 12 j=1,n
          gg=gg+g(j)**2
!         dgg=dgg+xi(j)**2
          dgg=dgg+(xi(j)+g(j))*xi(j)
 12             continue
        if(gg.eq.0.)return
        gam=dgg/gg
        do 13 j=1,n
          g(j)=-xi(j)
          h(j)=g(j)+gam*h(j)
          xi(j)=h(j)
 13             continue
 14                 continue
      write(*,*) 'frprmn maximum iterations exceeded'
      stop
      return
      END



      SUBROUTINE MNBRAK(ax,bx,cx,fa,fb,fc,func)
      REAL(8) ax,bx,cx,fa,fb,fc,func,GOLD,GLIMIT,TINY
      EXTERNAL func
      PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.d-20)

      REAL(8) dum,fu,q,r,u,ulim

      fa=func(ax)
      fb=func(bx)

       if(fb.gt.fa)then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
      cx=bx+GOLD*(bx-ax)
      fc=func(cx)
 1         if(fb.ge.fc)then
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),TINY),q-r))
        ulim=bx+GLIMIT*(cx-bx)
        if((bx-u)*(u-cx).gt.0.)then
          fu=func(u)
          if(fu.lt.fc)then
            ax=bx
            fa=fb
            bx=u
            fb=fu
            return
          else if(fu.gt.fb)then
            cx=u
            fc=fu
            return
          endif
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        else if((cx-u)*(u-ulim).gt.0.)then
          fu=func(u)
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
            fu=func(u)
          endif
        else if((u-ulim)*(ulim-cx).ge.0.)then
          u=ulim
          fu=func(u)
        else
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        endif
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        goto 1
      endif
      return
      END


      SUBROUTINE LINMIN(p,xi,n,fret)
      use nr_common_data
      INTEGER n,NMAX
      REAL(8) fret,p(n),xi(n),TOL
      PARAMETER (NMAX=60000,TOL=1.d-4)
!     USES brent,f1dim,mnbrak
      INTEGER j,ncom
      REAL(8) ax,bx,fa,fb,fx,xmin,xx,pcom(NMAX),xicom(NMAX),brent,dbrent
      COMMON /f1com/ pcom,xicom,ncom
      EXTERNAL f1dim,df1dim
      REAL(8) f1dim,df1dim

      ncom=n
      do 11 j=1,n
        pcom(j)=p(j)
        xicom(j)=xi(j)
 11         continue
      ax=0.
!      xx=1.
!------------------------------------------------
!     NOTE CHANGED SO THAT IT DOESN'T TAKE BIG STEPS
!--------------------------------------------------
!      xx = 1.e-9
!      xx = 1.0d0 

      xx = linmin_param

      call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
      fret=dbrent(ax,xx,bx,f1dim,df1dim,TOL,xmin)

      do 12 j=1,n
        xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
 12         continue
      return
      END


      subroutine sincos(angle,sine,cosine)
      use nr_common_data
      implicit none
      integer ix,xsign
      integer ssign,csign
      real*8 angle,sine,cosine
      real*8 x,y,z,sx,cx,sz,cz
      real*8 table(0:91)
      data table /&
     &   0.00000000000000000000d+0, 1.74524064372835128194d-2,&
     &   3.48994967025009716460d-2, 5.23359562429438327221d-2,&
     &   6.97564737441253007760d-2, 8.71557427476581735581d-2,&
     &   1.04528463267653471400d-1, 1.21869343405147481113d-1,&
     &   1.39173100960065444112d-1, 1.56434465040230869010d-1,&
     &   1.73648177666930348852d-1, 1.90808995376544812405d-1,&
     &   2.07911690817759337102d-1, 2.24951054343864998051d-1,&
     &   2.41921895599667722560d-1, 2.58819045102520762349d-1,&
     &   2.75637355816999185650d-1, 2.92371704722736728097d-1,&
     &   3.09016994374947424102d-1, 3.25568154457156668714d-1,&
     &   3.42020143325668733044d-1, 3.58367949545300273484d-1,&
     &   3.74606593415912035415d-1, 3.90731128489273755062d-1,&
     &   4.06736643075800207754d-1, 4.22618261740699436187d-1,&
     &   4.38371146789077417453d-1, 4.53990499739546791560d-1,&
     &   4.69471562785890775959d-1, 4.84809620246337029075d-1,&
     &   5.00000000000000000000d-1, 5.15038074910054210082d-1,&
     &   5.29919264233204954047d-1, 5.44639035015027082224d-1,&
     &   5.59192903470746830160d-1, 5.73576436351046096108d-1,&
     &   5.87785252292473129169d-1, 6.01815023152048279918d-1,&
     &   6.15661475325658279669d-1, 6.29320391049837452706d-1,&
     &   6.42787609686539326323d-1, 6.56059028990507284782d-1,&
     &   6.69130606358858213826d-1, 6.81998360062498500442d-1,&
     &   6.94658370458997286656d-1, 7.07106781186547524401d-1,&
     &   7.19339800338651139356d-1, 7.31353701619170483288d-1,&
     &   7.43144825477394235015d-1, 7.54709580222771997943d-1,&
     &   7.66044443118978035202d-1, 7.77145961456970879980d-1,&
     &   7.88010753606721956694d-1, 7.98635510047292846284d-1,&
     &   8.09016994374947424102d-1, 8.19152044288991789684d-1,&
     &   8.29037572555041692006d-1, 8.38670567945424029638d-1,&
     &   8.48048096156425970386d-1, 8.57167300702112287465d-1,&
     &   8.66025403784438646764d-1, 8.74619707139395800285d-1,&
     &   8.82947592858926942032d-1, 8.91006524188367862360d-1,&
     &   8.98794046299166992782d-1, 9.06307787036649963243d-1,&
     &   9.13545457642600895502d-1, 9.20504853452440327397d-1,&
     &   9.27183854566787400806d-1, 9.33580426497201748990d-1,&
     &   9.39692620785908384054d-1, 9.45518575599316810348d-1,&
     &   9.51056516295153572116d-1, 9.56304755963035481339d-1,&
     &   9.61261695938318861916d-1, 9.65925826289068286750d-1,&
     &   9.70295726275996472306d-1, 9.74370064785235228540d-1,&
     &   9.78147600733805637929d-1, 9.81627183447663953497d-1,&
     &   9.84807753012208059367d-1, 9.87688340595137726190d-1,&
     &   9.90268068741570315084d-1, 9.92546151641322034980d-1,&
     &   9.94521895368273336923d-1, 9.96194698091745532295d-1,&
     &   9.97564050259824247613d-1, 9.98629534754573873784d-1,&
     &   9.99390827019095730006d-1, 9.99847695156391239157d-1,&
     &   1.00000000000000000000d+0, 9.99847695156391239157d-1 /

      if(abs(angle).gt.1.d+11) then 
!     subroutine crashes if angle > 1.d+11, so if this happens, just calculate the 
!     functions the old way, and return
         sine = sin(angle/pifac)
         cosine = cos(angle/pifac)
         return
      endif


!
!
!     get the angle value in the range from 0 to 360 degrees
!
      x = angle
      xsign = 1
      if (x .lt. 0.0d0) then
          xsign = -1
          x = -x
      end if
      x = x - 360.0d0*int(x/360.0d0)
!
!     find nearest integer and residual of the angle value
!
      ix = nint(x)
      z = x - dble(ix)
      y = z * z
!
!     look up the sine and cosine of the nearest integer
!
      if (ix .le. 180) then
         ssign = 1
         csign = 1
      else
         ssign = -1
         csign = -1
         ix = ix - 180
      end if
      if (ix .gt. 90) then
         csign = -csign
         ix = 180 - ix
      end if
      sx = table(ix)
      cx = table(90-ix)
      if (ssign .lt. 0)  sx = -sx
      if (csign .lt. 0)  cx = -cx
!
!     find sine and cosine of residual angle to 5 decimal places
!
!     sz = 1.74531263774940077459d-2 * z
!     cz = 1.0d0 - 1.52307909153324666207d-4 * y
!
!     find sine and cosine of residual angle to 11 decimal places
!
!     sz = (-8.86092781698004819918d-7 * y&
!    &        + 1.74532925198378577601d-2) * z
!     cz = 1.0d0 - (-3.86631403698859047896d-9 * y&
!    &                 + 1.52308709893047593702d-4) * y
!
!     find sine and cosine of residual angle to 17 decimal places
!
      sz = ((1.34959795251974073996d-11 * y&
     &        - 8.86096155697856783296d-7) * y&
     &           + 1.74532925199432957214d-2) * z
      cz = 1.0d0 - ((3.92582397764340914444d-14 * y&
     &                  - 3.86632385155548605680d-9) * y&
     &                     + 1.52308709893354299569d-4) * y
!
!     combine the tabulated and calculated parts by trigonometry
!
      sine = sx*cz + cx*sz
      if (xsign .lt. 0)  sine = -sine
      cosine = cx*cz - sx*sz
      return
      end

      SUBROUTINE spline(x, y, n, yp1, ypn, y2)
!   use nrtype
!
! Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e.
! y(i)=f(x(i)), with x(1)<x(2)<...<x(n), and given values yp1 and ypn for
! the first derivative of the interpolating function at points 1 and n,
! respectively, this routine returns an array y2(1:n) of length n which
! contains the second derivatives of the interpolating function at the
! tabulated points x(i).  If yp1 and/or ypn are equal to 1.e30 or larger,
! the routine is signaled to set the corresponding boundary condition for a
! natural spline with zero second derivative on that boundary.
! Parameter: nmax is the largest anticipiated value of n
! (adopted from Numerical Recipes in FORTRAN 77)
!
      INTEGER, PARAMETER :: DP=KIND(1.0D0)
      INTEGER:: n
      INTEGER, PARAMETER:: nmax=500
      REAL(8):: yp1, ypn, x(n), y(n), y2(n)
      INTEGER:: i, k
      REAL(8):: p, qn, sig, un, u(nmax)
      real(8) :: dyda,dydb,dydx,dydax,dydbx
      
      if (yp1.gt..99e30) then
         y2(1)=0.
         u(1)=0.
      else
         y2(1)=-0.5
         u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif

      do i=2, n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/&
        & (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo

      if (ypn.gt..99e30) then
         qn=0.
         un=0.
      else
         qn=0.5
         un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif

      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)

      do k=n-1, 1, -1
         y2(k)=y2(k)*y2(k+1)+u(k)
      enddo

      return
      END


      subroutine Splint(xa,ya,y2a,n,x,y,dydx)
!
!     Spline interpolation  adapted from Numerical Recipes.
!     Give a value y at x. input xa(n) is equally distributed.
!     re-tested on 20/6/97 by Lu. 
!     note that x must be smaller xa(n)!!! otherwise set y = 0.

      implicit double precision (a-h,o-z)
      integer n, klo, k
      dimension xa(n),ya(n),y2a(n)
      real(8) :: sixi
      parameter (sixi = 0.166666666666667d0)

      if (x .gt. xa(n)) then
         y = 0.0d0
         return
         endif

         klo=1
         khi=n
 1       if (khi-klo.gt.1) then
            k=(khi+klo)/2
            if(xa(k).gt.x)then
               khi=k
            else
               klo=k
            endif
            goto 1
         endif
         h=xa(khi)-xa(klo)
         if (h.eq.0.) then
            write(*,*) 'bad xa input.'
            stop
         endif
         hi = 1.0d0/h
         h2sixi = (h**2)*sixi

         a=(xa(khi)-x)*hi
         b=(x-xa(klo))*hi
         y=a*ya(klo)+b*ya(khi)+ &
     &        ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*h2sixi

         dyda = ya(klo) + ((3.0d0 * a**2 - 1.0d0) * y2a(klo))*h2sixi
         dydb = ya(khi) + ((3.0d0 * b**2 - 1.0d0) * y2a(khi))*h2sixi
         dydx = hi*(-dyda + dydb)

         return
         end

      SUBROUTINE mat3inverse (A, AINV, OK_FLAG)

!***********************************************************************************************************************************
!  M33INV  -  Compute the inverse of a 3x3 matrix.
!
!  A       = input 3x3 matrix to be inverted
!  AINV    = output 3x3 inverse of matrix A
!  OK_FLAG = (output) .TRUE. if the input matrix could be inverted, and .FALSE. if the input matrix is singular.
!***********************************************************************************************************************************


      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A
      DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: AINV
      LOGICAL, INTENT(OUT) :: OK_FLAG

      DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
      DOUBLE PRECISION :: DET
      DOUBLE PRECISION, DIMENSION(3,3) :: COFACTOR


      DET =   A(1,1)*A(2,2)*A(3,3)  &
            - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3)  &
            + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2)  &
            - A(1,3)*A(2,2)*A(3,1)

      IF (ABS(DET) .LE. EPS) THEN
         AINV = 0.0D0
         OK_FLAG = .FALSE.
         RETURN
      END IF

      COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

      AINV = TRANSPOSE(COFACTOR) / DET

      OK_FLAG = .TRUE.

      RETURN

      END SUBROUTINE mat3inverse

      subroutine mat3multvec(matrix,vin,vout)
      real(8), dimension(3,3) :: matrix
      real(8), dimension(3) :: vin,vout
      integer i,j,k
      real(8) :: sum

      do i = 1,3
         sum = 0.0d0
         do j = 1,3
            sum = sum + matrix(i,j)*vin(j)
         end do 
         vout(i) = sum
      end do 
      end subroutine mat3multvec

      end module nr_mod

      FUNCTION f1dim(x)
      INTEGER NMAX
      REAL(8) f1dim,func,x
      PARAMETER (NMAX=60000)
!     USES func
      INTEGER j,ncom
      REAL(8) pcom(NMAX),xicom(NMAX),xt(NMAX)
      COMMON /f1com/ pcom,xicom,ncom
      do 11 j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
 11         continue
      f1dim=func(xt)

      return
      END

      FUNCTION df1dim(x)
      INTEGER NMAX
      REAL(8) df1dim,x
      PARAMETER (NMAX=60000)
!     USES dfunc
      INTEGER j,ncom
      REAL(8) df(NMAX),pcom(NMAX),xicom(NMAX),xt(NMAX)
      COMMON /f1com/ pcom,xicom,ncom
      do 11 j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
 11         continue
      call dfunc(xt,df)
      df1dim=0.
      do 12 j=1,ncom
        df1dim=df1dim+df(j)*xicom(j)
 12         continue
      return
      END

      FUNCTION gammln(xx)
      REAL(8) gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,&
     &24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,&
     &-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
 11         continue
      gammln=tmp+log(stp*ser/x)
      return
      END

      FUNCTION dbrent(ax,bx,cx,f,df,tol,xmin)
      INTEGER ITMAX
      REAL(8) dbrent,ax,bx,cx,tol,xmin,df,f,ZEPS
      EXTERNAL df,f
      PARAMETER (ITMAX=100,ZEPS=1.0e-10)
      INTEGER iter
      REAL(8) a,b,d,d1,d2,du,dv,dw,dx,e,fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,&
     &     v,w,x,xm
      LOGICAL ok1,ok2
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.
      fx=f(x)
      fv=fx
      fw=fx
      dx=df(x)
      dv=dx
      dw=dx
      do 11 iter=1,ITMAX
        xm=0.5*(a+b)
        tol1=tol*abs(x)+ZEPS
        tol2=2.*tol1
        if(abs(x-xm).le.(tol2-.5*(b-a))) goto 3
        if(abs(e).gt.tol1) then
          d1=2.*(b-a)
          d2=d1
          if(dw.ne.dx) d1=(w-x)*dx/(dx-dw)
          if(dv.ne.dx) d2=(v-x)*dx/(dx-dv)
          u1=x+d1
          u2=x+d2
          ok1=((a-u1)*(u1-b).gt.0.).and.(dx*d1.le.0.)
          ok2=((a-u2)*(u2-b).gt.0.).and.(dx*d2.le.0.)
          olde=e
          e=d
          if(.not.(ok1.or.ok2))then
            goto 1
          else if (ok1.and.ok2)then
            if(abs(d1).lt.abs(d2))then
              d=d1
            else
              d=d2
            endif
          else if (ok1)then
            d=d1
          else
            d=d2
          endif
          if(abs(d).gt.abs(0.5*olde))goto 1
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
        endif
 1             if(dx.ge.0.) then
          e=a-x
        else
          e=b-x
        endif
        d=0.5*e
 2             if(abs(d).ge.tol1) then
          u=x+d
          fu=f(u)
        else
          u=x+sign(tol1,d)
          fu=f(u)
          if(fu.gt.fx)goto 3
        endif
        du=df(u)
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          dv=dw
          w=x
          fw=fx
          dw=dx
          x=u
          fx=fu
          dx=du
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            dv=dw
            w=u
            fw=fu
            dw=du
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
            dv=du
          endif
        endif
 11         continue
            write(*,*) 'dbrent exceeded maximum iterations'
 3         xmin=x
      dbrent=fx
      return
      END


! TRUE if x1*x2 negative
integer Function RootBracketed(x1,x2)
  real*8 x1,x2 
  integer resultat
  if ((x1 > 0.and.x2 > 0).or.(x1 < 0.and.x2 < 0)) then 
    resultat = 0
  else
    resultat = 1
  endif
  RootBracketed = resultat
end

! returns the minimum of two real numbers
real*8 Function Minimum(x1,x2) 
  real*8 x1,x2,resultat
  if (x1 < x2) then
    resultat = x1
  else 
    resultat = x2
  endif
  Minimum = resultat
end
