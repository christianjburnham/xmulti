      implicit none
      include "make_prob.common"
      real(8) :: beta
      real(8) :: trans_mat(0:3,128,128)

      call init_random_seed

      open(10,file = 'molecular_geometry.mannitol')
      open(11,file = 'model.mannitol')
      open(20,file = 'prob.dat')

      beta = 1.0d0 

      call read_conformers()
      call calculate_iso_msd()
!      call test
      call calculate_probs(beta,trans_mat)
      call walk(trans_mat)

      end


      subroutine test
      implicit none
      include "make_prob.common"      
      real(8) :: energy0
      real(8), dimension(3) :: grad0
      real(8) :: delta
      integer i
      
      write(*,*) 'test'

      delta = 1.d-9

      phi = 2.0d0 
      theta = 3.0d0 
      psi = 1.0d0


      call get_energy
      energy0 = energy
      grad0 = grad

      do i = 1,3
         if(i.eq.1) phi = phi + delta
         if(i.eq.2) theta = theta + delta
         if(i.eq.3) psi = psi + delta

         call get_energy

         write(*,*) (energy - energy0)/delta, grad0(i)

         if(i.eq.1) phi = phi - delta
         if(i.eq.2) theta = theta - delta
         if(i.eq.3) psi = psi - delta
      end do 



      end subroutine test


      subroutine get_energy()
      implicit none
      include "make_prob.common"
      integer i,ivec
      real(8) :: xdif,ydif,zdif
      real(8) :: gphi_x,gphi_y,gphi_z
      real(8) :: gtheta_x,gtheta_y,gtheta_z
      real(8) :: gpsi_x,gpsi_y,gpsi_z
      real(8), dimension(3,3,3) :: drotmat
      real(8), dimension(3) :: vec1,vec2
      real(8), dimension(3,3) :: mat,gg
      real(8) :: totmass

!     first thing we need to do is to find the rotated coordinates
      call rotate()
      call get_drotmat(phi,theta,psi,drotmat)
      
      energy = 0.0d0 
      grad = 0.0d0 
      totmass = 0.0d0 

      do i = 1,natoms
         xdif = ss2(i,1) - rr1(i,1)
         ydif = ss2(i,2) - rr1(i,2)
         zdif = ss2(i,3) - rr1(i,3)
         energy = energy + atmass(i) * (xdif**2 + ydif**2 + zdif**2)
         totmass = totmass + atmass(i)

         vec1(1) = rr2(i,1)
         vec1(2) = rr2(i,2) 
         vec1(3) = rr2(i,3)

         do ivec = 1,3
            mat(:,:) = drotmat(:,:,ivec)
            call mat3multvec(mat,vec1,vec2)
            
            gg(:,ivec) = vec2(:)
         end do 

         grad(1) = grad(1) + 2.0d0 * atmass(i) * (xdif * gg(1,1) + ydif * gg(2,1) + zdif * gg(3,1))
         grad(2) = grad(2) + 2.0d0 * atmass(i) * (xdif * gg(1,2) + ydif * gg(2,2) + zdif * gg(3,2))
         grad(3) = grad(3) + 2.0d0 * atmass(i) * (xdif * gg(1,3) + ydif * gg(2,3) + zdif * gg(3,3))

      end do 

      energy = energy/totmass
      grad = grad / totmass
      
      end subroutine get_energy

      subroutine minimize()
      implicit none
      include "make_prob.common"
      real(8) :: xvec(3)
      real(8) :: ftol,fret
      integer iter,n

      xvec(1) = phi
      xvec(2) = theta
      xvec(3) = psi

      n = 3
      ftol = 0.0d0

      call frprmn(xvec,n,ftol,iter,fret)

      end subroutine minimize

      real(8) function func(xvec)
      implicit none
      include "make_prob.common"
      real(8), dimension(3) :: xvec

      phi = xvec(1)
      theta = xvec(2)
      psi = xvec(3)

      call get_energy()
      func = energy

      end function func

      subroutine dfunc(xvec,gvec)
      implicit none
      include "make_prob.common"
      real(8) :: xvec(3),gvec(3)
      gvec = grad
      end subroutine dfunc


      subroutine rotate()
      implicit none
      include "make_prob.common"
      real(8), dimension(3,3) :: rotmat
      real(8), dimension(3) :: vec1,vec2
      integer i

      call get_rotmat(phi,theta,psi,rotmat)

      do i = 1,natoms
         vec1(:) = rr2(i,:)
         call mat3multvec(rotmat,vec1,vec2)
         ss2(i,:) = vec2(:)
      end do 

      end subroutine rotate


      subroutine get_rotmat(phi,theta,psi,rotmat)
!     generates the rotation matrix from the euler angles

      implicit none
      real(8), dimension(3,3) :: rotmat
      real(8) :: phi,theta,psi
      real(8) :: s_phi,c_phi,s_theta,c_theta,s_psi,c_psi

      s_phi = dsin(phi)
      c_phi = dcos(phi)

      s_theta = dsin(theta)
      c_theta = dcos(theta)

      s_psi = dsin(psi)
      c_psi = dcos(psi)

      rotmat(1,1) = c_phi * c_psi - s_phi * c_theta * s_psi
      rotmat(1,2) = s_phi * c_psi + c_phi * c_theta * s_psi
      rotmat(1,3) = s_theta * s_psi
      
      rotmat(2,1) = -c_phi * s_psi - s_phi * c_theta * c_psi
      rotmat(2,2) = -s_phi * s_psi + c_phi * c_theta * c_psi
      rotmat(2,3) = s_theta * c_psi

      rotmat(3,1) = s_phi * s_theta
      rotmat(3,2) = -c_phi * s_theta
      rotmat(3,3) = c_theta

      end subroutine get_rotmat

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

      subroutine get_drotmat(phi,theta,psi,drotmat)
!     generates the derivatives of the rotation matrix from the euler angles

      implicit none
      real(8), dimension(3,3,3) :: drotmat
      real(8) :: phi,theta,psi
      real(8) :: s_phi,c_phi,s_theta,c_theta,s_psi,c_psi

      s_phi = dsin(phi)
      c_phi = dcos(phi)

      s_theta = dsin(theta)
      c_theta = dcos(theta)

      s_psi = dsin(psi)
      c_psi = dcos(psi)
      

!     derivatives wrt phi

      drotmat(1,1,1) = -s_phi * c_psi - c_phi * c_theta * s_psi
      drotmat(1,2,1) = c_phi * c_psi - s_phi * c_theta * s_psi
      drotmat(1,3,1) = 0.0d0
      
      drotmat(2,1,1) = s_phi * s_psi - c_phi * c_theta * c_psi   
      drotmat(2,2,1) = -c_phi * s_psi - s_phi * c_theta *c_psi    
      drotmat(2,3,1) = 0.0d0 

      drotmat(3,1,1) = c_phi * s_theta
      drotmat(3,2,1) = s_phi * s_theta
      drotmat(3,3,1) = 0.0d0 

!     derivatives wrt theta

      drotmat(1,1,2) = s_phi * s_theta * s_psi
      drotmat(1,2,2) = - c_phi * s_theta * s_psi
      drotmat(1,3,2) = c_theta * s_psi
      
      drotmat(2,1,2) = s_phi * s_theta * c_psi
      drotmat(2,2,2) = - c_phi * s_theta *c_psi
      drotmat(2,3,2) = c_theta * c_psi

      drotmat(3,1,2) = s_phi * c_theta
      drotmat(3,2,2) = -c_phi * c_theta
      drotmat(3,3,2) = -s_theta

!     derivatives wrt psi

      drotmat(1,1,3) = - c_phi * s_psi - s_phi * c_theta * c_psi
      drotmat(1,2,3) = - s_phi * s_psi + c_phi * c_theta * c_psi
      drotmat(1,3,3) = s_theta * c_psi
      
      drotmat(2,1,3) = -c_phi * c_psi + s_phi * c_theta * s_psi
      drotmat(2,2,3) = -s_phi * c_psi - c_phi * c_theta * s_psi
      drotmat(2,3,3) = -s_theta * s_psi

      drotmat(3,1,3) = 0.0d0 
      drotmat(3,2,3) = 0.0d0 
      drotmat(3,3,3) = 0.0d0 

      end subroutine get_drotmat


      SUBROUTINE FRPRMN(p,n,ftol,iter,fret)
      INTEGER iter,n,NMAX,ITMAX
      REAL(8) fret,ftol,p(n),EPS,func
      EXTERNAL func
      PARAMETER (NMAX=60000,ITMAX=2000000,EPS=1.d-10)
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
!           write(*,*) nsteps,fp
!           write(90,*) nsteps,fp
!           call flush(90)
!           call print_coordinates(89)
!           call print_rigid_coordinates(88)

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
      xx = 1.e-9
!      xx = 1.0d0 

      call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
      fret=dbrent(ax,xx,bx,f1dim,df1dim,TOL,xmin)

      do 12 j=1,n
        xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
 12         continue
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

      subroutine read_conformers()
      implicit none
      include "make_prob.common"
      integer iso,ios,i,jj,iatom
      character(len = 10) molname,aa,bb

!     read in the masses
      read(11,*) aa,bb,natoms
      do iatom = 1,natoms
         read(11,*) atname(iatom),jj,atmass(iatom)
      end do 

      iso = 0 
      do while(.true.)
         read(10,*,iostat = ios) molname,natoms
         if(ios.lt.0) exit
         iso = iso + 1
         do i = 1,natoms
            read(10,*) atname(i),rr(iso,i,1),rr(iso,i,2),rr(iso,i,3)
         end do 
      end do 
      
      nconformers = iso

      write(20,*) molname,nconformers,'gg'

      end subroutine read_conformers

      subroutine calculate_iso_msd()
      implicit none
      include "make_prob.common"
      integer i,iso1,iso2,m,mmax
      real(8) :: energy0
      real(8) :: rand,phibest,thetabest,psibest,ebest,pi

      pi = 4.d0 * datan(1.d0)

      do iso1 = 1,nconformers
         do iso2 = iso1+1,nconformers
            
            do i = 1,natoms
               rr1(i,:) = rr(iso1,i,:)
               rr2(i,:) = rr(iso2,i,:)
            end do 
            
            call get_energy
            energy0 = energy

            phi = 0.0d0 
            theta = 0.0d0 
            psi = 0.0d0 
            ebest = 1.d+20
            mmax = 8

            do m = 1,mmax
               if(m.gt.1) then 
                  call random_number(rand)
                  phi = rand * 2.0d0 * pi
                  call random_number(rand)
                  theta = rand * 2.0d0 * pi
                  call random_number(rand)
                  psi = rand * 2.0d0 * pi
               endif
               call minimize

               if(energy.lt.ebest) then
                  ebest = energy
                  phibest = phi
                  thetabest = theta
                  psibest = psi
               endif
            end do 
            phi = phibest
            theta = thetabest
            psi = psibest
            energy = ebest

            write(*,*) iso1,iso2,energy
            iso_mat(iso1,iso2,0) = energy

            iso_mat(iso1,iso2,1) = phi
            iso_mat(iso1,iso2,2) = theta
            iso_mat(iso1,iso2,3) = psi

         end do 
      end do 


      end subroutine calculate_iso_msd

      subroutine get_probs(beta,iso,prob,angle)
      implicit none
      include "make_prob.common"
      integer iso,i,i1,i2,n,nmax
      real(8), dimension(100) :: msd,prob
      real(8), dimension(100,3) :: angle
      real(8) :: sum,beta,aa

      sum = 0.0d0 

      n = 0 
      do i = 1,nconformers
         i1 = min(i,iso)
         i2 = max(i,iso)
         n = n + 1
         if(i1.ne.i2) then 
            msd(n) = iso_mat(i1,i2,0)
            sum = sum + dexp(-beta * msd(n))

            angle(n,1) = iso_mat(i1,i2,1)
            angle(n,2) = iso_mat(i1,i2,2)
            angle(n,3) = iso_mat(i1,i2,3)
         else
            msd(n) = 0.d0
            angle(n,1) = 0.0d0 
            angle(n,2) = 0.0d0 
            angle(n,3) = 0.0d0 
         endif
      end do 

!     calculate the probabilities 

      do n = 1,nconformers
         if(n.ne.iso) then 
         prob(n) = dexp(-beta * msd(n)) / sum
         else
            prob(n) = 0.0d0 
         endif
      end do 

      end subroutine get_probs


      subroutine init_random_seed()
      integer :: i, n, clock
      integer, dimension(:), allocatable :: seed

      call random_seed(size = n)
      allocate(seed(n))

      call system_clock(count=clock)

      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      call random_seed(put = seed)

      deallocate(seed)
      end subroutine init_random_seed

      subroutine print_conformer(iso,rotmat)
      implicit none
      include "make_prob.common"
      integer iso,i
      real(8) :: x,y,z,xs,ys,zs
      real(8), dimension(3) :: vec1,vec2
      real(8), dimension(3,3) :: rotmat

      write(66,*) natoms
      write(66,*) 
      do i = 1,natoms
         x = rr(iso,i,1)
         y = rr(iso,i,2)
         z = rr(iso,i,3)
         vec1 = (/x,y,z/)
         call mat3multvec(rotmat,vec1,vec2)
         xs = vec2(1)
         ys = vec2(2)
         zs = vec2(3) 

         write(66,*) atname(i),xs,ys,zs
      end do 

      end subroutine print_conformer


      subroutine calculate_probs(beta,trans_mat)
      implicit none
      include "make_prob.common"
      real(8) :: beta
      integer iso,i
      real(8), dimension(100) :: prob
      real(8), dimension(100,3) :: angle
      real(8) :: trans_mat(0:3,128,128)
      real(8) :: phi1,theta1,psi1
      real(8),dimension(3,3) :: rotmat

      write(*,*) 'CALC PROBS'

      do iso = 1,nconformers

         call get_probs(beta,iso,prob,angle)

         do i = 1,nconformers

            if(i.lt.iso) then 
               phi1 = angle(i,1)
               theta1 = angle(i,2)
               psi1 = angle(i,3)
               call get_rotmat(phi1,theta1,psi1,rotmat)
               rotmat = transpose(rotmat)
               call get_euler_from_rotmat(rotmat,phi1,theta1,psi1)
               angle(i,1) = phi1
               angle(i,2) = theta1
               angle(i,3) = psi1
            endif

            trans_mat(0,iso,i) = prob(i)
            trans_mat(1,iso,i) = angle(i,1)
            trans_mat(2,iso,i) = angle(i,2)
            trans_mat(3,iso,i) = angle(i,3)
            write(*,*) 'HELLO'
            write(20,*) iso,i,prob(i),angle(i,:)
         end do 
      end do 

      end subroutine calculate_probs

      subroutine get_euler_from_rotmat(rotmat,phi,theta,psi)
      implicit none
      real(8), dimension(3,3) :: rotmat
      real(8) :: phi,theta,psi,stheta,spsi

      phi = datan2(rotmat(3,1),-rotmat(3,2))
      psi = datan2(rotmat(1,3),rotmat(2,3))

      spsi = dsin(psi)

      if(dabs(spsi).ge.1.d-5) then 
         stheta = rotmat(1,3) / spsi
      else
         stheta = rotmat(2,3) / dcos(psi)
      endif
      theta = datan2(stheta,rotmat(3,3))

      end subroutine get_euler_from_rotmat

      subroutine walk(trans_mat)
      implicit none
      include "make_prob.common"
      integer conformer,n,nmax
      real(8), dimension(3,3) :: rotmat,rotmat1
      real(8) :: sum
      real(8) :: phi1,theta1,psi1
      real(8) :: random
      real(8) :: trans_mat(0:3,128,128)
      integer, dimension(2048) :: hist
      integer k,new_conformer,i

      hist = 0
      conformer = 1
      nmax = 512

      rotmat = 0.0d0 
      rotmat(1,1) = 1.0d0 
      rotmat(2,2) = 1.0d0 
      rotmat(3,3) = 1.0d0 

      do n = 1,nmax
         sum = 0.0d0 
         call random_number(random)
         do k = 1,nconformers
            sum = sum + trans_mat(0,conformer,k)
            if(sum.gt.random) exit
         end do 
         new_conformer = k
         hist(new_conformer) = hist(new_conformer) + 1

         phi1 = trans_mat(1,conformer,new_conformer)
         theta1 = trans_mat(2,conformer,new_conformer)
         psi1 = trans_mat(3,conformer,new_conformer)

         call get_rotmat(phi1,theta1,psi1,rotmat1)
         rotmat = matmul(rotmat1,rotmat)

         conformer = new_conformer
         write(*,*) 'conformer = ',conformer
         call print_conformer(conformer,rotmat)
      end do 

      write(*,*) 'HIST'
      do i = 1,nconformers
         write(*,*) i,hist(i)
      end do 

      end subroutine walk
