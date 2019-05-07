      integer iostat
      real(8), dimension(3,3) :: cvec
      real(8) :: amag,bmag,cmag,alpha,beta,gamma
      character(len = 8) name,type
      real(8), dimension(6) :: cparam
      real(8) :: ax,ay,az,bx,by,bz,cx,cy,cz
      open(10,file = 'best.xyz')
      open(11,file = 'super.xyz')

      rcut = 8.0d0 

      n = 0 
      cvec = 0.0d0 

      do  while(.true.)
         read(10,*,iostat = ios) natoms
         if(ios.lt.0) exit

         read(10,*) type,cparam,uuu
         do i = 1,natoms
            read(10,*) name,x,y,z
         end do
      end do 
      rewind(10)

      call get_cvec(type,cparam,cvec)

      ax = cvec(1,1) 
      ay = cvec(2,1) 
      az = cvec(3,1)

      bx = cvec(1,2)
      by = cvec(2,2)
      bz = cvec(3,2)

      cx = cvec(1,3)
      cy = cvec(2,3)
      cz = cvec(3,3)

      vol = ax * (by * cz - bz * cy) &
     &     +ay * (bz * cx - bx * cz) & 
     &     +az * (bx * cy - by * cx)

      axbx = ay*bz - az*by
      axby = az*bx - ax*bz
      axbz = ax*by - ay*bx

      axbmag = sqrt(axbx**2 + axby**2 + axbz**2)

      bxcx = by*cz - bz*cy
      bxcy = bz*cx - bx*cz
      bxcz = bx*cy - by*cx

      bxcmag = sqrt(bxcx**2 + bxcy**2 + bxcz**2)

      cxax = cy*az - cz*ay
      cxay = cz*ax - cx*az
      cxaz = cx*ay - cy*ax

      cxamag = sqrt(cxax**2 + cxay**2 + cxaz**2)
      
      jcell1max = int(2.0d0 * rcut * bxcmag / vol) + 1
      jcell2max = int(2.0d0 * rcut * cxamag / vol) + 1
      jcell3max = int(2.0d0 * rcut * axbmag / vol) + 1

!      jcell1max = 3
!      jcell2max = 3
!      jcell3max = 3

      nrep =  jcell1max*jcell2max*jcell3max
      write(*,*) jcell1max,jcell2max,jcell3max
      write(*,*)

      do while(.true.)
         read(10,*,iostat = ios) natoms
         if(ios.lt.0) exit
         n = n + 1
         write(*,*) n,natoms
         write(11,*) natoms*nrep
         read(10,*) type,cparam,uuu

         call get_cvec(type,cparam,cvec)
         call get_angle_cvec(cvec,amag,bmag,cmag,alpha,beta,gamma)

         ax = cvec(1,1) 
         ay = cvec(2,1) 
         az = cvec(3,1)
         
         bx = cvec(1,2)
         by = cvec(2,2)
         bz = cvec(3,2)

         cx = cvec(1,3)
         cy = cvec(2,3)
         cz = cvec(3,3)

         write(11,*) 'ANGLE',amag,bmag,cmag,alpha,beta,gamma,uuu,jcell1max,jcell2max,jcell3max
         do i = 1,natoms
            read(10,*) name,x,y,z
            do j1 = 0,jcell1max-1
               do j2 = 0,jcell2max-1
                  do j3 = 0,jcell3max-1
                     xx =  x + cvec(1,1) * j1 + cvec(1,2) * j2 + cvec(1,3) * j3
                     yy =  y + cvec(2,1) * j1 + cvec(2,2) * j2 + cvec(2,3) * j3
                     zz =  z + cvec(3,1) * j1 + cvec(3,2) * j2 + cvec(3,3) * j3
                     write(11,*) name,xx,yy,zz
                  end do 
               end do 
            end do 
      end do 
         if(iostat.lt.0) exit
      end do 
      end



      subroutine get_cartesian_cvec(amag,bmag,cmag,alpha,beta,gamma,cvec)
      implicit none
      real(8) :: alpha,beta,gamma,amag,bmag,cmag
      real(8) :: alpha_rad,beta_rad,gamma_rad
      real(8), dimension(6) :: param
      real(8) :: angfac,pi
      real(8), dimension(3,3) :: cvec

      pi = 4.0d0 * datan(1.d0)
      angfac = 180.0d0 / pi

      alpha_rad = alpha / angfac
      beta_rad = beta / angfac
      gamma_rad = gamma / angfac

      cvec = 0.0d0

      cvec(1,1) = amag
      cvec(1,2) = bmag * dcos(gamma_rad)
      cvec(2,2) = dsqrt(bmag**2 - cvec(1,2)**2)
      cvec(1,3) = cmag * dcos(beta_rad)
      cvec(2,3) = (bmag * cmag * dcos(alpha_rad) - cvec(1,2) * cvec(1,3))/cvec(2,2)
      cvec(3,3) = dsqrt(cmag**2 - cvec(1,3)**2 - cvec(2,3)**2)

      end subroutine get_cartesian_cvec

      subroutine get_cvec(type,cparam,cvec)
      real(8) :: amag,bmag,cmag,alpha,beta,gamma
      real(8), dimension(3,3) :: cvec
      real(8), dimension(6) :: cparam
      character(len = 8) type

      if(type.eq.'ANGLE') then
         amag = cparam(1)
         bmag = cparam(2)
         cmag = cparam(3)
         alpha = cparam(4)
         beta = cparam(5)
         gamma = cparam(6)
         call get_cartesian_cvec(amag,bmag,cmag,alpha,beta,gamma,cvec)
      else if(type.eq.'XYZ') then
         cvec(1,1) = cparam(1)
         cvec(1,2) = cparam(2)
         cvec(2,2) = cparam(3)
         cvec(1,3) = cparam(4)
         cvec(2,3) = cparam(5)
         cvec(3,3) = cparam(6)
      endif

      end subroutine get_cvec


      subroutine get_angle_cvec(cvec,amag,bmag,cmag,alpha,beta,gamma)
      implicit none
      real(8) :: angfac
      real(8) :: ax,ay,az,bx,by,bz,cx,cy,cz
      real(8) :: amag,bmag,cmag,alpha,beta,gamma
      real(8), dimension(6) :: param,param2
      real(8), dimension(3,3) :: cvec
      real(8) :: pi

      pi = 4.00d0 * datan(1.0d0) 

      angfac = 180.0d0 / pi

      ax = cvec(1,1)
      ay = 0.0d0
      az = 0.0d0

      bx = cvec(1,2)
      by = cvec(2,2)
      bz = 0.0d0

      cx = cvec(1,3)
      cy = cvec(2,3)
      cz = cvec(3,3)

      amag = dsqrt(ax**2 + ay**2 + az**2)
      bmag = dsqrt(bx**2 + by**2 + bz**2)
      cmag = dsqrt(cx**2 + cy**2 + cz**2)

      alpha = dacos((bx*cx + by*cy + bz*cz)/(bmag*cmag)) * angfac
      beta = dacos((cx*ax + cy*ay + cz*az)/(cmag*amag)) * angfac
      gamma = dacos((ax*bx + ay*by + az*bz)/(amag*bmag)) * angfac

      end subroutine get_angle_cvec

