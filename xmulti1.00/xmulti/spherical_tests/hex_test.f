      implicit none
      real(8) a_xxxx,a_xxxy,a_xxxz,a_xxyy,a_xxyz,a_xxzz,a_xyyy,&
     & a_xyyz,a_xyzz,a_xzzz,a_yyyy,a_yyyz,a_yyzz,a_yzzz,a_zzzz
      real(8) b_xxxx,b_xxxy,b_xxxz,b_xxyy,b_xxyz,b_xxzz,b_xyyy,&
     & b_xyyz,b_xyzz,b_xzzz,b_yyyy,b_yyyz,b_yyzz,b_yzzz,b_zzzz
      real(8), dimension(9) :: a_solid4,b_solid4
      real(8) :: sum1,sum2
      real(8) :: rmin,rmax
      real(8) :: random
      integer i

      write(*,*) '*********************************************'
      write(*,*) '**  HEX TEST'
      write(*,*) '**'
      write(*,*) '**  short test showing conversion of cartesians'
      write(*,*) '**  to spherical harmonics for hexadecapoles.'
      write(*,*) '**'
      write(*,*) '** Christian J. Burnham, University &
     &College Dublin, 2019'
      write(*,*) '** @christianjburnham@gmail.com'
      write(*,*) '*********************************************'

      call init_random_seed()

      rmin = -1.0d0 
      rmax = 1.0d0 

      a_xxxx = random(rmin,rmax)
      a_xxxy = random(rmin,rmax)
      a_xxxz = random(rmin,rmax)

      a_xxyy = random(rmin,rmax) 
      a_xxyz = random(rmin,rmax)
      a_xxzz = random(rmin,rmax)

      a_xyyy = random(rmin,rmax)
      a_xyyz = random(rmin,rmax)
      a_xyzz = random(rmin,rmax)

      a_xzzz = random(rmin,rmax)
      a_yyyy = random(rmin,rmax)
      a_yyyz = random(rmin,rmax)

      a_yyzz = random(rmin,rmax)
      a_yzzz = random(rmin,rmax)
      a_zzzz = random(rmin,rmax)

      b_xxxx = random(rmin,rmax)
      b_xxxy = random(rmin,rmax)
      b_xxxz = random(rmin,rmax)

      b_xxyy = random(rmin,rmax)
      b_xxyz = random(rmin,rmax)
      b_xxzz = random(rmin,rmax)

      b_xyyy = random(rmin,rmax)
      b_xyyz = random(rmin,rmax)
      b_xyzz = random(rmin,rmax)

      b_xzzz = random(rmin,rmax)
      b_yyyy = random(rmin,rmax)
      b_yyyz = random(rmin,rmax)

      b_yyzz = random(rmin,rmax)
      b_yzzz = random(rmin,rmax)
      b_zzzz = random(rmin,rmax)

      write(*,*) 'initial cartesians'
      write(*,*) 'a'
      write(*,17) a_xxxx,a_xxxy,a_xxxz&
     &           ,a_xxyy,a_xxyz,a_xxzz&
     &           ,a_xyyy,a_xyyz,a_xyzz&
     &           ,a_xzzz,a_yyyy,a_yyyz&
     &           ,a_yyzz,a_yzzz,a_zzzz
      write(*,*) 'b'
      write(*,17) b_xxxx,b_xxxy,b_xxxz&
     &           ,b_xxyy,b_xxyz,b_xxzz&
     &           ,b_xyyy,b_xyyz,b_xyzz&
     &           ,b_xzzz,b_yyyy,b_yyyz&
     &           ,b_yyzz,b_yzzz,b_zzzz


      call convert_hex_to_spherical(&     
     & a_xxxx,a_xxxy,a_xxxz,a_xxyy,a_xxyz,a_xxzz,a_xyyy,&
     & a_xyyz,a_xyzz,a_xzzz,a_yyyy,a_yyyz,a_yyzz,a_yzzz,a_zzzz,&
     & a_solid4)
      call convert_hex_to_cartesian(a_solid4,&
     & a_xxxx,a_xxxy,a_xxxz,a_xxyy,a_xxyz,a_xxzz,a_xyyy,&
     & a_xyyz,a_xyzz,a_xzzz,a_yyyy,a_yyyz,a_yyzz,a_yzzz,a_zzzz)
      call convert_hex_to_spherical(&     
     & b_xxxx,b_xxxy,b_xxxz,b_xxyy,b_xxyz,b_xxzz,b_xyyy,&
     & b_xyyz,b_xyzz,b_xzzz,b_yyyy,b_yyyz,b_yyzz,b_yzzz,b_zzzz,&
     & b_solid4)
      call convert_hex_to_cartesian(b_solid4,&
     & b_xxxx,b_xxxy,b_xxxz,b_xxyy,b_xxyz,b_xxzz,b_xyyy,&
     & b_xyyz,b_xyzz,b_xzzz,b_yyyy,b_yyyz,b_yyzz,b_yzzz,b_zzzz)

      write(*,*)
      write(*,*) 'detraced cartesians'
      write(*,*) 'a'
      write(*,17) a_xxxx,a_xxxy,a_xxxz&
     &           ,a_xxyy,a_xxyz,a_xxzz&
     &           ,a_xyyy,a_xyyz,a_xyzz&
     &           ,a_xzzz,a_yyyy,a_yyyz&
     &           ,a_yyzz,a_yzzz,a_zzzz
      write(*,*) 'b'
      write(*,17) b_xxxx,b_xxxy,b_xxxz&
     &           ,b_xxyy,b_xxyz,b_xxzz&
     &           ,b_xyyy,b_xyyz,b_xyzz&
     &           ,b_xzzz,b_yyyy,b_yyyz&
     &           ,b_yyzz,b_yzzz,b_zzzz

      call convert_hex_to_spherical(&     
     & a_xxxx,a_xxxy,a_xxxz,a_xxyy,a_xxyz,a_xxzz,a_xyyy,&
     & a_xyyz,a_xyzz,a_xzzz,a_yyyy,a_yyyz,a_yyzz,a_yzzz,a_zzzz,&
     & a_solid4)
      call convert_hex_to_cartesian(a_solid4,&
     & a_xxxx,a_xxxy,a_xxxz,a_xxyy,a_xxyz,a_xxzz,a_xyyy,&
     & a_xyyz,a_xyzz,a_xzzz,a_yyyy,a_yyyz,a_yyzz,a_yzzz,a_zzzz)
      call convert_hex_to_spherical(&     
     & b_xxxx,b_xxxy,b_xxxz,b_xxyy,b_xxyz,b_xxzz,b_xyyy,&
     & b_xyyz,b_xyzz,b_xzzz,b_yyyy,b_yyyz,b_yyzz,b_yzzz,b_zzzz,&
     & b_solid4)
      call convert_hex_to_cartesian(b_solid4,&
     & b_xxxx,b_xxxy,b_xxxz,b_xxyy,b_xxyz,b_xxzz,b_xyyy,&
     & b_xyyz,b_xyzz,b_xzzz,b_yyyy,b_yyyz,b_yyzz,b_yzzz,b_zzzz)

      write(*,*) 
      write(*,*) 'check: applying the detracing operator a 2nd time.'
      write(*,*) 'a'
      write(*,17) a_xxxx,a_xxxy,a_xxxz&
     &           ,a_xxyy,a_xxyz,a_xxzz&
     &           ,a_xyyy,a_xyyz,a_xyzz&
     &           ,a_xzzz,a_yyyy,a_yyyz&
     &           ,a_yyzz,a_yzzz,a_zzzz
      write(*,*) 'b'
      write(*,17) b_xxxx,b_xxxy,b_xxxz&
     &           ,b_xxyy,b_xxyz,b_xxzz&
     &           ,b_xyyy,b_xyyz,b_xyzz&
     &           ,b_xzzz,b_yyyy,b_yyyz&
     &           ,b_yyzz,b_yzzz,b_zzzz

      write(*,*) 
      write(*,*) 'check: traces'
      write(*,*) 'trace(a)_xx',a_xxxx+a_xxyy+a_xxzz
      write(*,*) 'trace(a)_yy',a_xxyy+a_yyyy+a_yyzz
      write(*,*) 'trace(a)_zz',a_xxzz+a_yyzz+a_zzzz
      write(*,*) 'trace(a)_xy',a_xxxy+a_xyyy+a_xyzz
      write(*,*) 'trace(a)_xz',a_xxxz+a_xyyz+a_xzzz
      write(*,*) 'trace(a)_yz',a_xxyz+a_yyyz+a_yzzz

!     calculate the inner product in two ways.

      sum1 = a_xxxx * b_xxxx & 
     &      +a_xxxy * b_xxxy * 4.0d0& 
     &      +a_xxxz * b_xxxz * 4.0d0& 

     &      +a_xxyy * b_xxyy * 6.0d0& 
     &      +a_xxyz * b_xxyz * 12.0d0& 
     &      +a_xxzz * b_xxzz * 6.0d0& 

     &      +a_xyyy * b_xyyy * 4.0d0& 
     &      +a_xyyz * b_xyyz * 12.0d0& 
     &      +a_xyzz * b_xyzz * 12.0d0&

     &      +a_xzzz * b_xzzz * 4.0d0&
     &      +a_yyyy * b_yyyy &
     &      +a_yyyz * b_yyyz * 4.0d0&

     &      +a_yyzz * b_yyzz *6.0d0&
     &      +a_yzzz * b_yzzz *4.0d0&
     &      +a_zzzz * b_zzzz 


      sum2 = 0.0d0 
      do i = 1,9
         sum2 = sum2 + a_solid4(i)*b_solid4(i)
      end do 

      write(*,*)
      write(*,*) 'inner products'
      write(*,*) 'calculated in cartesians          = ',sum1
      write(*,*) 'calculated in spherical harmonics = ',sum2

 17   format('xxxx',f12.6,'  xxxy',f12.6,'  xxxz',f12.6&
     &      ,'  xxyy',f12.6,'  xxyz',f12.6,'  xxzz',f12.6&
     &      ,'  xyyy',f12.6,'  xyyz',f12.6,'  xyzz',f12.6&
     &      ,'  xzzz',f12.6,'  yyyy',f12.6,'  yyyz',f12.6&
     &      ,'  yyzz',f12.6,'  yzzz',f12.6,'  zzzz',f12.6)      

      end

      subroutine convert_hex_to_cartesian(solid4,&
     & xxxx,xxxy,xxxz,xxyy,xxyz,xxzz,xyyy,&
     & xyyz,xyzz,xzzz,yyyy,yyyz,yyzz,yzzz,zzzz)
      implicit none
      real(8) xxxx,xxxy,xxxz,xxyy,xxyz,xxzz,xyyy,&
     & xyyz,xyzz,xzzz,yyyy,yyyz,yyzz,yzzz,zzzz
      real(8), dimension(9) :: solid4
      real(8) :: s

      s = dsqrt(8.0d0 / 35.0d0)

      xxxx = s * ((3.0d0 / 8.0d0) * solid4(1) - (dsqrt(5.d0)/4.0d0) * solid4(4) + (dsqrt(35.0d0)/8.0d0) * solid4(8))
      xxxy = s*(1.0d0 / 8.0d0) * (-dsqrt(5.0d0) * solid4(5) + dsqrt(35.0d0) * solid4(9))
      xxyy = s*((1.0d0 / 8.0d0) * solid4(1) - (dsqrt(35.0d0)/8.0d0) * solid4(8))
      xyyy = -s*(1.0d0 / 8.0d0) * (dsqrt(5.0d0) * solid4(5) + dsqrt(35.0d0) * solid4(9))
      yyyy = s*((3.0d0 / 8.0d0) * solid4(1) + (dsqrt(5.0d0)/4.0d0) * solid4(4) + (dsqrt(35.0d0)/8.0d0) * solid4(8))
      xxxz = s*(1.0d0 / 16.0d0) * (-3.0d0 * dsqrt(10.0d0) * solid4(2) + dsqrt(70.0d0) * solid4(6))
      xxyz = s*(1.0d0 / 16.0d0) * (-dsqrt(10.0d0) * solid4(3) + dsqrt(70.0d0) * solid4(7))
      xyyz = -s*(1.0d0 / 16.0d0) * (dsqrt(10.0d0) * solid4(2) + dsqrt(70.0d0) * solid4(6))
      yyyz = -s*(1.0d0 / 16.0d0) * (3.0d0 * dsqrt(10.0d0) * solid4(3) + dsqrt(70.0d0) * solid4(7))      

      xxzz = s*(-0.5d0 * solid4(1) + (dsqrt(5.0d0)/4.0d0) * solid4(4))

      xyzz = s*(dsqrt(5.0d0)/4.0d0) * solid4(5)
      yyzz = s*(-0.5d0 * solid4(1) - (dsqrt(5.0d0)/4.0d0 ) * solid4(4))
      xzzz = s*dsqrt(5.0d0 / 8.0d0) * solid4(2)
      yzzz = s*dsqrt(5.0d0 / 8.0d0) * solid4(3)
      zzzz = s*solid4(1)
      
      end subroutine convert_hex_to_cartesian

      subroutine convert_hex_to_spherical(&
     & xxxx,xxxy,xxxz,xxyy,xxyz,xxzz,xyyy,&
     & xyyz,xyzz,xzzz,yyyy,yyyz,yyzz,yzzz,zzzz,&
     & solid4)
      implicit none
      real(8) xxxx,xxxy,xxxz,xxyy,xxyz,xxzz,xyyy,&
     & xyyz,xyzz,xzzz,yyyy,yyyz,yyzz,yzzz,zzzz
      real(8), dimension(9) :: solid4
      real(8) xxrr,xyrr,xzrr,yyrr,yzrr,zzrr,rrrr
      real(8) :: s

      s = dsqrt(8.0d0/35.0d0) 

      xxrr = xxxx + xxyy + xxzz
      yyrr = xxyy + yyyy + yyzz
      zzrr = xxzz + yyzz + zzzz

      rrrr = xxrr + yyrr + zzrr

!     40
      solid4(1) = s*(1.0d0/8.0d0) * (35.0d0 * zzzz - 30.0d0 * zzrr + 3.0d0 * rrrr)
!     41c
      solid4(2) = s*dsqrt(5.0d0/8.0d0) * (7.0d0 * xzzz - 3.0d0 * xzrr)
!     41s
      solid4(3) = s*dsqrt(5.0d0/8.0d0) * (7.0d0 * yzzz - 3.0d0 * yzrr)
!     42c
      solid4(4) = s*dsqrt(5.0d0/16.0d0) * (-xxrr + yyrr + 7.0d0 * (xxzz - yyzz))
!     42s
      solid4(5) = s*dsqrt(5.0d0/4.0d0) * (7.0d0 * xyzz - xyrr)
!     43c
      solid4(6) = s*dsqrt(35.0d0/8.0d0) * (xxxz - 3.0d0 * xyyz)
!     43s
      solid4(7) = s*dsqrt(35.0d0/8.0d0) * (3.0d0 * xxyz - yyyz)
!     44c
      solid4(8) = s*dsqrt(35.0d0/64.0d0) * (xxxx - 6.0d0 * xxyy + yyyy)
!     44s
      solid4(9) = s*dsqrt(35.0d0/4.0d0) * (xxxy - xyyy)

      end subroutine convert_hex_to_spherical


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

      function random(rmin,rmax)
      implicit none
      real(8) :: rand,rmin,rmax
      real(8) :: random
      call random_number(rand)
      
      random = rmin + (rmax - rmin) * rand

      end function random
