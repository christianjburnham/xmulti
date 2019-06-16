      implicit none
      real(8) :: a_xxx,a_xxy,a_xyy,a_yyy,a_xxz,a_xyz,a_yyz,a_xzz,a_yzz,a_zzz
      real(8) :: b_xxx,b_xxy,b_xyy,b_yyy,b_xxz,b_xyz,b_yyz,b_xzz,b_yzz,b_zzz
      real(8), dimension(7) :: a_solid3,b_solid3
      real(8) :: sum1,sum2
      real(8) :: rmin,rmax
      real(8) :: random
      integer i

      write(*,*) '*********************************************'
      write(*,*) '**  OCT TEST'
      write(*,*) '**'
      write(*,*) '**  short test showing conversion of cartesians'
      write(*,*) '**  to spherical harmonics for octopoles.'
      write(*,*) '**'
      write(*,*) '** Christian J. Burnham, University &
     &College Dublin, 2019'
      write(*,*) '** @christianjburnham@gmail.com'
      write(*,*) '*********************************************'

      call init_random_seed()

      rmin = -1.0d0 
      rmax = 1.0d0 

      a_xxx = random(rmin,rmax)
      a_xxy = random(rmin,rmax)
      a_xyy = random(rmin,rmax)
      a_yyy = random(rmin,rmax)
      a_xxz = random(rmin,rmax)
      a_xyz = random(rmin,rmax)
      a_yyz = random(rmin,rmax)
      a_xzz = random(rmin,rmax)
      a_yzz = random(rmin,rmax)
      a_zzz = random(rmin,rmax)

      b_xxx = random(rmin,rmax)
      b_xxy = random(rmin,rmax)
      b_xyy = random(rmin,rmax)
      b_yyy = random(rmin,rmax)
      b_xxz = random(rmin,rmax)
      b_xyz = random(rmin,rmax)
      b_yyz = random(rmin,rmax)
      b_xzz = random(rmin,rmax)
      b_yzz = random(rmin,rmax)
      b_zzz = random(rmin,rmax)

      write(*,*) 'initial cartesians'
      write(*,*) 'a'
      write(*,17) a_xxx,a_xxy,a_xyy,a_yyy,a_xxz,a_xyz,a_yyz,a_xzz,a_yzz,a_zzz
      write(*,*) 'b'
      write(*,17) b_xxx,b_xxy,b_xyy,b_yyy,b_xxz,b_xyz,b_yyz,b_xzz,b_yzz,b_zzz

      call convert_oct_to_spherical(a_xxx,a_xxy,a_xyy,a_yyy,a_xxz,a_xyz,&
     &     a_yyz,a_xzz,a_yzz,a_zzz,a_solid3)
      call convert_oct_to_cartesian(a_solid3,&
     &           a_xxx,a_xxy,a_xyy,a_yyy,a_xxz,a_xyz,a_yyz,a_xzz,a_yzz,a_zzz)
      call convert_oct_to_spherical(b_xxx,b_xxy,b_xyy,b_yyy,b_xxz,b_xyz,&
     &     b_yyz,b_xzz,b_yzz,b_zzz,b_solid3)
      call convert_oct_to_cartesian(b_solid3,&
     &           b_xxx,b_xxy,b_xyy,b_yyy,b_xxz,b_xyz,b_yyz,b_xzz,b_yzz,b_zzz)

      write(*,*)
      write(*,*) 'detraced cartesians'
      write(*,*) 'a'
      write(*,17) a_xxx,a_xxy,a_xyy,a_yyy,a_xxz,a_xyz,a_yyz,a_xzz,a_yzz,a_zzz
      write(*,*) 'b'
      write(*,17) b_xxx,b_xxy,b_xyy,b_yyy,b_xxz,b_xyz,b_yyz,b_xzz,b_yzz,b_zzz

      call convert_oct_to_spherical(a_xxx,a_xxy,a_xyy,a_yyy,a_xxz,a_xyz,&
     &     a_yyz,a_xzz,a_yzz,a_zzz,a_solid3)
      call convert_oct_to_cartesian(a_solid3,&
     &           a_xxx,a_xxy,a_xyy,a_yyy,a_xxz,a_xyz,a_yyz,a_xzz,a_yzz,a_zzz)
      call convert_oct_to_spherical(b_xxx,b_xxy,b_xyy,b_yyy,b_xxz,b_xyz,&
     &     b_yyz,b_xzz,b_yzz,b_zzz,b_solid3)
      call convert_oct_to_cartesian(b_solid3,&
     &           b_xxx,b_xxy,b_xyy,b_yyy,b_xxz,b_xyz,b_yyz,b_xzz,b_yzz,b_zzz)

      write(*,*)
      write(*,*) 'check: applying the detracing operator a 2nd time.'
      write(*,*) 'a'
      write(*,17) a_xxx,a_xxy,a_xyy,a_yyy,a_xxz,a_xyz,a_yyz,a_xzz,a_yzz,a_zzz
      write(*,*) 'b'
      write(*,17) b_xxx,b_xxy,b_xyy,b_yyy,b_xxz,b_xyz,b_yyz,b_xzz,b_yzz,b_zzz

      write(*,*) 
      write(*,*) 'check: traces'

      write(*,*) 'trace(a)_x',a_xxx + a_xyy + a_xzz
      write(*,*) 'trace(a)_y',a_xxy + a_yyy + a_yzz
      write(*,*) 'trace(a)_z',a_xxz + a_yyz + a_zzz
      write(*,*) 'trace(b)_x',b_xxx + b_xyy + b_xzz
      write(*,*) 'trace(b)_y',b_xxy + b_yyy + b_yzz
      write(*,*) 'trace(b)_z',b_xxz + b_yyz + b_zzz

!     calculate the inner product in two ways.
      
      sum1 = a_xxx * b_xxx &
     &     + a_xxy * b_xxy * 3.0d0 &
     &     + a_xyy * b_xyy * 3.0d0 &
     &     + a_yyy * b_yyy & 
     &     + a_xxz * b_xxz * 3.0d0 &
     &     + a_xyz * b_xyz * 6.0d0 & 
     &     + a_yyz * b_yyz * 3.d00 &
     &     + a_xzz * b_xzz * 3.0d0 & 
     &     + a_yzz * b_yzz * 3.0d0 &
     &     + a_zzz * b_zzz

      sum2 = 0.0d0 
      do i = 1,7
         sum2 = sum2 + a_solid3(i)*b_solid3(i)
      end do 

      write(*,*)
      write(*,*) 'inner products'
      write(*,*) 'calculated in cartesians          = ',sum1
      write(*,*) 'calculated in spherical harmonics = ',sum2


 17   format('xxx',f12.6,'  xxy',f12.6,'  xyy',f12.6&
     &      ,'  yyy',f12.6,'  xxz',f12.6,'  xyz',f12.6&
     &      ,'  yyz',f12.6,'  xzz',f12.6,'  yzz',f12.6,'  zzz',f12.6)


      end

      subroutine convert_oct_to_cartesian(solid3,&
     &           xxx,xxy,xyy,yyy,xxz,xyz,yyz,xzz,yzz,zzz)
      implicit none
      real(8) :: xxx,xxy,xyy,yyy,xxz,xyz,yyz,xzz,yzz,zzz
      real(8), dimension(7) :: solid3
      real(8) :: fac1 = dsqrt(15.0d0)/10.0d0, fac2 = dsqrt(15.0d0)/30.0d0, fac3 = dsqrt(10.0d0)/10.0d0 & 
     &     ,fac4 = dsqrt(6.0d0)/6.0d0, fac5 = dsqrt(15.0d0)/15.0d0, fac6 = dsqrt(2.0d0 / 5.0d0)

!     gives traceless octopole moments in Cartesian form

      xxx = 0.5d0 * solid3(6) - fac1 * solid3(2)
      xxy = 0.5d0 * solid3(7) - fac2 * solid3(3)
      xyy = -0.5d0 * solid3(6) - fac2  * solid3(2)
      yyy = - 0.5d0 * solid3(7) - fac1 * solid3(3)
      xxz = fac4 * solid3(4) - fac3 * solid3(1)
      xyz = fac4 * solid3(5)
      yyz = -fac4 * solid3(4) - fac3 * solid3(1)
      xzz = 2.0d0 * fac5 * solid3(2)
      yzz = 2.0d0 * fac5 * solid3(3)
      zzz = fac6 * solid3(1)

      end subroutine convert_oct_to_cartesian

      subroutine convert_oct_to_spherical(xxx,xxy,xyy,yyy,xxz,xyz,yyz,xzz,yzz,zzz,&
     &     solid3)
      implicit none
      real(8) :: xxx,xxy,xyy,yyy,xxz,xyz,yyz,xzz,yzz,zzz,xrr,yrr,zrr
      real(8), dimension(7) :: solid3
      real(8) :: fac1 = dsqrt(10.0d0)/10.0d0, fac2 = dsqrt(15.0d0)/10.0d0, fac3 = dsqrt(3.0d0/2.0d0) & 
     &     , fac4 = dsqrt(6.0d0)

!     converts octopoles to spherical form
!     factor of 6 because oct_internal = (1/3!) sum_i C_i r_i,alpha r_i,beta r_i,gamma

      xrr = xxx + xyy + xzz
      yrr = xxy + yyy + yzz
      zrr = xxz + yyz + zzz

!     30
      solid3(1) = fac1 * (5 * zzz - 3.0d0 * zrr)
!     31c
      solid3(2) = fac2 * (5.0d0 * xzz - xrr)
!     31s
      solid3(3) = fac2 * (5.0d0 * yzz - yrr)
!     32c
      solid3(4) = fac3 * (xxz - yyz)
!     32s
      solid3(5) = fac4 * xyz
!     33c
      solid3(6) = 0.5d0 * (xxx - 3.0d0 * xyy)
!     33s
      solid3(7) = 0.5d0 * (3.0d0 * xxy - yyy)
      end subroutine convert_oct_to_spherical

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
