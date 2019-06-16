      implicit none
      real(8) :: a_xx,a_xy,a_xz,a_yy,a_yz,a_zz
      real(8) :: b_xx,b_xy,b_xz,b_yy,b_yz,b_zz
      real(8), dimension(5) :: a_solid2,b_solid2
      real(8) :: sum1,sum2
      real(8) :: rmin,rmax
      real(8) :: random
      integer i

      write(*,*) '*********************************************'
      write(*,*) '**  QUAD TEST'
      write(*,*) '**'
      write(*,*) '**  short test showing conversion of cartesians'
      write(*,*) '**  to spherical harmonics for quadrupoles.'
      write(*,*) '**'
      write(*,*) '** Christian J. Burnham, University &
     &College Dublin, 2019'
      write(*,*) '** @christianjburnham@gmail.com'
      write(*,*) '*********************************************'

      call init_random_seed()

      rmin = -1.0d0 
      rmax = 1.0d0 

      a_xx = random(rmin,rmax)
      a_xy = random(rmin,rmax)
      a_xz = random(rmin,rmax)
      a_yy = random(rmin,rmax)
      a_yz = random(rmin,rmax)
      a_zz = random(rmin,rmax)

      b_xx = random(rmin,rmax)
      b_xy = random(rmin,rmax)
      b_xz = random(rmin,rmax)
      b_yy = random(rmin,rmax)
      b_yz = random(rmin,rmax)
      b_zz = random(rmin,rmax)

      write(*,*) 'initial cartesians'
      write(*,*) 'a'
      write(*,17) a_xx,a_xy,a_xz,a_yy,a_yz,a_zz
      write(*,*) 'b'
      write(*,17) b_xx,b_xy,b_xz,b_yy,b_yz,b_zz

      call convert_quad_to_spherical(a_xx,a_xy,a_xz,a_yy,a_yz,a_zz,a_solid2)
      call convert_quad_to_cartesian(a_solid2,a_xx,a_xy,a_xz,a_yy,a_yz,a_zz)

      call convert_quad_to_spherical(b_xx,b_xy,b_xz,b_yy,b_yz,b_zz,b_solid2)
      call convert_quad_to_cartesian(b_solid2,b_xx,b_xy,b_xz,b_yy,b_yz,b_zz)

      write(*,*)
      write(*,*) 'detraced cartesians'
      write(*,*) 'a'
      write(*,17) a_xx,a_xy,a_xz,a_yy,a_yz,a_zz
      write(*,*) 'b'
      write(*,17) b_xx,b_xy,b_xz,b_yy,b_yz,b_zz

      call convert_quad_to_spherical(a_xx,a_xy,a_xz,a_yy,a_yz,a_zz,a_solid2)
      call convert_quad_to_cartesian(a_solid2,a_xx,a_xy,a_xz,a_yy,a_yz,a_zz)

      call convert_quad_to_spherical(b_xx,b_xy,b_xz,b_yy,b_yz,b_zz,b_solid2)
      call convert_quad_to_cartesian(b_solid2,b_xx,b_xy,b_xz,b_yy,b_yz,b_zz)

      write(*,*)
      write(*,*) 'check: applying the detracing operator a 2nd time.'
      write(*,*) 'a'
      write(*,17) a_xx,a_xy,a_xz,a_yy,a_yz,a_zz
      write(*,*) 'b'
      write(*,17) b_xx,b_xy,b_xz,b_yy,b_yz,b_zz

      write(*,*)
      write(*,*) 'check: trace'
      write(*,*) 'trace(a)',a_xx + a_yy + a_zz
      write(*,*) 'trace(b)',b_xx + b_yy + b_zz

!     calculate the inner product in two ways.

      sum1 = a_xx * b_xx + a_yy * b_yy + a_zz * b_zz &
     &     + 2.0d0 * (a_xy * b_xy + a_xz * b_xz + a_yz * b_yz)

      sum2 = 0.0d0 
      do i = 1,5
         sum2 = sum2 + a_solid2(i)*b_solid2(i)
      end do 

      write(*,*)
      write(*,*) 'inner products'
      write(*,*) 'calculated in cartesians          = ',sum1
      write(*,*) 'calculated in spherical harmonics = ',sum2

 17   format('xx',f12.6,'  xy',f12.6,'  xz',f12.6,'  yy',f12.6,'  yz',f12.6,'  zz',f12.6)

      end

      subroutine convert_quad_to_spherical(xx,xy,xz,yy,yz,zz,solid2)
      implicit none
      real(8) :: xx,xy,xz,yy,yz,zz,rr
      real(8), dimension(5) :: solid2
      real(8) :: fac1 = dsqrt(6.0d0)/6.0d0, fac2 = dsqrt(2.0d0), fac3 = dsqrt(2.0d0)/2.0d0

!     converts quadrupoles to spherical form

      rr = xx + yy + zz

!     20
      solid2(1) = fac1 * (3.0d0 * zz - rr)
!     21c
      solid2(2) = fac2 * xz
!     21s
      solid2(3) = fac2 * yz
!     22c
      solid2(4) = fac3 * (xx - yy)
!     22s
      solid2(5) = fac2 * xy

      end subroutine convert_quad_to_spherical

      subroutine convert_quad_to_cartesian(solid2,xx,xy,xz,yy,yz,zz)
      implicit none
      real(8) :: xx,xy,xz,yy,yz,zz
      real(8), dimension(5) :: solid2
      real(8) :: fac1 = dsqrt(6.0d0)/6.0d0, fac2 = dsqrt(2.0d0)/2.0d0, fac3 = dsqrt(6.0d0)/3.0d0

!     gives traceless quadrupole moments in Cartesian form

      xx = -fac1 * solid2(1) + fac2 * solid2(4)
      yy = -fac1 * solid2(1) - fac2 * solid2(4)
      zz = solid2(1) *fac3

      xy = fac2 * solid2(5)
      xz = fac2 * solid2(2)
      yz = fac2 * solid2(3)

      end subroutine convert_quad_to_cartesian


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
