!     RANDOMEULER
!     generates random euler angles
!     Christian J. Burnham, UCD, 2018

      implicit none
      real(8) :: random
      real(8) :: phi,theta,psi
      integer i
      call init_random_seed()
      

      call get_random_euler(phi,theta,psi,random)
      write(*,*) phi,theta,psi
      end

      subroutine get_random_euler(phi,theta,psi,random)
      implicit none
      real(8) :: random
      real(8) :: x1,y1,z1,x2,y2,z2,x3,y3,z3
      real(8) :: dot,r1,r2
      real(8) :: rotmat(3,3)
      real(8) :: phi,theta,psi,stheta

      call make_rand_vec(random,x1,y1,z1)
      call make_rand_vec(random,x2,y2,z2)

!     project the r2 vector such that it's normal to r1, and in the 12 plane

      dot = x1*x2 + y1*y2 + z1*z2

      x2 = x2 - dot * x1
      y2 = y2 - dot * y1
      z2 = z2 - dot * z1

      r1 = dsqrt(x1**2 + y1**2 + z1**2)
      x1 = x1 / r1
      y1 = y1 / r1
      z1 = z1 / r1

      r2 = dsqrt(x2**2 + y2**2 + z2**2)
      x2 = x2 / r2
      y2 = y2 / r2
      z2 = z2 / r2

!     take cross product to get the components of r3
      x3 = y1*z2 - z1*y2
      y3 = z1*x2 - x1*z2
      z3 = x1*y2 - y1*x2

!     now that we have r1,r2,r3, generate the rotation matrix
      rotmat(1,1) = x1
      rotmat(1,2) = y1
      rotmat(1,3) = z1

      rotmat(2,1) = x2
      rotmat(2,2) = y2
      rotmat(2,3) = z2

      rotmat(3,1) = x3
      rotmat(3,2) = y3
      rotmat(3,3) = z3

!      write(*,*)
!      write(*,*) 'ROTMAT'
!      write(*,*) rotmat(1,1),rotmat(1,2),rotmat(1,3)
!      write(*,*) rotmat(2,1),rotmat(2,2),rotmat(2,3)
!      write(*,*) rotmat(3,1),rotmat(3,2),rotmat(3,3)

      phi = datan2(rotmat(3,1),-rotmat(3,2))
      psi = datan2(rotmat(1,3),rotmat(2,3))
      stheta = rotmat(1,3) / dsin(psi)
      theta = datan2(stheta,rotmat(3,3))

!     check

!      call get_rotmat(phi,theta,psi,rotmat)

!      write(*,*)
!      write(*,*) 'ROTMAT'
!      write(*,*) rotmat(1,1),rotmat(1,2),rotmat(1,3)
!      write(*,*) rotmat(2,1),rotmat(2,2),rotmat(2,3)
!      write(*,*) rotmat(3,1),rotmat(3,2),rotmat(3,3)


      end subroutine get_random_euler


      subroutine make_rand_vec(random,x,y,z)
      implicit none
      real(8) :: random
      real(8) :: pi
      real(8) :: theta
      real(8) :: x,y,z
      pi = 4.0d0 * atan(1.0d0)

      call random_number(random)
      theta = random * 2.0d0 * pi
      call random_number(random)
      z = (random-0.5d0)*2.0d0 

      x = dsqrt(1.0d0-z**2)*dcos(theta)
      y = dsqrt(1.0d0-z**2)*dsin(theta)

      end subroutine make_rand_vec

      subroutine init_random_seed()
      integer :: i, n, clock
      integer, dimension(:), allocatable :: seed

      call random_seed(size = n)
      allocate(seed(n))

      call system_clock(count=clock)

      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      call random_seed(put = seed)

      deallocate(seed)
      end subroutine

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
