      implicit none
      real(8) :: random
      real(8) :: rr(2000,3)
      real(8) :: phi,theta,psi
      character(8),dimension(30) :: molname
      integer, dimension(30) :: nmol
      real(8) :: xdif,ydif,zdif,x,y,z,rdis,rmin
      real(8) :: pi
      integer :: nreject,nmax,nmoltotal,n,imoltype,imol,naccept
      real(8) :: cx,cy,cz
      logical :: accept
      character(len = 8) name
      logical periodic
      integer nmoltypes
      integer ios,nitems
      character(len=2048) :: buffer
      character(len = 2048), dimension(64) :: textlist

!     default
      periodic = .false.
      rmin = 3.0d0 

      open(10,file = 'makestruct.control')

      imoltype = 0
      ios = 0
      
      do while(ios == 0)
         read(10,'(A)',iostat = ios) buffer
         if(ios.lt.0) exit
         call split_text_to_list(buffer,textlist,nitems)      
         if(textlist(1).eq.'PERIODIC') then
            periodic = .true.
         else if(textlist(1).eq.'CVEC') then 
            read(textlist(2),*) cx
            read(textlist(3),*) cy
            read(textlist(4),*) cz
         else if(textlist(1).eq.'RMIN') then 
            read(textlist(2),*) rmin
         else 
            imoltype = imoltype + 1
            if(textlist(1).ne.'   ') then 
               molname(imoltype) = textlist(1)
               read(textlist(2),*) nmol(imoltype)
            endif
         endif
      end do 
      nmoltypes = imoltype-1

      nmoltotal = 0
      do imoltype = 1,nmoltypes
         nmoltotal = nmoltotal + nmol(imoltype)
      end do  

      call init_random_seed()

      pi = 4.0d0 * atan(1.0) 
      nmax = 0 

      write(*,*) nmoltotal
      if(periodic) then
         write(*,*) 'XYZ',cx,0,cy,0,0,cz
         write(*,*) 'RIGID FRAC'
      else
         write(*,*) 
         write(*,*) 'RIGID XYZ'
      endif
      do imoltype = 1,nmoltypes
         do imol = 1,nmol(imoltype)
         
         accept = .false.
         nreject = 0
         do while(.not.accept)
            call random_number(random)
            x = random * cx
            call random_number(random)
            y = random * cy
            call random_number(random)
            z = random * cz

            naccept = 0 

            do n = 1,nmax
               xdif = rr(n,1) - x
               ydif = rr(n,2) - y
               zdif = rr(n,3) - z

               if(periodic) then
                  xdif = xdif - cx * anint(xdif / cx)
                  ydif = ydif - cy * anint(ydif / cy)
                  zdif = zdif - cz * anint(zdif / cz)
               endif
               rdis = sqrt(xdif**2 + ydif**2 + zdif**2) 

               if(rdis.gt.rmin) naccept = naccept + 1
            end do 


            if(naccept.eq.nmax) accept = .true.

            if(.not.accept) then
            nreject = nreject + 1
            if(nreject.gt.1e+6) then
               write(*,*) "ERROR: COULDN'T FIT IN MOLECULES. TRY MAKING BOX BIGGER."
               stop
            endif
            endif
         end do 

         nmax = nmax + 1
         rr(nmax,1) = x
         rr(nmax,2) = y
         rr(nmax,3) = z
         
         call get_random_euler(phi,theta,psi,random)
         
         name = molname(imoltype)

         if(periodic) then 
            write(*,*) name,x/cx,y/cy,z/cz,phi,theta,psi
         else
            write(*,*) name,x,y,z,phi,theta,psi
         endif
      end do 
      end do 

      end




      subroutine split_text_to_list(buffer,textlist,n)
      character(len=2048) buffer,text1,text2,string
      character(len = 2048), dimension(64) :: textlist
      character(len=2048) :: justifyl
      integer i,n

      string = justifyl(buffer)
      i = 0 
      do while(.true.)
         call split_text(string,text1,text2)
         i = i + 1
         textlist(i) = text1
         if(text2.eq.'  ') exit
         string = text2
      end do 
      n = i 
      end

      subroutine split_text(buffer,text1,text2)
      implicit none
      character(len=2048) :: buffer,text2,text1
      integer pos,first_tab,first_space
      character(len=2048) :: justifyl

! Find the first instance of whitespace.  Split label and data.

      first_space = scan(buffer,' ')
      first_tab = scan(buffer,char(9))
      if(first_tab.eq.0) first_tab = 512
      if(first_space.eq.0) first_space = 512

      pos = min(first_tab,first_space)
      text1 = buffer(1:pos)
      text2 = buffer(pos+1:)

      text1 = justifyl(text1)
      text2 = justifyl(text2)

      end subroutine split_text

     FUNCTION JUSTIFYL(STRING) 
      IMPLICIT NONE 
      CHARACTER(LEN = *), INTENT(IN) :: STRING 
      CHARACTER(LEN = LEN(STRING)) JUSTIFYL 
      INTEGER I, L 
      L=LEN(STRING) 
      JUSTIFYL = REPEAT(' ',L)  ! Clear out any rubbish (try leaving this out) 
      I=1 
      DO WHILE ((STRING(I:I) == ' ' .OR. STRING(I:I) == CHAR(9)) &
     &          .AND. I < L)              ! Look for first non-(space or tab) character 
        I=I+1 
      END DO 
      JUSTIFYL(1:L-I+1) = STRING(I:L)       ! Shift left 
      JUSTIFYL(L-I+2:L) = REPEAT(' ',I-1)   ! Replace end with spaces 
      END FUNCTION JUSTIFYL


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
